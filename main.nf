#!/usr/bin/env nextflow

nextflow.enable.dsl=2

CONFIG = params.config
ADAPTERS = params.adapters
PRIMERS = params.primers
BARCODES = params.barcodes
FASTA = params.fasta
GTF = params.gtf


def parse_config(config_path) {
    Channel.fromPath(config_path)
        .splitCsv(sep: '\t', header: true)
        .map { row -> [row.library, row.readgroup, file(row.bam)] }
}


process segment {

    container "docker://porchard/skera:1.4.0"
    memory { 7.GB * task.attempt }
    time '24h'
    maxRetries 3
    errorStrategy {task.attempt <= maxRetries ? 'retry' : 'ignore'}
    cache 'lenient'
    tag "${library}"
    cpus 10

    input:
    tuple val(library), val(readgroup), path(bam), path(adapters)

    output:
    tuple val(library), val(readgroup), path("${library}.skera.bam"), emit: bam
    tuple val(library), val(readgroup), path("${library}.*")

    """
    skera split --log-level INFO --num-threads $task.cpus $bam $adapters ${library}.skera.bam
    """

}


process remove_primers {

    container 'docker://porchard/isoseq:4.3.0'
    memory '4 GB'
    cpus 10
    time '24h'

    input:
    tuple val(library), val(readgroup), path(bam), path(primers)

    output:
    tuple val(library), val(readgroup), path("${library}.fl.5p--3p.bam")

    """
    lima --num-threads 10 $bam $primers ${library}.fl.bam --isoseq
    """

}


process tag {

    container 'docker://porchard/isoseq:4.3.0'
    memory '4 GB'
    cpus 10
    time '5h'

    input:
    tuple val(library), val(readgroup), path(bam)

    output:
    tuple val(library), val(readgroup), path("${library}.flt.bam")

    """
    isoseq tag --num-threads 10 $bam ${library}.flt.bam --design T-12U-16B
    """

}


// isoseq tag extracts barcodes as the reverse complement of the 10x whitelist.
// RC them so that barcodes match the standard whitelist, ensuring consistency with matched ATAC data.
process rc_tags {

    container 'library://porchard/default/general:20220107'
    memory '4 GB'
    time '5h'
    tag "${library}"

    input:
    tuple val(library), val(readgroup), path(bam)

    output:
    tuple val(library), val(readgroup), path("${library}.rc.bam")

    """
    #!/usr/bin/env python3
    import pysam

    COMP = str.maketrans('ACGTNacgtn', 'TGCANtgcan')

    def rc(seq):
        return seq.translate(COMP)[::-1]

    with pysam.AlignmentFile('${bam}', check_sq=False) as bam_in, \
         pysam.AlignmentFile('${library}.rc.bam', 'wb', header=bam_in.header) as bam_out:
        for read in bam_in:
            if read.has_tag('XC'):
                read.set_tag('XC', rc(read.get_tag('XC')))
            if read.has_tag('XM'):
                read.set_tag('XM', rc(read.get_tag('XM')))
            bam_out.write(read)
    """

}


process refine {

    container 'docker://porchard/isoseq:4.3.0'
    cpus 10
    time '24h'
    memory '5 GB'

    input:
    tuple val(library), val(readgroup), path(bam), path(primers)

    output:
    tuple val(library), val(readgroup), path("${library}.${readgroup}.fltnc.bam")

    """
    isoseq refine --num-threads 10 $bam $primers ${library}.${readgroup}.fltnc.bam --require-polya
    """

}


process merge_readgroups {
    
    container 'library://porchard/default/general:20220107'
    cpus 10
    time '24h'
    memory '10 GB'

    input:
    tuple val(library), val(readgroups), path(bams)

    output:
    tuple val(library), path("${library}.bam")

    script:
    if (bams.size() == 1) {
        """
        ln -s ${bams[0]} ${library}.bam
        """
    } else {
        """
        samtools merge -@ ${task.cpus-1} -o ${library}.bam ${bams.join(' ')}
        """
    }

}


process correct_barcodes {

    container 'docker://porchard/isoseq:4.3.0'
    cpus 10
    time '24h'
    memory '20 GB'
    label 'largemem'

    input:
    tuple val(library), path(bam), path(barcodes)

    output:
    tuple val(library), path("${library}.corrected.bam")

    """
    isoseq correct --num-threads 10 --barcodes $barcodes $bam ${library}.corrected.bam
    """

}


// isoseq correct leaves uncorrected barcodes in CB for reads that fail correction (gp=0).
// Reset those to '-' so they don't get treated as valid cells downstream.
process reset_failed_barcodes {

    container 'library://porchard/default/general:20220107'
    memory '4 GB'
    time '5h'
    tag "${library}"

    input:
    tuple val(library), path(bam)

    output:
    tuple val(library), path("${library}.reset.bam")

    """
    #!/usr/bin/env python3
    import pysam

    with pysam.AlignmentFile('${bam}', check_sq=False) as bam_in, \
         pysam.AlignmentFile('${library}.reset.bam', 'wb', header=bam_in.header) as bam_out:
        for read in bam_in:
            if read.has_tag('gp') and read.get_tag('gp') == 0:
                read.set_tag('CB', '-')
            bam_out.write(read)
    """

}


process make_index {

    container 'docker://porchard/isoseq:4.3.0'
    cpus 5
    time '24h'
    memory '30 GB'
    label 'largemem'

    input:
    path(fasta)

    output:
    path('genome.mmi')

    """
    pbmm2 index --preset ISOSEQ --num-threads 5 $fasta genome.mmi
    """

}


process align {

    container 'docker://porchard/isoseq:4.3.0'
    cpus 20
    time '48h'
    memory '20 GB'

    input:
    tuple val(library), path(bam), path(index)

    output:
    tuple val(library), path("${library}.bam")

    """
    pbmm2 align --preset ISOSEQ --unmapped --sort --num-threads 20 --bam-index BAI $index $bam ${library}.bam
    """

}


process index_bam {

    container "library://porchard/default/general:20220107"
    memory '5 GB'
    time '24h'
    cache 'lenient'
    tag "${library}"

    input:
    tuple val(library), path(bam)

    output:
    tuple val(library), path(bam), path("${bam}.bai")

    """
    samtools index $bam
    """

}


process make_gtfdb {

    memory '10 GB'
    container 'docker://porchard/isoquant:20240911'
    time '24h'

    input:
    tuple val(x), path(gtf)

    output:
    path("${x}.db")

    """
    #!/usr/bin/env python

    import re
    import gffutils

    def fix_duplicate_ids(input_gtf, output_gtf):
        with open(input_gtf) as f_in, open(output_gtf, 'w') as f_out:
            for line in f_in:
                if line.startswith('#'):
                    f_out.write(line)
                    continue
                fields = line.rstrip('\\n').split('\\t')
                if fields[2] == 'transcript':
                    attrs = fields[8]
                    gene_id = re.search(r'gene_id "([^"]+)"', attrs)
                    transcript_id = re.search(r'transcript_id "([^"]+)"', attrs)
                    if gene_id and transcript_id and gene_id.group(1) == transcript_id.group(1):
                        new_tid = transcript_id.group(1) + '_transcript'
                        attrs = attrs.replace(
                            'transcript_id "' + transcript_id.group(1) + '"',
                            'transcript_id "' + new_tid + '"'
                        )
                        fields[8] = attrs
                f_out.write('\\t'.join(fields) + '\\n')

    fix_duplicate_ids('${gtf}', 'fixed.gtf')
    db = gffutils.create_db('fixed.gtf', '${x}.db', disable_infer_genes=True, disable_infer_transcripts=True)
    """

}


process isoquant {

    cpus 20
    memory { 40.GB * task.attempt }
    publishDir "${params.results}/isoquant"
    container 'docker://porchard/isoquant:20240911'
    time '24h'
    label 'largemem'
    maxRetries 1
    errorStrategy 'retry'

    input:
    tuple val(library), path(bam), path(bam_index), path(fasta), path(gtfdb)

    output:
    tuple val(library), path("${library}/${library}.read_assignments.tsv.gz")

    """
    # HOME is reset to current working directory because isoquant tries to create a directory in and write in $HOME, which results in an error when run in a read-only singularity container
    mkdir -p genedb
    export HOME=. && isoquant.py --output . --bam $bam --data_type pacbio --reference $fasta --no_model_construction --complete_genedb --genedb $gtfdb --stranded none --labels $library --bam_tags CR,UR --threads ${task.cpus - 1} --prefix $library --genedb_output genedb/
    """

}


process assign_reads_to_transcripts {

    memory '100 GB'
    container 'library://porchard/default/general:20220107'
    time '5h'
    label 'largemem'

    input:
    tuple val(library), path(isoquant)

    output:
    tuple val(library), path("${library}.transcript-assignments.txt")

    """
    isoquant-read-assignments-to-transcript-assignment.py $isoquant > ${library}.transcript-assignments.txt
    """
    
}


process assign_reads_to_genes {

    memory '7 GB'
    container 'library://porchard/default/general:20220107'
    time '10h'

    input:
    tuple val(library), path(bam), path(bam_index), path(gtf)

    output:
    tuple val(library), path(bam), path("${library}.gene-assignments.txt")

    """
    assign-reads-to-genes.py --strandedness forward --gtf $gtf --bam $bam > ${library}.gene-assignments.txt
    """

}


process add_tags {

    memory '30 GB'
    container 'library://porchard/default/general:20220107'
    time '10h'
    label 'largemem'

    input:
    tuple val(library), path('in.bam'), path(gene_assignments), path(transcript_assignments), path(gtf)

    output:
    tuple val(library), path("${library}.bam")

    """
    add-tags.py --gtf $gtf --transcript-assignments $transcript_assignments --gene-assignments $gene_assignments --bam in.bam --output-bam ${library}.bam
    """

}

process sort_by_cb {

    memory '10 GB'
    container 'library://porchard/default/general:20220107'
    time '10h'
    cpus 10

    input:
    tuple val(library), path('in.bam')

    output:
    tuple val(library), path("${library}.bam")

    """
    samtools sort -@ ${task.cpus-1} -t CB -o ${library}.bam in.bam
    """

}

process correct_umis {

    memory '7 GB'
    time '24h'
    container 'docker://ontresearch/wf-single-cell:sha0fcdf10929fbef2d426bb985e16b81153a88c6f4'
    maxRetries 1
    errorStrategy 'retry'

    input:
    tuple val(library), path('in.bam')

    output:
    tuple val(library), path("${library}.bam")

    """
    correct-umis.py in.bam ${library}.bam
    """

}


process sort_bam {

    memory '10 GB'
    publishDir "${params.results}/bam"
    container 'library://porchard/default/general:20220107'
    time '10h'
    cpus 10

    input:
    tuple val(library), path('in.bam')

    output:
    tuple val(library), path("${library}.bam")

    """
    samtools sort -@ ${task.cpus-1} -o ${library}.bam in.bam
    """

}

process make_count_matrices {

    memory { 50.GB * task.attempt }
    publishDir "${params.results}/count-matrices"
    time '24h'
    label 'largemem'
    container 'library://porchard/default/general:20220107'
    maxRetries 2
    errorStrategy 'retry'

    input:
    tuple val(library), path(bam), path(gtf)

    output:
    tuple val(library), path("${library}.genes.matrix.mtx"), path("${library}.genes.features.tsv"), path("${library}.genes.barcodes.tsv"), emit: gene_matrices
    tuple val(library), path("${library}.transcripts.matrix.mtx"), path("${library}.transcripts.features.tsv"), path("${library}.transcripts.barcodes.tsv"), emit: transcript_matrices

    """
    make-count-matrices.py --bam $bam --gtf $gtf --prefix ${library}.
    """

}

process calculate_qc_metrics {

    memory { 100.GB * task.attempt }
    publishDir "${params.results}/qc"
    time '24h'
    label 'largemem'
    container 'library://porchard/default/general:20220107'
    maxRetries 1
    errorStrategy 'retry'

    input:
    tuple val(library), path(bam)

    output:
    tuple val(library), path("${library}.qc.txt")

    """
    calculate-qc-metrics.py $bam > ${library}.qc.txt
    """

}


process trim_for_scafe {

    publishDir "${params.results}/preprocess-for-scafe"
    container 'docker://porchard/general:20220406125608'
    time '10h'
    memory '5 GB'

    input:
    tuple val(library), path(bam)

    output:
    tuple val(library), path("${library}.trimmed.bam"), path("${library}.trimmed.bam.bai")

    """
    trim-for-scafe.py --input-bam $bam --output-bam ${library}.trimmed.unsorted.bam --trim-to 100 --max-softclipping 1
    samtools sort -m 3G -o ${library}.trimmed.bam ${library}.trimmed.unsorted.bam
    samtools index ${library}.trimmed.bam
    """

}


process filter_5prime_ends_for_scafe {

    publishDir "${params.results}/preprocess-for-scafe"
    container 'docker://porchard/general:20220406125608'
    memory '75 GB'
    time '10h'
    label 'largemem'

    input:
    tuple val(library), path(bam), path(bam_index)

    output:
    tuple val(library), path("${library}.trimmed-and-filtered.bam"), path("${library}.trimmed-and-filtered.bam.bai")
    path("*.png")

    """
    filter-bam-to-most-supported-5prime-ends.py --bam-in $bam --bam-out ${library}.trimmed-and-filtered.bam --prefix ${library}.
    samtools index ${library}.trimmed-and-filtered.bam
    """

}


process cellbender {

    memory '40 GB'
    publishDir "${params.results}/cellbender"
    container 'docker://porchard/cellbender:0.3.0'
    time '24h'
    label 'gpu'

    input:
    tuple val(library), path('matrix.mtx'), path('genes.tsv'), path('barcodes.tsv'), val(feature_type)

    output:
    path("${library}*")
    path("${library}*.h5"), emit: h5_files

    """
    cellbender remove-background --cuda --epochs 150 --fpr 0.01 0.05 0.1 --input . --output ./${library}.${feature_type}.cellbender.h5
    cp .command.log ${library}.${feature_type}.log
    """

}


process make_chrom_sizes {

    container 'library://porchard/default/general:20220107'
    memory '1 GB'
    time '1h'

    input:
    path(fasta)

    output:
    path('chrom.sizes')

    """
    samtools faidx $fasta
    cut -f1,2 ${fasta}.fai > chrom.sizes
    """

}


process make_bedgraphs {

    memory '20 GB'
    time '24h'
    container 'library://porchard/default/general:20220107'
    tag "${library}"

    input:
    tuple val(library), path(bam)

    output:
    path("${library}.*.bdg")

    """
    make-scaled-and-stranded-bedgraphs.py --rev-strand-negative --scale-per-reads 1000000 $bam ${library}
    """

}


process make_bigwigs {

    publishDir "${params.results}/bigwig"
    memory '20 GB'
    time '24h'
    container 'library://porchard/default/general:20220107'

    input:
    path(chrom_sizes)
    each path(bedgraph)

    output:
    path("${prefix}.bw")

    script:
    prefix = bedgraph.getName().replaceAll('.bdg', '')

    """
    LC_COLLATE=C sort -k1,1 -k2n,2 $bedgraph | grep -w -P -e 'chr[\\dMXY]+' > sorted.bedgraph
    bedClip sorted.bedgraph $chrom_sizes clipped.bedgraph
    bedGraphToBigWig clipped.bedgraph $chrom_sizes ${prefix}.bw
    rm sorted.bedgraph clipped.bedgraph
    """

}


workflow {

    bams = parse_config(CONFIG)
    adapters = Channel.fromPath(ADAPTERS)
    primers = Channel.fromPath(PRIMERS)
    barcodes = Channel.fromPath(BARCODES)
    fasta = Channel.fromPath(FASTA)
    gtf = Channel.fromPath(GTF)

    // make some index / reference files
    gtfdb = make_gtfdb(gtf.map({it -> [it.getName().tokenize('.')[0], it]}))
    chrom_sizes = make_chrom_sizes(fasta)
    mm_index = make_index(fasta)

    // main pipeline
    sreads = segment(bams.combine(adapters)).bam
    tagged = remove_primers(sreads.combine(primers)) | tag | rc_tags
    refined = refine(tagged.combine(primers))
    merged = merge_readgroups(refined.groupTuple(by: 0))
    corrected = correct_barcodes(merged.combine(barcodes)) | reset_failed_barcodes

    aligned = align(corrected.combine(mm_index)) | index_bam

    transcript_assignments = aligned.combine(fasta).combine(gtfdb) | isoquant | assign_reads_to_transcripts
    gene_assignments = aligned.combine(gtf) | assign_reads_to_genes

    processed_bams = gene_assignments.combine(transcript_assignments, by: 0).combine(gtf) | add_tags | sort_by_cb | correct_umis | sort_bam
    count_matrices = make_count_matrices(processed_bams.combine(gtf))
    calculate_qc_metrics(processed_bams)

    count_matrices.gene_matrices.map({it -> it + ['genes']}).mix(count_matrices.transcript_matrices.map({it -> it + ['transcripts']})) | cellbender

    make_bigwigs(chrom_sizes, make_bedgraphs(processed_bams).flatten())

    trim_for_scafe(processed_bams) | filter_5prime_ends_for_scafe
    
}
