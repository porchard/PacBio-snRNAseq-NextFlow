#!/usr/bin/env python
# coding: utf-8

import pysam
import pandas as pd
import csv
import logging
import argparse

logging.basicConfig(level=logging.DEBUG, format='%(asctime)s - %(levelname)s: %(message)s')


def parse_attribute(attribute_series: pd.Series, attribute_name: str) -> pd.Series:
    """
    Parse the attributes column of a (GENCODE/RefSeq) GTF file.

    >>> attributes = pd.Series(['gene_id "ENSG00000223972.5"; transcript_id "ENSG00000223972.5";', 'gene_id "ENSG00000227232.5"; transcript_id "ENSG00000227232.5";'])
    >>> parse_attribute(attributes, 'gene_id').to_list()
    ['ENSG00000223972.5', 'ENSG00000227232.5']
    """
    if not isinstance(attribute_series, pd.Series):
        raise TypeError('attribute_series must be a pandas Series')
    if not isinstance(attribute_name, str):
        raise TypeError('attribute_name must be a string')
    
    return attribute_series.str.extract(f'{attribute_name} "(.*?)"')[0]


def gtf_to_df(gtf: str, parse_attributes: list=None) -> pd.DataFrame:
    """
    Load a GTF file into a dataframe, parsing requested attributes.

    >>> df = gtf_to_df('test/gencode.v30.GRCh38.ERCC.genes.collapsed_only.gtf.gz', parse_attributes=['gene_id'])
    >>> list(df.columns)
    ['chrom', 'source', 'feature', 'start', 'end', 'score', 'strand', 'frame', 'attributes', 'gene_id']
    >>> df.gene_id.to_list()[0]
    'ENSG00000223972.5'

    """
    df = pd.read_csv(gtf, sep='\t', header=None, names=['chrom', 'source', 'feature', 'start', 'end', 'score', 'strand', 'frame', 'attributes'], comment='#')
    if parse_attributes is not None:
        for a in parse_attributes:
            df[a] = parse_attribute(df.attributes, a)
    return df


def load_transcript_assignments(f):
    transcript_assignments = dict()
    with open(f, 'r') as fh:
        reader = csv.DictReader(fh, delimiter='\t')
        for line in reader:
            if line['isoform_id'] != '-':
                assert(line['read_id'] not in transcript_assignments)
                transcript_assignments[line['read_id']] = line['isoform_id']
    return transcript_assignments


# load up the read_id --> gene, transcript assignments from isoquant
# GTF = '/scratch/scjp_root/scjp0/porchard/PacBio-pipeline-development/data/ref-data/genes/genes.gtf'
# TRANSCRIPT_ASSIGNMENTS = '/scratch/scjp_root/scjp0/porchard/PacBio-pipeline-development/work/post-align/results/assign_reads_to_transcripts/15245-AH-1.transcript-assignments.txt'
# GENE_ASSIGNMENTS = '/scratch/scjp_root/scjp0/porchard/PacBio-pipeline-development/work/post-align/results/assign_reads_to_genes/15245-AH-1.gene-assignments.txt'
# BAM_IN = '/scratch/scjp_root/scjp0/porchard/PacBio-pipeline-development/work/post-align/data/15245-AH-1.bam'
# BAM_OUT = 'test.bam'

parser = argparse.ArgumentParser(description='Add gene and transcript tags to a BAM file')
parser.add_argument('--gtf', required=True, help='GTF file used for gene and transcript definitions')
parser.add_argument('--transcript-assignments', required=True, help='Output of assign_reads_to_transcripts.py')
parser.add_argument('--gene-assignments', required=True, help='Output of assign_reads_to_genes.py')
parser.add_argument('--bam', required=True, help='Input BAM file')
parser.add_argument('--output-bam', required=True, help='Output BAM file with gene and transcript tags added')
args = parser.parse_args()

GTF = args.gtf
TRANSCRIPT_ASSIGNMENTS = args.transcript_assignments
GENE_ASSIGNMENTS = args.gene_assignments
BAM_IN = args.bam
BAM_OUT = args.output_bam



gtf_df = gtf_to_df(GTF, parse_attributes=['gene_id', 'transcript_id'])
gtf_df = gtf_df[gtf_df.feature=='transcript']
transcript_id_to_gene_id = dict(zip(gtf_df.transcript_id, gtf_df.gene_id))

transcript_assignments = load_transcript_assignments(TRANSCRIPT_ASSIGNMENTS)

with pysam.AlignmentFile(BAM_IN, 'rb') as bam_in:
    with pysam.AlignmentFile(BAM_OUT, 'wb', template=bam_in) as bam_out:
        with open(GENE_ASSIGNMENTS, 'r') as gene_assignments:
            for count, read in enumerate(bam_in.fetch(until_eof=True)):
                if count % 1000000 == 0:
                    logging.info('Processed {:,} reads'.format(count))
                read_id, gene_id, overlaps_exon = gene_assignments.readline().rstrip().split('\t')
                assert(read.query_name == read_id)
                assert(not read.has_tag('GX'))
                assert(not read.has_tag('OE'))
                assert(not read.has_tag('TX'))
                read.set_tag('GX', gene_id)
                read.set_tag('OE', overlaps_exon)
                if not read.is_secondary and not read.is_supplementary:
                    if read.query_name in transcript_assignments and transcript_id_to_gene_id[transcript_assignments[read.query_name]] == gene_id:
                        read.set_tag('TX', transcript_assignments[read.query_name])
                bam_out.write(read)
