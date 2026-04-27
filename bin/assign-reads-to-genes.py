#!/usr/bin/env python
# coding: utf-8


import argparse
import logging

import pysam
import pandas as pd

logging.basicConfig(level=logging.DEBUG, format='%(asctime)s - %(levelname)s: %(message)s')


def parse_attribute(attribute_series: pd.Series, attribute_name: str) -> pd.Series:
    """
    Parse the attributes column of a (GENCODE/RefSeq) GTF file.

    Input:
    * a [str]: the attributes element (column 9 of the GTF file)
    * regex [str]: a regular expression that will be iteratively applied to the attribute string to capture attribute key, val pairs. Default should work for GENCODE/RefSeq
    """
    if not isinstance(attribute_series, pd.Series):
        raise TypeError('attribute_series must be a pandas Series')
    if not isinstance(attribute_name, str):
        raise TypeError('attribute_name must be a string')

    return attribute_series.str.extract(f'{attribute_name} "(.*?)"')


def gtf_to_df(gtf: str, parse_attributes: list=None) -> pd.DataFrame:
    df = pd.read_csv(gtf, sep='\t', low_memory=False, header=None, names=['chrom', 'source', 'feature', 'start', 'end', 'score', 'strand', 'frame', 'attributes'], comment='#')
    if parse_attributes is not None:
        for a in parse_attributes:
            df[a] = parse_attribute(df.attributes, a)
    return df


def block_overlaps_block(block_1, block_2):
    block_1_start, block_1_end = block_1
    block_2_start, block_2_end = block_2
    # if one block ends before the other starts, they don't overlap
    if block_1_end < block_2_start or block_2_end < block_1_start:
        return False
    else:
        return True


def blocks_overlap(read_blocks, gene_blocks):
    for rb in read_blocks:
        for gb in gene_blocks:
            if block_overlaps_block(rb, gb):
                return True
    return False


# use the BAM file gene assignment (which correlates strongly with starsolo's GeneFull)
# load up the GTF file (exons only, and full gene body)
# and for each read, see if the blocks overlap with (1) full gene body and (2) exons
# then output a table with all of that

parser = argparse.ArgumentParser()
parser.add_argument('--strandedness', required=False, default='unstranded', help='Strandedness (unstranded, forward, or reverse)')
parser.add_argument('--gtf', required=True, help='GTF file with gene annotations')
parser.add_argument('--bam', required=True, help='BAM file with reads')
# args = parser.parse_args(['--gtf', '/scratch/scjp_root/scjp0/porchard/2023-HSM-ONT/data/ref-data/genes/genes.gtf', '--bam', '/scratch/scjp_root/scjp0/porchard/2023-HSM-ONT/work/fix-swaps/9266-VD-1/tagged.bam'])
args = parser.parse_args()

assert(args.strandedness in ['unstranded', 'forward', 'reverse'])


# BAM should be sorted
# need to sort genes by start position
# then can sweep along the chromosome


gtf_df = gtf_to_df(args.gtf, parse_attributes=['gene_id'])

gene_bodies = {} # chrom --> [(gene_id, start, end, strand), ...]
gene_bodies = {chrom: [(row['gene_id'], row['start'], row['end'], row['strand']) for i, row in df.sort_values(['start', 'end']).iterrows()] for chrom, df in gtf_df[gtf_df.feature=='gene'].groupby('chrom')}

gene_to_gene_body = {}
for chrom in gene_bodies:
    for gene_id, start, end, strand in gene_bodies[chrom]:
        gene_to_gene_body[gene_id] = (start, end, strand)


exons = {} # chrom --> gene_id --> [(start, end), (start, end), ...]
for (chrom, gene_id), df in gtf_df[gtf_df.feature=='exon'].groupby(['chrom', 'gene_id']):
    if chrom not in exons:
        exons[chrom] = {}
    df_sorted = df.sort_values(['start', 'end'])
    exons[chrom][gene_id] = list(zip(df_sorted.start, df_sorted.end))


gene_body_index = {chrom: 0 for chrom in gene_bodies}

with pysam.AlignmentFile(args.bam, 'rb') as bam:
    count = 0
    for read in bam.fetch(until_eof=True):
        count += 1
        if count % 1000000 == 0:
            logging.info('Processed {:,} reads'.format(count))
        # output, for each read:
        # read_name, gene_name, overlaps_exon
        read_gene_assignment = '-'
        read_overlaps_gene_exon = False
        #if read.is_secondary or read.is_supplementary:
        #    continue

        chrom = read.reference_name
        strand = '-' if read.is_reverse else '+'
        if chrom not in gene_bodies:
            print(f'{read.query_name}\t{read_gene_assignment}\t{read_overlaps_gene_exon}')
            continue
        while (gene_body_index[chrom] + 1) < len(gene_bodies[chrom]) and gene_bodies[chrom][gene_body_index[chrom]][2] < read.reference_start:
            gene_body_index[chrom] += 1
        if gene_bodies[chrom][gene_body_index[chrom]][2] < read.reference_start:
            # no more genes left on this chromosome
            print(f'{read.query_name}\t{read_gene_assignment}\t{read_overlaps_gene_exon}')
            continue
        # get all genes that could possibly overlap with this read
        genes = []
        for gene in gene_bodies[chrom][gene_body_index[chrom]:]:
            if gene[1] > read.reference_end:
                break
            genes.append(gene[0])
        
        # so now we have our list of gene candidates
        # check for overlap
        # if there is overlap with > 1 gene, favor exon overlap
        overlaps_gene_body = []
        overlaps_exons = []

        for gene in genes:
            gene_body_block = gene_to_gene_body[gene][0:2] # [gene start, gene end]
            exon_blocks = exons[chrom][gene]
            matches_strand = args.strandedness == 'unstranded' or (args.strandedness == 'forward' and strand == gene_to_gene_body[gene][2]) or (args.strandedness == 'reverse' and strand != gene_to_gene_body[gene][2])
            overlaps_gene_body.append(matches_strand and blocks_overlap(read.get_blocks(), [gene_body_block]))
            overlaps_exons.append(matches_strand and blocks_overlap(read.get_blocks(), exon_blocks))
        
        if sum(overlaps_gene_body) == 0:
            pass
        elif sum(overlaps_gene_body) == 1:
            read_gene_assignment = genes[overlaps_gene_body.index(True)]
            read_overlaps_gene_exon = overlaps_exons[overlaps_gene_body.index(True)]
        elif sum(overlaps_gene_body) > 1:
            if sum(overlaps_exons) == 1:
                read_gene_assignment = genes[overlaps_exons.index(True)]
                read_overlaps_gene_exon = True
        
        print(f'{read.query_name}\t{read_gene_assignment}\t{read_overlaps_gene_exon}')


logging.info('Finished processing reads')
