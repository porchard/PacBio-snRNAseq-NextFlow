#!/usr/bin/env python
# coding: utf-8

import argparse
import logging

import pysam
import pandas as pd

logging.basicConfig(level=logging.DEBUG, format='%(asctime)s - %(levelname)s: %(message)s')

    
def write_mm(prefix, df, features_col, barcodes_col, count_col, all_features=None, all_barcodes=None):
    features = df[features_col].unique() if all_features is None else all_features
    barcodes = df[barcodes_col].unique() if all_barcodes is None else all_barcodes
    if all_features is not None:
        assert(df[features_col].isin(all_features).all())
    if all_barcodes is not None:
        assert(df[barcodes_col].isin(all_barcodes).all())
    feature_to_feature_index = {feature: i for i, feature in enumerate(features, 1)}
    barcode_to_barcode_index = {barcode: i for i, barcode in enumerate(barcodes, 1)}

    df['feature_index'] = df[features_col].map(feature_to_feature_index)
    df['barcode_index'] = df[barcodes_col].map(barcode_to_barcode_index)

    with open(f'{prefix}barcodes.tsv', 'w') as fh:
        for b in barcodes:
            fh.write(b + '\n')

    with open(f'{prefix}features.tsv', 'w') as fh:
        for f in features:
            fh.write(f + '\n')

    with open(f'{prefix}matrix.mtx', 'w') as fh:
        # generate the header
        header = ['%%MatrixMarket matrix coordinate integer general', '%', ' '.join([str(len(features)), str(len(barcodes)), str(len(df))])]
        fh.write('\n'.join(header) + '\n')

        # write the final file
        for feature_index, barcode_index, count in zip(df['feature_index'], df['barcode_index'], df[count_col]):
            fh.write(f'{feature_index} {barcode_index} {count}\n')

    return True


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


parser = argparse.ArgumentParser()
parser.add_argument('--gtf', required=True)
parser.add_argument('--bam', required=True)
parser.add_argument('--prefix', required=True)
args = parser.parse_args()


gtf_df = gtf_to_df(args.gtf)
genes = gtf_df[gtf_df.feature=='gene']
transcripts = gtf_df[gtf_df.feature=='transcript']
genes['gene_id'] = parse_attribute(genes.attributes, 'gene_id')
genes['gene_name'] = parse_attribute(genes.attributes, 'gene_name')
transcripts['gene_id'] = parse_attribute(transcripts.attributes, 'gene_id')
transcripts['transcript_id'] = parse_attribute(transcripts.attributes, 'transcript_id')

gene_id_to_gene_name = dict(zip(genes.gene_id, genes.gene_name))
gene_id_to_feature = {gene_id: '{}\t{}\tGene Expression'.format(gene_id, gene_id_to_gene_name[gene_id]) for gene_id in gene_id_to_gene_name}

transcript_counts = {} # CB -> feature -> set(UMIs)
gene_counts = {} # CB -> feature -> set(UMIs)

with pysam.AlignmentFile(args.bam, 'rb') as bam:
    for count, read in enumerate(bam.fetch(until_eof=True)):
        if count % 1000000 == 0:
            logging.info('Processed {:,} reads'.format(count))
        if read.is_secondary or read.is_supplementary:
            continue
        if read.has_tag('CB') and read.get_tag('CB') != '-' and read.has_tag('UB') and read.get_tag('UB') != '-' and read.has_tag('GX') and read.get_tag('GX') != '-':
            if read.get_tag('CB') not in gene_counts:
                gene_counts[read.get_tag('CB')] = {}
            if read.get_tag('GX') not in gene_counts[read.get_tag('CB')]:
                gene_counts[read.get_tag('CB')][read.get_tag('GX')] = set()
            gene_counts[read.get_tag('CB')][read.get_tag('GX')].add(read.get_tag('UB'))
            
            if read.has_tag('TX') and read.get_tag('TX') != '-':
                if read.get_tag('CB') not in transcript_counts:
                    transcript_counts[read.get_tag('CB')] = {}
                if read.get_tag('TX') not in transcript_counts[read.get_tag('CB')]:
                    transcript_counts[read.get_tag('CB')][read.get_tag('TX')] = set()
                transcript_counts[read.get_tag('CB')][read.get_tag('TX')].add(read.get_tag('UB'))


logging.info('Writing gene count matrix')
gene_counts_df = []

for cb, cb_dict in gene_counts.items():
    for feature, umis in cb_dict.items():
        gene_counts_df.append([cb, gene_id_to_feature[feature], len(umis)])

gene_counts_df = pd.DataFrame(gene_counts_df, columns=['CB', 'feature', 'umis'])
write_mm(args.prefix + 'genes.', gene_counts_df, 'feature', 'CB', 'umis', all_features=list(gene_id_to_feature.values()))


logging.info('Writing transcript count matrix')
transcript_counts_df = []

for cb, cb_dict in transcript_counts.items():
    for feature, umis in cb_dict.items():
        transcript_counts_df.append([cb, feature, len(umis)])

transcript_counts_df = pd.DataFrame(transcript_counts_df, columns=['CB', 'feature', 'umis'])
write_mm(args.prefix + 'transcripts.', transcript_counts_df, 'feature', 'CB', 'umis', all_features=transcripts.transcript_id.unique())

logging.info('Done.')