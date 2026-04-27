#!/usr/bin/env python
# coding: utf-8

import argparse
import logging

import pysam
import pandas as pd
import numpy as np

logging.basicConfig(level=logging.DEBUG, format='%(asctime)s - %(levelname)s: %(message)s')


def mean(d):
    """Return the mean of values represented as {value: count}.

    >>> mean({1: 2, 3: 1})
    1.6666666666666667
    >>> mean({5: 10})
    5.0
    >>> mean({0: 1, 10: 1})
    5.0
    """
    total = sum(v * c for v, c in d.items())
    n = sum(d.values())
    return total / n


def median(d):
    """Return the median of values represented as {value: count}.

    >>> median({1: 2, 3: 1})
    1
    >>> median({1: 1, 2: 1, 3: 1, 4: 1})
    2.5
    >>> median({5: 10})
    5
    >>> median({1: 3, 2: 1, 3: 1})
    1
    """
    sorted_values = sorted(d.keys())
    n = sum(d.values())
    cumulative = 0
    mid = n / 2

    if n % 2 == 1:
        target = (n + 1) // 2
        for v in sorted_values:
            cumulative += d[v]
            if cumulative >= target:
                return v
    else:
        lower = None
        for v in sorted_values:
            cumulative += d[v]
            if lower is None and cumulative >= mid:
                lower = v
            if cumulative >= mid + 1:
                return (lower + v) / 2 if lower != v else v


class Cell:

    def __init__(self, barcode):
        self.barcode = barcode
        self.supplementary_alignments = 0
        self.secondary_alignments = 0
        self.total_reads = 0
        # everything from here on down is primary alignments only
        self.umis = set() # (feature, UMI)
        self.primary_alignments = 0
        self.mapped_and_overlaps_exon = 0
        self.mapped_and_does_not_overlap_exon = 0
        self.mapped = 0
        self.uniquely_mapped = 0
        self.chromosome_read_counts = dict() # chrom -> count
        self.fragment_length_counts = dict() # fragment length -> count
        self.mapq = dict() # mapq -> count
        self.assigned_to_gene = 0

    def record_alignment(self, read):
        self.total_reads += 1
        
        if read.is_secondary or read.is_supplementary:
            if read.is_secondary:
                self.secondary_alignments += 1
            if read.is_supplementary:
                self.supplementary_alignments += 1
            return None
        
        self.primary_alignments += 1
        
        if read.has_tag('UB') and read.get_tag('UB') != '-' and read.has_tag('GX') and read.get_tag('GX') != '-':
            self.umis.add((read.get_tag('UB'), read.get_tag('GX')))
        
        if not read.is_unmapped:
            self.mapped += 1
            assert(read.get_tag('OE') in ['True', 'False'])
            if read.get_tag('OE') == 'True':
                self.mapped_and_overlaps_exon += 1
            else:
                self.mapped_and_does_not_overlap_exon += 1
        
        if read.mapping_quality > 0:
            self.uniquely_mapped += 1
            
        chrom = read.reference_name
        if chrom not in self.chromosome_read_counts:
            self.chromosome_read_counts[chrom] = 0
        self.chromosome_read_counts[chrom] += 1
        
        read_length = read.infer_query_length()
        if read_length is not None:
            if read_length not in self.fragment_length_counts:
                self.fragment_length_counts[read_length] = 0
            self.fragment_length_counts[read_length] += 1
        
        if read.mapping_quality not in self.mapq:
            self.mapq[read.mapping_quality] = 0
        self.mapq[read.mapping_quality] += 1
        
        if read.has_tag('GX') and read.get_tag('GX') != '-':
            self.assigned_to_gene += 1
        
        return None


    def gather_metrics(self):
        metrics = dict()
        metrics['barcode'] = self.barcode
        metrics['total_reads'] = self.total_reads
        metrics['secondary_alignments'] = self.secondary_alignments
        metrics['supplementary_alignments'] = self.supplementary_alignments
        metrics['primary_alignments'] = self.primary_alignments
        metrics['umis'] = len(self.umis)
        metrics['fraction_exonic'] = self.mapped_and_overlaps_exon / self.mapped if self.mapped > 0 else 0
        metrics['mapped_primary_alignments'] = self.mapped
        metrics['fraction_primary_alignments_mapped'] = self.mapped / self.primary_alignments if self.primary_alignments > 0 else 0
        metrics['unmapped_primary_alignments'] = self.primary_alignments - self.mapped
        metrics['fraction_primary_alignments_unmapped'] = (self.primary_alignments - self.mapped) / self.primary_alignments if self.primary_alignments > 0 else 0
        metrics['uniquely_mapped_primary_alignments'] = self.uniquely_mapped
        metrics['fraction_primary_alignments_uniquely_mapped'] = self.uniquely_mapped / self.primary_alignments if self.primary_alignments > 0 else 0
        metrics['fraction_mitochondrial'] = self.chromosome_read_counts['chrM'] / self.primary_alignments if 'chrM' in self.chromosome_read_counts else 0
        metrics['median_fl'] = median(self.fragment_length_counts) if len(self.fragment_length_counts) > 0 else 0
        metrics['mean_fl'] = mean(self.fragment_length_counts) if len(self.fragment_length_counts) > 0 else 0
        metrics['primary_alignments_assigned_to_gene'] = self.assigned_to_gene
        metrics['fraction_primary_alignments_assigned_to_gene'] = self.assigned_to_gene / self.primary_alignments if self.primary_alignments > 0 else 0
        metrics['primary_alignments_with_mapq_0'] = self.mapq.get(0, 0)
        metrics['fraction_primary_alignments_mapq_0'] = self.mapq.get(0, 0) / self.primary_alignments if self.primary_alignments > 0 else 0
        metrics['mean_mapq'] = mean(self.mapq)
        metrics['median_mapq'] = median(self.mapq)
        return metrics


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
parser.add_argument('bam')
#a = ['/scratch/scjp_root/scjp0/porchard/PacBio-pipeline-development/work/post-align/work/0c/c519c02cd4d17fa554051bc91b23c0/15245-AH-1.bam']
#args = parser.parse_args(a)
args = parser.parse_args()

# QC metrics
cells = dict()
no_cell_tag = 0
total_reads = 0

with pysam.AlignmentFile(args.bam, 'rb') as bam:
    for read in bam.fetch(until_eof=True):
        total_reads += 1
        if total_reads % 1000000 == 0:
            logging.info('Processed {:,} reads'.format(total_reads))                
        barcode = read.get_tag('CB')
        if barcode not in cells:
            cells[barcode] = Cell(barcode)
        cells[barcode].record_alignment(read)


logging.info('Processed {:,} reads'.format(total_reads))

logging.info('Writing QC metrics')


print_metrics = ['barcode', 'total_reads', 'secondary_alignments', 'supplementary_alignments', 'primary_alignments', 'umis', 'fraction_exonic', 'mapped_primary_alignments', 'fraction_primary_alignments_mapped', 'unmapped_primary_alignments', 'fraction_primary_alignments_unmapped', 'uniquely_mapped_primary_alignments', 'fraction_primary_alignments_uniquely_mapped', 'fraction_mitochondrial', 'median_fl', 'mean_fl', 'primary_alignments_assigned_to_gene', 'fraction_primary_alignments_assigned_to_gene', 'primary_alignments_with_mapq_0', 'fraction_primary_alignments_mapq_0', 'mean_mapq', 'median_mapq']
print('\t'.join(print_metrics))

for cell in cells.values():
    metrics = cell.gather_metrics()
    to_print = [str(metrics[i]) for i in print_metrics]
    print('\t'.join(to_print))

logging.info('Done')