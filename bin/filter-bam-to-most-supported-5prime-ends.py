#!/usr/bin/env python
# coding: utf-8

import logging
import argparse

import pysam
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns


parser = argparse.ArgumentParser()
parser.add_argument('--bam-in', required=True, help='Input BAM file with CB and UB tags')
parser.add_argument('--bam-out', required=True, help='Output BAM file with filtered reads')
parser.add_argument('--prefix', default='test.', help='Prefix for output plots')
args = parser.parse_args()


# CB --> UMI --> (chrom, gene, 5' end)
# BAM_IN = '/scratch/scjp_root/scjp0/porchard/2023-HSM-ONT/work/fix-swaps-bc23eb7-inc-skipped-basecalling/9266-VD-1/9266-VD-1.tagged.bam'
# BAM_IN = '/nfs/turbo/umms-scjp-pfizer/human-multiome/peter-processing/2023-HSM/work/gex5/gex5/results/starsolo/9267-VD-1-hg38/9267-VD-1-hg38.Aligned.sortedByCoord.out.bam'
# BAM_OUT = 'out.bam'
# PREFIX = 'test.'

BAM_IN = args.bam_in
BAM_OUT = args.bam_out
PREFIX = args.prefix

logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(levelname)s: %(message)s')

logging.info("Determining most likely 5' ends")

counts = dict()
# ends = []

with pysam.AlignmentFile(BAM_IN, 'rb') as f:
    total_reads = 0
    for read in f.fetch(until_eof=True):
        total_reads += 1
        if total_reads % 1000000 == 0:
            logging.info('Processed {:,} reads'.format(total_reads))
        if read.is_unmapped or read.is_secondary or read.is_supplementary or read.is_read2:
            continue
        if read.has_tag('CB') and read.has_tag('UB'):
            chrom = read.reference_name
            five_prime_end = read.reference_start if read.is_forward else read.reference_end
            strand = '+' if read.is_forward else '-'
            CB = read.get_tag('CB')
            UMI = read.get_tag('UB')
            assert read.has_tag('GX') or read.has_tag('GN'), 'Read is missing gene name tag (GX/GN)?'
            gene_tag = 'GX' if read.has_tag('GX') else 'GN'
            gene = read.get_tag(gene_tag)
            gene_or_region = gene if gene != '-' else 'r' + str(int(five_prime_end // 1e6))
            if CB == '-' or UMI == '-':
                continue
            # ends.append([CB, UMI, chrom, five_prime_end, strand, gene, gene_or_region])
            key = (CB, UMI, chrom, strand, gene_or_region)
            if key not in counts:
                counts[key] = dict()
            if five_prime_end not in counts[key]:
                counts[key][five_prime_end] = 0
            counts[key][five_prime_end] += 1

logging.info('Read {:,} total reads; found {:,} unique (CB, UMI, chrom, strand, (gene or region)) combinations'.format(total_reads, len(counts)))

logging.info('Identifying most likely 5\' end for each (CB, UMI, chrom, strand, (gene or region)) combination')

selected_ends = dict()

for key, c in counts.items():
    strand = key[3]
    assert(strand in ['+', '-'])
    max_value = max(c.values())
    keys_with_max_value = list(sorted([k for k, v in c.items() if v == max_value]))
    # if '+' strand, then the smallest value is the most upstream
    selected_ends[key] = keys_with_max_value[0] if strand == '+' else keys_with_max_value[-1]  


logging.info("Filtering to reads with most likely 5' ends")


with pysam.AlignmentFile(BAM_IN, 'rb') as f:
    with pysam.AlignmentFile(BAM_OUT, 'wb', header=f.header) as out:
        total_reads = 0
        keep = 0
        drop_bad_end = 0
        drop_other = 0
        for read in f.fetch(until_eof=True):
            total_reads += 1
            if total_reads % 1000000 == 0:
                logging.info('Processed {:,} reads'.format(total_reads))
            if read.is_unmapped or read.is_secondary or read.is_supplementary or read.is_read2:
                drop_other += 1
                continue
            if read.has_tag('CB') and read.has_tag('UB'):
                chrom = read.reference_name
                five_prime_end = read.reference_start if read.is_forward else read.reference_end
                strand = '+' if read.is_forward else '-'
                CB = read.get_tag('CB')
                UMI = read.get_tag('UB')
                assert read.has_tag('GX') or read.has_tag('GN'), 'Read is missing gene name tag (GX/GN)?'
                gene_tag = 'GX' if read.has_tag('GX') else 'GN'
                gene = read.get_tag(gene_tag)
                gene_or_region = gene if gene != '-' else 'r' + str(int(five_prime_end // 1e6))
                if CB == '-' or UMI == '-':
                    drop_other += 1
                    continue
                key = (CB, UMI, chrom, strand, gene_or_region)
                assert(key in selected_ends)
                if selected_ends[key] == five_prime_end:
                    out.write(read)
                    keep += 1
                else:
                    drop_bad_end += 1

assert(drop_other + drop_bad_end + keep == total_reads)
logging.info("Processed {:,} total reads; kept {:,}, dropped {:,} for having likely incorrect 5' end, dropped {:,} for other reasons".format(total_reads, keep, drop_bad_end, drop_other))


logging.info("Generating summary plots")

unique_end_counts = [['-'.join(k), len(v), max(v.keys()) - min(v.keys())] for k, v in counts.items()]
unique_end_counts = pd.DataFrame(unique_end_counts, columns=['key', 'unique_ends', 'range'])
unique_end_counts = unique_end_counts.sort_values('range')

fig, ax = plt.subplots()

sns.histplot(ax=ax, data=unique_end_counts, x='unique_ends', discrete=True)
ax.set_xlabel("Unique 5' ends")
ax.set_ylabel('Read sets (CB+UMI+chrom+strand+(gene or region))')

fig.savefig(f'{PREFIX}end-counts.png', dpi=300, facecolor='white', bbox_inches='tight')


range_bins = [0, 5, 10, 50, 100, 500, 1000, 5000]
bin_names = ['0', '1-5', '6-10', '11-50', '51-100', '101-500', '501-1000', '1001-5000', '5001+']
unique_end_counts['range_bin'] = np.digitize(unique_end_counts['range'], range_bins, right=True)
unique_end_counts['range_bin'] = pd.Categorical(unique_end_counts['range_bin'].map(lambda x: bin_names[x]), categories=bin_names, ordered=True)


fig, ax = plt.subplots()

sns.histplot(ax=ax, data=unique_end_counts, x='range_bin', stat='count')
for t in ax.get_xticklabels():
    t.set(rotation=45)
ax.set_xlabel("Distance between most upstream and downstream 5' ends (bps)")
ax.set_ylabel('Read sets (CB+UMI+chrom+strand+(gene or region))')

fig.savefig(f'{PREFIX}end-distance.png', dpi=300, facecolor='white', bbox_inches='tight')


logging.info('Done')