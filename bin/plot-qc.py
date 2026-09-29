#!/usr/bin/env python
# coding: utf-8


import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import numpy as np
import upsetplot
import argparse


def barcode_rank_plot(metrics, ax):
    df = metrics.sort_values('rna_umis', ascending=False)
    df['barcode_rank'] = range(1, len(df) + 1)
    sns.scatterplot(x='barcode_rank', y='rna_umis', data=df, ax=ax, edgecolor=None, alpha=0.1, s=5)
    ax.set_ylabel('UMIs')
    ax.set_xscale('log')
    ax.set_yscale('log')
    ax.set_xlabel('Barcode rank')
    return ax


def scatter(x, y, metrics, ax):
    log_scale = ['rna_umis']
    labels = {
        'rna_umis': 'UMIs (RNA)',
        'rna_genes': 'Genes detected (RNA)',
        'rna_fraction_mitochondrial': 'Fraction mito. (RNA)',
        'rna_fraction_exonic': 'Exon/full-gene-body ratio (RNA)',
        'rna_fraction_primary_alignments_mapq_0': 'Frac. primary alignments w/ mapq=0 (RNA)',
        'rna_fraction_primary_alignments_assigned_to_gene': 'Frac. assigned to gene (RNA)'
    }
    if x not in labels:
        labels[x] = x
    if y not in labels:
        labels[y] = y
    kwargs = {'edgecolor': None, 'alpha': 0.1, 's': 3}
    sns.scatterplot(x=x, y=y, data=metrics, ax=ax, **kwargs)
    ax.set_xlabel(labels[x])
    ax.set_ylabel(labels[y])
    if x in log_scale:
        ax.set_xscale('log')
    if y in log_scale:
        ax.set_yscale('log')
    return x


parser = argparse.ArgumentParser()
parser.add_argument('--rna-metrics', required=True)
parser.add_argument('--prefix', required=True)
args = parser.parse_args()

RNA_METRICS = args.rna_metrics
PREFIX = args.prefix

rna_metrics = pd.read_csv(RNA_METRICS, sep='\t')
rna_metrics = rna_metrics[rna_metrics.barcode!='-']

metrics = rna_metrics.set_index('barcode').rename(columns=lambda x: 'rna_' + x).reset_index()

# Plot QC metrics
fig, axs = plt.subplots(ncols=5, nrows=1, figsize=(5*4, 1*4))

ax = axs[0]
barcode_rank_plot(metrics, ax)

ax = axs[1]
scatter('rna_umis', 'rna_fraction_mitochondrial', metrics, ax)

ax = axs[2]
scatter('rna_umis', 'rna_fraction_exonic', metrics, ax)

ax = axs[3]
scatter('rna_umis', 'rna_fraction_primary_alignments_mapq_0', metrics, ax)

ax = axs[4]
scatter('rna_umis', 'rna_fraction_primary_alignments_assigned_to_gene', metrics, ax)

fig.tight_layout()
fig.savefig(f'{PREFIX}qc.png', dpi=300, facecolor='white', bbox_inches='tight')
