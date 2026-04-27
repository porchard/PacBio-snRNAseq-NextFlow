#!/usr/bin/env python
# coding: utf-8

# In[7]:


import sys
import pandas as pd
import logging

logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(levelname)s: %(message)s')

ISOQUANT = sys.argv[1]

logging.info(f'Loading {ISOQUANT}')
df = pd.read_csv(ISOQUANT, sep='\t', comment='#', header=None, names=['read_id', 'chrom', 'strand', 'isoform_id', 'gene_id', 'assignment_type', 'assignment_events', 'exons', 'additional_info'])


logging.info('Calculating transcript counts')
KEEP_TRANSCRIPT_ASSIGNMENT_TYPES = ['unique', 'unique_minor_difference', 'inconsistent', 'inconsistent_non_intronic']
transcript_assignments = df.loc[df.assignment_type.isin(KEEP_TRANSCRIPT_ASSIGNMENT_TYPES),['read_id', 'isoform_id', 'assignment_type']].drop_duplicates()
transcript_assignments.assignment_type = pd.Categorical(transcript_assignments.assignment_type, categories=KEEP_TRANSCRIPT_ASSIGNMENT_TYPES, ordered=True)
assert(transcript_assignments.read_id.value_counts().max() == 1)


missing = pd.DataFrame({'read_id': df[~df.read_id.isin(transcript_assignments.read_id)].read_id.unique()})
missing['isoform_id'] = '-'
missing['assignment_type'] = '-'

logging.info('{:,} reads ({}%) successfully assigned to transcripts'.format(len(transcript_assignments), round(100*len(transcript_assignments)/(len(transcript_assignments) + len(missing)), 2)))


transcript_assignments = pd.concat([transcript_assignments, missing])
transcript_assignments.to_csv(sys.stdout, sep='\t', index=False)