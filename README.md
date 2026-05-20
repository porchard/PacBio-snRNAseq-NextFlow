# PacBio Kinnex snRNA-seq Pipeline

Nextflow DSL2 pipeline for preprocessing PacBio Kinnex (MAS-Seq) single-nucleus RNA-seq data. Starting from unmapped HiFi BAM files -- where each read is an array of ~16 cDNA molecules -- the pipeline produces tagged BAMs, gene and transcript count matrices, per-cell QC metrics, CellBender-corrected counts, and strand-specific bigwigs.

## Requirements

- [Nextflow](https://www.nextflow.io/) (DSL2)
- Singularity
- GPU (for CellBender ambient RNA removal)

## Input

The pipeline is configured via command-line parameters. The primary input is a **tab-separated config file** with three columns, including a header line:

| Column | Description |
|--------|-------------|
| `library` | Library name (used to name output files) |
| `readgroup` | Read group identifier |
| `bam` | Path to an unmapped HiFi BAM file |

Example config file:

```
library	readgroup	bam
sample1	rg1	/path/to/sample1_rg1.bam
sample1	rg2	/path/to/sample1_rg2.bam
sample2	rg1	/path/to/sample2.bam
```

Multiple readgroups per library are supported -- they are merged after the initial processing steps.

## Reference files

| Parameter | Description |
|-----------|-------------|
| `--adapters` | MAS-Seq adapter FASTA (see `data/masseq-adapters/`) |
| `--primers` | Primer FASTA (see `data/primers/`) |
| `--barcodes` | Cell barcode allowlist |
| `--fasta` | Reference genome FASTA |
| `--gtf` | Gene annotation GTF |

## Running

```bash
nextflow run -resume \
  --config config.tsv \
  --adapters data/masseq-adapters/mas16_primers.fasta \
  --primers data/primers/primers.10x3p.fa \
  --barcodes /path/to/barcodes.txt \
  --fasta /path/to/genome.fa \
  --gtf /path/to/genome.gtf \
  --results results \
  main.nf
```

There is no `nextflow.config`; all configuration is via command-line parameters.

## Output

```
results/
├── bam/                    # Final tagged, coordinate-sorted BAMs
├── count-matrices/         # Gene and transcript MatrixMarket files (.matrix.mtx, .features.tsv, .barcodes.tsv)
├── qc/                     # Per-cell QC metrics (TSV)
├── isoquant/               # IsoQuant read assignments
├── cellbender/             # CellBender-corrected matrices (.h5)
├── preprocess-for-scafe/   # BAM files preprocessed for SCAFE
└── bigwig/                 # Strand-specific bigwigs
```

## Bundled data files

### MAS-Seq adapters

MAS-Seq adapters in `data/masseq-adapters/` were downloaded from:

```
https://downloads.pacbcloud.com/public/dataset/MAS-Seq/REF-MAS_adapters/MAS-Seq_Adapter_v1/
```

### Primers

Primer FASTA files in `data/primers/` contain the 5' and 3' primer sequences used by `lima` to identify and remove technical sequences from IsoSeq reads. These are derived from the [Iso-Seq single-cell CLI workflow documentation](https://isoseq.how/umi/cli-workflow.html):

| File | 10x Kit |
|------|---------|
| `primers.10x3p.fa` | 10x Chromium 3' (e.g., v3.1) |
| `primers.10x5p.fa` | 10x Chromium 5' Gene Expression |

Choose the primer file matching the 10x library chemistry used.
