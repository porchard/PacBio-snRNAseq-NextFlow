#!/usr/bin/env python

from umi_tools import UMIClusterer
import pysam
import argparse
import logging

logging.basicConfig(level=logging.DEBUG, format='%(asctime)s - %(levelname)s: %(message)s')


def read_to_group(read):
    if read.get_tag('GX') != '-':
        return read.get_tag('GX')
    else:
        if read.is_unmapped:
            return None
        else:
            return '{}:{}'.format(read.reference_name, read.pos // 1e6)


def reads_to_group_umi_counts(reads, uncorrected_umi_tag):
    group_to_umi_counts = dict()
    for read in reads:
        group = read_to_group(read)
        umi = read.get_tag(uncorrected_umi_tag)
        if group is None:
            continue
        if group not in group_to_umi_counts:
            group_to_umi_counts[group] = dict()
        if umi not in group_to_umi_counts[group]:
            group_to_umi_counts[group][umi] = 0
        group_to_umi_counts[group][umi] += 1
    return group_to_umi_counts


def correct_group_umis(umis, clusterer):
    clustered_umis = clusterer({k.encode(): v for k, v in umis.items()}, threshold=1)
    corrections = dict()
    for x in clustered_umis:
        for i in x:
            assert(i not in corrections)
            corrections[i] = x[0]
    return {k.decode(): v.decode() for k, v in corrections.items()}


def correct_umis(reads, uncorrected_umi_tag, corrected_umi_tag):
    group_umi_counts = reads_to_group_umi_counts(reads, uncorrected_umi_tag)
    clusterer = UMIClusterer(cluster_method="directional")
    corrections = {group: correct_group_umis(umis, clusterer) for group, umis in group_umi_counts.items()}
    for r in reads:
        group = read_to_group(r)
        if group is None:
            r.set_tag(corrected_umi_tag, '-')
        else:
            r.set_tag(corrected_umi_tag, corrections[group][r.get_tag(uncorrected_umi_tag)])
    return reads


def get_barcode(read, barcode_tag):
    """Return the cell barcode, or None if missing/empty."""
    try:
        val = read.get_tag(barcode_tag)
    except KeyError:
        return None
    if val in ('', '-'):
        return None
    return val


class BamReader:
    """Iterate over a BAM file grouped by cell barcode.

    Each call to next() returns (barcode, reads) where barcode is the CB
    value (or None for reads without a valid barcode) and reads is the list
    of consecutive reads sharing that barcode.
    """

    def __init__(self, bam_path, barcode_tag='CB'):
        self._bam = pysam.AlignmentFile(bam_path, 'rb')
        self._barcode_tag = barcode_tag
        self._iter = self._bam.fetch(until_eof=True)
        self._next_read = None
        self._exhausted = False
        self._advance()

    def _advance(self):
        try:
            self._next_read = next(self._iter)
        except StopIteration:
            self._next_read = None
            self._exhausted = True

    @property
    def header(self):
        return self._bam.header

    def __iter__(self):
        return self

    def __next__(self):
        if self._exhausted:
            raise StopIteration
        barcode = get_barcode(self._next_read, self._barcode_tag)
        reads = [self._next_read]
        self._advance()
        while not self._exhausted and get_barcode(self._next_read, self._barcode_tag) == barcode:
            reads.append(self._next_read)
            self._advance()
        return (barcode, reads)

    def close(self):
        self._bam.close()

    def __enter__(self):
        return self

    def __exit__(self, *args):
        self.close()


def process_bam(bam_in_path, bam_out_path, barcode_tag, uncorrected_umi_tag, corrected_umi_tag):
    total_reads = 0
    total_barcodes = 0

    with BamReader(bam_in_path, barcode_tag) as reader:
        with pysam.AlignmentFile(bam_out_path, 'wb', header=reader.header) as bam_out:
            for barcode, reads in reader:
                total_reads += len(reads)
                if barcode is None:
                    logging.info('Skipping {} reads with no barcode'.format(len(reads)))
                    for r in reads:
                        r.set_tag(corrected_umi_tag, '-')
                else:
                    logging.info('Correcting UMIs for CB {}'.format(barcode))
                    correct_umis(reads, uncorrected_umi_tag, corrected_umi_tag)
                    total_barcodes += 1
                    if total_barcodes % 1000 == 0:
                        logging.info('Processed {:,} barcodes ({:,} reads)'.format(total_barcodes, total_reads))
                for r in reads:
                    bam_out.write(r)

    logging.info('Done. Processed {:,} reads across {:,} barcodes.'.format(total_reads, total_barcodes))


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Correct UMIs in a BAM file using directional clustering. Input BAM must be sorted by cell barcode tag.')
    parser.add_argument('bam_in', help='Input BAM file (sorted by cell barcode tag)')
    parser.add_argument('bam_out', help='Output BAM file')
    parser.add_argument('--barcode-tag', default='CB', help='BAM tag for cell barcode (default: CB)')
    parser.add_argument('--uncorrected-umi-tag', default='XM', help='BAM tag for uncorrected UMI (default: XM)')
    parser.add_argument('--corrected-umi-tag', default='UB', help='BAM tag for corrected UMI (default: UB)')
    args = parser.parse_args()

    process_bam(args.bam_in, args.bam_out, args.barcode_tag, args.uncorrected_umi_tag, args.corrected_umi_tag)
