#!/usr/bin/env python3
"""Make histograms of mate locations split by SNP allele.
See --help for how to use.

Takes reads which overlap an RSID and splits them by ref/alt allele.
Takes mates which map to the same chromosome and makes histograms.

Mates of ref/alt reads are plotted either upwards or downwards.
Mates located on a different chromosome are dropped entirely.

Mostly written by ChatGPT.
"""

import argparse  # Command-line argument parsing
import math # Floor and ceiling
import os  # File operations

import pysam  # Reading BAM files
import matplotlib.pyplot as plt # Make plots
import numpy as np # Simple list operations

from typing import Dict, List, Tuple

MBP = 1_000_000
"""Mega base pair (1 million bp)."""
HET_GENO = [(0, 1), (1, 0)]
"""Possible phased heterozygous genotypes."""
REF_COLOR = '#E1BE6A'
"""Color used in ref histograms."""
ALT_COLOR = '#40B0A6'
"""Color used in alt histograms."""

def get_snp_info(vcf_path: str, chrom: str, 
                 pos: int) -> Tuple[str, str, str, Dict[str, bool]]:
    """Get SNP info from a VCF file.

    Parameters
    ----------
    vcf_path : str
        Path to the VCF file.
    chrom : str
        Chromosome identifier.
    pos : int
        1-based bp position of the SNP on the chromosome.

    Returns
    -------
    Tuple[str, str, str, Dict[str, bool]]
        A tuple (SNP ID, ref allele, alt allele, is_het).
        `is_het` is whether each sample's genotype is heterozygous

    Raises
    ------
    ValueError
        If the SNP is not found at the specified position.
    """

    with pysam.VariantFile(vcf_path) as vcf:
        for record in vcf.fetch(chrom, pos - 1, pos):
            is_het = {name: s['GT'] in HET_GENO for name, s
                      in record.samples.items()}
            return (record.id, record.ref, record.alts[0], is_het)

    raise ValueError(f"SNP not found at {chrom}:{pos}")

def get_allele_mate_pos(bam_path: str, chrom: str, snp_pos: int, 
                        ref: str, alt: str) -> Dict[str, List[float]]:
    """Get mate positions for reads supporting ref/alt alleles.

    Ignores reads with unmapped mates and reads that do not cover the SNP.
    BAM file must be position-sorted and indexed for fast random access.

    Parameters
    ----------
    bam_path : str
        Path to the BAM file.
    chrom : str
        Chromosome of the SNP.
    snp_pos : int
        Position of the SNP in bp.
    ref : str
        Reference allele of the SNP.
    alt : str
        Alternate allele of the SNP.

    Returns
    -------
    Dict[str, List[float]]
        Two lists of bp positions: mates of ref/alt-supporting reads.
    """

    ref_mates = []
    alt_mates = []
    with pysam.AlignmentFile(bam_path, 'rb') as bam_file:
        for read in bam_file.fetch(chrom, snp_pos - 1, snp_pos):
            # We want mate pairs with high-quality mapping to the same chr
            if (read.is_secondary or read.mate_is_unmapped
                or read.reference_name != read.next_reference_name):
                continue
            
            ref_positions = read.get_reference_positions(full_length=True)
            # Handle deletions of the SNP position
            if snp_pos - 1 not in ref_positions: continue
            
            base = read.query_sequence[ref_positions.index(snp_pos - 1)]
            # Convert positions to Mbp for readability
            if base == ref:
                ref_mates.append(read.next_reference_start)
            elif base == alt:
                alt_mates.append(read.next_reference_start)
    
    return {'ref': sorted(ref_mates), 'alt': sorted(alt_mates)}

def plot_read_positions(rsid: str, snp_pos: float,
                        mates: Dict[str, Dict[str, List[float]]]) -> None:
    """Plot up/down histograms of mate positions.

    Reference mate location histogram faces upwards,
    alternate, downwards. Both normalized by number
    of read mates within 0.5Mbp. Line drawn at SNP pos.
    
    Parameters
    ----------
    rsid : str
        SNP RSID for plot title
    snp_pos : float
        Position of the SNP in bp
    mates : Dict[str, Dict[str, List[float]]]
        Dictionary of sample to two lists of bp positions: 
        mates of ref/alt-supporting reads.
    """

    # Use Mbp for readability
    ref_pos = np.concatenate([pos['ref'] for pos in mates.values()]) / MBP
    alt_pos = np.concatenate([pos['alt'] for pos in mates.values()]) / MBP

    if len(ref_pos) == 0:
        print('No ref read mates to plot.')
        return
    if len(alt_pos) == 0:
        print('No alt read mates to plot.')
        return

    fig, (ax1, ax2) = plt.subplots(nrows=2, sharex=True, figsize=(8, 4))

    # Use a common range for both plots for direct comparability
    pos_min = 5 * math.floor(min(ref_pos.min(), alt_pos.min()) / 5)
    pos_max = 5 * math.ceil(max(ref_pos.max(), alt_pos.max()) / 5)
    pos_range = (pos_min, pos_max)
    bins = int((pos_max - pos_min) / 5)

    ax1.hist(ref_pos, bins=bins, range=pos_range, color=REF_COLOR, label='Ref')
    # Alt histogram is plotted upside-down
    ax2.hist(alt_pos, bins=bins, range=pos_range, color=ALT_COLOR, label='Alt')
    ax2.invert_yaxis()
    # Vertical line at SNP position
    ax2.axvline(snp_pos / MBP, color='black', linestyle='--', label=f'{rsid}')

    fig.legend(fontsize=12, loc=(0.7, 0.2))
    fig.suptitle(f'Mate pairs of reads overlapping {rsid}', fontsize=15)
    fig.text(0.01, 0.5, 'Read counts', rotation='vertical',
             verticalalignment='center', fontsize=12)
    ax2.set_xlabel('Position (Mbp)', fontsize=12)

    fig.tight_layout()
    # Allow Y axis label room & remove gap between subplots
    fig.subplots_adjust(left=0.075, hspace=0)
    fig.savefig(f'plots/{rsid}_mate_locations.png')

if __name__ == '__main__':
    # Set up argparse
    parser = argparse.ArgumentParser(
        description='Makes visual representations of allele biased segments.',
        formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument('vcf', help='Path to the SNP VCF file')
    parser.add_argument('bam_dir', help='Directory containing BAM files')
    parser.add_argument('snp', help='SNP of interest, in chr:pos format')

    parser.add_argument('--bam-ext', help='Extension for BAMs (<sample>.ext)',
                        required=False, default='.pos.sorted.bam')
    args = parser.parse_args()

    # Basic validity checks
    if not os.path.isfile(args.vcf):
        raise FileNotFoundError(f'VCF {args.vcf} does not exist')
    if not os.path.isdir(args.bam_dir):
        raise FileNotFoundError(f'BAM dir {args.bam_dir} is not a directory')
    
    # Enforce SNP format
    try:
        chrom, snp_pos = args.snp.split(':')
        snp_pos = int(snp_pos)
    except ValueError:
        raise ValueError(f'SNP {args.snp} not in chr:pos format')

    # Check files in the given directory
    bams = {file.split('.')[0]: f'{args.bam_dir}/{file}'
            for file in os.listdir(args.bam_dir) if file.endswith(args.bam_ext)}
    if not bams:
        raise FileNotFoundError(f'No {args.bam_dir}/*.{args.bam_ext}')

    # Main logic
    rsid, ref, alt, is_het = get_snp_info(args.vcf, chrom, snp_pos)

    mates = {sample: get_allele_mate_pos(file, chrom, snp_pos, ref, alt)
             for sample, file in bams.items() if is_het[sample]}
    plot_read_positions(rsid, snp_pos, mates)