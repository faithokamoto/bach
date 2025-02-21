#!/usr/bin/env python3
"""Simulates Hi-C reads for input to bach.
See --help for documentation.

Generates per-sample BAM files with read pairs
where at least one mate overlaps a variant given in a VCF file,
and following sample genotypes in the VCF.

Depth and variants with allele bias are parameters.

Example run:

python scripts/simulate_reads.py --depth 30 \
    --biased NC_000020.11:1182052 --biased NC_000020.11:1182049 \
    test/sequence.fasta test/snps.vcf test/bams

Written by ChatGPT.
"""

import argparse  # Parse command-line arguments
import os  # Handle file paths
import random  # Generate random insert sizes and mate positions
import pysam  # Read/write BAM and VCF files
from collections import defaultdict  # Organize reads by sample

def parse_fasta(fasta_path: str) -> dict:
    """Reads a FASTA file and returns a dictionary of sequences."""
    sequences = {}
    with pysam.FastaFile(fasta_path) as fasta:
        for ref in fasta.references:
            sequences[ref] = fasta.fetch(ref)
    return sequences

def validate_vcf(vcf_path: str, reference: dict):
    """Validates the VCF file against the reference genome."""
    with pysam.VariantFile(vcf_path) as vcf:
        for record in vcf:
            ref_base = reference[record.chrom][record.pos - 1]
            if record.ref != ref_base:
                raise ValueError(f"VCF reference mismatch at {record.chrom}:{record.pos} (VCF: {record.ref}, FASTA: {ref_base})")

def assign_random_genotype():
    """Randomly assigns a genotype as homozygous or heterozygous."""
    return (0, 0) if random.random() < 0.5 else (0, 1)

def get_mate_offset(biased: bool) -> int:
    """Generates a mate offset where most reads are near, but some are far apart."""
    if biased:
        if random.random() < 0.5:
            return random.randint(5_000_000, 20_000_000) * (-1 if random.random() < 0.5 else 1)
    if random.random() < 0.95:
        return random.randint(500, 5000)  # Most mates are nearby
    else:
        return random.randint(1_000_000, 50_000_000) * (-1 if random.random() < 0.5 else 1)  # Some mates are far

def simulate_reads(vcf_path: str, reference: dict, depth: int, biased_variants: set) -> dict:
    """Simulates Hi-C reads per sample based on VCF variants."""
    reads = defaultdict(list)
    with pysam.VariantFile(vcf_path) as vcf:
        for record in vcf:
            for sample in record.samples:
                genotype = record.samples[sample]['GT']
                if None in genotype:
                    genotype = assign_random_genotype()
                
                alleles = [record.ref] + list(record.alts)
                for _ in range(depth):
                    selected_allele = random.choice(alleles) if 1 in genotype else alleles[genotype[0]]
                    read_start = record.pos - 75  # Assume 150bp reads
                    mate_offset = get_mate_offset(record.pos in biased_variants)
                    if mate_offset < 0: mate_offset = 0
                    if mate_offset >= len(reference[record.chrom]):
                        mate_offset = len(reference[record.chrom]) - 1
                    read_seq = list(reference[record.chrom][read_start:read_start+150])
                    if selected_allele != record.ref:
                        read_seq[75] = selected_allele  # Mutate SNP position
                    mate_seq = reference[record.chrom][read_start + mate_offset: read_start + mate_offset + 150]
                    reads[sample].append((record.chrom, read_start, ''.join(read_seq), read_start + mate_offset, mate_seq))
    return reads

def write_bam(output_dir: str, reads: dict, reference: dict):
    """Writes simulated reads to sorted BAM files and indexes them."""
    os.makedirs(output_dir, exist_ok=True)
    header = {
        "HD": {"VN": "1.6"},
        "SQ": [{"SN": chrom, "LN": len(seq)} for chrom, seq in reference.items()]
    }
    
    for sample, read_data in reads.items():
        bam_path = os.path.join(output_dir, f"{sample}.pos.sorted.bam")
        with pysam.AlignmentFile(bam_path, "wb", header=header) as bam_file:
            for chrom, start, seq, mate_start, mate_seq in read_data:
                a = pysam.AlignedSegment(bam_file.header)
                a.query_name = f"read_{start}"
                a.reference_id = list(reference.keys()).index(chrom)
                a.reference_start = start
                a.seq = seq
                a.flag = 99  # Properly paired read
                a.next_reference_id = a.reference_id
                a.next_reference_start = mate_start
                a.cigar = [(0, len(seq))]  # Proper CIGAR string
                bam_file.write(a)
                
                mate = pysam.AlignedSegment(bam_file.header)
                mate.query_name = f"read_{start}/2"
                mate.reference_id = a.reference_id
                mate.reference_start = mate_start
                mate.seq = mate_seq
                mate.flag = 147  # Mate properly paired
                mate.next_reference_id = a.reference_id
                mate.next_reference_start = start
                mate.cigar = [(0, len(mate_seq))]
                bam_file.write(mate)
        pysam.sort("-o", bam_path, bam_path)
        pysam.index(bam_path)

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("fasta", help="Reference FASTA file")
    parser.add_argument("vcf", help="VCF file with variants")
    parser.add_argument("output_dir", help="Directory to store BAM files")
    parser.add_argument("--depth", type=int, default=30, help="Read depth per SNP")
    parser.add_argument("--biased", action='append', default=[], help="List of biased variants (chr:pos)")
    args = parser.parse_args()
    
    reference = parse_fasta(args.fasta)
    validate_vcf(args.vcf, reference)
    biased_variants = set(args.biased)
    reads = simulate_reads(args.vcf, reference, args.depth, biased_variants)
    write_bam(args.output_dir, reads, reference)
