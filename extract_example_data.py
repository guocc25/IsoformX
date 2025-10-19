#!/usr/bin/env python3
"""Extract first N reads from BAM files using pysam."""

import pysam
import sys
import os

def extract_reads(input_bam, output_bam, num_reads=6000):
    """Extract first N reads from BAM file."""
    print(f"Extracting first {num_reads} reads from {input_bam}...")
    
    with pysam.AlignmentFile(input_bam, "rb") as infile:
        with pysam.AlignmentFile(output_bam, "wb", template=infile) as outfile:
            count = 0
            for read in infile:
                if count >= num_reads:
                    break
                outfile.write(read)
                count += 1
            
            print(f"Extracted {count} reads to {output_bam}")

def main():
    # Extract from elongating BAM
    input_elongating = "../input_data/13_tagged_elongating.bam"
    output_elongating = "example_data/13_elongating_6k.bam"
    
    if os.path.exists(input_elongating):
        extract_reads(input_elongating, output_elongating, 6000)
    else:
        print(f"Error: {input_elongating} not found")
    
    # Extract from polyadenylation BAM
    input_polya = "../input_data/13_tagged_polyadenylation.bam"
    output_polya = "example_data/13_polyadenylation_6k.bam"
    
    if os.path.exists(input_polya):
        extract_reads(input_polya, output_polya, 6000)
    else:
        print(f"Error: {input_polya} not found")

if __name__ == "__main__":
    main()
