"""Isoform proportion analysis across samples."""

from typing import List, Dict, Tuple, Optional
import pandas as pd
import pysam
import numpy as np
from pathlib import Path
from collections import defaultdict
import gffutils
from dataclasses import dataclass
import re


@dataclass
class IsoformInfo:
    """Information about an isoform from GTF."""
    isoform_id: str
    gene_id: str
    chromosome: str
    strand: str
    start: int
    end: int
    exons: List[Tuple[int, int]]


def load_isoforms_from_gtf(gtf_file: str) -> Dict[str, IsoformInfo]:
    """
    Load isoform information from GTF file.
    
    Args:
        gtf_file: Path to GTF file
    
    Returns:
        Dictionary mapping isoform_id to IsoformInfo
    """
    print(f"Loading isoforms from {gtf_file}...")
    
    # Create database from GTF file
    db = gffutils.create_db(gtf_file, ':memory:', force=True, keep_order=True)
    
    isoforms = {}
    
    # Process each transcript
    for transcript in db.features_of_type('transcript'):
        # Get exons for this transcript
        exons = []
        for exon in db.children(transcript, featuretype='exon', order_by='start'):
            exons.append((exon.start - 1, exon.end))  # Convert to 0-based coordinates
        
        # Sort exons by start position
        exons.sort(key=lambda x: x[0])
        
        # Extract attributes
        isoform_id = transcript.attributes.get('transcript_id', [transcript.id])[0]
        
        # Try to get gene_id from various sources
        gene_id = transcript.attributes.get('gene_id', ['unknown'])[0]
        if gene_id == 'unknown' or gene_id == 'None':
            # Try to get from Parent attribute
            parent = transcript.attributes.get('Parent', ['unknown'])[0]
            if parent != 'unknown':
                gene_id = parent
            else:
                # Try to get reference transcript information
                ref_transcript = transcript.attributes.get('reference_transcript', ['unknown'])[0]
                if ref_transcript != 'unknown':
                    gene_id = ref_transcript
        
        # Create isoform info
        isoform_info = IsoformInfo(
            isoform_id=isoform_id,
            gene_id=gene_id,
            chromosome=transcript.chrom,
            strand=transcript.strand,
            start=transcript.start - 1,  # Convert to 0-based
            end=transcript.end,
            exons=exons
        )
        
        isoforms[isoform_id] = isoform_info
    
    print(f"Loaded {len(isoforms)} isoforms")
    return isoforms


def count_reads_for_isoform(bam_file: str, isoform: IsoformInfo, min_overlap: float = 0.5) -> int:
    """
    Count reads that overlap with an isoform.
    
    Args:
        bam_file: Path to BAM file
        isoform: Isoform information
        min_overlap: Minimum overlap ratio required
    
    Returns:
        Number of reads overlapping the isoform
    """
    try:
        bam = pysam.AlignmentFile(bam_file, "rb")
        count = 0
        
        # Get reads overlapping the isoform region
        for read in bam.fetch(isoform.chromosome, isoform.start, isoform.end):
            if read.is_unmapped or read.is_duplicate or read.is_qcfail:
                continue
            
            # Calculate overlap with isoform
            read_start = read.reference_start
            read_end = read.reference_end
            
            # Calculate overlap length
            overlap_start = max(read_start, isoform.start)
            overlap_end = min(read_end, isoform.end)
            overlap_length = max(0, overlap_end - overlap_start)
            
            # Calculate overlap ratio
            read_length = read_end - read_start
            isoform_length = isoform.end - isoform.start
            
            if read_length > 0 and isoform_length > 0:
                overlap_ratio = overlap_length / min(read_length, isoform_length)
                
                if overlap_ratio >= min_overlap:
                    count += 1
        
        bam.close()
        return count
        
    except Exception as e:
        print(f"Error processing {bam_file}: {e}")
        return 0


def analyze_isoform_proportions(
    bam_files: List[str],
    isoform_gtf: str,
    min_reads: int = 1,
    min_overlap: float = 0.5,
    verbose: bool = False
) -> pd.DataFrame:
    """
    Analyze isoform proportions across samples.
    
    Args:
        bam_files: List of BAM file paths
        isoform_gtf: Path to isoform GTF file
        min_reads: Minimum reads per isoform to include
        min_overlap: Minimum overlap ratio for read-isoform matching
        verbose: Verbose output
    
    Returns:
        DataFrame with isoform proportions
    """
    if verbose:
        print(f"Analyzing {len(bam_files)} BAM files...")
    
    # Load isoform information
    isoforms = load_isoforms_from_gtf(isoform_gtf)
    
    # Initialize results
    results = []
    sample_names = []
    
    # Extract sample names from BAM files
    for bam_file in bam_files:
        sample_name = Path(bam_file).stem
        sample_names.append(sample_name)
    
    if verbose:
        print(f"Sample names: {sample_names}")
    
    # Count reads for each isoform in each sample
    for isoform_id, isoform_info in isoforms.items():
        if verbose and len(results) % 100 == 0:
            print(f"Processing isoform {len(results) + 1}/{len(isoforms)}: {isoform_id}")
        
        # Count reads in each sample
        sample_counts = {}
        total_reads = 0
        
        for i, bam_file in enumerate(bam_files):
            count = count_reads_for_isoform(bam_file, isoform_info, min_overlap)
            sample_counts[sample_names[i]] = count
            total_reads += count
        
        # Only include isoforms with sufficient reads
        if total_reads >= min_reads:
            # Calculate proportions
            sample_proportions = {}
            for sample_name in sample_names:
                count = sample_counts[sample_name]
                proportion = count / total_reads if total_reads > 0 else 0
                sample_proportions[f"{sample_name}_reads"] = count
                sample_proportions[f"{sample_name}_proportion"] = proportion
            
            # Create result row
            result_row = {
                'isoform_id': isoform_id,
                'gene_id': isoform_info.gene_id,
                'chromosome': isoform_info.chromosome,
                'strand': isoform_info.strand,
                'start': isoform_info.start,
                'end': isoform_info.end,
                'total_reads': total_reads,
                'samples_with_reads': sum(1 for count in sample_counts.values() if count > 0),
                'max_reads_per_sample': max(sample_counts.values()),
                'min_reads_per_sample': min(sample_counts.values()),
                'mean_reads_per_sample': np.mean(list(sample_counts.values())),
                'std_reads_per_sample': np.std(list(sample_counts.values())),
            }
            
            # Add sample-specific data
            result_row.update(sample_proportions)
            
            results.append(result_row)
    
    # Create DataFrame
    df = pd.DataFrame(results)
    
    if verbose:
        print(f"Found {len(df)} isoforms with at least {min_reads} reads")
    
    # Only sort if we have results
    if len(df) > 0:
        # Natural chromosome sort, then by start position
        def natural_sort_key(text: str):
            parts = re.split(r'(\d+)', str(text))
            return tuple(int(p) if p.isdigit() else p.lower() for p in parts)
        df['__chrom_key'] = df['chromosome'].map(natural_sort_key)
        df = df.sort_values(['__chrom_key', 'start', 'isoform_id'], ascending=[True, True, True])
        df = df.drop(columns=['__chrom_key'])
    else:
        if verbose:
            print("No isoforms found with sufficient reads")
    
    return df


def create_proportion_summary(df: pd.DataFrame) -> Dict[str, any]:
    """
    Create summary statistics for isoform proportions.
    
    Args:
        df: DataFrame with isoform proportions
    
    Returns:
        Dictionary with summary statistics
    """
    summary = {
        'total_isoforms': len(df),
        'total_reads': df['total_reads'].sum(),
        'isoforms_with_multiple_samples': len(df[df['samples_with_reads'] > 1]),
        'isoforms_single_sample': len(df[df['samples_with_reads'] == 1]),
        'mean_reads_per_isoform': df['total_reads'].mean(),
        'median_reads_per_isoform': df['total_reads'].median(),
        'max_reads_per_isoform': df['total_reads'].max(),
        'min_reads_per_isoform': df['total_reads'].min(),
    }
    
    return summary
