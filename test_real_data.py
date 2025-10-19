#!/usr/bin/env python3
"""Test script to run IsoformX on real BAM data."""

import sys
import os
sys.path.insert(0, os.path.join(os.path.dirname(__file__), 'src'))

from isoformx.io import load_bam_data
from isoformx.cluster import cluster_read_signatures, ClusterParams
from isoformx.consensus import build_consensus_isoforms, ConsensusParams
from isoformx.score import filter_by_confidence, ScoringParams
from isoformx.io import save_isoforms

def test_real_data(bam_file, output_file):
    """Test the pipeline with real BAM data."""
    print(f"=== Testing IsoformX with {bam_file} ===")
    
    # Load reads from BAM file
    print("Loading reads from BAM file...")
    read_signatures = load_bam_data(bam_file, min_length=1000)
    print(f"Loaded {len(read_signatures)} reads")
    
    if len(read_signatures) == 0:
        print("No reads found! Check BAM file format and read length filter.")
        return
    
    # Print some statistics
    print(f"Read length range: {min(r.read_length for r in read_signatures)} - {max(r.read_length for r in read_signatures)}")
    print(f"Chromosomes: {set(r.chromosome for r in read_signatures)}")
    print(f"Strands: {set(r.strand for r in read_signatures)}")
    
    # Use relaxed clustering parameters for real data
    cluster_params = ClusterParams(
        max_distance=0.3,  # Moderate distance
        min_reads=3,       # Require at least 3 reads
        splice_site_tolerance=20,
        position_tolerance=50
    )
    
    print(f"Clustering with max_distance={cluster_params.max_distance}, min_reads={cluster_params.min_reads}")
    
    # Cluster reads
    clusters = cluster_read_signatures(read_signatures, cluster_params)
    print(f"Found {len(clusters)} clusters")
    
    for i, cluster in enumerate(clusters):
        print(f"  Cluster {i}: {len(cluster.read_signatures)} reads")
        if cluster.centroid:
            print(f"    Centroid: {cluster.centroid.chromosome}:{cluster.centroid.start}-{cluster.centroid.end}")
    
    # Build consensus isoforms
    consensus_params = ConsensusParams(
        min_coverage=0.2,
        min_support_reads=3,
        splice_site_threshold=0.5,
        quality_threshold=0.2
    )
    
    print(f"Building consensus with min_coverage={consensus_params.min_coverage}")
    
    isoforms = build_consensus_isoforms(clusters, consensus_params)
    print(f"Built {len(isoforms)} consensus isoforms")
    
    # Filter by confidence
    scoring_params = ScoringParams(
        min_reads=3,
        min_coverage=0.2,
        min_quality=0.2
    )
    
    isoforms = filter_by_confidence(isoforms, min_confidence=0.3, params=scoring_params)
    print(f"Retained {len(isoforms)} high-confidence isoforms")
    
    # Print details of found isoforms
    for i, isoform in enumerate(isoforms):
        print(f"  Isoform {i}: {isoform.isoform_id}")
        print(f"    Gene: {isoform.gene_id}")
        print(f"    Position: {isoform.chromosome}:{isoform.start}-{isoform.end}")
        print(f"    Support reads: {isoform.support_reads}")
        print(f"    Confidence: {isoform.confidence_score:.3f}")
        print(f"    Coverage: {isoform.coverage:.3f}")
        print(f"    Exons: {isoform.exon_count}")
        print(f"    Splice sites: {len(isoform.splice_sites)}")
        print()
    
    # Save results
    if isoforms:
        save_isoforms(isoforms, output_file)
        print(f"Saved {len(isoforms)} isoforms to {output_file}")
    else:
        print("No isoforms found to save.")

if __name__ == "__main__":
    if len(sys.argv) != 3:
        print("Usage: python test_real_data.py <bam_file> <output_file>")
        sys.exit(1)
    
    bam_file = sys.argv[1]
    output_file = sys.argv[2]
    
    if not os.path.exists(bam_file):
        print(f"Error: BAM file {bam_file} not found")
        sys.exit(1)
    
    test_real_data(bam_file, output_file)
