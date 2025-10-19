#!/usr/bin/env python3
"""Test script to debug the isoform pipeline."""

import sys
import os
sys.path.insert(0, os.path.join(os.path.dirname(__file__), 'src'))

from isoformx.io import load_mock_data, MockDataParams
from isoformx.cluster import cluster_read_signatures, ClusterParams
from isoformx.consensus import build_consensus_isoforms, ConsensusParams
from isoformx.score import filter_by_confidence, ScoringParams

def test_pipeline():
    """Test the pipeline with debug output."""
    print("=== Testing IsoformX Pipeline ===")
    
    # Generate mock data with relaxed parameters
    params = MockDataParams(n_reads=50, n_genes=10)
    read_signatures, reference_transcripts = load_mock_data(params)
    
    print(f"Generated {len(read_signatures)} reads and {len(reference_transcripts)} reference transcripts")
    
    # Use very relaxed clustering parameters
    cluster_params = ClusterParams(
        max_distance=0.8,  # Very loose
        min_reads=2,
        splice_site_tolerance=50,
        position_tolerance=100
    )
    
    print(f"Clustering with max_distance={cluster_params.max_distance}, min_reads={cluster_params.min_reads}")
    
    # Cluster reads
    clusters = cluster_read_signatures(read_signatures, cluster_params)
    print(f"Found {len(clusters)} clusters")
    
    for i, cluster in enumerate(clusters):
        print(f"  Cluster {i}: {len(cluster.read_signatures)} reads")
    
    # Build consensus isoforms with relaxed parameters
    consensus_params = ConsensusParams(
        min_coverage=0.1,  # Very low threshold
        min_support_reads=2,
        splice_site_threshold=0.3,  # Low threshold
        quality_threshold=0.1
    )
    
    print(f"Building consensus with min_coverage={consensus_params.min_coverage}")
    
    isoforms = build_consensus_isoforms(clusters, consensus_params)
    print(f"Built {len(isoforms)} consensus isoforms")
    
    # Filter by confidence with low threshold
    scoring_params = ScoringParams(
        min_reads=2,
        min_coverage=0.1,
        min_quality=0.1
    )
    
    isoforms = filter_by_confidence(isoforms, min_confidence=0.1, params=scoring_params)
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
        print()

if __name__ == "__main__":
    test_pipeline()
