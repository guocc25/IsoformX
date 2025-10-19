"""Unit tests for IsoformX."""

import pytest
import numpy as np
from typing import List
from src.isoformx.models import ReadSignature, Isoform, SpliceSite, Cluster
from src.isoformx.utils import (
    weighted_median,
    jaccard_distance,
    splice_site_distance,
    compute_confidence_score
)
from src.isoformx.cluster import cluster_read_signatures, ClusterParams
from src.isoformx.consensus import build_consensus_isoforms, ConsensusParams


class TestUtils:
    """Test utility functions."""
    
    def test_weighted_median(self):
        """Test weighted median computation."""
        values = [1, 2, 3, 4, 5]
        weights = [1, 1, 1, 1, 1]
        
        result = weighted_median(values, weights)
        assert result == 3.0
        
        # Test with different weights
        weights = [1, 2, 1, 1, 1]
        result = weighted_median(values, weights)
        assert result == 2.0
    
    def test_jaccard_distance(self):
        """Test Jaccard distance computation."""
        set1 = {1, 2, 3}
        set2 = {2, 3, 4}
        
        result = jaccard_distance(set1, set2)
        expected = 1.0 - (2 / 4)  # 2 common, 4 total
        assert abs(result - expected) < 1e-6
        
        # Test empty sets
        result = jaccard_distance(set(), set())
        assert result == 0.0
    
    def test_splice_site_distance(self):
        """Test splice site distance computation."""
        sites1 = [100, 200, 300]
        sites2 = [105, 205, 305]  # Within tolerance
        
        result = splice_site_distance(sites1, sites2, tolerance=10)
        assert result == 0.0  # All sites match
        
        sites2 = [100, 200, 400]  # One site different
        result = splice_site_distance(sites1, sites2, tolerance=10)
        assert result == 1.0 / 3.0  # 1 out of 3 sites match
    
    def test_compute_confidence_score(self):
        """Test confidence score computation."""
        support_reads = 5
        coverage = 0.8
        quality_scores = [0.9, 0.8, 0.7, 0.6, 0.5]
        splice_consensus = 0.9
        
        result = compute_confidence_score(
            support_reads, coverage, quality_scores, splice_consensus
        )
        
        assert 0.0 <= result <= 1.0
        assert result > 0.5  # Should be high confidence


class TestModels:
    """Test data models."""
    
    def test_read_signature(self):
        """Test ReadSignature model."""
        splice_sites = [
            SpliceSite(position=100, strand="+", donor_acceptor="donor"),
            SpliceSite(position=200, strand="+", donor_acceptor="acceptor")
        ]
        
        read = ReadSignature(
            read_id="test_read",
            transcript_id="test_tx",
            gene_id="test_gene",
            chromosome="chr1",
            strand="+",
            start=50,
            end=250,
            splice_sites=splice_sites,
            coverage=0.8,
            quality_score=0.9
        )
        
        assert read.read_id == "test_read"
        assert read.exon_count == 3  # 2 splice sites + 1
        assert read.intron_count == 2
        
        exons = read.get_exon_bounds()
        assert len(exons) == 3
        assert exons[0] == (50, 100)
        assert exons[1] == (101, 200)
        assert exons[2] == (201, 250)
    
    def test_isoform(self):
        """Test Isoform model."""
        splice_sites = [
            SpliceSite(position=100, strand="+", donor_acceptor="donor"),
            SpliceSite(position=200, strand="+", donor_acceptor="acceptor")
        ]
        
        isoform = Isoform(
            isoform_id="test_isoform",
            gene_id="test_gene",
            chromosome="chr1",
            strand="+",
            start=50,
            end=250,
            splice_sites=splice_sites
        )
        
        assert isoform.isoform_id == "test_isoform"
        assert isoform.exon_count == 3
        assert isoform.intron_count == 2
        
        exons = isoform.get_exon_bounds()
        assert len(exons) == 3


class TestClustering:
    """Test clustering functionality."""
    
    def test_cluster_reads(self):
        """Test read clustering."""
        # Create test reads
        reads = []
        for i in range(10):
            read = ReadSignature(
                read_id=f"read_{i}",
                transcript_id=f"tx_{i}",
                gene_id="test_gene",
                chromosome="chr1",
                strand="+",
                start=100 + i * 10,
                end=200 + i * 10,
                coverage=0.8,
                quality_score=0.9
            )
            reads.append(read)
        
        # Cluster reads
        params = ClusterParams(max_distance=0.1, min_reads=2)
        clusters = cluster_read_signatures(reads, params)
        
        assert len(clusters) > 0
        assert all(len(cluster.read_signatures) >= params.min_reads for cluster in clusters)
    
    def test_cluster_params(self):
        """Test cluster parameters."""
        params = ClusterParams()
        
        assert params.max_distance == 0.2
        assert params.min_reads == 2
        assert params.splice_site_tolerance == 10


class TestConsensus:
    """Test consensus building."""
    
    def test_build_consensus(self):
        """Test consensus building."""
        # Create test cluster
        reads = []
        for i in range(5):
            read = ReadSignature(
                read_id=f"read_{i}",
                transcript_id=f"tx_{i}",
                gene_id="test_gene",
                chromosome="chr1",
                strand="+",
                start=100 + i,
                end=200 + i,
                coverage=0.8,
                quality_score=0.9
            )
            reads.append(read)
        
        cluster = Cluster(
            cluster_id="test_cluster",
            read_signatures=reads
        )
        
        # Build consensus
        params = ConsensusParams()
        isoforms = build_consensus_isoforms([cluster], params)
        
        assert len(isoforms) > 0
        assert all(isoform.confidence_score >= 0.0 for isoform in isoforms)
    
    def test_consensus_params(self):
        """Test consensus parameters."""
        params = ConsensusParams()
        
        assert params.min_coverage == 0.5
        assert params.min_support_reads == 2
        assert params.splice_site_threshold == 0.7


class TestIntegration:
    """Test integration scenarios."""
    
    def test_mock_data_pipeline(self):
        """Test complete pipeline with mock data."""
        from src.isoformx.io import load_mock_data, MockDataParams
        
        # Generate mock data
        params = MockDataParams(n_reads=100, n_genes=10)
        read_signatures, reference_transcripts = load_mock_data(params)
        
        assert len(read_signatures) == 100
        assert len(reference_transcripts) == 30  # 10 genes * 3 transcripts
        
        # Run clustering
        clusters = cluster_read_signatures(read_signatures)
        assert len(clusters) > 0
        
        # Build consensus
        isoforms = build_consensus_isoforms(clusters)
        assert len(isoforms) > 0
        
        # Check that isoforms have reasonable properties
        for isoform in isoforms:
            assert isoform.confidence_score >= 0.0
            assert isoform.confidence_score <= 1.0
            assert isoform.support_reads > 0
            assert isoform.coverage >= 0.0
    
    def test_edge_cases(self):
        """Test edge cases."""
        # Empty input
        clusters = cluster_read_signatures([])
        assert len(clusters) == 0
        
        # Single read
        read = ReadSignature(
            read_id="single_read",
            transcript_id="single_tx",
            gene_id="single_gene",
            chromosome="chr1",
            strand="+",
            start=100,
            end=200,
            coverage=0.8,
            quality_score=0.9
        )
        
        clusters = cluster_read_signatures([read])
        assert len(clusters) == 0  # Should be filtered out due to min_reads=2


if __name__ == "__main__":
    pytest.main([__file__])
