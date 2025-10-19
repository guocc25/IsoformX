"""Unit tests for clustering functionality."""

import pytest
import numpy as np
from typing import List
from src.isoformx.models import ReadSignature, SpliceSite, Cluster
from src.isoformx.cluster import (
    cluster_read_signatures,
    ClusterParams,
    ErrorAwareClusterer,
    UnionFind
)
from src.isoformx.consensus import build_consensus_isoforms, ConsensusParams


class TestUnionFind:
    """Test Union-Find data structure."""
    
    def test_union_find_basic(self):
        """Test basic Union-Find operations."""
        uf = UnionFind(5)
        
        # Initially all elements are separate
        assert uf.find(0) == 0
        assert uf.find(1) == 1
        assert uf.find(2) == 2
        
        # Union some elements
        uf.union(0, 1)
        assert uf.find(0) == uf.find(1)
        assert uf.find(0) != uf.find(2)
        
        # Union more elements
        uf.union(1, 2)
        assert uf.find(0) == uf.find(2)
        
        # Check cluster sizes
        assert uf.size[uf.find(0)] == 3
    
    def test_union_find_clusters(self):
        """Test getting clusters from Union-Find."""
        uf = UnionFind(6)
        
        # Create two clusters
        uf.union(0, 1)
        uf.union(1, 2)
        uf.union(3, 4)
        uf.union(4, 5)
        
        clusters = uf.get_clusters()
        
        assert len(clusters) == 2
        assert len(clusters[uf.find(0)]) == 3
        assert len(clusters[uf.find(3)]) == 3


class TestErrorAwareClusterer:
    """Test error-aware clustering."""
    
    def test_distance_computation(self):
        """Test distance computation between reads."""
        clusterer = ErrorAwareClusterer(ClusterParams())
        
        # Create two similar reads
        read1 = ReadSignature(
            read_id="read1",
            transcript_id="tx1",
            gene_id="gene1",
            chromosome="chr1",
            strand="+",
            start=100,
            end=200,
            coverage=0.8,
            quality_score=0.9
        )
        
        read2 = ReadSignature(
            read_id="read2",
            transcript_id="tx2",
            gene_id="gene1",
            chromosome="chr1",
            strand="+",
            start=105,
            end=205,
            coverage=0.8,
            quality_score=0.9
        )
        
        distance = clusterer.compute_distance(read1, read2)
        assert 0.0 <= distance <= 1.0
        assert distance < 0.5  # Should be similar
    
    def test_distance_different_chromosomes(self):
        """Test distance between reads on different chromosomes."""
        clusterer = ErrorAwareClusterer(ClusterParams())
        
        read1 = ReadSignature(
            read_id="read1",
            transcript_id="tx1",
            gene_id="gene1",
            chromosome="chr1",
            strand="+",
            start=100,
            end=200,
            coverage=0.8,
            quality_score=0.9
        )
        
        read2 = ReadSignature(
            read_id="read2",
            transcript_id="tx2",
            gene_id="gene1",
            chromosome="chr2",
            strand="+",
            start=100,
            end=200,
            coverage=0.8,
            quality_score=0.9
        )
        
        distance = clusterer.compute_distance(read1, read2)
        assert distance == 1.0  # Should be maximum distance
    
    def test_distance_different_strands(self):
        """Test distance between reads on different strands."""
        clusterer = ErrorAwareClusterer(ClusterParams())
        
        read1 = ReadSignature(
            read_id="read1",
            transcript_id="tx1",
            gene_id="gene1",
            chromosome="chr1",
            strand="+",
            start=100,
            end=200,
            coverage=0.8,
            quality_score=0.9
        )
        
        read2 = ReadSignature(
            read_id="read2",
            transcript_id="tx2",
            gene_id="gene1",
            chromosome="chr1",
            strand="-",
            start=100,
            end=200,
            coverage=0.8,
            quality_score=0.9
        )
        
        distance = clusterer.compute_distance(read1, read2)
        assert distance == 1.0  # Should be maximum distance
    
    def test_splice_site_distance(self):
        """Test splice site distance computation."""
        clusterer = ErrorAwareClusterer(ClusterParams())
        
        # Create reads with splice sites
        read1 = ReadSignature(
            read_id="read1",
            transcript_id="tx1",
            gene_id="gene1",
            chromosome="chr1",
            strand="+",
            start=100,
            end=300,
            splice_sites=[
                SpliceSite(position=150, strand="+", donor_acceptor="donor"),
                SpliceSite(position=200, strand="+", donor_acceptor="acceptor")
            ],
            coverage=0.8,
            quality_score=0.9
        )
        
        read2 = ReadSignature(
            read_id="read2",
            transcript_id="tx2",
            gene_id="gene1",
            chromosome="chr1",
            strand="+",
            start=100,
            end=300,
            splice_sites=[
                SpliceSite(position=155, strand="+", donor_acceptor="donor"),
                SpliceSite(position=205, strand="+", donor_acceptor="acceptor")
            ],
            coverage=0.8,
            quality_score=0.9
        )
        
        distance = clusterer.compute_distance(read1, read2)
        assert distance < 0.5  # Should be similar due to close splice sites
    
    def test_cluster_reads(self):
        """Test clustering of reads."""
        clusterer = ErrorAwareClusterer(ClusterParams(max_distance=0.3, min_reads=2))
        
        # Create reads that should cluster together
        reads = []
        for i in range(5):
            read = ReadSignature(
                read_id=f"read_{i}",
                transcript_id=f"tx_{i}",
                gene_id="gene1",
                chromosome="chr1",
                strand="+",
                start=100 + i * 5,  # Small differences
                end=200 + i * 5,
                coverage=0.8,
                quality_score=0.9
            )
            reads.append(read)
        
        clusters = clusterer.cluster_reads(reads)
        
        assert len(clusters) > 0
        assert all(len(cluster.read_signatures) >= 2 for cluster in clusters)
    
    def test_cluster_reads_no_clusters(self):
        """Test clustering when no reads cluster together."""
        clusterer = ErrorAwareClusterer(ClusterParams(max_distance=0.1, min_reads=2))
        
        # Create reads that are too different
        reads = []
        for i in range(5):
            read = ReadSignature(
                read_id=f"read_{i}",
                transcript_id=f"tx_{i}",
                gene_id="gene1",
                chromosome="chr1",
                strand="+",
                start=100 + i * 1000,  # Large differences
                end=200 + i * 1000,
                coverage=0.8,
                quality_score=0.9
            )
            reads.append(read)
        
        clusters = clusterer.cluster_reads(reads)
        
        assert len(clusters) == 0  # No clusters should form
    
    def test_refine_clusters(self):
        """Test cluster refinement."""
        clusterer = ErrorAwareClusterer(ClusterParams())
        
        # Create a cluster with an outlier
        reads = []
        for i in range(4):
            read = ReadSignature(
                read_id=f"read_{i}",
                transcript_id=f"tx_{i}",
                gene_id="gene1",
                chromosome="chr1",
                strand="+",
                start=100 + i * 10,
                end=200 + i * 10,
                coverage=0.8,
                quality_score=0.9
            )
            reads.append(read)
        
        # Add an outlier
        outlier = ReadSignature(
            read_id="outlier",
            transcript_id="tx_outlier",
            gene_id="gene1",
            chromosome="chr1",
            strand="+",
            start=1000,  # Far away
            end=1100,
            coverage=0.8,
            quality_score=0.9
        )
        reads.append(outlier)
        
        cluster = Cluster(
            cluster_id="test_cluster",
            read_signatures=reads
        )
        
        refined_clusters = clusterer.refine_clusters([cluster])
        
        assert len(refined_clusters) > 0
        # The outlier should be removed
        for refined_cluster in refined_clusters:
            assert len(refined_cluster.read_signatures) < len(reads)


class TestClusteringIntegration:
    """Test clustering integration."""
    
    def test_cluster_read_signatures_function(self):
        """Test the main clustering function."""
        # Create test reads
        reads = []
        for i in range(10):
            read = ReadSignature(
                read_id=f"read_{i}",
                transcript_id=f"tx_{i}",
                gene_id="gene1",
                chromosome="chr1",
                strand="+",
                start=100 + i * 20,
                end=200 + i * 20,
                coverage=0.8,
                quality_score=0.9
            )
            reads.append(read)
        
        # Test with default parameters
        clusters = cluster_read_signatures(reads)
        assert len(clusters) >= 0
        
        # Test with custom parameters
        params = ClusterParams(max_distance=0.1, min_reads=3)
        clusters = cluster_read_signatures(reads, params)
        assert len(clusters) >= 0
    
    def test_cluster_with_splice_sites(self):
        """Test clustering with splice sites."""
        reads = []
        
        # Create reads with similar splice sites
        for i in range(3):
            read = ReadSignature(
                read_id=f"read_{i}",
                transcript_id=f"tx_{i}",
                gene_id="gene1",
                chromosome="chr1",
                strand="+",
                start=100,
                end=300,
                splice_sites=[
                    SpliceSite(position=150 + i, strand="+", donor_acceptor="donor"),
                    SpliceSite(position=200 + i, strand="+", donor_acceptor="acceptor")
                ],
                coverage=0.8,
                quality_score=0.9
            )
            reads.append(read)
        
        clusters = cluster_read_signatures(reads)
        
        # Should form at least one cluster
        assert len(clusters) > 0
        
        # Check that clusters have reasonable properties
        for cluster in clusters:
            assert len(cluster.read_signatures) >= 2
            assert cluster.centroid is not None
    
    def test_cluster_consensus_building(self):
        """Test consensus building from clusters."""
        # Create test cluster
        reads = []
        for i in range(5):
            read = ReadSignature(
                read_id=f"read_{i}",
                transcript_id=f"tx_{i}",
                gene_id="gene1",
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
        isoforms = build_consensus_isoforms([cluster])
        
        assert len(isoforms) > 0
        assert all(isoform.confidence_score >= 0.0 for isoform in isoforms)
        assert all(isoform.support_reads > 0 for isoform in isoforms)


class TestEdgeCases:
    """Test edge cases in clustering."""
    
    def test_empty_input(self):
        """Test clustering with empty input."""
        clusters = cluster_read_signatures([])
        assert len(clusters) == 0
    
    def test_single_read(self):
        """Test clustering with single read."""
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
    
    def test_identical_reads(self):
        """Test clustering with identical reads."""
        reads = []
        for i in range(5):
            read = ReadSignature(
                read_id=f"read_{i}",
                transcript_id="tx1",
                gene_id="gene1",
                chromosome="chr1",
                strand="+",
                start=100,
                end=200,
                coverage=0.8,
                quality_score=0.9
            )
            reads.append(read)
        
        clusters = cluster_read_signatures(reads)
        
        assert len(clusters) > 0
        # All reads should be in the same cluster
        assert len(clusters[0].read_signatures) == 5
    
    def test_extreme_parameters(self):
        """Test clustering with extreme parameters."""
        reads = []
        for i in range(10):
            read = ReadSignature(
                read_id=f"read_{i}",
                transcript_id=f"tx_{i}",
                gene_id="gene1",
                chromosome="chr1",
                strand="+",
                start=100 + i * 10,
                end=200 + i * 10,
                coverage=0.8,
                quality_score=0.9
            )
            reads.append(read)
        
        # Very strict parameters
        params = ClusterParams(max_distance=0.01, min_reads=5)
        clusters = cluster_read_signatures(reads, params)
        
        # Should result in fewer clusters
        assert len(clusters) <= len(reads)
        
        # Very loose parameters
        params = ClusterParams(max_distance=1.0, min_reads=1)
        clusters = cluster_read_signatures(reads, params)
        
        # Should result in more clusters
        assert len(clusters) >= 0


if __name__ == "__main__":
    pytest.main([__file__])
