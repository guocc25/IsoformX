"""Error-aware clustering using Union-Find algorithm."""

from typing import List, Dict, Set, Tuple, Optional
import numpy as np
from dataclasses import dataclass
from .models import ReadSignature, Cluster, Isoform
from .utils import (
    splice_site_distance,
    overlap_coefficient,
    jaccard_distance,
    weighted_median
)


@dataclass
class ClusterParams:
    """Parameters for clustering algorithm."""
    max_distance: float = 0.2
    min_reads: int = 2
    splice_site_tolerance: int = 10
    position_tolerance: int = 50
    coverage_weight: float = 0.3
    splice_weight: float = 0.4
    position_weight: float = 0.3


class UnionFind:
    """Union-Find data structure for clustering."""
    
    def __init__(self, n: int):
        self.parent = list(range(n))
        self.rank = [0] * n
        self.size = [1] * n
    
    def find(self, x: int) -> int:
        """Find root of x with path compression."""
        if self.parent[x] != x:
            self.parent[x] = self.find(self.parent[x])
        return self.parent[x]
    
    def union(self, x: int, y: int) -> bool:
        """Union sets containing x and y."""
        px, py = self.find(x), self.find(y)
        if px == py:
            return False
        
        # Union by rank
        if self.rank[px] < self.rank[py]:
            px, py = py, px
        
        self.parent[py] = px
        self.size[px] += self.size[py]
        
        if self.rank[px] == self.rank[py]:
            self.rank[px] += 1
        
        return True
    
    def get_clusters(self) -> Dict[int, List[int]]:
        """Get all clusters as dict mapping root to members."""
        clusters = {}
        for i in range(len(self.parent)):
            root = self.find(i)
            if root not in clusters:
                clusters[root] = []
            clusters[root].append(i)
        return clusters


class ErrorAwareClusterer:
    """Error-aware clustering for read signatures."""
    
    def __init__(self, params: ClusterParams):
        self.params = params
    
    def compute_distance(self, sig1: ReadSignature, sig2: ReadSignature) -> float:
        """
        Compute error-aware distance between two read signatures.
        
        Args:
            sig1: First read signature
            sig2: Second read signature
        
        Returns:
            Distance score (0-1, lower is more similar)
        """
        # Check basic compatibility
        if (sig1.chromosome != sig2.chromosome or 
            sig1.strand != sig2.strand):
            return 1.0
        
        # Position distance
        pos_dist = self._position_distance(sig1, sig2)
        
        # Splice site distance
        splice_dist = self._splice_site_distance(sig1, sig2)
        
        # Coverage similarity
        coverage_sim = self._coverage_similarity(sig1, sig2)
        
        # Weighted combination
        distance = (
            self.params.position_weight * pos_dist +
            self.params.splice_weight * splice_dist +
            (1 - coverage_sim) * self.params.coverage_weight
        )
        
        return min(distance, 1.0)
    
    def _position_distance(self, sig1: ReadSignature, sig2: ReadSignature) -> float:
        """Compute normalized position distance."""
        # Start position distance
        start_dist = abs(sig1.start - sig2.start)
        
        # End position distance
        end_dist = abs(sig1.end - sig2.end)
        
        # Length difference
        len1 = sig1.end - sig1.start
        len2 = sig2.end - sig2.start
        len_dist = abs(len1 - len2)
        
        # Normalize by transcript length
        max_len = max(len1, len2)
        if max_len == 0:
            return 0.0
        
        normalized_dist = (start_dist + end_dist + len_dist) / (3 * max_len)
        return min(normalized_dist, 1.0)
    
    def _splice_site_distance(self, sig1: ReadSignature, sig2: ReadSignature) -> float:
        """Compute splice site distance."""
        sites1 = [s.position for s in sig1.splice_sites]
        sites2 = [s.position for s in sig2.splice_sites]
        
        return splice_site_distance(sites1, sites2, self.params.splice_site_tolerance)
    
    def _coverage_similarity(self, sig1: ReadSignature, sig2: ReadSignature) -> float:
        """Compute coverage similarity."""
        # Simple coverage overlap
        overlap_start = max(sig1.start, sig2.start)
        overlap_end = min(sig1.end, sig2.end)
        
        if overlap_start >= overlap_end:
            return 0.0
        
        overlap_length = overlap_end - overlap_start + 1
        total_length = max(sig1.end - sig1.start, sig2.end - sig2.start) + 1
        
        return overlap_length / total_length
    
    def cluster_reads(self, read_signatures: List[ReadSignature]) -> List[Cluster]:
        """
        Cluster read signatures using error-aware distance.
        
        Args:
            read_signatures: List of read signatures to cluster
        
        Returns:
            List of clusters
        """
        if not read_signatures:
            return []
        
        n = len(read_signatures)
        uf = UnionFind(n)
        
        # Compute pairwise distances and union similar reads
        for i in range(n):
            for j in range(i + 1, n):
                distance = self.compute_distance(read_signatures[i], read_signatures[j])
                
                if distance <= self.params.max_distance:
                    uf.union(i, j)
        
        # Extract clusters
        cluster_dict = uf.get_clusters()
        clusters = []
        
        for cluster_id, indices in cluster_dict.items():
            if len(indices) >= self.params.min_reads:
                cluster_reads = [read_signatures[i] for i in indices]
                cluster = Cluster(
                    cluster_id=f"cluster_{cluster_id}",
                    read_signatures=cluster_reads
                )
                clusters.append(cluster)
        
        return clusters
    
    def refine_clusters(self, clusters: List[Cluster]) -> List[Cluster]:
        """
        Refine clusters by removing outliers and splitting if needed.
        
        Args:
            clusters: Initial clusters
        
        Returns:
            Refined clusters
        """
        refined_clusters = []
        
        for cluster in clusters:
            if len(cluster.read_signatures) < self.params.min_reads:
                continue
            
            # Remove outliers
            filtered_reads = self._remove_outliers(cluster.read_signatures)
            
            if len(filtered_reads) >= self.params.min_reads:
                refined_cluster = Cluster(
                    cluster_id=cluster.cluster_id,
                    read_signatures=filtered_reads
                )
                refined_clusters.append(refined_cluster)
        
        return refined_clusters
    
    def _remove_outliers(self, reads: List[ReadSignature]) -> List[ReadSignature]:
        """Remove outliers from a cluster."""
        if len(reads) <= 2:
            return reads
        
        # Compute distances to centroid
        centroid = self._compute_centroid(reads)
        distances = []
        
        for read in reads:
            dist = self.compute_distance(read, centroid)
            distances.append(dist)
        
        # Remove reads with distance > 2 * median distance
        median_dist = np.median(distances)
        threshold = 2 * median_dist
        
        filtered_reads = []
        for read, dist in zip(reads, distances):
            if dist <= threshold:
                filtered_reads.append(read)
        
        return filtered_reads
    
    def _compute_centroid(self, reads: List[ReadSignature]) -> ReadSignature:
        """Compute centroid read signature."""
        if not reads:
            raise ValueError("Cannot compute centroid for empty list")
        
        if len(reads) == 1:
            return reads[0]
        
        # Use weighted median for positions
        starts = [r.start for r in reads]
        ends = [r.end for r in reads]
        weights = [r.coverage for r in reads]
        
        centroid_start = int(weighted_median(starts, weights))
        centroid_end = int(weighted_median(ends, weights))
        
        # Consensus splice sites
        consensus_sites = self._consensus_splice_sites(reads)
        
        return ReadSignature(
            read_id=f"centroid_{len(reads)}",
            transcript_id=None,
            gene_id=reads[0].gene_id,
            chromosome=reads[0].chromosome,
            strand=reads[0].strand,
            start=centroid_start,
            end=centroid_end,
            splice_sites=consensus_sites,
            coverage=np.mean([r.coverage for r in reads]),
            quality_score=np.mean([r.quality_score for r in reads])
        )
    
    def _consensus_splice_sites(self, reads: List[ReadSignature]) -> List:
        """Compute consensus splice sites from reads."""
        from .models import SpliceSite
        
        if not reads:
            return []
        
        # Collect all splice sites
        all_sites = []
        for read in reads:
            all_sites.extend(read.splice_sites)
        
        if not all_sites:
            return []
        
        # Group sites by proximity
        consensus_sites = []
        tolerance = self.params.splice_site_tolerance
        
        # Sort by position
        all_sites.sort(key=lambda x: x.position)
        
        current_group = [all_sites[0]]
        for site in all_sites[1:]:
            if abs(site.position - current_group[-1].position) <= tolerance:
                current_group.append(site)
            else:
                # Process current group
                consensus_sites.append(self._consensus_site(current_group))
                current_group = [site]
        
        # Process final group
        consensus_sites.append(self._consensus_site(current_group))
        
        return consensus_sites
    
    def _consensus_site(self, sites: List) -> object:
        """Compute consensus splice site from a group."""
        from .models import SpliceSite
        
        if not sites:
            raise ValueError("Cannot compute consensus from empty site group")
        
        if len(sites) == 1:
            return sites[0]
        
        # Use median position
        positions = [s.position for s in sites]
        median_pos = int(np.median(positions))
        
        # Use most common strand and type
        strands = [s.strand for s in sites]
        types = [s.donor_acceptor for s in sites]
        
        consensus_strand = max(set(strands), key=strands.count)
        consensus_type = max(set(types), key=types.count)
        
        # Average confidence
        avg_confidence = np.mean([s.confidence for s in sites])
        
        return SpliceSite(
            position=median_pos,
            strand=consensus_strand,
            donor_acceptor=consensus_type,
            confidence=avg_confidence
        )


def cluster_read_signatures(
    read_signatures: List[ReadSignature],
    params: Optional[ClusterParams] = None
) -> List[Cluster]:
    """
    Main function to cluster read signatures.
    
    Args:
        read_signatures: List of read signatures to cluster
        params: Clustering parameters (uses defaults if None)
    
    Returns:
        List of clusters
    """
    if params is None:
        params = ClusterParams()
    
    clusterer = ErrorAwareClusterer(params)
    
    # Initial clustering
    clusters = clusterer.cluster_reads(read_signatures)
    
    # Refine clusters
    refined_clusters = clusterer.refine_clusters(clusters)
    
    return refined_clusters
