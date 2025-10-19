"""Consensus building for isoform identification."""

from typing import List, Optional, Tuple, Dict
import numpy as np
from dataclasses import dataclass
from .models import Cluster, Isoform, ReadSignature, SpliceSite
from .utils import weighted_median, compute_confidence_score


@dataclass
class ConsensusParams:
    """Parameters for consensus building."""
    min_coverage: float = 0.5
    min_support_reads: int = 2
    splice_site_threshold: float = 0.7
    position_tolerance: int = 10
    quality_threshold: float = 0.3


class ConsensusBuilder:
    """Builds consensus isoforms from clusters."""
    
    def __init__(self, params: ConsensusParams):
        self.params = params
    
    def build_consensus(self, cluster: Cluster) -> Optional[Isoform]:
        """
        Build consensus isoform from a cluster.
        
        Args:
            cluster: Cluster of read signatures
        
        Returns:
            Consensus isoform or None if insufficient support
        """
        if len(cluster.read_signatures) < self.params.min_support_reads:
            return None
        
        # Compute consensus boundaries
        start, end = self._consensus_boundaries(cluster.read_signatures)
        
        # Compute consensus splice sites
        consensus_sites = self._consensus_splice_sites(cluster.read_signatures)
        
        # Create consensus isoform
        isoform = Isoform(
            isoform_id=f"isoform_{cluster.cluster_id}",
            gene_id=cluster.read_signatures[0].gene_id,
            chromosome=cluster.read_signatures[0].chromosome,
            strand=cluster.read_signatures[0].strand,
            start=start,
            end=end,
            splice_sites=consensus_sites,
            read_signatures=cluster.read_signatures
        )
        
        # Compute confidence score
        confidence = self._compute_isoform_confidence(isoform)
        isoform.confidence_score = confidence
        
        return isoform
    
    def _consensus_boundaries(self, reads: List[ReadSignature]) -> Tuple[int, int]:
        """Compute consensus start and end positions."""
        starts = [r.start for r in reads]
        ends = [r.end for r in reads]
        weights = [r.coverage for r in reads]
        
        # Use weighted median for robustness
        consensus_start = int(weighted_median(starts, weights))
        consensus_end = int(weighted_median(ends, weights))
        
        return consensus_start, consensus_end
    
    def _consensus_splice_sites(self, reads: List[ReadSignature]) -> List[SpliceSite]:
        """Compute consensus splice sites from reads."""
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
        tolerance = self.params.position_tolerance
        
        # Sort by position
        all_sites.sort(key=lambda x: x.position)
        
        current_group = [all_sites[0]]
        for site in all_sites[1:]:
            if abs(site.position - current_group[-1].position) <= tolerance:
                current_group.append(site)
            else:
                # Process current group
                consensus_site = self._consensus_site_group(current_group)
                if consensus_site:
                    consensus_sites.append(consensus_site)
                current_group = [site]
        
        # Process final group
        consensus_site = self._consensus_site_group(current_group)
        if consensus_site:
            consensus_sites.append(consensus_site)
        
        return consensus_sites
    
    def _consensus_site_group(self, sites: List[SpliceSite]) -> Optional[SpliceSite]:
        """Compute consensus splice site from a group."""
        if not sites:
            return None
        
        # Check if group meets threshold
        support_ratio = len(sites) / len(sites)
        if support_ratio < self.params.splice_site_threshold:
            return None
        
        # Use median position
        positions = [s.position for s in sites]
        median_pos = int(np.median(positions))
        
        # Use most common strand and type
        strands = [s.strand for s in sites]
        types = [s.donor_acceptor for s in sites]
        
        consensus_strand = max(set(strands), key=strands.count)
        consensus_type = max(set(types), key=types.count)
        
        # Average confidence weighted by support
        weights = [s.confidence for s in sites]
        avg_confidence = np.average(weights)
        
        return SpliceSite(
            position=median_pos,
            strand=consensus_strand,
            donor_acceptor=consensus_type,
            confidence=avg_confidence
        )
    
    def _compute_isoform_confidence(self, isoform: Isoform) -> float:
        """Compute confidence score for isoform."""
        if not isoform.read_signatures:
            return 0.0
        
        # Extract metrics
        support_reads = len(isoform.read_signatures)
        coverage = isoform.coverage
        quality_scores = [r.quality_score for r in isoform.read_signatures]
        
        # Splice site consensus
        splice_consensus = self._compute_splice_consensus(isoform)
        
        return compute_confidence_score(
            support_reads=support_reads,
            coverage=coverage,
            quality_scores=quality_scores,
            splice_site_consensus=splice_consensus
        )
    
    def _compute_splice_consensus(self, isoform: Isoform) -> float:
        """Compute splice site consensus score."""
        if not isoform.splice_sites:
            return 1.0  # Mono-exon isoforms
        
        consensus_scores = []
        for site in isoform.splice_sites:
            # Count supporting reads for this site
            supporting_reads = 0
            for read in isoform.read_signatures:
                for read_site in read.splice_sites:
                    if (abs(read_site.position - site.position) <= self.params.position_tolerance and
                        read_site.donor_acceptor == site.donor_acceptor):
                        supporting_reads += 1
                        break
            
            # Compute consensus score
            consensus_score = supporting_reads / len(isoform.read_signatures)
            consensus_scores.append(consensus_score)
        
        return np.mean(consensus_scores) if consensus_scores else 0.0


class MultiLocusConsensusBuilder:
    """Builds consensus isoforms across multiple loci."""
    
    def __init__(self, params: ConsensusParams):
        self.params = params
        self.builder = ConsensusBuilder(params)
    
    def build_consensus_isoforms(self, clusters: List[Cluster]) -> List[Isoform]:
        """
        Build consensus isoforms from clusters.
        
        Args:
            clusters: List of clusters
        
        Returns:
            List of consensus isoforms
        """
        isoforms = []
        
        for cluster in clusters:
            isoform = self.builder.build_consensus(cluster)
            if isoform and self._meets_quality_threshold(isoform):
                isoforms.append(isoform)
        
        return isoforms
    
    def _meets_quality_threshold(self, isoform: Isoform) -> bool:
        """Check if isoform meets quality thresholds."""
        return (
            isoform.confidence_score >= self.params.quality_threshold and
            isoform.coverage >= self.params.min_coverage and
            isoform.support_reads >= self.params.min_support_reads
        )
    
    def merge_similar_isoforms(self, isoforms: List[Isoform]) -> List[Isoform]:
        """
        Merge similar isoforms to reduce redundancy.
        
        Args:
            isoforms: List of isoforms
        
        Returns:
            List of merged isoforms
        """
        if not isoforms:
            return []
        
        # Group isoforms by gene
        gene_groups = {}
        for isoform in isoforms:
            gene_id = isoform.gene_id
            if gene_id not in gene_groups:
                gene_groups[gene_id] = []
            gene_groups[gene_id].append(isoform)
        
        merged_isoforms = []
        
        for gene_id, gene_isoforms in gene_groups.items():
            merged = self._merge_gene_isoforms(gene_isoforms)
            merged_isoforms.extend(merged)
        
        return merged_isoforms
    
    def _merge_gene_isoforms(self, isoforms: List[Isoform]) -> List[Isoform]:
        """Merge isoforms within a gene."""
        if len(isoforms) <= 1:
            return isoforms
        
        # Sort by confidence score
        isoforms.sort(key=lambda x: x.confidence_score, reverse=True)
        
        merged = []
        used = set()
        
        for i, isoform in enumerate(isoforms):
            if i in used:
                continue
            
            # Find similar isoforms to merge
            similar_indices = [i]
            for j in range(i + 1, len(isoforms)):
                if j in used:
                    continue
                
                if self._are_similar(isoform, isoforms[j]):
                    similar_indices.append(j)
                    used.add(j)
            
            # Merge similar isoforms
            if len(similar_indices) > 1:
                merged_isoform = self._merge_isoform_group([isoforms[k] for k in similar_indices])
                merged.append(merged_isoform)
            else:
                merged.append(isoform)
            
            used.add(i)
        
        return merged
    
    def _are_similar(self, isoform1: Isoform, isoform2: Isoform) -> bool:
        """Check if two isoforms are similar enough to merge."""
        # Check basic compatibility
        if (isoform1.chromosome != isoform2.chromosome or
            isoform1.strand != isoform2.strand):
            return False
        
        # Check position overlap
        overlap_start = max(isoform1.start, isoform2.start)
        overlap_end = min(isoform1.end, isoform2.end)
        
        if overlap_start >= overlap_end:
            return False
        
        overlap_length = overlap_end - overlap_start + 1
        total_length = max(isoform1.end - isoform1.start, isoform2.end - isoform2.start) + 1
        overlap_ratio = overlap_length / total_length
        
        if overlap_ratio < 0.8:
            return False
        
        # Check splice site similarity
        sites1 = [s.position for s in isoform1.splice_sites]
        sites2 = [s.position for s in isoform2.splice_sites]
        
        if len(sites1) != len(sites2):
            return False
        
        # Check if splice sites are within tolerance
        for site1, site2 in zip(sites1, sites2):
            if abs(site1 - site2) > self.params.position_tolerance:
                return False
        
        return True
    
    def _merge_isoform_group(self, isoforms: List[Isoform]) -> Isoform:
        """Merge a group of similar isoforms."""
        if not isoforms:
            raise ValueError("Cannot merge empty isoform group")
        
        if len(isoforms) == 1:
            return isoforms[0]
        
        # Use the highest confidence isoform as base
        base_isoform = max(isoforms, key=lambda x: x.confidence_score)
        
        # Merge read signatures
        all_reads = []
        for isoform in isoforms:
            all_reads.extend(isoform.read_signatures)
        
        # Remove duplicates
        unique_reads = []
        seen_ids = set()
        for read in all_reads:
            if read.read_id not in seen_ids:
                unique_reads.append(read)
                seen_ids.add(read.read_id)
        
        # Create merged isoform
        merged_isoform = Isoform(
            isoform_id=base_isoform.isoform_id,
            gene_id=base_isoform.gene_id,
            chromosome=base_isoform.chromosome,
            strand=base_isoform.strand,
            start=base_isoform.start,
            end=base_isoform.end,
            splice_sites=base_isoform.splice_sites,
            read_signatures=unique_reads,
            consensus_sequence=base_isoform.consensus_sequence,
            isoform_type=base_isoform.isoform_type
        )
        
        # Recompute confidence score
        merged_isoform.confidence_score = self.builder._compute_isoform_confidence(merged_isoform)
        
        return merged_isoform


def build_consensus_isoforms(
    clusters: List[Cluster],
    params: Optional[ConsensusParams] = None
) -> List[Isoform]:
    """
    Main function to build consensus isoforms from clusters.
    
    Args:
        clusters: List of clusters
        params: Consensus parameters (uses defaults if None)
    
    Returns:
        List of consensus isoforms
    """
    if params is None:
        params = ConsensusParams()
    
    builder = MultiLocusConsensusBuilder(params)
    
    # Build consensus isoforms
    isoforms = builder.build_consensus_isoforms(clusters)
    
    # Merge similar isoforms
    merged_isoforms = builder.merge_similar_isoforms(isoforms)
    
    return merged_isoforms
