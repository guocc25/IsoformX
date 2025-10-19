"""Data models for isoform identification."""

from dataclasses import dataclass, field
from typing import List, Optional, Tuple, Dict, Any
from enum import Enum


class IsoformType(Enum):
    """Isoform classification types."""
    FSM = "FSM"  # Full Splice Match
    ISM = "ISM"  # Incomplete Splice Match
    NIC = "NIC"  # Novel In Catalog
    NNC = "NNC"  # Novel Not in Catalog
    FSM_MONO = "FSM_MONO"  # FSM mono-exon
    ISM_MONO = "ISM_MONO"  # ISM mono-exon
    NIC_MONO = "NIC_MONO"  # NIC mono-exon
    NNC_MONO = "NNC_MONO"  # NNC mono-exon
    ANTISENSE = "ANTISENSE"  # Antisense
    INTERGENIC = "INTERGENIC"  # Intergenic


@dataclass
class SpliceSite:
    """Represents a splice site."""
    position: int
    strand: str
    donor_acceptor: str  # "donor" or "acceptor"
    confidence: float = 1.0


@dataclass
class ReadSignature:
    """Represents the signature of a long read for clustering."""
    read_id: str
    transcript_id: Optional[str]
    gene_id: Optional[str]
    chromosome: str
    strand: str
    start: int
    end: int
    splice_sites: List[SpliceSite] = field(default_factory=list)
    tss: Optional[int] = None  # Transcription start site
    tes: Optional[int] = None   # Transcription end site
    coverage: float = 1.0
    quality_score: float = 0.0
    read_length: int = 0
    alignment_score: float = 0.0
    
    def __post_init__(self):
        """Sort splice sites by position after initialization."""
        self.splice_sites.sort(key=lambda x: x.position)
    
    @property
    def exon_count(self) -> int:
        """Number of exons (splice sites + 1)."""
        return len(self.splice_sites) + 1
    
    @property
    def intron_count(self) -> int:
        """Number of introns."""
        return len(self.splice_sites)
    
    def get_exon_bounds(self) -> List[Tuple[int, int]]:
        """Get exon boundaries as (start, end) tuples."""
        if not self.splice_sites:
            return [(self.start, self.end)]
        
        exons = []
        current_start = self.start
        
        for site in self.splice_sites:
            if site.donor_acceptor == "donor":
                exons.append((current_start, site.position))
                current_start = site.position + 1
            elif site.donor_acceptor == "acceptor":
                current_start = site.position
        
        # Add final exon
        exons.append((current_start, self.end))
        return exons


@dataclass
class Isoform:
    """Represents an identified isoform."""
    isoform_id: str
    gene_id: str
    chromosome: str
    strand: str
    start: int
    end: int
    splice_sites: List[SpliceSite] = field(default_factory=list)
    read_signatures: List[ReadSignature] = field(default_factory=list)
    consensus_sequence: Optional[str] = None
    isoform_type: Optional[IsoformType] = None
    confidence_score: float = 0.0
    support_reads: int = 0
    coverage: float = 0.0
    annotations: Dict[str, Any] = field(default_factory=dict)
    
    def __post_init__(self):
        """Sort splice sites and update derived fields."""
        self.splice_sites.sort(key=lambda x: x.position)
        self.support_reads = len(self.read_signatures)
        if self.read_signatures:
            self.coverage = sum(rs.coverage for rs in self.read_signatures) / len(self.read_signatures)
    
    @property
    def exon_count(self) -> int:
        """Number of exons."""
        return len(self.splice_sites) + 1
    
    @property
    def intron_count(self) -> int:
        """Number of introns."""
        return len(self.splice_sites)
    
    def get_exon_bounds(self) -> List[Tuple[int, int]]:
        """Get exon boundaries as (start, end) tuples."""
        if not self.splice_sites:
            return [(self.start, self.end)]
        
        exons = []
        current_start = self.start
        
        for site in self.splice_sites:
            if site.donor_acceptor == "donor":
                exons.append((current_start, site.position))
                current_start = site.position + 1
            elif site.donor_acceptor == "acceptor":
                current_start = site.position
        
        # Add final exon
        exons.append((current_start, self.end))
        return exons
    
    def add_read_signature(self, read_sig: ReadSignature) -> None:
        """Add a supporting read signature."""
        self.read_signatures.append(read_sig)
        self.support_reads = len(self.read_signatures)
        if self.read_signatures:
            self.coverage = sum(rs.coverage for rs in self.read_signatures) / len(self.read_signatures)


@dataclass
class Cluster:
    """Represents a cluster of read signatures."""
    cluster_id: str
    read_signatures: List[ReadSignature]
    centroid: Optional[ReadSignature] = None
    consensus_isoform: Optional[Isoform] = None
    
    def __post_init__(self):
        """Update centroid after initialization."""
        if self.read_signatures and not self.centroid:
            self.centroid = self._compute_centroid()
    
    def _compute_centroid(self) -> ReadSignature:
        """Compute centroid read signature."""
        if not self.read_signatures:
            raise ValueError("Cannot compute centroid for empty cluster")
        
        if len(self.read_signatures) == 1:
            return self.read_signatures[0]
        
        # Use weighted median for positions
        starts = [rs.start for rs in self.read_signatures]
        ends = [rs.end for rs in self.read_signatures]
        
        centroid_start = self._weighted_median(starts)
        centroid_end = self._weighted_median(ends)
        
        # Create centroid with consensus splice sites
        consensus_sites = self._consensus_splice_sites()
        
        return ReadSignature(
            read_id=f"centroid_{self.cluster_id}",
            transcript_id=None,
            gene_id=self.read_signatures[0].gene_id,
            chromosome=self.read_signatures[0].chromosome,
            strand=self.read_signatures[0].strand,
            start=centroid_start,
            end=centroid_end,
            splice_sites=consensus_sites,
            coverage=sum(rs.coverage for rs in self.read_signatures) / len(self.read_signatures),
            quality_score=sum(rs.quality_score for rs in self.read_signatures) / len(self.read_signatures)
        )
    
    def _weighted_median(self, values: List[int]) -> int:
        """Compute weighted median of values."""
        if not values:
            return 0
        
        sorted_values = sorted(values)
        n = len(sorted_values)
        
        if n % 2 == 1:
            return sorted_values[n // 2]
        else:
            return (sorted_values[n // 2 - 1] + sorted_values[n // 2]) // 2
    
    def _consensus_splice_sites(self) -> List[SpliceSite]:
        """Compute consensus splice sites from cluster."""
        if not self.read_signatures:
            return []
        
        # Collect all splice sites
        all_sites = []
        for rs in self.read_signatures:
            all_sites.extend(rs.splice_sites)
        
        if not all_sites:
            return []
        
        # Group sites by proximity and compute consensus
        consensus_sites = []
        tolerance = 10  # Base pair tolerance for grouping
        
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
    
    def _consensus_site(self, sites: List[SpliceSite]) -> SpliceSite:
        """Compute consensus splice site from a group."""
        if not sites:
            raise ValueError("Cannot compute consensus from empty site group")
        
        if len(sites) == 1:
            return sites[0]
        
        # Use median position
        positions = [s.position for s in sites]
        median_pos = self._weighted_median(positions)
        
        # Use most common strand and type
        strands = [s.strand for s in sites]
        types = [s.donor_acceptor for s in sites]
        
        consensus_strand = max(set(strands), key=strands.count)
        consensus_type = max(set(types), key=types.count)
        
        # Average confidence
        avg_confidence = sum(s.confidence for s in sites) / len(sites)
        
        return SpliceSite(
            position=median_pos,
            strand=consensus_strand,
            donor_acceptor=consensus_type,
            confidence=avg_confidence
        )
