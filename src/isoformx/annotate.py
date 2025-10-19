"""SQANTI-style annotation for isoforms."""

from typing import List, Dict, Optional, Tuple, Set
import numpy as np
import gffutils
from dataclasses import dataclass
from .models import Isoform, IsoformType, ReadSignature
from .utils import compute_exon_overlap, merge_intervals


@dataclass
class ReferenceTranscript:
    """Reference transcript annotation."""
    transcript_id: str
    gene_id: str
    chromosome: str
    strand: str
    start: int
    end: int
    exons: List[Tuple[int, int]]
    cds_start: Optional[int] = None
    cds_end: Optional[int] = None
    utr5_start: Optional[int] = None
    utr5_end: Optional[int] = None
    utr3_start: Optional[int] = None
    utr3_end: Optional[int] = None


@dataclass
class AnnotationResult:
    """Result of isoform annotation."""
    isoform_type: IsoformType
    reference_transcript: Optional[ReferenceTranscript]
    overlap_ratio: float
    splice_site_matches: int
    total_splice_sites: int
    structural_category: str
    subcategory: str
    annotations: Dict[str, any]


class SQANTIAnnotator:
    """SQANTI-style annotation for isoforms."""
    
    def __init__(self, reference_transcripts: List[ReferenceTranscript]):
        self.reference_transcripts = reference_transcripts
        self._build_index()
    
    def _build_index(self):
        """Build index for efficient lookup."""
        self.transcript_index = {}
        self.gene_index = {}
        self.chromosome_index = {}
        
        for transcript in self.reference_transcripts:
            # Transcript ID index
            self.transcript_index[transcript.transcript_id] = transcript
            
            # Gene ID index
            if transcript.gene_id not in self.gene_index:
                self.gene_index[transcript.gene_id] = []
            self.gene_index[transcript.gene_id].append(transcript)
            
            # Chromosome index
            if transcript.chromosome not in self.chromosome_index:
                self.chromosome_index[transcript.chromosome] = []
            self.chromosome_index[transcript.chromosome].append(transcript)

    def _overlapping_gene_for_isoform(self, isoform: Isoform) -> Optional[str]:
        """Return the gene_id with the strongest locus overlap on the same chromosome/strand.
        We define gene locus as the union of all transcript spans for that gene.
        """
        candidates = self.chromosome_index.get(isoform.chromosome, [])
        if not candidates:
            return None
        # Aggregate per gene span (min start, max end) on same strand
        gene_spans: Dict[str, Tuple[int, int]] = {}
        for t in candidates:
            if t.strand != isoform.strand:
                continue
            start, end = gene_spans.get(t.gene_id, (t.start, t.end))
            start = min(start, t.start)
            end = max(end, t.end)
            gene_spans[t.gene_id] = (start, end)
        # Pick gene with maximum overlap length with isoform span
        best_gene = None
        best_overlap = 0
        for gene_id, (gstart, gend) in gene_spans.items():
            ovl = max(0, min(isoform.end, gend) - max(isoform.start, gstart))
            if ovl > best_overlap:
                best_overlap = ovl
                best_gene = gene_id
        return best_gene if best_overlap > 0 else None
    
    def annotate_isoform(self, isoform: Isoform) -> AnnotationResult:
        """
        Annotate an isoform using SQANTI-style classification.
        
        Args:
            isoform: Isoform to annotate
        
        Returns:
            Annotation result
        """
        # Find best matching reference transcript
        best_match = self._find_best_match(isoform)
        
        if best_match is None:
            return self._annotate_novel(isoform)
        
        # Classify based on match
        return self._classify_match(isoform, best_match)
    
    def _find_best_match(self, isoform: Isoform) -> Optional[ReferenceTranscript]:
        """Find best matching reference transcript."""
        # Get candidates from same chromosome
        candidates = self.chromosome_index.get(isoform.chromosome, [])
        
        if not candidates:
            return None
        
        best_match = None
        best_score = 0.0
        
        for candidate in candidates:
            # Check strand compatibility
            if candidate.strand != isoform.strand:
                continue
            
            # Compute overlap score
            score = self._compute_match_score(isoform, candidate)
            
            if score > best_score:
                best_score = score
                best_match = candidate
        
        return best_match if best_score > 0.5 else None
    
    def _compute_match_score(self, isoform: Isoform, reference: ReferenceTranscript) -> float:
        """Compute match score between isoform and reference."""
        # Position overlap
        isoform_exons = isoform.get_exon_bounds()
        reference_exons = reference.exons
        
        overlap_ratio = compute_exon_overlap(isoform_exons, reference_exons)
        
        # Splice site matches
        splice_matches = self._count_splice_matches(isoform, reference)
        total_splice_sites = max(len(isoform.splice_sites), len(reference.exons) - 1)
        
        splice_score = splice_matches / total_splice_sites if total_splice_sites > 0 else 1.0
        
        # Combined score
        return 0.7 * overlap_ratio + 0.3 * splice_score
    
    def _count_splice_matches(self, isoform: Isoform, reference: ReferenceTranscript) -> int:
        """Count matching splice sites."""
        isoform_sites = [s.position for s in isoform.splice_sites]
        reference_sites = []
        
        # Extract splice sites from reference exons
        for i in range(len(reference.exons) - 1):
            reference_sites.append(reference.exons[i][1])  # Donor
            reference_sites.append(reference.exons[i + 1][0])  # Acceptor
        
        # Count matches within tolerance
        tolerance = 10
        matches = 0
        
        for iso_site in isoform_sites:
            for ref_site in reference_sites:
                if abs(iso_site - ref_site) <= tolerance:
                    matches += 1
                    break
        
        return matches
    
    def _annotate_novel(self, isoform: Isoform) -> AnnotationResult:
        """Annotate novel isoform (no reference match)."""
        # Check if it's intergenic
        if self._is_intergenic(isoform):
            return AnnotationResult(
                isoform_type=IsoformType.INTERGENIC,
                reference_transcript=None,
                overlap_ratio=0.0,
                splice_site_matches=0,
                total_splice_sites=len(isoform.splice_sites),
                structural_category="intergenic",
                subcategory="intergenic",
                annotations={"novel": True, "intergenic": True}
            )
        
        # Check if it's antisense
        if self._is_antisense(isoform):
            return AnnotationResult(
                isoform_type=IsoformType.ANTISENSE,
                reference_transcript=None,
                overlap_ratio=0.0,
                splice_site_matches=0,
                total_splice_sites=len(isoform.splice_sites),
                structural_category="antisense",
                subcategory="antisense",
                annotations={"novel": True, "antisense": True}
            )
        
        # Novel in catalog (NIC) or novel not in catalog (NNC)
        overlapping_gene = self._overlapping_gene_for_isoform(isoform)
        if self._is_nic(isoform):
            return AnnotationResult(
                isoform_type=IsoformType.NIC,
                reference_transcript=None,
                overlap_ratio=0.0,
                splice_site_matches=0,
                total_splice_sites=len(isoform.splice_sites),
                structural_category="NIC",
                subcategory="NIC",
                annotations={
                    "novel": True,
                    "nic": True,
                    **({"reference_gene": overlapping_gene, "novel_for_gene": True} if overlapping_gene else {})
                }
            )
        else:
            # Novel transcript for an overlapped gene if any
            if overlapping_gene:
                return AnnotationResult(
                    isoform_type=IsoformType.NNC,
                    reference_transcript=None,
                    overlap_ratio=0.0,
                    splice_site_matches=0,
                    total_splice_sites=len(isoform.splice_sites),
                    structural_category="novel_gene",
                    subcategory="NNC",
                    annotations={"novel": True, "nnc": True, "reference_gene": overlapping_gene, "novel_for_gene": True}
                )
            return AnnotationResult(
                isoform_type=IsoformType.NNC,
                reference_transcript=None,
                overlap_ratio=0.0,
                splice_site_matches=0,
                total_splice_sites=len(isoform.splice_sites),
                structural_category="NNC",
                subcategory="NNC",
                annotations={"novel": True, "nnc": True}
            )
    
    def _is_intergenic(self, isoform: Isoform) -> bool:
        """Check if isoform is intergenic."""
        candidates = self.chromosome_index.get(isoform.chromosome, [])
        
        for candidate in candidates:
            # Check if isoform overlaps with any gene
            if (isoform.start <= candidate.end and isoform.end >= candidate.start):
                return False
        
        return True
    
    def _is_antisense(self, isoform: Isoform) -> bool:
        """Check if isoform is antisense."""
        candidates = self.chromosome_index.get(isoform.chromosome, [])
        
        for candidate in candidates:
            # Check if isoform overlaps with gene on opposite strand
            if (candidate.strand != isoform.strand and
                isoform.start <= candidate.end and isoform.end >= candidate.start):
                return True
        
        return False
    
    def _is_nic(self, isoform: Isoform) -> bool:
        """Check if isoform is novel in catalog (NIC)."""
        # NIC: novel combination of known splice sites
        candidates = self.chromosome_index.get(isoform.chromosome, [])
        
        for candidate in candidates:
            if candidate.strand != isoform.strand:
                continue
            
            # Check if all splice sites are known
            isoform_sites = [s.position for s in isoform.splice_sites]
            reference_sites = []
            
            for i in range(len(candidate.exons) - 1):
                reference_sites.append(candidate.exons[i][1])
                reference_sites.append(candidate.exons[i + 1][0])
            
            # Check if all isoform splice sites are known
            tolerance = 10
            known_sites = 0
            
            for iso_site in isoform_sites:
                for ref_site in reference_sites:
                    if abs(iso_site - ref_site) <= tolerance:
                        known_sites += 1
                        break
            
            if known_sites == len(isoform_sites):
                return True
        
        return False
    
    def _classify_match(self, isoform: Isoform, reference: ReferenceTranscript) -> AnnotationResult:
        """Classify isoform based on reference match."""
        # Compute overlap ratio
        isoform_exons = isoform.get_exon_bounds()
        reference_exons = reference.exons
        
        overlap_ratio = compute_exon_overlap(isoform_exons, reference_exons)
        
        # Count splice site matches
        splice_matches = self._count_splice_matches(isoform, reference)
        total_splice_sites = max(len(isoform.splice_sites), len(reference.exons) - 1)
        
        # Determine structural category
        if overlap_ratio >= 0.95 and splice_matches == total_splice_sites:
            # Full splice match
            if len(isoform.splice_sites) == 0:
                isoform_type = IsoformType.FSM_MONO
                category = "FSM"
                subcategory = "mono-exon"
            else:
                isoform_type = IsoformType.FSM
                category = "FSM"
                subcategory = "multi-exon"
        elif splice_matches == total_splice_sites:
            # Incomplete splice match
            if len(isoform.splice_sites) == 0:
                isoform_type = IsoformType.ISM_MONO
                category = "ISM"
                subcategory = "mono-exon"
            else:
                isoform_type = IsoformType.ISM
                category = "ISM"
                subcategory = "multi-exon"
        else:
            # Novel in catalog
            if len(isoform.splice_sites) == 0:
                isoform_type = IsoformType.NIC_MONO
                category = "NIC"
                subcategory = "mono-exon"
            else:
                isoform_type = IsoformType.NIC
                category = "NIC"
                subcategory = "multi-exon"
        
        return AnnotationResult(
            isoform_type=isoform_type,
            reference_transcript=reference,
            overlap_ratio=overlap_ratio,
            splice_site_matches=splice_matches,
            total_splice_sites=total_splice_sites,
            structural_category=category,
            subcategory=subcategory,
            annotations={
                "reference_transcript": reference.transcript_id,
                "reference_gene": reference.gene_id,
                "overlap_ratio": overlap_ratio,
                "splice_matches": splice_matches,
                "total_splice_sites": total_splice_sites
            }
        )


def load_reference_transcripts(gff_file: str) -> List[ReferenceTranscript]:
    """
    Load reference transcripts from GFF/GTF file.
    
    Args:
        gff_file: Path to GFF/GTF file
    
    Returns:
        List of reference transcripts
    """
    print(f"Loading reference transcripts from {gff_file}...")
    
    # Create database from GFF file
    db = gffutils.create_db(gff_file, ':memory:', force=True, keep_order=True)
    
    transcripts = []
    
    # Process each transcript
    for transcript in db.features_of_type('transcript'):
        # Get gene information
        gene = None
        try:
            gene = db[transcript.attributes.get('gene_id', [transcript.attributes.get('Parent', [''])[0]])[0]]
        except:
            # If no gene found, create a dummy gene
            gene_id = transcript.attributes.get('gene_id', [transcript.attributes.get('Parent', ['unknown'])[0]])[0]
        
        # Get exons for this transcript
        exons = []
        for exon in db.children(transcript, featuretype='exon', order_by='start'):
            exons.append((exon.start - 1, exon.end))  # Convert to 0-based coordinates
        
        # Sort exons by start position
        exons.sort(key=lambda x: x[0])
        
        # Extract attributes
        transcript_id = transcript.attributes.get('transcript_id', [transcript.id])[0]
        gene_id = transcript.attributes.get('gene_id', [transcript.attributes.get('Parent', ['unknown'])[0]])[0]
        
        # Create reference transcript
        ref_transcript = ReferenceTranscript(
            transcript_id=transcript_id,
            gene_id=gene_id,
            chromosome=transcript.chrom,
            strand=transcript.strand,
            start=transcript.start - 1,  # Convert to 0-based
            end=transcript.end,
            exons=exons
        )
        
        transcripts.append(ref_transcript)
    
    print(f"Loaded {len(transcripts)} reference transcripts")
    return transcripts


def annotate_isoforms(
    isoforms: List[Isoform],
    reference_transcripts: List[ReferenceTranscript]
) -> List[AnnotationResult]:
    """
    Annotate a list of isoforms.
    
    Args:
        isoforms: List of isoforms to annotate
        reference_transcripts: List of reference transcripts
    
    Returns:
        List of annotation results
    """
    annotator = SQANTIAnnotator(reference_transcripts)
    
    results = []
    for isoform in isoforms:
        result = annotator.annotate_isoform(isoform)
        results.append(result)
    
    return results


def annotate_isoforms_with_reference(
    isoforms: List[Isoform],
    reference_file: str
) -> List[AnnotationResult]:
    """
    Annotate isoforms using reference GFF/GTF file.
    
    Args:
        isoforms: List of isoforms to annotate
        reference_file: Path to reference GFF/GTF file
    
    Returns:
        List of annotation results
    """
    # Load reference transcripts
    reference_transcripts = load_reference_transcripts(reference_file)
    
    # Annotate isoforms
    return annotate_isoforms(isoforms, reference_transcripts)
