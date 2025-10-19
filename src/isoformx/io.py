"""I/O operations for isoform identification."""

from typing import List, Dict, Optional, Tuple, Iterator
import pysam
import pandas as pd
import numpy as np
from pathlib import Path
from dataclasses import dataclass
from .models import ReadSignature, Isoform, SpliceSite
from .annotate import ReferenceTranscript


@dataclass
class MockDataParams:
    """Parameters for generating mock data."""
    n_reads: int = 1000
    n_genes: int = 50
    n_transcripts_per_gene: int = 3
    chromosome_length: int = 1000000
    read_length_mean: int = 2000
    read_length_std: int = 500
    error_rate: float = 0.05


class MockDataGenerator:
    """Generate mock long-read data for testing."""
    
    def __init__(self, params: MockDataParams):
        self.params = params
        np.random.seed(42)  # For reproducible results
    
    def generate_mock_reads(self) -> List[ReadSignature]:
        """Generate mock read signatures."""
        reads = []
        
        for i in range(self.params.n_reads):
            # Generate random gene and transcript
            gene_id = f"GENE_{i % self.params.n_genes:03d}"
            transcript_id = f"TX_{i % self.params.n_transcripts_per_gene:03d}"
            
            # Generate random position
            start = np.random.randint(0, self.params.chromosome_length - 1000)
            read_length = int(np.random.normal(self.params.read_length_mean, self.params.read_length_std))
            end = start + read_length
            
            # Generate splice sites
            splice_sites = self._generate_splice_sites(start, end)
            
            # Generate quality metrics
            coverage = np.random.uniform(0.3, 1.0)
            quality_score = np.random.uniform(0.2, 1.0)
            alignment_score = np.random.uniform(0.5, 1.0)
            
            read = ReadSignature(
                read_id=f"read_{i:06d}",
                transcript_id=transcript_id,
                gene_id=gene_id,
                chromosome="chr1",
                strand="+" if np.random.random() > 0.5 else "-",
                start=start,
                end=end,
                splice_sites=splice_sites,
                coverage=coverage,
                quality_score=quality_score,
                read_length=read_length,
                alignment_score=alignment_score
            )
            
            reads.append(read)
        
        return reads
    
    def _generate_splice_sites(self, start: int, end: int) -> List[SpliceSite]:
        """Generate random splice sites."""
        splice_sites = []
        
        # Random number of splice sites (0-5)
        n_sites = np.random.randint(0, 6)
        
        if n_sites == 0:
            return splice_sites
        
        # Generate positions
        positions = np.random.randint(start + 100, end - 100, size=n_sites)
        positions = np.sort(positions)
        
        for i, pos in enumerate(positions):
            # Alternate between donor and acceptor
            donor_acceptor = "donor" if i % 2 == 0 else "acceptor"
            
            site = SpliceSite(
                position=int(pos),
                strand="+",
                donor_acceptor=donor_acceptor,
                confidence=np.random.uniform(0.7, 1.0)
            )
            
            splice_sites.append(site)
        
        return splice_sites
    
    def generate_mock_reference(self) -> List[ReferenceTranscript]:
        """Generate mock reference transcripts."""
        transcripts = []
        
        for gene_idx in range(self.params.n_genes):
            gene_id = f"GENE_{gene_idx:03d}"
            
            for tx_idx in range(self.params.n_transcripts_per_gene):
                transcript_id = f"TX_{tx_idx:03d}"
                
                # Generate random position
                start = np.random.randint(0, self.params.chromosome_length - 5000)
                end = start + np.random.randint(1000, 5000)
                
                # Generate exons
                exons = self._generate_exons(start, end)
                
                transcript = ReferenceTranscript(
                    transcript_id=transcript_id,
                    gene_id=gene_id,
                    chromosome="chr1",
                    strand="+",
                    start=start,
                    end=end,
                    exons=exons
                )
                
                transcripts.append(transcript)
        
        return transcripts
    
    def _generate_exons(self, start: int, end: int) -> List[Tuple[int, int]]:
        """Generate random exons."""
        exons = []
        
        # Random number of exons (1-5)
        n_exons = np.random.randint(1, 6)
        
        if n_exons == 1:
            return [(start, end)]
        
        # Generate exon boundaries
        boundaries = np.random.randint(start + 100, end - 100, size=n_exons - 1)
        boundaries = np.sort(boundaries)
        
        # Create exons
        current_start = start
        for boundary in boundaries:
            exons.append((current_start, boundary))
            current_start = boundary + 1
        
        exons.append((current_start, end))
        
        return exons


class BAMReader:
    """Read long-read data from BAM files."""
    
    def __init__(self, bam_path: str):
        self.bam_path = bam_path
        self.bamfile = pysam.AlignmentFile(bam_path, "rb")
    
    def read_signatures(self, min_length: int = 1000) -> List[ReadSignature]:
        """
        Extract read signatures from BAM file.
        
        Args:
            min_length: Minimum read length to include
        
        Returns:
            List of read signatures
        """
        signatures = []
        
        for read in self.bamfile:
            if read.is_unmapped or read.is_secondary or read.is_supplementary:
                continue
            
            if read.query_length < min_length:
                continue
            
            # Extract basic information
            read_id = read.query_name
            chromosome = self.bamfile.get_reference_name(read.reference_id)
            strand = "-" if read.is_reverse else "+"
            start = read.reference_start
            end = read.reference_end
            
            # Extract splice sites from CIGAR
            splice_sites = self._extract_splice_sites(read)
            
            # Extract quality metrics
            coverage = self._compute_coverage(read)
            quality_score = self._compute_quality_score(read)
            alignment_score = self._compute_alignment_score(read)
            
            signature = ReadSignature(
                read_id=read_id,
                transcript_id=None,
                gene_id=None,
                chromosome=chromosome,
                strand=strand,
                start=start,
                end=end,
                splice_sites=splice_sites,
                coverage=coverage,
                quality_score=quality_score,
                read_length=read.query_length,
                alignment_score=alignment_score
            )
            
            signatures.append(signature)
        
        return signatures
    
    def _extract_splice_sites(self, read) -> List[SpliceSite]:
        """Extract splice sites from CIGAR string."""
        splice_sites = []
        
        if not read.cigartuples:
            return splice_sites
        
        current_pos = read.reference_start
        
        for op, length in read.cigartuples:
            if op == 3:  # N (splice junction)
                # Donor site
                strand = "-" if read.is_reverse else "+"
                donor_site = SpliceSite(
                    position=current_pos - 1,
                    strand=strand,
                    donor_acceptor="donor",
                    confidence=1.0
                )
                splice_sites.append(donor_site)
                
                # Acceptor site
                acceptor_site = SpliceSite(
                    position=current_pos + length,
                    strand=strand,
                    donor_acceptor="acceptor",
                    confidence=1.0
                )
                splice_sites.append(acceptor_site)
            
            if op in [0, 2, 3, 7, 8]:  # M, D, N, =, X
                current_pos += length
        
        return splice_sites
    
    def _compute_coverage(self, read) -> float:
        """Compute coverage for read."""
        # Simple coverage based on read length and alignment
        if read.query_length == 0:
            return 0.0
        
        aligned_length = sum(length for op, length in read.cigartuples if op == 0)
        return aligned_length / read.query_length
    
    def _compute_quality_score(self, read) -> float:
        """Compute quality score for read."""
        if not read.query_qualities:
            return 0.0
        
        # Convert Phred scores to probabilities
        qual_scores = np.array(read.query_qualities)
        error_probs = 10 ** (-qual_scores / 10)
        avg_error_prob = np.mean(error_probs)
        
        # Convert back to quality score
        return 1.0 - avg_error_prob
    
    def _compute_alignment_score(self, read) -> float:
        """Compute alignment score for read."""
        if read.mapping_quality == 255:  # Unknown mapping quality
            return 0.5
        
        return read.mapping_quality / 60.0  # Normalize to 0-1
    
    def close(self):
        """Close BAM file."""
        self.bamfile.close()


class GFFWriter:
    """Write isoforms to GFF format."""
    
    def __init__(self, output_path: str):
        self.output_path = output_path
    
    def write_isoforms(self, isoforms: List[Isoform], annotations: List[Dict] = None):
        """
        Write isoforms to GFF file, sorted by chromosome and position.
        
        Args:
            isoforms: List of isoforms to write
            annotations: Optional list of annotation dictionaries
        """
        import re
        
        def natural_sort_key(text):
            """Natural sorting for chromosome names (works for any species)"""
            parts = re.split(r'(\d+)', text)
            return [int(part) if part.isdigit() else part.lower() for part in parts]
        
        # Group isoforms by chromosome and sort chromosomes
        chromosome_groups = {}
        for isoform in isoforms:
            chrom = isoform.chromosome
            if chrom not in chromosome_groups:
                chromosome_groups[chrom] = []
            chromosome_groups[chrom].append(isoform)
        
        # Sort chromosomes naturally
        sorted_chromosomes = sorted(chromosome_groups.keys(), key=natural_sort_key)
        
        # Collect all GFF lines in proper order
        gff_lines = []
        
        for chrom in sorted_chromosomes:
            # Sort isoforms within chromosome by start position
            chrom_isoforms = sorted(chromosome_groups[chrom], key=lambda x: x.start)
            
            for i, isoform in enumerate(chrom_isoforms):
                # Create transcript line
                attributes = [
                    f"ID={isoform.isoform_id}",
                    f"Name={isoform.isoform_id}",
                    f"gene_id={isoform.gene_id}",
                    f"confidence={isoform.confidence_score:.3f}",
                    f"support_reads={isoform.support_reads}",
                    f"coverage={isoform.coverage:.3f}"
                ]
                
                if annotations and i < len(annotations):
                    ann = annotations[i]
                    if "structural_category" in ann:
                        attributes.append(f"structural_category={ann['structural_category']}")
                    if "reference_transcript" in ann:
                        attributes.append(f"reference_transcript={ann['reference_transcript']}")
                
                transcript_line = f"{isoform.chromosome}\tIsoformX\ttranscript\t{isoform.start}\t{isoform.end}\t.\t{isoform.strand}\t.\t{';'.join(attributes)}"
                gff_lines.append(transcript_line)
                
                # Create exon lines immediately after transcript
                exons = isoform.get_exon_bounds()
                for j, (exon_start, exon_end) in enumerate(exons):
                    exon_id = f"{isoform.isoform_id}.exon{j+1}"
                    exon_attributes = [
                        f"ID={exon_id}",
                        f"Parent={isoform.isoform_id}"
                    ]
                    
                    exon_line = f"{isoform.chromosome}\tIsoformX\texon\t{exon_start}\t{exon_end}\t.\t{isoform.strand}\t.\t{';'.join(exon_attributes)}"
                    gff_lines.append(exon_line)
        
        # Write to file
        with open(self.output_path, 'w') as f:
            # Write header
            f.write("##gff-version 3\n")
            f.write(f"##source IsoformX\n")
            f.write(f"##date {pd.Timestamp.now().strftime('%Y-%m-%d')}\n")
            
            # Write sorted lines
            for line in gff_lines:
                f.write(line + "\n")


def load_mock_data(params: Optional[MockDataParams] = None) -> Tuple[List[ReadSignature], List[ReferenceTranscript]]:
    """
    Load mock data for testing.
    
    Args:
        params: Mock data parameters (uses defaults if None)
    
    Returns:
        Tuple of (read_signatures, reference_transcripts)
    """
    if params is None:
        params = MockDataParams()
    
    generator = MockDataGenerator(params)
    
    reads = generator.generate_mock_reads()
    reference = generator.generate_mock_reference()
    
    return reads, reference


def load_bam_data(bam_path: str, min_length: int = 1000) -> List[ReadSignature]:
    """
    Load read signatures from BAM file.
    
    Args:
        bam_path: Path to BAM file
        min_length: Minimum read length to include
    
    Returns:
        List of read signatures
    """
    reader = BAMReader(bam_path)
    try:
        signatures = reader.read_signatures(min_length)
        return signatures
    finally:
        reader.close()


def save_isoforms(isoforms: List[Isoform], output_path: str, annotations: List[Dict] = None):
    """
    Save isoforms to GFF file, sorted by chromosome and position.
    
    Args:
        isoforms: List of isoforms to save
        output_path: Output file path
        annotations: Optional list of annotation dictionaries
    """
    writer = GFFWriter(output_path)
    writer.write_isoforms(isoforms, annotations)


def load_isoforms_from_gff(gff_file: str) -> List[Isoform]:
    """
    Load isoforms from GFF file.
    
    Args:
        gff_file: Path to GFF file
    
    Returns:
        List of isoforms
    """
    isoforms = []
    current_isoform = None
    
    with open(gff_file, 'r') as f:
        for line in f:
            line = line.strip()
            
            # Skip comments and empty lines
            if line.startswith('#') or not line:
                continue
            
            fields = line.split('\t')
            if len(fields) < 9:
                continue
            
            feature_type = fields[2]
            
            if feature_type == 'transcript':
                # Parse transcript line
                chromosome = fields[0]
                start = int(fields[3]) - 1  # Convert to 0-based
                end = int(fields[4])
                strand = fields[6]
                
                # Parse attributes
                attributes = {}
                for attr in fields[8].split(';'):
                    if '=' in attr:
                        key, value = attr.split('=', 1)
                        attributes[key.strip()] = value.strip()
                
                isoform_id = attributes.get('ID', f'isoform_{len(isoforms)}')
                gene_id = attributes.get('gene_id', 'unknown')
                confidence = float(attributes.get('confidence', 0.5))
                support_reads = int(attributes.get('support_reads', 1))
                coverage = float(attributes.get('coverage', 0.5))
                
                # Create isoform
                current_isoform = Isoform(
                    isoform_id=isoform_id,
                    gene_id=gene_id,
                    chromosome=chromosome,
                    strand=strand,
                    start=start,
                    end=end,
                    confidence_score=confidence,
                    support_reads=support_reads,
                    coverage=coverage,
                    splice_sites=[]  # Will be populated from exons
                )
                
                isoforms.append(current_isoform)
            
            elif feature_type == 'exon' and current_isoform:
                # Parse exon line
                exon_start = int(fields[3]) - 1  # Convert to 0-based
                exon_end = int(fields[4])
                
                # Add exon to current isoform
                # Note: This is a simplified approach - in practice you'd want to
                # collect all exons first, then create splice sites
                pass
    
    return isoforms
