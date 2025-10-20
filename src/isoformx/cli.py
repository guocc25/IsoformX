"""Command-line interface for IsoformX."""

import click
import sys
from pathlib import Path
from typing import Optional
import logging
from datetime import datetime
from .cluster import cluster_read_signatures, ClusterParams
from .consensus import build_consensus_isoforms, ConsensusParams
from .annotate import annotate_isoforms, ReferenceTranscript
from .score import score_isoforms, filter_by_confidence, ScoringParams
from .io import load_mock_data, load_bam_data, save_isoforms, MockDataParams


@click.group()
@click.version_option(version="0.1.0")
def main():
    """IsoformX: Long-read isoform identification pipeline."""
    # Configure root logger once
    logging.basicConfig(
        level=logging.INFO,
        format='[%(asctime)s] %(levelname)s %(message)s',
        datefmt='%Y-%m-%d %H:%M:%S'
    )
    logging.info("IsoformX CLI initialized")


@main.command()
@click.option("--input", "-i", required=True, help="Input BAM file")
@click.option("--output", "-o", required=True, help="Output GFF file")
@click.option("--reference", "-r", help="Reference genome FASTA file")
@click.option("--annotation", "-a", help="Reference annotation GTF file")
@click.option("--min-reads", default=2, help="Minimum reads per cluster")
@click.option("--max-distance", default=0.2, help="Maximum distance for clustering")
@click.option("--min-coverage", default=0.5, help="Minimum coverage for consensus")
@click.option("--min-confidence", default=0.5, help="Minimum confidence score")
@click.option("--mock-data", is_flag=True, help="Use mock data for testing")
@click.option("--verbose", "-v", is_flag=True, help="Verbose output")
def cluster(input: str, output: str, reference: Optional[str], annotation: Optional[str],
           min_reads: int, max_distance: float, min_coverage: float, min_confidence: float,
           mock_data: bool, verbose: bool):
    """Cluster reads and identify isoforms."""
    
    start_time = datetime.now()
    logging.info(f"cluster started | input={input} output={output}")
    if verbose:
        click.echo("Starting IsoformX clustering...")
    
    try:
        # Load data
        if mock_data:
            if verbose:
                click.echo("Loading mock data...")
            read_signatures, reference_transcripts = load_mock_data()
        else:
            if not Path(input).exists():
                click.echo(f"Error: Input file {input} not found", err=True)
                sys.exit(1)
            
            if verbose:
                click.echo(f"Loading reads from {input}...")
            read_signatures = load_bam_data(input)
        
        if verbose:
            click.echo(f"Loaded {len(read_signatures)} reads")
        
        # Set up parameters
        cluster_params = ClusterParams(
            min_reads=min_reads,
            max_distance=max_distance
        )
        
        consensus_params = ConsensusParams(
            min_coverage=min_coverage,
            min_support_reads=min_reads
        )
        
        scoring_params = ScoringParams(
            min_reads=min_reads,
            min_coverage=min_coverage
        )
        
        # Cluster reads
        if verbose:
            click.echo("Clustering reads...")
        clusters = cluster_read_signatures(read_signatures, cluster_params)
        
        if verbose:
            click.echo(f"Found {len(clusters)} clusters")
        
        # Build consensus isoforms
        if verbose:
            click.echo("Building consensus isoforms...")
        isoforms = build_consensus_isoforms(clusters, consensus_params)
        
        if verbose:
            click.echo(f"Built {len(isoforms)} consensus isoforms")
        
        # Score isoforms
        if verbose:
            click.echo("Scoring isoforms...")
        isoforms = filter_by_confidence(isoforms, min_confidence, scoring_params)
        
        if verbose:
            click.echo(f"Retained {len(isoforms)} high-confidence isoforms")
        
        # Annotate isoforms (if reference provided)
        annotations = None
        # Note: Reference annotation not implemented in cluster command yet
        # This would require loading reference transcripts from a GFF file
        
        # Save results
        if verbose:
            click.echo(f"Saving results to {output}...")
        save_isoforms(isoforms, output, annotations)
        
        if verbose:
            click.echo("Done!")
        
        click.echo(f"Successfully identified {len(isoforms)} isoforms")
        elapsed = (datetime.now() - start_time).total_seconds()
        logging.info(f"cluster completed in {elapsed:.2f}s | output={output}")
        
    except Exception as e:
        click.echo(f"Error: {e}", err=True)
        if verbose:
            import traceback
            traceback.print_exc()
        logging.exception("cluster failed")
        sys.exit(1)


@main.command()
@click.option("--input", "-i", required=True, help="Input GFF file")
@click.option("--output", "-o", required=True, help="Output GFF file")
@click.option("--reference", "-r", required=True, help="Reference annotation GTF file")
@click.option("--verbose", "-v", is_flag=True, help="Verbose output")
def annotate(input: str, output: str, reference: str, verbose: bool):
    """Annotate isoforms using reference annotation."""
    
    start_time = datetime.now()
    logging.info(f"annotate started | input={input} reference={reference} output={output}")
    if verbose:
        click.echo("Starting annotation...")
    
    try:
        # Import annotation functions
        from .annotate import annotate_isoforms_with_reference
        from .io import load_isoforms_from_gff, save_isoforms
        
        # Load isoforms from GFF file
        if verbose:
            click.echo(f"Loading isoforms from {input}...")
        
        isoforms = load_isoforms_from_gff(input)
        
        if verbose:
            click.echo(f"Loaded {len(isoforms)} isoforms")
        
        # Annotate isoforms
        if verbose:
            click.echo(f"Annotating isoforms using reference {reference}...")
        
        annotation_results = annotate_isoforms_with_reference(isoforms, reference, verbose=verbose)
        
        # Convert to annotation dictionaries
        annotations = []
        for result in annotation_results:
            ann = {
                "structural_category": result.structural_category,
                "subcategory": result.subcategory,
                "overlap_ratio": result.overlap_ratio,
                "splice_matches": result.splice_site_matches,
                "total_splice_sites": result.total_splice_sites
            }
            if result.reference_transcript:
                ann["reference_transcript"] = result.reference_transcript.transcript_id
                ann["reference_gene"] = result.reference_transcript.gene_id
            annotations.append(ann)
        
        # Save annotated isoforms
        if verbose:
            click.echo(f"Saving annotated isoforms to {output}...")
        
        save_isoforms(isoforms, output, annotations)
        
        # Print summary
        category_counts = {}
        for result in annotation_results:
            category = result.structural_category
            category_counts[category] = category_counts.get(category, 0) + 1
        
        click.echo("Annotation completed!")
        click.echo("Structural categories:")
        for category, count in sorted(category_counts.items()):
            click.echo(f"  {category}: {count}")
        elapsed = (datetime.now() - start_time).total_seconds()
        logging.info(f"annotate completed in {elapsed:.2f}s | output={output}")
        
    except Exception as e:
        click.echo(f"Error: {e}", err=True)
        if verbose:
            import traceback
            traceback.print_exc()
        logging.exception("annotate failed")
        sys.exit(1)


@main.command()
@click.option("--bam-files", "-b", required=True, multiple=True, help="Input BAM files (can specify multiple)")
@click.option("--isoform-gtf", "-g", required=True, help="Isoform GTF file")
@click.option("--output", "-o", required=True, help="Output summary file")
@click.option("--min-reads", default=1, help="Minimum reads per isoform")
@click.option("--verbose", "-v", is_flag=True, help="Verbose output")
def proportion(bam_files, isoform_gtf, output, min_reads, verbose):
    """Summarize isoform proportions across samples."""
    
    start_time = datetime.now()
    logging.info(f"proportion started | n_bams={len(bam_files)} isoform_gtf={isoform_gtf} output={output}")
    if verbose:
        click.echo("Starting isoform proportion analysis...")
        click.echo(f"BAM files: {len(bam_files)}")
        click.echo(f"Isoform GTF: {isoform_gtf}")
    
    try:
        from .proportion import analyze_isoform_proportions
        
        # Analyze proportions
        results = analyze_isoform_proportions(
            bam_files=list(bam_files),
            isoform_gtf=isoform_gtf,
            min_reads=min_reads,
            verbose=verbose
        )
        
        # Save results
        if verbose:
            click.echo(f"Saving results to {output}...")
        
        results.to_csv(output, sep='\t', index=False)
        
        # Print summary
        click.echo("Analysis completed!")
        click.echo(f"Total isoforms analyzed: {len(results)}")
        click.echo(f"Results saved to: {output}")
        elapsed = (datetime.now() - start_time).total_seconds()
        logging.info(f"proportion completed in {elapsed:.2f}s | output={output}")
        
        if verbose:
            click.echo("\nTop 10 most abundant isoforms:")
            top_isoforms = results.nlargest(10, 'total_reads')
            for _, row in top_isoforms.iterrows():
                click.echo(f"  {row['isoform_id']}: {row['total_reads']} reads across {row['samples_with_reads']} samples")
        
    except Exception as e:
        click.echo(f"Error: {e}", err=True)
        if verbose:
            import traceback
            traceback.print_exc()
        logging.exception("proportion failed")
        sys.exit(1)


@main.command()
@click.option("--input", "-i", required=True, help="Input GFF file")
@click.option("--output", "-o", required=True, help="Output GFF file")
@click.option("--min-confidence", default=0.5, help="Minimum confidence score")
@click.option("--verbose", "-v", is_flag=True, help="Verbose output")
def filter(input: str, output: str, min_confidence: float, verbose: bool):
    """Filter isoforms by confidence score."""
    
    start_time = datetime.now()
    logging.info(f"filter started | input={input} output={output} min_confidence={min_confidence}")
    if verbose:
        click.echo("Starting filtering...")
    
    try:
        # TODO: Implement GFF reading and filtering
        click.echo("Filtering not yet implemented")
        elapsed = (datetime.now() - start_time).total_seconds()
        logging.info(f"filter completed in {elapsed:.2f}s | output={output}")
        
    except Exception as e:
        click.echo(f"Error: {e}", err=True)
        if verbose:
            import traceback
            traceback.print_exc()
        logging.exception("filter failed")
        sys.exit(1)


@main.command()
@click.option("--n-reads", default=1000, help="Number of mock reads")
@click.option("--n-genes", default=50, help="Number of mock genes")
@click.option("--output", "-o", required=True, help="Output GFF file")
@click.option("--verbose", "-v", is_flag=True, help="Verbose output")
def mock(n_reads: int, n_genes: int, output: str, verbose: bool):
    """Generate mock data and run pipeline."""
    
    start_time = datetime.now()
    logging.info(f"mock started | n_reads={n_reads} n_genes={n_genes} output={output}")
    if verbose:
        click.echo("Generating mock data...")
    
    try:
        # Generate mock data
        params = MockDataParams(n_reads=n_reads, n_genes=n_genes)
        read_signatures, reference_transcripts = load_mock_data(params)
        
        if verbose:
            click.echo(f"Generated {len(read_signatures)} reads and {len(reference_transcripts)} reference transcripts")
        
        # Run pipeline
        if verbose:
            click.echo("Running pipeline...")
        
        # Cluster reads with relaxed parameters
        cluster_params = ClusterParams(
            max_distance=0.8,  # Very loose
            min_reads=2,
            splice_site_tolerance=50,
            position_tolerance=100
        )
        clusters = cluster_read_signatures(read_signatures, cluster_params)
        
        # Build consensus isoforms with relaxed parameters
        consensus_params = ConsensusParams(
            min_coverage=0.1,  # Very low threshold
            min_support_reads=2,
            splice_site_threshold=0.3,  # Low threshold
            quality_threshold=0.1
        )
        isoforms = build_consensus_isoforms(clusters, consensus_params)
        
        # Filter by confidence with low threshold
        scoring_params = ScoringParams(
            min_reads=2,
            min_coverage=0.1,
            min_quality=0.1
        )
        isoforms = filter_by_confidence(isoforms, min_confidence=0.1, params=scoring_params)
        
        # Annotate isoforms
        annotation_results = annotate_isoforms(isoforms, reference_transcripts)
        
        # Convert to annotation dictionaries
        annotations = []
        for result in annotation_results:
            ann = {
                "structural_category": result.structural_category,
                "subcategory": result.subcategory,
                "overlap_ratio": result.overlap_ratio,
                "splice_matches": result.splice_site_matches,
                "total_splice_sites": result.total_splice_sites
            }
            if result.reference_transcript:
                ann["reference_transcript"] = result.reference_transcript.transcript_id
                ann["reference_gene"] = result.reference_transcript.gene_id
            annotations.append(ann)
        
        # Save results
        save_isoforms(isoforms, output, annotations)
        
        if verbose:
            click.echo("Done!")
        
        click.echo(f"Successfully identified {len(isoforms)} isoforms")
        elapsed = (datetime.now() - start_time).total_seconds()
        logging.info(f"mock completed in {elapsed:.2f}s | output={output}")
        
    except Exception as e:
        click.echo(f"Error: {e}", err=True)
        if verbose:
            import traceback
            traceback.print_exc()
        logging.exception("mock failed")
        sys.exit(1)


if __name__ == "__main__":
    main()
