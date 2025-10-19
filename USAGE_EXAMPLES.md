# IsoformX Usage Examples and Documentation

## Table of Contents
1. [Environment Setup](#environment-setup)
2. [Basic Usage](#basic-usage)
3. [Command Reference](#command-reference)
4. [Real Data Examples](#real-data-examples)
5. [Output Formats](#output-formats)
6. [Troubleshooting](#troubleshooting)
7. [Advanced Configuration](#advanced-configuration)

## Environment Setup

### Why You Need to Activate the Environment

The `isoformx` command is only available when the conda environment is properly activated. Here's why:

1. **Package Installation**: The `isoformx` package is installed in the conda environment
2. **PATH Configuration**: The environment's `bin` directory contains the executable
3. **Dependencies**: All required libraries are installed in the environment

### How to Activate the Environment

```bash
# Activate the isoformx conda environment
conda activate isoformx

# Verify it's working
isoformx --help
which isoformx
```

### If Activation Doesn't Work

If you get "command not found", try:

```bash
# Check available environments
conda info --envs

# If isoformx environment doesn't exist, recreate it
conda env create -f environment.yml

# Then activate and reinstall
conda activate isoformx
pip install -e . --force-reinstall
```

## Basic Usage

### 1. Test with Mock Data

```bash
# Generate mock data and run pipeline
isoformx mock --n-reads 100 --n-genes 10 --output test_isoforms.gff3 --verbose

# With more relaxed parameters for testing
isoformx mock --n-reads 200 --n-genes 20 --output test_isoforms.gff3 --verbose
```

### 2. Process Real BAM Files

```bash
# Basic clustering from BAM file
isoformx cluster --input ../input_data/21_tagged_polyadenylation.bam --output isoforms.gff3 --verbose

# With custom parameters
isoformx cluster \
    --input ../input_data/21_tagged_polyadenylation.bam \
    --output isoforms.gff3 \
    --min-reads 5 \
    --max-distance 0.2 \
    --min-confidence 0.5 \
    --verbose
```

### 3. Filter Results

```bash
# Filter by confidence score
isoformx filter --input isoforms.gff3 --min-confidence 0.6 --output high_conf_isoforms.gff3

# Filter by support reads
isoformx filter --input isoforms.gff3 --min-reads 10 --output well_supported_isoforms.gff3
```

### 4. Annotate Isoforms

```bash
# Annotate against reference (if you have reference GFF)
isoformx annotate \
    --input isoforms.gff3 \
    --reference /path/to/reference.gff3 \
    --output annotated_isoforms.gff3
```

## Command Reference

### `isoformx mock`
Generate mock data and run the complete pipeline.

**Options:**
- `--n-reads INT`: Number of reads to generate (default: 100)
- `--n-genes INT`: Number of genes to generate (default: 10)
- `--output FILE`: Output GFF3 file (default: mock_isoforms.gff3)
- `--verbose`: Show detailed progress

**Example:**
```bash
isoformx mock --n-reads 500 --n-genes 25 --output mock_results.gff3 --verbose
```

### `isoformx cluster`
Cluster reads and identify isoforms from BAM files.

**Options:**
- `--input FILE`: Input BAM file (required)
- `--output FILE`: Output GFF3 file (default: isoforms.gff3)
- `--min-reads INT`: Minimum reads per cluster (default: 3)
- `--max-distance FLOAT`: Maximum distance for clustering (default: 0.3)
- `--min-confidence FLOAT`: Minimum confidence score (default: 0.3)
- `--min-length INT`: Minimum read length (default: 1000)
- `--verbose`: Show detailed progress

**Example:**
```bash
isoformx cluster \
    --input sample.bam \
    --output sample_isoforms.gff3 \
    --min-reads 5 \
    --max-distance 0.2 \
    --min-confidence 0.5 \
    --verbose
```

### `isoformx filter`
Filter isoforms by various criteria.

**Options:**
- `--input FILE`: Input GFF3 file (required)
- `--output FILE`: Output GFF3 file (default: filtered_isoforms.gff3)
- `--min-confidence FLOAT`: Minimum confidence score
- `--min-reads INT`: Minimum support reads
- `--min-coverage FLOAT`: Minimum coverage
- `--min-exons INT`: Minimum number of exons

**Example:**
```bash
isoformx filter \
    --input isoforms.gff3 \
    --output high_quality_isoforms.gff3 \
    --min-confidence 0.7 \
    --min-reads 10 \
    --min-coverage 0.8
```

### `isoformx annotate`
Annotate isoforms using reference annotation.

**Options:**
- `--input FILE`: Input GFF3 file (required)
- `--reference FILE`: Reference GFF3 file (required)
- `--output FILE`: Output GFF3 file (default: annotated_isoforms.gff3)
- `--overlap-threshold FLOAT`: Overlap threshold for annotation (default: 0.5)

**Example:**
```bash
isoformx annotate \
    --input isoforms.gff3 \
    --reference reference_annotation.gff3 \
    --output annotated_isoforms.gff3 \
    --overlap-threshold 0.7
```

## Real Data Examples

### Processing Your BAM Files

Based on your data structure, here are examples for processing your files:

```bash
# Process elongating reads
isoformx cluster \
    --input ../input_data/13_tagged_elongating.bam \
    --output sample13e_isoforms.gff3 \
    --min-reads 3 \
    --verbose

# Process polyadenylation reads
isoformx cluster \
    --input ../input_data/13_tagged_polyadenylation.bam \
    --output sample13p_isoforms.gff3 \
    --min-reads 3 \
    --verbose

# Process multiple samples
for sample in 13 14 15 19 20 21; do
    echo "Processing sample ${sample}..."
    
    # Elongating reads
    isoformx cluster \
        --input ../input_data/${sample}_tagged_elongating.bam \
        --output sample${sample}e_isoforms.gff3 \
        --min-reads 3 \
        --verbose
    
    # Polyadenylation reads
    isoformx cluster \
        --input ../input_data/${sample}_tagged_polyadenylation.bam \
        --output sample${sample}p_isoforms.gff3 \
        --min-reads 3 \
        --verbose
done
```

### Batch Processing Script

Create a script `process_all_samples.sh`:

```bash
#!/bin/bash

# Activate environment
conda activate isoformx

# Process all samples
samples=(13 14 15 19 20 21)
conditions=("elongating" "polyadenylation")

for sample in "${samples[@]}"; do
    for condition in "${conditions[@]}"; do
        echo "Processing sample ${sample} ${condition}..."
        
        isoformx cluster \
            --input ../input_data/${sample}_tagged_${condition}.bam \
            --output sample${sample}_${condition}_isoforms.gff3 \
            --min-reads 3 \
            --min-confidence 0.4 \
            --verbose
        
        # Filter high-confidence isoforms
        isoformx filter \
            --input sample${sample}_${condition}_isoforms.gff3 \
            --output sample${sample}_${condition}_high_conf.gff3 \
            --min-confidence 0.6 \
            --min-reads 5
        
        # Annotate isoforms
        isoformx annotate \
            --input sample${sample}_${condition}_isoforms.gff3 \
            --output sample${sample}_${condition}_annotated.gff3 \
            --reference ~/workshop/genome/Arabidopsis/Jin/Araport11_GTF_genes_transposons.Jan2023.gtf \
            --verbose
    done
done

echo "All samples processed!"
```

Make it executable and run:
```bash
chmod +x process_all_samples.sh
./process_all_samples.sh
```

## Annotation

### Basic Annotation

Annotate isoforms using reference annotation:

```bash
isoformx annotate \
    --input isoforms.gff3 \
    --output annotated_isoforms.gff3 \
    --reference reference.gtf \
    --verbose
```

### Structural Categories

The annotation classifies isoforms into SQANTI-style categories:

- **FSM**: Full Splice Match - matches reference transcript exactly
- **ISM**: Incomplete Splice Match - partial match with reference
- **NIC**: Novel In Catalog - novel combination of known splice sites
- **NNC**: Novel Not in Catalog - contains novel splice sites
- **antisense**: Antisense transcripts (opposite strand)
- **intergenic**: Intergenic transcripts (no gene overlap)

### Annotation Output

The annotated GFF includes additional attributes:

```
Chr1    IsoformX    transcript    11588   13164   .    -    .    ID=isoform_1;structural_category=NIC;reference_transcript=AT1G01030.1;overlap_ratio=0.85;splice_matches=3;total_splice_sites=4
```

### Batch Annotation

Annotate multiple files:

```bash
for file in *.gff3; do
    isoformx annotate \
        --input "$file" \
        --output "${file%.gff3}_annotated.gff3" \
        --reference reference.gtf \
        --verbose
done
```

### Analyze Annotation Results

Count structural categories:

```bash
# Count categories
grep "structural_category=" annotated_isoforms.gff3 | \
    cut -d';' -f7 | cut -d'=' -f2 | sort | uniq -c

# Extract specific categories
grep "structural_category=FSM" annotated_isoforms.gff3 > fsm_isoforms.gff3
grep "structural_category=NIC" annotated_isoforms.gff3 > nic_isoforms.gff3
```

## Proportion Analysis

### Basic Proportion Analysis

Analyze isoform proportions across multiple samples:

```bash
isoformx proportion \
    -b sample1.bam \
    -b sample2.bam \
    -b sample3.bam \
    -g isoforms.gtf \
    -o proportions.tsv \
    --verbose
```

### Filtered Analysis

Filter isoforms by minimum read count:

```bash
isoformx proportion \
    -b sample1.bam \
    -b sample2.bam \
    -g isoforms.gtf \
    -o proportions_filtered.tsv \
    --min-reads 10 \
    --verbose
```

### Output Format

The proportion analysis generates a TSV file with columns:

- **Basic info**: `isoform_id`, `gene_id`, `chromosome`, `strand`, `start`, `end`
- **Statistics**: `total_reads`, `samples_with_reads`, `max_reads_per_sample`, `min_reads_per_sample`, `mean_reads_per_sample`, `std_reads_per_sample`
- **Per-sample data**: `{sample}_reads`, `{sample}_proportion`

### Analyze Results

```bash
# Count isoforms by sample presence
awk -F'\t' '$8 == 2' proportions.tsv | wc -l  # Both samples
awk -F'\t' '$8 == 1' proportions.tsv | wc -l  # Single sample

# Top abundant isoforms
sort -k7 -nr proportions.tsv | head -10

# Sample-specific isoforms
awk -F'\t' '$13 > 0 && $15 == 0' proportions.tsv > sample1_specific.tsv
```

## Output Formats

### GFF3 Format

The pipeline outputs standard GFF3 format with IsoformX-specific attributes:

```
##gff-version 3
##source IsoformX
##date 2025-10-18
Chr3    IsoformX        transcript      6191558 6194306 .       +       .       ID=isoform_cluster_9000;Name=isoform_cluster_9000;gene_id=None;confidence=0.781;support_reads=18;coverage=0.926
Chr3    IsoformX        exon    6191558 6191825 .       +       .       ID=isoform_cluster_9000.exon1;Parent=isoform_cluster_9000
```

**Key Attributes:**
- `ID`: Unique isoform identifier
- `confidence`: Confidence score (0-1)
- `support_reads`: Number of supporting reads
- `coverage`: Coverage score (0-1)
- `gene_id`: Associated gene (if known)

### Statistics Output

The pipeline provides statistics during processing:

```
=== Testing IsoformX with sample.bam ===
Loading reads from BAM file...
Loaded 2847 reads
Read length range: 1000 - 15000
Chromosomes: {'Chr1', 'Chr2', 'Chr3', 'Chr4', 'Chr5'}
Strands: {'+', '-'}
Clustering with max_distance=0.3, min_reads=3
Found 1250 clusters
Building consensus with min_coverage=0.2
Built 1250 consensus isoforms
Retained 606 high-confidence isoforms
Saved 606 isoforms to isoforms.gff3
```

## Troubleshooting

### Common Issues

1. **"command not found: isoformx"**
   ```bash
   # Solution: Activate conda environment
   conda activate isoformx
   ```

2. **"No reads found"**
   ```bash
   # Check BAM file format and read length filter
   isoformx cluster --input sample.bam --min-length 500 --verbose
   ```

3. **"No isoforms found"**
   ```bash
   # Use more relaxed parameters
   isoformx cluster --input sample.bam --min-reads 2 --max-distance 0.5 --verbose
   ```

4. **Memory issues with large files**
   ```bash
   # Process smaller chunks or use more restrictive filters
   isoformx cluster --input sample.bam --min-reads 5 --min-confidence 0.5
   ```

### Debug Mode

For detailed debugging, you can run the pipeline with verbose output:

```bash
isoformx cluster --input sample.bam --output debug.gff3 --verbose
```

This will show:
- Number of reads loaded
- Clustering progress
- Consensus building details
- Filtering statistics

## Advanced Configuration

### Custom Parameters

You can adjust the pipeline behavior by modifying parameters:

```bash
# Very strict clustering (high precision)
isoformx cluster \
    --input sample.bam \
    --output strict_isoforms.gff3 \
    --min-reads 10 \
    --max-distance 0.1 \
    --min-confidence 0.8

# Very permissive clustering (high recall)
isoformx cluster \
    --input sample.bam \
    --output permissive_isoforms.gff3 \
    --min-reads 2 \
    --max-distance 0.5 \
    --min-confidence 0.2
```

### Quality Control

```bash
# Generate quality statistics
isoformx cluster --input sample.bam --output qc_isoforms.gff3 --verbose

# Filter by multiple criteria
isoformx filter \
    --input qc_isoforms.gff3 \
    --output high_quality.gff3 \
    --min-confidence 0.6 \
    --min-reads 5 \
    --min-coverage 0.7 \
    --min-exons 2
```

### Integration with Other Tools

```bash
# Convert to BED format for visualization
grep "transcript" isoforms.gff3 | \
    awk '{print $1"\t"$4"\t"$5"\t"$9}' | \
    sed 's/ID=//g' | \
    sed 's/;.*//g' > isoforms.bed

# Extract sequences (if you have genome FASTA)
bedtools getfasta -fi genome.fa -bed isoforms.bed -fo isoforms.fa
```

## Performance Tips

1. **Use appropriate read length filters** based on your sequencing technology
2. **Adjust clustering parameters** based on your data quality
3. **Process samples in parallel** for multiple files
4. **Use filtering** to reduce downstream analysis load
5. **Monitor memory usage** for large BAM files

## Next Steps

After running the pipeline:

1. **Visualize results** in genome browsers (IGV, JBrowse)
2. **Compare with reference annotation** to identify novel isoforms
3. **Perform functional analysis** on identified isoforms
4. **Validate high-confidence candidates** with experimental methods
