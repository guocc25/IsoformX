# IsoformX: Long-read Isoform Identification Pipeline

A pipeline for identifying novel isoforms from long-read sequencing data using error-aware clustering and consensus building.

## Features

- **Error-aware clustering**: Uses Union-Find with distance metrics that account for sequencing errors
- **Consensus building**: Generates high-confidence isoform sequences from clustered reads
- **SQANTI-style annotation**: Categorizes isoforms (FSM, ISM, NIC, NNC, etc.)
- **Confidence scoring**: Provides reliability metrics for each identified isoform
- **Mock data support**: Includes synthetic datasets for testing and development

## Installation

### Using Conda (Recommended)

```bash
conda env create -f environment.yml
conda activate isoformx
pip install -e .
```

### Environment Activation

**Important**: You must activate the conda environment before using IsoformX:

```bash
# Activate the environment
conda activate isoformx

# Verify it's working
isoformx --help
which isoformx
```

If you get "command not found: isoformx", the environment is not properly activated. Make sure to:
1. Check available environments: `conda info --envs`
2. Activate the correct environment: `conda activate isoformx`
3. Reinstall if needed: `pip install -e . --force-reinstall`

## Quick Start

### Process All Your BAM Files

```bash
# Make sure environment is activated
conda activate isoformx

# Run the quick start script
./quick_start.sh
```

This will process all BAM files in `../input_data/` and save results in `results/` directory.

### Manual Processing

```bash
# Test with mock data first
isoformx mock --n-reads 100 --output test.gff3 --verbose

# Process a single BAM file
isoformx cluster --input ../input_data/21_tagged_polyadenylation.bam --output isoforms.gff3 --verbose

# Filter high-confidence isoforms
isoformx filter --input isoforms.gff3 --min-confidence 0.6 --output high_conf.gff3
```

## Usage

### Basic clustering

```bash
isoformx cluster input.bam --reference ref.fa --output isoforms.gff3
```

### With custom parameters

```bash
isoformx cluster input.bam \
    --reference ref.fa \
    --output isoforms.gff3 \
    --min-reads 3 \
    --max-distance 0.1 \
    --min-coverage 0.8
```

### Annotation

```bash
isoformx annotate isoforms.gff3 \
    --reference ref.fa \
    --annotation ref.gtf \
    --output annotated.gff3
```

## Pipeline Overview

1. **Read Processing**: Extract read signatures (splice sites, TSS/TES)
2. **Clustering**: Group reads by error-aware distance metrics
3. **Consensus**: Build consensus sequences for each cluster
4. **Annotation**: Categorize isoforms relative to reference
5. **Scoring**: Assign confidence scores

## File Structure

```
isoformx/
├── src/isoformx/
│   ├── models.py      # Data classes (ReadSig, Isoform)
│   ├── utils.py       # Utility functions
│   ├── cluster.py     # Error-aware clustering
│   ├── consensus.py   # Consensus building
│   ├── annotate.py    # SQANTI-style annotation
│   ├── score.py       # Confidence scoring
│   ├── io.py          # I/O operations
│   └── cli.py         # Command-line interface
├── tests/             # Unit tests
├── environment.yml    # Conda environment
└── pyproject.toml     # Package configuration
```

## Development

Run tests:

```bash
pytest tests/
```

Run with mock data:

```bash
isoformx cluster --mock-data --output test_output.gff3
```
## How to use the pipeline

## Activate the environment
conda activate isoformx

## Run with mock data (for testing)
isoformx mock --n-reads 100 --n-genes 10 --output test.gff3

## Run with real BAM data
isoformx cluster --input ../input_data/21_tagged_polyadenylation.bam --output isoforms.gff3 --verbose

## Filter isoforms by confidence
isoformx filter --input isoforms.gff3 --min-confidence 0.6 --output high_conf_isoforms.gff3

## Annotate isoforms (if you have reference annotation)
isoformx annotate --input isoforms.gff3 --reference ref.gff3 --output annotated_isoforms.gff3
