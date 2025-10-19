#!/bin/bash
# Example: Isoform Proportion Analysis
# This script demonstrates how to analyze isoform proportions across multiple samples

echo "=== IsoformX Proportion Analysis Example ==="
echo

# Activate environment
conda activate isoformx

# Set paths
ISOFORM_GTF="test_output/elongating_annotated.gff3"
BAM_FILES=(
    "example_data/13_elongating_6k.bam"
    "example_data/13_polyadenylation_6k.bam"
)

echo "Isoform GTF: $ISOFORM_GTF"
echo "BAM files:"
for bam in "${BAM_FILES[@]}"; do
    echo "  - $bam"
done
echo

# Check if files exist
missing_files=0
for bam in "${BAM_FILES[@]}"; do
    if [ ! -f "$bam" ]; then
        echo "❌ Missing BAM file: $bam"
        missing_files=1
    fi
done

if [ ! -f "$ISOFORM_GTF" ]; then
    echo "❌ Missing isoform GTF: $ISOFORM_GTF"
    missing_files=1
fi

if [ $missing_files -eq 1 ]; then
    echo "Please run the annotation example first to generate required files."
    exit 1
fi

echo "✅ All required files found"
echo

# Step 1: Basic proportion analysis
echo "=== Step 1: Basic Proportion Analysis ==="
isoformx proportion \
    -b "${BAM_FILES[0]}" \
    -b "${BAM_FILES[1]}" \
    -g "$ISOFORM_GTF" \
    -o proportion_basic.tsv \
    --verbose

if [ $? -eq 0 ]; then
    echo "✅ Basic analysis completed"
    
    # Show summary
    total_isoforms=$(tail -n +2 proportion_basic.tsv | wc -l)
    echo "   Total isoforms analyzed: $total_isoforms"
else
    echo "❌ Basic analysis failed"
    exit 1
fi

echo

# Step 2: Filtered analysis (minimum 10 reads)
echo "=== Step 2: Filtered Analysis (min 10 reads) ==="
isoformx proportion \
    -b "${BAM_FILES[0]}" \
    -b "${BAM_FILES[1]}" \
    -g "$ISOFORM_GTF" \
    -o proportion_filtered.tsv \
    --min-reads 10 \
    --verbose

if [ $? -eq 0 ]; then
    echo "✅ Filtered analysis completed"
    
    # Show summary
    filtered_isoforms=$(tail -n +2 proportion_filtered.tsv | wc -l)
    echo "   Filtered isoforms (≥10 reads): $filtered_isoforms"
else
    echo "❌ Filtered analysis failed"
    exit 1
fi

echo

# Step 3: Analyze results
echo "=== Step 3: Analyze Results ==="

echo "Top 5 most abundant isoforms:"
head -6 proportion_basic.tsv | tail -5 | while IFS=$'\t' read -r isoform_id gene_id chromosome strand start end total_reads samples_with_reads max_reads min_reads mean_reads std_reads elongating_reads elongating_prop polya_reads polya_prop; do
    echo "  $isoform_id: $total_reads total reads ($samples_with_reads samples)"
done

echo
echo "Isoforms with reads in both samples:"
both_samples=$(tail -n +2 proportion_basic.tsv | awk -F'\t' '$8 == 2' | wc -l)
echo "  Found $both_samples isoforms present in both samples"

echo
echo "Isoforms specific to elongating sample:"
elongating_only=$(tail -n +2 proportion_basic.tsv | awk -F'\t' '$13 > 0 && $15 == 0' | wc -l)
echo "  Found $elongating_only elongating-specific isoforms"

echo
echo "Isoforms specific to polyadenylation sample:"
polya_only=$(tail -n +2 proportion_basic.tsv | awk -F'\t' '$13 == 0 && $15 > 0' | wc -l)
echo "  Found $polya_only polyadenylation-specific isoforms"

echo

# Step 4: Create summary statistics
echo "=== Step 4: Summary Statistics ==="

# Calculate total reads per sample
total_elongating=$(tail -n +2 proportion_basic.tsv | awk -F'\t' '{sum += $13} END {print sum}')
total_polya=$(tail -n +2 proportion_basic.tsv | awk -F'\t' '{sum += $15} END {print sum}')

echo "Total reads per sample:"
echo "  Elongating: $total_elongating reads"
echo "  Polyadenylation: $total_polya reads"

# Calculate mean proportions
mean_elongating=$(tail -n +2 proportion_basic.tsv | awk -F'\t' '{sum += $14; count++} END {print sum/count}')
mean_polya=$(tail -n +2 proportion_basic.tsv | awk -F'\t' '{sum += $16; count++} END {print sum/count}')

echo
echo "Mean proportions:"
echo "  Elongating: $(printf "%.3f" $mean_elongating)"
echo "  Polyadenylation: $(printf "%.3f" $mean_polya)"

echo

# Step 5: Output files
echo "=== Step 5: Output Files ==="
echo "Files created:"
ls -lh proportion_*.tsv

echo
echo "File descriptions:"
echo "  proportion_basic.tsv: All isoforms with ≥1 reads"
echo "  proportion_filtered.tsv: Isoforms with ≥10 reads"
echo
echo "Column descriptions:"
echo "  isoform_id: Unique isoform identifier"
echo "  gene_id: Associated gene ID"
echo "  chromosome, strand, start, end: Genomic coordinates"
echo "  total_reads: Total reads across all samples"
echo "  samples_with_reads: Number of samples with reads"
echo "  max/min/mean/std_reads_per_sample: Read count statistics"
echo "  *_reads: Read count per sample"
echo "  *_proportion: Proportion of reads per sample"

echo
echo "=== Next Steps ==="
echo "1. Visualize proportions with R/Python"
echo "2. Identify differentially expressed isoforms"
echo "3. Compare with reference annotation"
echo "4. Validate novel isoforms experimentally"

echo
echo "Proportion analysis completed successfully! 🎉"
