#!/bin/bash
# Quick Start Script for IsoformX Pipeline
# This script processes all your BAM files automatically

echo "=== IsoformX Quick Start ==="
echo "Processing all BAM files in input_data/"

# Activate conda environment
echo "Activating conda environment..."
conda activate isoformx

# Check if environment is active
if ! command -v isoformx &> /dev/null; then
    echo "Error: isoformx command not found. Please activate the conda environment first:"
    echo "conda activate isoformx"
    exit 1
fi

# Create output directory
mkdir -p results

# Process all BAM files
echo "Processing BAM files..."

# Array of sample numbers
samples=(13 14 15 19 20 21)
#samples=(13)
conditions=("elongating" "polyadenylation")

for sample in "${samples[@]}"; do
    for condition in "${conditions[@]}"; do
        bam_file="../input_data/${sample}_tagged_${condition}.bam"
        
        if [ -f "$bam_file" ]; then
            echo "Processing: ${sample}_${condition}"
            
            # Run clustering
            isoformx cluster \
                --input "$bam_file" \
                --output "results/sample${sample}_${condition}_isoforms.gff3" \
				--annotation ~/workshop/genome/Arabidopsis/Jin/Athaliana_447_TAIR10.gtf \
                --min-reads 3 \
                --min-confidence 0.3 \
                --verbose
            
            # Filter high-confidence isoforms
            isoformx filter \
                --input "results/sample${sample}_${condition}_isoforms.gff3" \
                --output "results/sample${sample}_${condition}_high_conf.gff3" \
                --min-confidence 0.6
            
            # Annotate isoforms
            isoformx annotate \
                --input "results/sample${sample}_${condition}_isoforms.gff3" \
                --output "results/sample${sample}_${condition}_annotated.gff3" \
                --reference ~/workshop/genome/Arabidopsis/Jin/Araport11_GTF_genes_transposons.Jan2023.gtf \
                --verbose
            
            echo "Completed: ${sample}_${condition}"
        else
            echo "Warning: ${bam_file} not found, skipping..."
        fi
    done
done

echo ""
echo "=== Running Proportion Analysis ==="

# Collect all BAM files for proportion analysis
bam_files=()
for sample in "${samples[@]}"; do
    for condition in "${conditions[@]}"; do
        bam_file="../input_data/${sample}_tagged_${condition}.bam"
        if [ -f "$bam_file" ]; then
            bam_files+=("$bam_file")
        fi
    done
done

# Run proportion analysis if we have BAM files and annotated isoforms
if [ ${#bam_files[@]} -gt 0 ]; then
    echo "Running proportion analysis across ${#bam_files[@]} samples..."
    
    # Use the first annotated file as reference (they should all be similar)
    first_annotated=""
    for sample in "${samples[@]}"; do
        for condition in "${conditions[@]}"; do
            annotated_file="results/sample${sample}_${condition}_annotated.gff3"
            if [ -f "$annotated_file" ]; then
                first_annotated="$annotated_file"
                break 2
            fi
        done
    done
    
    if [ -n "$first_annotated" ]; then
        # Build proportion command with -b flags
        proportion_cmd="isoformx proportion"
        for bam_file in "${bam_files[@]}"; do
            proportion_cmd="$proportion_cmd -b $bam_file"
        done
        proportion_cmd="$proportion_cmd --isoform-gtf $first_annotated --output results/isoform_proportions.tsv --min-reads 5 --verbose"
        
        eval $proportion_cmd
        
        echo "Proportion analysis completed: results/isoform_proportions.tsv"
    else
        echo "Warning: No annotated files found for proportion analysis"
    fi
else
    echo "Warning: No BAM files found for proportion analysis"
fi

echo ""
echo "=== Processing Complete ==="
echo "Results saved in results/ directory:"
ls -la results/

echo ""
echo "Summary of high-confidence isoforms:"
for file in results/*_high_conf.gff3; do
    if [ -f "$file" ]; then
        count=$(grep -c "transcript" "$file" 2>/dev/null || echo "0")
        echo "$(basename $file): $count isoforms"
    fi
done

echo ""
echo "Summary of annotated isoforms:"
for file in results/*_annotated.gff3; do
    if [ -f "$file" ]; then
        count=$(grep -c "transcript" "$file" 2>/dev/null || echo "0")
        echo "$(basename $file): $count isoforms"
    fi
done

if [ -f "results/isoform_proportions.tsv" ]; then
    echo ""
    echo "Proportion analysis results:"
    total_isoforms=$(tail -n +2 results/isoform_proportions.tsv | wc -l)
    echo "Total isoforms analyzed: $total_isoforms"
    
    echo "Top 5 most abundant isoforms:"
    head -6 results/isoform_proportions.tsv | tail -5 | while IFS=$'\t' read -r isoform_id gene_id chromosome strand start end total_reads samples_with_reads max_reads min_reads mean_reads std_reads elongating_reads elongating_prop polya_reads polya_prop; do
        echo "  $isoform_id: $total_reads total reads ($samples_with_reads samples)"
    done
fi

echo ""
echo "Next steps:"
echo "1. Review results in results/ directory"
echo "2. Check annotated isoforms for structural categories"
echo "3. Analyze proportion results for differential expression"
echo "4. Visualize isoforms in genome browser (IGV, JBrowse)"
echo "5. Compare with reference annotation"
echo "6. Perform downstream analysis"
