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
            
            echo "Completed: ${sample}_${condition}"
        else
            echo "Warning: ${bam_file} not found, skipping..."
        fi
    done
done

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
echo "Next steps:"
echo "1. Review results in results/ directory"
echo "2. Visualize isoforms in genome browser"
echo "3. Compare with reference annotation"
echo "4. Perform downstream analysis"
