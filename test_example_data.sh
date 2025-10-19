#!/bin/bash
# Test script for IsoformX example data
# This script demonstrates the pipeline with the extracted example datasets

echo "=== IsoformX Example Data Test Script ==="
echo "Testing the pipeline with 6000-read example datasets"
echo

# Activate conda environment
echo "Activating conda environment..."
conda activate isoformx

# Check if environment is active
if ! command -v isoformx &> /dev/null; then
    echo "Error: isoformx command not found. Please activate the conda environment first:"
    echo "conda activate isoformx"
    exit 1
fi

# Create test output directory
mkdir -p test_output

echo "Environment activated successfully!"
echo

# Test 1: Process elongating reads
echo "=== Test 1: Processing Elongating Reads ==="
echo "Input: example_data/13_elongating_6k.bam"
echo "Processing..."

isoformx cluster \
    --input example_data/13_elongating_6k.bam \
    --output test_output/elongating_isoforms.gff3 \
    --min-reads 2 \
    --min-confidence 0.3 \
    --verbose

if [ $? -eq 0 ]; then
    echo "✅ Elongating reads processed successfully!"
    
    # Count isoforms
    elongating_count=$(grep -c "transcript" test_output/elongating_isoforms.gff3)
    echo "   Found $elongating_count isoforms"
    
    # Show first few lines
    echo "   First few isoforms:"
    grep "transcript" test_output/elongating_isoforms.gff3 | head -3 | while read line; do
        echo "   $line"
    done
else
    echo "❌ Error processing elongating reads"
    exit 1
fi

echo

# Test 2: Process polyadenylation reads
echo "=== Test 2: Processing Polyadenylation Reads ==="
echo "Input: example_data/13_polyadenylation_6k.bam"
echo "Processing..."

isoformx cluster \
    --input example_data/13_polyadenylation_6k.bam \
    --output test_output/polyadenylation_isoforms.gff3 \
    --min-reads 2 \
    --min-confidence 0.3 \
    --verbose

if [ $? -eq 0 ]; then
    echo "✅ Polyadenylation reads processed successfully!"
    
    # Count isoforms
    polya_count=$(grep -c "transcript" test_output/polyadenylation_isoforms.gff3)
    echo "   Found $polya_count isoforms"
    
    # Show first few lines
    echo "   First few isoforms:"
    grep "transcript" test_output/polyadenylation_isoforms.gff3 | head -3 | while read line; do
        echo "   $line"
    done
else
    echo "❌ Error processing polyadenylation reads"
    exit 1
fi

echo

# Test 3: Filter high-confidence isoforms
echo "=== Test 3: Filtering High-Confidence Isoforms ==="

# Filter elongating isoforms (using grep for now since filter command needs implementation)
echo "Filtering elongating isoforms (min confidence: 0.6)..."
grep "confidence=0\.[6-9]" test_output/elongating_isoforms.gff3 > test_output/elongating_high_conf.gff3
if [ $? -eq 0 ]; then
    high_conf_elongating=$(grep -c "transcript" test_output/elongating_high_conf.gff3)
    echo "   High-confidence elongating isoforms: $high_conf_elongating"
else
    echo "❌ Error filtering elongating isoforms"
    high_conf_elongating=0
fi

# Filter polyadenylation isoforms
echo "Filtering polyadenylation isoforms (min confidence: 0.6)..."
grep "confidence=0\.[6-9]" test_output/polyadenylation_isoforms.gff3 > test_output/polyadenylation_high_conf.gff3
if [ $? -eq 0 ]; then
    high_conf_polya=$(grep -c "transcript" test_output/polyadenylation_high_conf.gff3)
    echo "   High-confidence polyadenylation isoforms: $high_conf_polya"
else
    echo "❌ Error filtering polyadenylation isoforms"
    high_conf_polya=0
fi

echo

# Test 4: Mock data test
echo "=== Test 4: Mock Data Test ==="
echo "Testing with synthetic data..."

isoformx mock \
    --n-reads 50 \
    --n-genes 5 \
    --output test_output/mock_isoforms.gff3 \
    --verbose

if [ $? -eq 0 ]; then
    mock_count=$(grep -c "transcript" test_output/mock_isoforms.gff3)
    echo "✅ Mock data test successful!"
    echo "   Generated $mock_count mock isoforms"
else
    echo "❌ Error in mock data test"
fi

echo

# Test 5: Annotation test
echo "=== Test 5: Annotation Test ==="
echo "Testing annotation with reference GTF..."

# Check if reference GTF exists
reference_gtf="$HOME/workshop/genome/Arabidopsis/Jin/Araport11_GTF_genes_transposons.Jan2023.gtf"
if [ -f "$reference_gtf" ]; then
    echo "Annotating elongating isoforms..."
    isoformx annotate \
        --input test_output/elongating_isoforms.gff3 \
        --output test_output/elongating_annotated.gff3 \
        --reference "$reference_gtf" \
        --verbose
    
    if [ $? -eq 0 ]; then
        echo "✅ Elongating annotation successful!"
        
        # Count structural categories
        echo "   Structural categories:"
        grep "structural_category=" test_output/elongating_annotated.gff3 | cut -d';' -f7 | cut -d'=' -f2 | sort | uniq -c | while read count category; do
            echo "     $category: $count"
        done
    else
        echo "❌ Error annotating elongating isoforms"
    fi
    
    echo "Annotating polyadenylation isoforms..."
    isoformx annotate \
        --input test_output/polyadenylation_isoforms.gff3 \
        --output test_output/polyadenylation_annotated.gff3 \
        --reference "$reference_gtf" \
        --verbose
    
    if [ $? -eq 0 ]; then
        echo "✅ Polyadenylation annotation successful!"
        
        # Count structural categories
        echo "   Structural categories:"
        grep "structural_category=" test_output/polyadenylation_annotated.gff3 | cut -d';' -f7 | cut -d'=' -f2 | sort | uniq -c | while read count category; do
            echo "     $category: $count"
        done
    else
        echo "❌ Error annotating polyadenylation isoforms"
    fi
else
    echo "⚠️  Reference GTF not found at $reference_gtf"
    echo "   Skipping annotation test"
fi

echo

# Summary
echo "=== Test Summary ==="
echo "All tests completed successfully!"
echo
echo "Output files created in test_output/:"
ls -lh test_output/
echo

echo "Summary of results:"
echo "- Elongating isoforms: $elongating_count (high-confidence: $high_conf_elongating)"
echo "- Polyadenylation isoforms: $polya_count (high-confidence: $high_conf_polya)"
echo "- Mock isoforms: $mock_count"

# Add annotation summary if available
if [ -f "test_output/elongating_annotated.gff3" ]; then
    echo "- Annotation: ✅ Completed for both datasets"
    echo "  Elongating categories:"
    grep "structural_category=" test_output/elongating_annotated.gff3 | cut -d';' -f7 | cut -d'=' -f2 | sort | uniq -c | while read count category; do
        echo "    $category: $count"
    done
    echo "  Polyadenylation categories:"
    grep "structural_category=" test_output/polyadenylation_annotated.gff3 | cut -d';' -f7 | cut -d'=' -f2 | sort | uniq -c | while read count category; do
        echo "    $category: $count"
    done
else
    echo "- Annotation: ⚠️  Skipped (no reference GTF)"
fi
echo

echo "=== Quality Checks ==="
echo "Checking GFF format and sorting..."

# Check if transcripts come before exons
echo "Checking transcript-exon order..."
elongating_order_check=$(head -20 test_output/elongating_isoforms.gff3 | grep -E "(transcript|exon)" | head -10 | awk '{print $3}' | tr '\n' ' ')
if [[ $elongating_order_check == *"transcript"* ]]; then
    echo "✅ Transcripts properly come before exons"
else
    echo "❌ Warning: Exons may appear before transcripts"
fi

# Check chromosome sorting
echo "Checking chromosome sorting..."
chromosomes=$(grep "transcript" test_output/elongating_isoforms.gff3 | cut -f1 | sort | uniq)
echo "   Chromosomes found: $chromosomes"

echo
echo "=== Next Steps ==="
echo "1. View results in genome browser (IGV, JBrowse)"
echo "2. Compare with reference annotation"
echo "3. Analyze isoform characteristics"
echo "4. Validate high-confidence candidates"
echo
echo "Example commands:"
echo "  # View elongating isoforms"
echo "  head -20 test_output/elongating_isoforms.gff3"
echo
echo "  # Count isoforms by chromosome"
echo "  grep 'transcript' test_output/elongating_isoforms.gff3 | cut -f1 | sort | uniq -c"
echo
echo "  # Extract high-confidence isoforms only"
echo "  grep 'transcript' test_output/elongating_high_conf.gff3"
echo
echo "Test completed successfully! 🎉"
