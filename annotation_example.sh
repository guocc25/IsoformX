#!/bin/bash
# Example: Annotation workflow with IsoformX
# This script demonstrates the complete workflow from BAM to annotated isoforms

echo "=== IsoformX Annotation Workflow Example ==="
echo

# Activate environment
conda activate isoformx

# Set paths
REFERENCE_GTF="$HOME/workshop/genome/Arabidopsis/Jin/Araport11_GTF_genes_transposons.Jan2023.gtf"
INPUT_BAM="example_data/13_elongating_6k.bam"
OUTPUT_PREFIX="annotation_example"

echo "Input BAM: $INPUT_BAM"
echo "Reference GTF: $REFERENCE_GTF"
echo

# Step 1: Identify isoforms
echo "=== Step 1: Identify Isoforms ==="
isoformx cluster \
    --input "$INPUT_BAM" \
    --output "${OUTPUT_PREFIX}_isoforms.gff3" \
    --min-reads 2 \
    --min-confidence 0.3 \
	--max-distance 0.5 \
    --verbose

if [ $? -eq 0 ]; then
    isoform_count=$(grep -c "transcript" "${OUTPUT_PREFIX}_isoforms.gff3")
    echo "✅ Identified $isoform_count isoforms"
else
    echo "❌ Error identifying isoforms"
    exit 1
fi

echo

# Step 2: Annotate isoforms
echo "=== Step 2: Annotate Isoforms ==="
isoformx annotate \
    --input "${OUTPUT_PREFIX}_isoforms.gff3" \
    --output "${OUTPUT_PREFIX}_annotated.gff3" \
    --reference "$REFERENCE_GTF" \
    --verbose

if [ $? -eq 0 ]; then
    echo "✅ Annotation completed"
else
    echo "❌ Error during annotation"
    exit 1
fi

echo

# Step 3: Analyze results
echo "=== Step 3: Analyze Results ==="

echo "Structural categories:"
grep "structural_category=" "${OUTPUT_PREFIX}_annotated.gff3" | \
    cut -d';' -f7 | cut -d'=' -f2 | sort | uniq -c | \
    while read count category; do
        echo "  $category: $count"
    done

echo
echo "High-confidence isoforms (confidence ≥ 0.6):"
high_conf_count=$(grep "confidence=0\.[6-9]" "${OUTPUT_PREFIX}_annotated.gff3" | grep -c "transcript")
echo "  Found $high_conf_count high-confidence isoforms"

echo
echo "Novel isoforms (NIC + NNC):"
novel_count=$(grep -E "structural_category=(NIC|NNC)" "${OUTPUT_PREFIX}_annotated.gff3" | grep -c "transcript")
echo "  Found $novel_count novel isoforms"

echo
echo "Known isoforms (FSM + ISM):"
known_count=$(grep -E "structural_category=(FSM|ISM)" "${OUTPUT_PREFIX}_annotated.gff3" | grep -c "transcript")
echo "  Found $known_count known isoforms"

echo

# Step 4: Extract specific categories
echo "=== Step 4: Extract Specific Categories ==="

echo "Extracting high-confidence novel isoforms..."
grep "confidence=0\.[6-9]" "${OUTPUT_PREFIX}_annotated.gff3" | \
    grep -E "structural_category=(NIC|NNC)" > "${OUTPUT_PREFIX}_novel_high_conf.gff3"

novel_high_conf_count=$(grep -c "transcript" "${OUTPUT_PREFIX}_novel_high_conf.gff3")
echo "  Extracted $novel_high_conf_count high-confidence novel isoforms"

echo "Extracting antisense isoforms..."
grep "structural_category=antisense" "${OUTPUT_PREFIX}_annotated.gff3" > "${OUTPUT_PREFIX}_antisense.gff3"

antisense_count=$(grep -c "transcript" "${OUTPUT_PREFIX}_antisense.gff3")
echo "  Extracted $antisense_count antisense isoforms"

echo

# Step 5: Summary
echo "=== Summary ==="
echo "Files created:"
ls -lh "${OUTPUT_PREFIX}"*.gff3

echo
echo "Workflow completed successfully! 🎉"
echo
echo "Next steps:"
echo "1. View results in genome browser (IGV, JBrowse)"
echo "2. Validate novel isoforms experimentally"
echo "3. Compare with other samples"
echo "4. Analyze functional implications"
