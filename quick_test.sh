#!/bin/bash
# Quick test script for IsoformX example data
# Simple version for quick validation

echo "=== IsoformX Quick Test ==="

# Activate environment
conda activate isoformx

# Test elongating data
echo "Testing elongating data..."
isoformx cluster --input example_data/13_elongating_6k.bam --output quick_test_elongating.gff3 --verbose
elongating_count=$(grep -c "transcript" quick_test_elongating.gff3)
echo "Found $elongating_count elongating isoforms"

# Test polyadenylation data  
echo "Testing polyadenylation data..."
isoformx cluster --input example_data/13_polyadenylation_6k.bam --output quick_test_polyadenylation.gff3 --verbose
polya_count=$(grep -c "transcript" quick_test_polyadenylation.gff3)
echo "Found $polya_count polyadenylation isoforms"

# Test mock data
echo "Testing mock data..."
isoformx mock --n-reads 20 --n-genes 3 --output quick_test_mock.gff3 --verbose
mock_count=$(grep -c "transcript" quick_test_mock.gff3)
echo "Found $mock_count mock isoforms"

echo
echo "=== Results ==="
echo "Elongating: $elongating_count isoforms"
echo "Polyadenylation: $polya_count isoforms" 
echo "Mock: $mock_count isoforms"
echo
echo "Files created:"
ls -lh quick_test_*.gff3

echo
echo "Quick test completed! ✅"
