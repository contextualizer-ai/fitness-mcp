#!/usr/bin/env python3
"""
Process fit_t.tab file to create gene-condition pairs where |value| > 4
"""
import csv
import sys
from pathlib import Path


def process_fit_t_file(input_path, output_path, threshold=4):
    """
    Process fit_t.tab file and create pairs for values where |x| > threshold
    
    Args:
        input_path: Path to input fit_t.tab file
        output_path: Path to output pairs file
        threshold: Absolute value threshold (default: 4)
    """
    pairs = []
    
    with open(input_path, 'r', encoding='utf-8') as infile:
        reader = csv.reader(infile, delimiter='\t')
        
        # Read header
        header = next(reader)
        
        # Get condition columns (skip first 3 columns: locusId, sysName, desc)
        condition_columns = header[3:]
        
        # Process data rows
        for row in reader:
            gene_id = row[0]  # locusId
            
            # Check each condition value
            for i, value in enumerate(row[3:]):
                try:
                    float_value = float(value)
                    if abs(float_value) > threshold:
                        condition_id = condition_columns[i]
                        pairs.append((gene_id, condition_id, float_value))
                except (ValueError, IndexError):
                    # Skip non-numeric values or missing data
                    continue
    
    # Write pairs to output file
    with open(output_path, 'w', encoding='utf-8') as outfile:
        writer = csv.writer(outfile, delimiter='\t')
        # Write header
        writer.writerow(['gene_id', 'condition_id', 'value'])
        # Write pairs
        for gene_id, condition_id, value in pairs:
            writer.writerow([gene_id, condition_id, value])
    
    print(f"Processed {len(pairs)} gene-condition pairs with |value| > {threshold}")
    print(f"Output written to: {output_path}")


def main():
    # Default paths and threshold
    project_root = Path(__file__).parent.parent.parent
    input_file = project_root / 'data' / 'fit_t.tab'
    threshold = 2
    output_file = project_root / 'data' / f'fit_t_pairs_threshold_{threshold}_long.tab'
    
    # Check if input file exists
    if not input_file.exists():
        print(f"Error: Input file {input_file} not found")
        sys.exit(1)
    
    # Process the file
    process_fit_t_file(input_file, output_file, threshold)


if __name__ == '__main__':
    main()