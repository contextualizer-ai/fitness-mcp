#!/usr/bin/env python3
"""
Process fit_t.tab file to create gene-condition pairs where |value| > 4
"""

import csv
import sys
from pathlib import Path
from typing import List, Tuple, Union


def generate_significant_fitness_pairs(
    fitness_file: Union[str, Path], pairs_file: Union[str, Path], threshold: float = 2.0
) -> None:
    """
    Generate gene-condition pairs with significant fitness effects from fitness data.

    Extracts gene-condition pairs where the absolute fitness value exceeds the threshold,
    creating a filtered dataset for network analysis and significant effect identification.

    Args:
        fitness_file: Path to input fitness data file (fit_t.tab)
        pairs_file: Path to output significant pairs file
        threshold: Minimum absolute fitness value to be considered significant (default: 2.0)
    """
    pairs: List[Tuple[str, str, float]] = []

    with open(fitness_file, "r", encoding="utf-8") as infile:
        reader = csv.reader(infile, delimiter="\t")

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
    with open(pairs_file, "w", encoding="utf-8") as outfile:
        writer = csv.writer(outfile, delimiter="\t")
        # Write header
        writer.writerow(["gene_id", "condition_id", "value"])
        # Write pairs
        for gene_id, condition_id, float_value in pairs:
            writer.writerow([gene_id, condition_id, str(float_value)])

    print(
        f"Generated {len(pairs)} significant gene-condition pairs with |value| > {threshold}"
    )
    print(f"Pairs file written to: {pairs_file}")


def main() -> None:
    # Default paths and threshold
    project_root = Path(__file__).parent.parent.parent
    input_file = project_root / "data" / "fit_t.tab"
    threshold = 4
    output_file = project_root / "data" / f"fit_t_pairs_threshold_{threshold}_long.tab"

    # Check if input file exists
    if not input_file.exists():
        print(f"Error: Input file {input_file} not found")
        sys.exit(1)

    # Generate significant pairs
    generate_significant_fitness_pairs(input_file, output_file, threshold)


if __name__ == "__main__":
    main()
