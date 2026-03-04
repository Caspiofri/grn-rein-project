#!/usr/bin/env python3
"""
L-Value Comparison: RE:IN Output vs Ground Truth

Compares the L-values inferred by the RE:IN engine against the ground truth
L-values derived from the original .aeon Boolean network.

Inputs:
  1. Ground truth JSON (from ground_truth_extractor.py) — contains truth tables
  2. RE:IN L-values CSV (from the RE:IN engine) — columns: solution_index, gene, L_value

The comparison:
  - Converts ground truth truth tables to L-values (binary encoding)
  - For each RE:IN solution, compares gene-by-gene against ground truth
  - Reports: match/mismatch per gene, overall accuracy, per-solution accuracy

Usage:
  python3 compare_l_values.py <ground_truth.json> <rein_l_values.csv>

Output:
  Prints comparison report to stdout.
  Optionally writes detailed CSV with --csv <output.csv>
"""

import csv
import json
import sys
from pathlib import Path
from typing import Dict, List, Optional, Tuple

# Import the L-value extraction logic from sibling module
try:
    from ground_truth_l_values import extract_l_values, truth_table_to_l_value
except ImportError:
    # If running standalone, add parent to path
    sys.path.insert(0, str(Path(__file__).resolve().parent))
    from ground_truth_l_values import extract_l_values, truth_table_to_l_value


def load_rein_csv(csv_path: Path) -> Dict[int, Dict[str, int]]:
    """Load RE:IN L-values from CSV.

    Returns:
        Dict mapping solution_index -> {gene_name: L_value}
    """
    solutions: Dict[int, Dict[str, int]] = {}

    with open(csv_path, newline="", encoding="utf-8") as f:
        reader = csv.DictReader(f)
        if reader.fieldnames is None:
            raise ValueError(f"Empty CSV: {csv_path}")

        # Detect column names (flexible)
        fields = set(reader.fieldnames)
        has_solution_idx = "solution_index" in fields

        for row in reader:
            sol_idx = int(row["solution_index"]) if has_solution_idx else 0
            gene = row["gene"].strip()
            l_value = int(row["L_value"])

            if sol_idx not in solutions:
                solutions[sol_idx] = {}
            solutions[sol_idx][gene] = l_value

    return solutions


def compare(
    ground_truth_path: Path,
    rein_csv_path: Path,
    output_csv_path: Optional[Path] = None,
):
    """Run full comparison and print report."""

    # --- Load ground truth ---
    with open(ground_truth_path, encoding="utf-8") as f:
        gt_data = json.load(f)

    gt_results = extract_l_values(gt_data)
    gt_l_values: Dict[str, Optional[int]] = {}
    gt_regulators: Dict[str, int] = {}
    for r in gt_results:
        gt_l_values[r["gene"]] = r["l_value"]
        gt_regulators[r["gene"]] = r["num_regulators"]

    # --- Load RE:IN output ---
    rein_solutions = load_rein_csv(rein_csv_path)

    model_name = gt_data.get("source_file", ground_truth_path.stem)

    # --- Header ---
    print("=" * 70)
    print(f"L-VALUE COMPARISON REPORT")
    print(f"=" * 70)
    print(f"Model:        {model_name}")
    print(f"Ground truth: {ground_truth_path}")
    print(f"RE:IN output: {rein_csv_path}")
    print(f"RE:IN solutions: {len(rein_solutions)}")
    print()

    # --- Per-solution comparison ---
    all_rows = []  # for CSV output
    best_accuracy = 0.0
    best_solution = -1

    for sol_idx in sorted(rein_solutions):
        rein_genes = rein_solutions[sol_idx]

        print(f"--- Solution {sol_idx} ---")
        matches = 0
        mismatches = 0
        skipped = 0
        total_compared = 0

        for gene in sorted(gt_l_values.keys()):
            gt_l = gt_l_values[gene]
            rein_l = rein_genes.get(gene)
            n_regs = gt_regulators.get(gene, 0)

            if gt_l is None:
                # Input node — skip comparison
                status = "INPUT"
                skipped += 1
                row = {
                    "solution_index": sol_idx,
                    "gene": gene,
                    "gt_l_value": "",
                    "rein_l_value": rein_l if rein_l is not None else "",
                    "num_regulators": n_regs,
                    "status": "INPUT",
                }
            elif rein_l is None:
                status = "MISSING"
                mismatches += 1
                total_compared += 1
                row = {
                    "solution_index": sol_idx,
                    "gene": gene,
                    "gt_l_value": gt_l,
                    "rein_l_value": "",
                    "num_regulators": n_regs,
                    "status": "MISSING_IN_REIN",
                }
                print(f"  {gene:30s}  GT={gt_l:5d}  REIN=???  MISSING")
            elif gt_l == rein_l:
                status = "MATCH"
                matches += 1
                total_compared += 1
                row = {
                    "solution_index": sol_idx,
                    "gene": gene,
                    "gt_l_value": gt_l,
                    "rein_l_value": rein_l,
                    "num_regulators": n_regs,
                    "status": "MATCH",
                }
                print(f"  {gene:30s}  GT={gt_l:5d}  REIN={rein_l:5d}  MATCH")
            else:
                status = "MISMATCH"
                mismatches += 1
                total_compared += 1
                row = {
                    "solution_index": sol_idx,
                    "gene": gene,
                    "gt_l_value": gt_l,
                    "rein_l_value": rein_l,
                    "num_regulators": n_regs,
                    "status": "MISMATCH",
                }
                print(f"  {gene:30s}  GT={gt_l:5d}  REIN={rein_l:5d}  MISMATCH")

            all_rows.append(row)

        # Extra genes in RE:IN not in ground truth
        for gene in sorted(rein_genes.keys()):
            if gene not in gt_l_values:
                all_rows.append({
                    "solution_index": sol_idx,
                    "gene": gene,
                    "gt_l_value": "",
                    "rein_l_value": rein_genes[gene],
                    "num_regulators": "?",
                    "status": "EXTRA_IN_REIN",
                })

        accuracy = (matches / total_compared * 100) if total_compared > 0 else 0.0
        print(f"\n  Matches: {matches}/{total_compared}  Accuracy: {accuracy:.1f}%")
        if skipped:
            print(f"  (Skipped {skipped} input nodes)")
        print()

        if accuracy > best_accuracy:
            best_accuracy = accuracy
            best_solution = sol_idx

    # --- Summary ---
    print("=" * 70)
    print("SUMMARY")
    print("=" * 70)
    if len(rein_solutions) > 1:
        print(f"Best solution: #{best_solution} with {best_accuracy:.1f}% accuracy")
    elif len(rein_solutions) == 1:
        print(f"Single solution accuracy: {best_accuracy:.1f}%")
    else:
        print("No RE:IN solutions found.")

    # --- Write CSV ---
    if output_csv_path and all_rows:
        fieldnames = ["solution_index", "gene", "gt_l_value", "rein_l_value",
                       "num_regulators", "status"]
        with open(output_csv_path, "w", newline="", encoding="utf-8") as f:
            writer = csv.DictWriter(f, fieldnames=fieldnames)
            writer.writeheader()
            writer.writerows(all_rows)
        print(f"\nDetailed CSV written to: {output_csv_path}")


def main():
    if len(sys.argv) < 3:
        print("Usage: python3 compare_l_values.py <ground_truth.json> <rein_l_values.csv> [--csv output.csv]")
        sys.exit(1)

    gt_path = Path(sys.argv[1])
    rein_path = Path(sys.argv[2])

    csv_out = None
    if "--csv" in sys.argv:
        idx = sys.argv.index("--csv")
        if idx + 1 < len(sys.argv):
            csv_out = Path(sys.argv[idx + 1])

    if not gt_path.exists():
        print(f"Error: ground truth not found: {gt_path}")
        sys.exit(1)
    if not rein_path.exists():
        print(f"Error: RE:IN CSV not found: {rein_path}")
        sys.exit(1)

    compare(gt_path, rein_path, csv_out)


if __name__ == "__main__":
    main()
