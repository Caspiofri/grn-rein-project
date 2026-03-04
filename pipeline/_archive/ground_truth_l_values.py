#!/usr/bin/env python3
"""
Ground Truth L-Value Extractor

Reads a ground truth JSON (produced by ground_truth_extractor.py) and converts
each gene's truth table into an integer L-value using binary encoding.

Encoding convention:
  For a gene with k regulators [r1, r2, ..., rk] sorted alphabetically,
  the truth table has 2^k rows. The rows are ordered by binary counting
  of the regulator values (all 0s first, incrementing). The L-value is the
  integer whose binary representation is the sequence of output bits:

    row 0 (all regulators=0) -> bit 0 (least significant)
    row 1                    -> bit 1
    ...
    row 2^k - 1 (all=1)     -> bit 2^k - 1 (most significant)

  Example: truth table [1, 1, 1, 0] -> L = 0b0111 = 7

Input nodes (no update function) get L_value = None in JSON, skipped in CSV.

Usage:
  python3 ground_truth_l_values.py <ground_truth.json> [--csv output.csv] [--json output.json]

If no output flags given, prints summary to stdout and writes CSV to
<input_stem>_l_values.csv alongside the input file.
"""

import json
import csv
import re
import sys
from pathlib import Path
from typing import Dict, List, Optional, Tuple


def parse_truth_table_key(key_str: str) -> List[Tuple[str, int]]:
    """Parse a truth table key like 'v_Pelle=0, v_Slmb=1' into sorted pairs."""
    pairs = []
    for part in key_str.split(","):
        part = part.strip()
        name, val = part.rsplit("=", 1)
        pairs.append((name.strip(), int(val)))
    return sorted(pairs, key=lambda x: x[0])


def truth_table_to_l_value(truth_table: Dict[str, int], regulators: List[str]) -> int:
    """Convert a truth table dict to an integer L-value.

    Args:
        truth_table: Dict mapping 'r1=0, r2=1' -> 0 or 1
        regulators: Sorted list of regulator names

    Returns:
        Integer L-value (binary encoding of the truth table outputs).
    """
    if not truth_table or not regulators:
        return 0

    k = len(regulators)
    l_value = 0

    # Iterate over all 2^k combinations in binary counting order
    for row_idx in range(2 ** k):
        # Build the assignment for this row
        assignment = {}
        for bit_pos, reg in enumerate(regulators):
            assignment[reg] = (row_idx >> bit_pos) & 1

        # Build the key string to look up in the truth table
        key_str = ", ".join(f"{reg}={assignment[reg]}" for reg in regulators)

        output = truth_table.get(key_str)
        if output is None:
            # Try alternative key formats
            for tt_key, tt_val in truth_table.items():
                parsed = parse_truth_table_key(tt_key)
                parsed_dict = dict(parsed)
                if parsed_dict == assignment:
                    output = tt_val
                    break

        if output is None:
            print(f"  WARNING: missing truth table entry for {key_str}", file=sys.stderr)
            output = 0

        if output:
            l_value |= (1 << row_idx)

    return l_value


def extract_l_values(ground_truth: Dict) -> List[Dict]:
    """Extract L-values for all genes from a ground truth JSON.

    Args:
        ground_truth: Parsed ground truth JSON dict.

    Returns:
        List of dicts with keys: gene, l_value, num_regulators, is_input, truth_table_size
    """
    results = []

    for gene_info in ground_truth.get("genes", []):
        name = gene_info["name"]
        is_input = gene_info.get("is_input", False)

        if is_input:
            results.append({
                "gene": name,
                "l_value": None,
                "num_regulators": 0,
                "is_input": True,
                "truth_table_size": 0,
            })
            continue

        regulators = sorted(r["name"] for r in gene_info.get("regulators", []))
        truth_table = gene_info.get("truth_table", {})

        l_value = truth_table_to_l_value(truth_table, regulators)

        results.append({
            "gene": name,
            "l_value": l_value,
            "num_regulators": len(regulators),
            "is_input": False,
            "truth_table_size": len(truth_table),
        })

    return results


def write_csv(results: List[Dict], csv_path: Path):
    """Write L-values to CSV (matching l_values_extracted.csv format)."""
    with open(csv_path, "w", newline="", encoding="utf-8") as f:
        writer = csv.writer(f)
        writer.writerow(["gene", "l_value", "num_regulators", "is_input"])
        for r in results:
            if r["is_input"]:
                writer.writerow([r["gene"], "", r["num_regulators"], "TRUE"])
            else:
                writer.writerow([r["gene"], r["l_value"], r["num_regulators"], "FALSE"])


def write_json(results: List[Dict], json_path: Path):
    """Write L-values to JSON."""
    with open(json_path, "w", encoding="utf-8") as f:
        json.dump(results, f, indent=2)


def main():
    if len(sys.argv) < 2:
        print("Usage: python3 ground_truth_l_values.py <ground_truth.json> [--csv out.csv] [--json out.json]")
        sys.exit(1)

    gt_path = Path(sys.argv[1])
    if not gt_path.exists():
        print(f"Error: file not found: {gt_path}")
        sys.exit(1)

    # Parse optional output flags
    csv_path = None
    json_path = None
    args = sys.argv[2:]
    i = 0
    while i < len(args):
        if args[i] == "--csv" and i + 1 < len(args):
            csv_path = Path(args[i + 1])
            i += 2
        elif args[i] == "--json" and i + 1 < len(args):
            json_path = Path(args[i + 1])
            i += 2
        else:
            i += 1

    # Default CSV output alongside input file
    if csv_path is None and json_path is None:
        csv_path = gt_path.parent / f"{gt_path.stem}_l_values.csv"

    # Load and process
    with open(gt_path, encoding="utf-8") as f:
        ground_truth = json.load(f)

    results = extract_l_values(ground_truth)

    # Write outputs
    if csv_path:
        write_csv(results, csv_path)
        print(f"Wrote CSV: {csv_path}")

    if json_path:
        write_json(results, json_path)
        print(f"Wrote JSON: {json_path}")

    # Print summary
    model = ground_truth.get("source_file", gt_path.name)
    non_input = [r for r in results if not r["is_input"]]
    inputs = [r for r in results if r["is_input"]]

    print(f"\nModel: {model}")
    print(f"Total genes: {len(results)} ({len(non_input)} with update functions, {len(inputs)} inputs)")
    print(f"\nGround truth L-values:")
    for r in results:
        if r["is_input"]:
            print(f"  {r['gene']:30s}  INPUT (no L-value)")
        else:
            print(f"  {r['gene']:30s}  L={r['l_value']:5d}  ({r['num_regulators']} regulators, {r['truth_table_size']} rows)")


if __name__ == "__main__":
    main()
