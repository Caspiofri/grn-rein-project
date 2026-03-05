#!/usr/bin/env python3
"""
Experiment 2 — Phase 1: Prepare RE:IN input files for 20 models.

Steps:
  1. Copy selected .aeon files into models/experiment_2_stress/
  2. Compute structural properties → e2_model_properties.csv
  3. Generate .rein + _ground_truth.json per model (K=15, N=10, seed=42)
  4. Validate all outputs

Usage:
  python3 pipeline/prepare_e2.py
"""

import csv
import json
import os
import re
import shutil
import sys
from pathlib import Path

# Add pipeline dir to path for generate_experiment imports
sys.path.insert(0, str(Path(__file__).resolve().parent))
from generate_experiment import parse_aeon, process_model

ROOT = Path(__file__).resolve().parent.parent
MODELS_DIR = ROOT / "models"
E2_MODELS_DIR = MODELS_DIR / "experiment_2_stress"
E2_REIN_DIR = ROOT / "e2_rein_files"

# Fixed experiment parameters
K = 15
N = 10
SEED = 42

# ── Selected 20 models ──
# Format: (aeon_filename_stem, bin, bio_name)
SELECTED = [
    # Bin A: 5-10 vars (must include 110_9_0 and 031_9_0 as anchors)
    ("110_9_0", "A", "Asymmetric Cell Division B (C. crescentus)"),
    ("031_9_0", "A", "Cell Cycle Transcription (Yeast)"),
    ("007_5_0", "A", "Cortical Area Development"),
    ("088_6_0", "A", "miR-9 Neurogenesis"),
    ("097_8_2", "A", "Drosophila Wings AP"),
    # Bin B: 11-20 vars (mix sparse and dense connectivity)
    ("177_11_0", "B", "Myeloid Progenitors"),
    ("058_14_0", "B", "Arabidopsis Thaliana Cell Cycle"),
    ("015_14_2", "B", "Neurotransmitter Signaling Pathway"),
    ("026_18_0", "B", "Budding Yeast Cell Cycle 2009"),
    ("086_18_2", "B", "Tumour Cell Invasion & Migration (Reduced)"),
    # Bin C: 21-30 vars
    ("108_23_2", "C", "Geroconversion"),
    ("279_22_2", "C", "EMT Hedgehog Signaling"),
    ("005_28_0", "C", "FA-BRCA Pathway"),
    ("044_25_1", "C", "Trichostrongylus Retortaeformis"),
    ("065_30_2", "C", "Tumour Cell Invasion & Migration"),
    # Bin D: 31-40 vars (expect majority timeouts)
    ("043_33_0", "D", "Bordetella Bronchiseptica"),
    ("013_32_2", "D", "Cholesterol Regulatory Pathway"),
    ("073_31_3", "D", "Lymphoid & Myeloid Cell Specification"),
    ("134_35_3", "D", "Rheumatoid Arthritis"),
    ("032_37_3", "D", "T-Cell Signalling 2006"),
]


def compute_structural_properties(aeon_path: Path) -> dict:
    """Parse an .aeon file and compute structural properties."""
    n_genes = 0
    n_edges = 0
    gene_regulators = {}  # gene -> {pos: set, neg: set}

    with open(aeon_path) as f:
        for line in f:
            line = line.split("//")[0].strip()
            if not line:
                continue

            # Count variable declarations (update functions)
            if line.startswith("$"):
                match = re.match(r"^\$(\S+):\s*(.+)$", line)
                if match:
                    n_genes += 1
                    gene = match.group(1)
                    if gene not in gene_regulators:
                        gene_regulators[gene] = {"positive": set(), "negative": set()}
                continue

            # Count and classify edges
            clean = line.replace(" ", "")
            edge_patterns = [
                ("-??", "both"), ("->?", "positive"), ("?->", "positive"),
                ("-|?", "negative"), ("?-|", "negative"), ("-?", "both"),
                ("->", "positive"), ("-|", "negative"),
            ]
            for symbol, kind in edge_patterns:
                if symbol in clean:
                    parts = clean.split(symbol, 1)
                    if len(parts) == 2:
                        src, dst = parts[0], parts[1].rstrip(";").strip()
                        if kind == "both":
                            n_edges += 2
                            if dst in gene_regulators:
                                gene_regulators[dst]["positive"].add(src)
                                gene_regulators[dst]["negative"].add(src)
                            else:
                                gene_regulators[dst] = {"positive": {src}, "negative": {src}}
                        else:
                            n_edges += 1
                            if dst not in gene_regulators:
                                gene_regulators[dst] = {"positive": set(), "negative": set()}
                            gene_regulators[dst][kind].add(src)
                    break

    # Compute pct_unbalanced
    # Unbalanced = gene has ≥1 regulator where ALL regulators carry the same sign
    genes_with_regs = 0
    unbalanced = 0
    for gene, regs in gene_regulators.items():
        all_regs = regs["positive"] | regs["negative"]
        if not all_regs:
            continue
        genes_with_regs += 1
        has_pos = len(regs["positive"]) > 0
        has_neg = len(regs["negative"]) > 0
        if not (has_pos and has_neg):
            unbalanced += 1

    pct_unbalanced = (unbalanced / genes_with_regs * 100) if genes_with_regs > 0 else 0

    return {
        "n_genes": n_genes,
        "n_edges": n_edges,
        "pct_unbalanced": round(pct_unbalanced, 1),
    }


def main():
    # ── Step 1: Copy .aeon files to experiment_2_stress/ ──
    print("=" * 70)
    print("EXPERIMENT 2 — Phase 1: Preparing RE:IN Input Files")
    print("=" * 70)

    E2_MODELS_DIR.mkdir(parents=True, exist_ok=True)
    E2_REIN_DIR.mkdir(parents=True, exist_ok=True)

    print(f"\nStep 1: Copying {len(SELECTED)} models to {E2_MODELS_DIR}")
    for stem, bin_label, bio_name in SELECTED:
        src = MODELS_DIR / f"{stem}.aeon"
        dst = E2_MODELS_DIR / f"{stem}.aeon"
        if not src.exists():
            print(f"  ERROR: {src} not found!")
            sys.exit(1)
        shutil.copy2(src, dst)
        print(f"  [{bin_label}] {stem} — {bio_name}")

    # ── Step 2: Compute structural properties ──
    print(f"\nStep 2: Computing structural properties")
    props_rows = []
    for stem, bin_label, bio_name in SELECTED:
        aeon_path = E2_MODELS_DIR / f"{stem}.aeon"
        props = compute_structural_properties(aeon_path)
        row = {
            "model_id": stem,
            "bin": bin_label,
            "n_genes": props["n_genes"],
            "n_edges": props["n_edges"],
            "pct_unbalanced": props["pct_unbalanced"],
            "bio_name": bio_name,
        }
        props_rows.append(row)
        print(f"  {stem:>10}  bin={bin_label}  genes={props['n_genes']:>2}  "
              f"edges={props['n_edges']:>3}  unbal={props['pct_unbalanced']:>5}%")

    props_path = ROOT / "e2_model_properties.csv"
    with open(props_path, "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=["model_id", "bin", "n_genes", "n_edges",
                                           "pct_unbalanced", "bio_name"])
        w.writeheader()
        w.writerows(props_rows)
    print(f"\n  Saved: {props_path}")

    # ── Step 3: Generate .rein + ground_truth.json ──
    print(f"\nStep 3: Generating .rein + ground_truth.json (K={K}, N={N}, seed={SEED})")
    for i, (stem, bin_label, bio_name) in enumerate(SELECTED):
        aeon_path = E2_MODELS_DIR / f"{stem}.aeon"
        model_seed = SEED + i
        process_model(aeon_path, E2_REIN_DIR, N=N, K=K, seed=model_seed)

    # ── Step 4: Validate ──
    print(f"\nStep 4: Validation")
    print("-" * 50)
    all_ok = True
    for stem, bin_label, bio_name in SELECTED:
        model_dir = E2_REIN_DIR / stem
        rein = model_dir / f"{stem}.rein"
        gt = model_dir / f"{stem}_ground_truth.json"

        errors = []
        if not rein.exists():
            errors.append("missing .rein")
        if not gt.exists():
            errors.append("missing ground_truth.json")

        if not errors:
            content = rein.read_text()
            exp_count = content.count("#Experiment")
            if exp_count < N:
                errors.append(f"only {exp_count} #Experiment (expected >= {N})")

            gt_data = json.loads(gt.read_text())
            if len(gt_data.get("genes", [])) == 0:
                errors.append("empty genes list")

        if errors:
            print(f"  FAIL: {stem} — {', '.join(errors)}")
            all_ok = False
        else:
            gene_count = gt_data["total_genes"]
            func_count = gt_data["total_update_functions"]
            input_count = gt_data["total_inputs"]
            print(f"  OK: {stem} ({gene_count} genes, {func_count} functions, "
                  f"{input_count} inputs, {exp_count} experiments)")

    print()
    if all_ok:
        print("ALL 20 MODELS VALIDATED SUCCESSFULLY")
    else:
        print("SOME MODELS FAILED VALIDATION")
        sys.exit(1)

    print(f"\nDeliverables:")
    print(f"  {E2_REIN_DIR}/  — 20 subdirectories with .rein + _ground_truth.json")
    print(f"  {props_path}    — structural properties CSV")


if __name__ == "__main__":
    main()
