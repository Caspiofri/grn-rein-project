# GRN-REIN Project

Algorithms for Synthesis of Gene Regulatory Networks using RE:IN.

Final-year Conputer Engineering project evaluating whether the RE:IN inference engine can
reconstruct Boolean network logic from (initial state -> final state)
observations.

## Quick Start

See `VM_ACCESS_GUIDE.md` to connect to the engine VM.
See `GRN_EXPERIMENT_PLAN_v2.md` for the full experiment strategy.

## Repository Structure

```
grn-rein-project/
├── pipeline/                          # Active pipeline scripts
│   ├── generate_experiment.py         # Generates ground truth JSON + .rein files from .aeon models
│   ├── compare_results.py             # Compares RE:IN L-value output against ground truth truth tables
│   └── _archive/                      # Superseded scripts (kept for reference)
│       ├── aeon_to_rein.py            # Old AEON-to-REIN converter
│       ├── rein_add_experiments.py    # Old experiment injector
│       ├── ground_truth_extractor.py  # Old GT extractor
│       ├── ground_truth_l_values.py   # Old binary-encoding GT L-value extractor
│       ├── compare_l_values.py        # Old binary-encoding comparison
│       ├── rein_regulatory_conditions.py  # Old R0-R19 matching attempt
│       └── vm_notebook_l_value_extraction.py  # F# extraction code reference
│
├── models/                            # Boolean network models in AEON format
│   ├── *.aeon                         # Full library (285 models, named ID_VARS_INPUTS)
│   ├── experiment_1_sanity/           # 5 small models for E1 sanity check
│   ├── experiment_2_stress/           # 20 models for E2 stress test
│   └── experiment_3_observation/      # 2 models for E3 observation sufficiency grid
│       ├── 110_9_0/110_9_0.aeon       #   C. crescentus cell cycle (best E1 performer)
│       └── 031_9_0/031_9_0.aeon       #   Yeast cell cycle (worst E1 performer)
│
├── results/
│   ├── experiment_1_N5/               # E1 sanity check (N=5, K=10) — completed
│   │   ├── <model_id>/
│   │   │   ├── <model_id>_ground_truth.json   # Pre-computed truth tables per gene
│   │   │   ├── <model_id>.rein                # Input file for RE:IN engine
│   │   │   └── <model_id>_l_values_extracted.csv  # RE:IN output L-values
│   │   └── experiment_1_comparison.csv         # Comparison results (30.6% recovery)
│   │
│   ├── experiment_1_N20/              # E1 re-run with more observations (N=20, K=15) — completed
│   │   ├── <model_id>/               #   (same structure as above)
│   │   └── experiment_1_N20_comparison.csv     # Comparison results (27.8% recovery)
│   │
│   └── experiment_3/                  # E3 observation sufficiency grid — awaiting RE:IN runs
│       ├── e3_run_manifest.csv        # Master tracking file (30 runs, all status=ready)
│       ├── 110_9_0/
│       │   └── K{5,10,15}_N{1,3,5,10,20}/     # 15 (K,N) combinations
│       │       ├── 110_9_0_ground_truth.json
│       │       └── 110_9_0.rein
│       └── 031_9_0/
│           └── K{5,10,15}_N{1,3,5,10,20}/     # 15 (K,N) combinations
│               ├── 031_9_0_ground_truth.json
│               └── 031_9_0.rein
│
├── reports/
│   └── model_screening_report.csv     # Initial model screening results
│
├── GRN_EXPERIMENT_PLAN_v2.md          # Current experiment plan
├── GRN_EXPERIMENT_PLAN_old.md         # Previous experiment plan (superseded)
└── VM_ACCESS_GUIDE.md                 # How to connect to the RE:IN VM
```

## Pipeline Usage

### Generate experiment files

```bash
python3 pipeline/generate_experiment.py <aeon_dir> <output_dir> --N <N> --K <K> --seed <seed>
```

Reads `.aeon` model files, simulates ground truth trajectories, and produces:
- `*_ground_truth.json` — truth tables for every gene + full trajectories
- `*.rein` — input file for the RE:IN engine (regulations + observations)

### Compare RE:IN output against ground truth

```bash
python3 pipeline/compare_results.py [results_dir]
```

For each gene, takes RE:IN's output L-value, computes the truth table that
L-value produces using RE:IN's R0-R19 regulatory condition formulas, and
compares it to the pre-computed ground truth truth table. Outputs a per-gene
match table and a recovery ratio.

## Experiments

| Experiment | Question | Models | Status |
|------------|----------|--------|--------|
| E1 Sanity Check | Does the basic pipeline work? | 5 small (5-9 genes, 0 inputs) | Completed (N=5 and N=20 runs) |
| E2 Stress Test | How does RE:IN scale? | 20 models (9-32 genes, 0-14 inputs) | Models selected, not yet run |
| E3 Observation Sufficiency | Does more data improve recovery? | 2 models x 15 (K,N) combos = 30 runs | .rein files generated, awaiting RE:IN |

## Key Findings (E1)

- N=5, K=10 recovery: **30.6%** (11/36 genes)
- N=20, K=15 recovery: **27.8%** (10/36 genes)
- RE:IN is mechanically correct but needs sufficient observations to
  distinguish the correct regulatory condition from degenerate alternatives
- Pure-activator and pure-repressor genes are especially vulnerable to
  degenerate L-values (e.g., `NoRepressors` is vacuously True when there
  are no repressors, making the gene always-ON)
