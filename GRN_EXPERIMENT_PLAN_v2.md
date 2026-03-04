# GRN-REIN Project — Experiment Strategy & Guidelines
> **Version 2.0** — Updated after E1 sanity check revealed degenerate L-value problem  
> Keep this at the root of the GitHub repo alongside `README.md`.

---

## 🏗️ Overall Architecture

```
┌─────────────────┐     git push      ┌──────────────────────┐
│  Raspberry Pi   │ ─────────────────▶│     GitHub Repo      │
│  - pipelines    │                   │  - .aeon files       │
│  - ground truth │ ◀─────────────────│  - .rein files       │
│  - analysis     │     git pull      │  - results JSONs     │
└─────────────────┘                   └──────────────────────┘
                                                 │
                                            git pull
                                                 ▼
┌──────────────────────────────────────────────────────┐
│              Work Computer (Windows / VM)            │
│  - RE:IN engine                                      │
│  - pulls .rein files, pushes result JSONs back       │
└──────────────────────────────────────────────────────┘
```

---

## 📐 The Three Evaluation Metrics

Before describing the experiments, it is critical to understand the three-layer
evaluation framework. All experiments report these three metrics consistently.

### Metric 1 — SAT Rate
**Question:** Did RE:IN find any solution at all?

- **How measured:** Check if the output JSON contains at least one L-value assignment.
- **SAT** = engine returned solutions.
- **UNSAT/Timeout** = engine failed to synthesize any network.
- **Why it matters:** This is the baseline. If the engine cannot find a solution,
  nothing else is measurable.

### Metric 2 — Functional Accuracy (the real correctness test)
**Question:** Does RE:IN's synthesized network *behave* identically to the original?

- **How measured:** Generate 200 random initial states the engine never saw during
  synthesis. Simulate both the ground-truth AEON network and RE:IN's synthesized
  network for K synchronous steps. Count the percentage of runs where final states
  match exactly.
- **Why 200 fresh states:** RE:IN only saw N×(init,final) pairs during synthesis.
  Fresh states test generalization — whether the synthesized logic is genuinely
  correct or just satisfied the training constraints by coincidence.
- **Why this is the real test:** Two networks can have different L-values but
  identical behavior on the specific states RE:IN was shown. Functional accuracy
  tests the *global* behavior of the network on the full state space.

### Metric 3 — Gene-Level Truth Table Match
**Question:** For each gene, did RE:IN recover the exact same regulatory logic?

- **How measured:** Compute the truth table for each gene using RE:IN's L-value
  and the gene's regulator topology. Compare to ground truth truth table.
- **Important split — Balanced vs. Unbalanced genes:**
  - **Balanced genes** (have both activators AND repressors): truth table match
    is a clean, unambiguous test.
  - **Unbalanced genes** (only activators OR only repressors): subject to
    degenerate L-value problem (see below). Report separately.

---

## ⚠️ The Degenerate L-Value Problem (Key Finding from E1)

During Experiment 1, a systematic failure pattern was identified and confirmed
across multiple real biological models.

### What Happens

For genes that have *only repressors* or *only activators*, certain R-conditions
contain sub-expressions that become vacuously True or False, producing a truth
table that is observationally valid but biologically wrong:

- **Pure-repressor gene:** `NoActivators` is always True (empty activator set),
  so `¬NoActivators` is always False. Any R-condition containing `¬NoActivators`
  degenerates — one half of the OR/AND is dead.
- **Pure-activator gene:** `NoRepressors` is always True (no repressors), so
  any R-condition containing `NoRepressors` is trivially satisfied → gene becomes
  always ON regardless of its activators.

### Concrete Examples from E1

**`v_Coup_fti`** — model 007_5_0 (mammalian cortex, Bhatt et al. 2010):
- Regulators: repressors Fgf8 and Sp8; zero activators.
- Ground truth: L=9 (NoRepressors) → truth table [1,0,0,0]
  → active only when BOTH repressors absent (strict inhibition by either repressor)
- RE:IN returned: L=17 (¬NoActivators ∨ ¬AllRepressors) → collapses to ¬AllRepressors
  → truth table [1,1,1,0] → active unless both repressors present simultaneously
- Biological error: in the real cortex, one repressor alone silences Coup-tfi.
  RE:IN's version requires both repressors to act together.

**`v_HCM1`** — model 031_9_0 (yeast cell cycle):
- Regulators: activators MBF and SBF; zero repressors.
- Ground truth: L=4 (AllActivators) → truth table [0,0,0,1]
  → active only when both activators present
- RE:IN returned: L=10 (NR ∨ (AA ∧ ¬AR)) → NR is always True → always ON
  → truth table [1,1,1,1]
- Consequence: HCM1 is constitutively active, cascading errors through the
  yeast cell cycle network and producing 0% functional accuracy.

### Why RE:IN Picks These

The degenerate L-value satisfies all N×(init→final) training observations because
the distinguishing regulator states — where one repressor is present but the other
is absent — rarely appear in the *final states* of K-step trajectories. The network
dynamics cause co-repressors to synchronize quickly, so RE:IN sees no evidence
to distinguish L=9 from L=17.

**This is a finding about RE:IN's fundamental limitations with unbalanced network
topologies — not a bug in your pipeline.**

### Formal Definition: Degenerate L-value

A gene G with L-value L is classified as **degenerate** if any of the following hold:
1. G has no repressors AND the R-condition for L contains `NoRepressors` (vacuously True)
2. G has no activators AND the R-condition for L contains `¬NoActivators` (vacuously False)
3. The resulting effective truth table differs from the ground truth truth table

---

## 🗓️ Execution Order

| Step | Experiment | Status       | Why this order                                        |
|------|------------|--------------|-------------------------------------------------------|
| 1    | **E1** Sanity Check      | ✅ Done     | Confirmed pipeline; revealed degenerate L-value issue |
| 2    | **E3** Observation Grid  | **Next**    | 30 fast runs; answers whether more data helps         |
| 3    | **E2** Stress Test       | Last        | Most expensive; do after methodology is validated     |

---

## 🔬 Experiment 1: Sanity Check ✅ COMPLETED

### Results Summary (N=20, K=15)

| Model   | Biological System              | SAT | Functional Acc. | Genes Recovered |
|---------|--------------------------------|-----|-----------------|-----------------|
| 007_5_0 | Mammalian cortex (5 genes)     | ✅  | 23.5%           | 0/5             |
| 031_9_0 | Yeast cell cycle (9 genes)     | ✅  | 0.0%            | 3/9             |
| 088_6_0 | Zebrafish neural dev. (6 genes)| ✅  | 0.0%            | 0/6             |
| 110_9_0 | C. crescentus cell cycle (9)   | ✅  | 25.0%           | 7/9             |
| 158_7_0 | —                              | ❌  | N/A (Timeout)   | —               |

**Overall: 4/5 SAT · 12.1% average functional accuracy · 10/29 genes recovered**

### Key Findings

1. RE:IN reliably finds *some* solution (80% SAT on these models).
2. Functional accuracy is low overall (12.1%) — synthesized networks diverge
   from ground truth when tested on unseen initial states.
3. The accuracy gap is directly caused by degenerate L-values on unbalanced genes.
4. Model 110_9_0 is the clear outlier: 7/9 correct, 25% functional accuracy.
   Its network topology has fewer unbalanced genes, making it more favorable for
   RE:IN. This motivates using it as the reference model for E3.
5. 158_7_0 timed out — this is a complexity limit, not a logic error.

---

## 🔬 Experiment 3: Observation Sufficiency — *"Can More Data Fix the Degeneracy?"*

### Goal

The most scientifically original experiment in this project. It directly tests
the central question raised by E1: is the degenerate L-value problem caused by
*insufficient observations*, or is it a structural property of the network topology
that no amount of data can resolve?

**Two competing hypotheses:**
- **H1 (Optimistic):** With enough (K, N), RE:IN will eventually encounter the
  distinguishing states and converge on the correct L-value.
- **H2 (Pessimistic):** The degenerate L-values satisfy all possible (init→final)
  constraints for these topologies. The distinguishing states never appear in final
  states regardless of K or N. More data cannot help.

### Model Selection

Run the full grid on **two reference models** chosen from E1:

| Reference | Model    | E1 Performance | Role                                   |
|-----------|----------|----------------|----------------------------------------|
| Primary   | 110_9_0  | 25% func. acc. | Best E1 performer; shows improvement ceiling |
| Secondary | 031_9_0  | 0% func. acc.  | Worst performer; tests whether hard cases can recover |

Holding the model constant means the only variable is observation quantity —
giving a clean, unconfounded result.

### Protocol

Run the full (K, N) grid on each reference model:

|           | N=1 | N=3 | N=5 | N=10 | N=20 |
|-----------|-----|-----|-----|------|------|
| **K=5**   | run | run | run | run  | run  |
| **K=10**  | run | run | run | run  | run  |
| **K=15**  | run | run | run | run  | run  |

**15 runs × 2 models = 30 runs total.**

For each run record:

| Column              | Description                                           |
|---------------------|-------------------------------------------------------|
| model_id            | e.g. 110_9_0                                          |
| K                   | Steps per trajectory                                  |
| N                   | Number of experiments injected                        |
| sat_status          | SAT / UNSAT / Timeout                                 |
| functional_accuracy | % of 200 fresh trajectories matching ground truth     |
| gene_recovery       | % genes with exact truth table match                  |
| stability_ratio     | % genes with variance=0 across all solutions          |
| degenerate_rate     | % genes with a degenerate L-value                     |
| n_solutions         | Number of concrete networks found by SMT              |
| runtime_s           | Wall-clock seconds                                    |

### Interpreting Results

| Pattern                                              | Conclusion                                             |
|------------------------------------------------------|--------------------------------------------------------|
| Functional accuracy increases with N and K           | H1 confirmed — more data helps                         |
| Accuracy plateaus below 100% regardless of (K, N)   | H2 confirmed — topology is the barrier                 |
| Stability increases but accuracy stays low           | RE:IN confidently converges on the *wrong* answer      |
| Only balanced genes improve; unbalanced stay wrong   | Degeneracy is topology-driven, data cannot fix it      |

### Output for Project Book

1. **Two heatmaps** (K × N grid): cell value = functional accuracy %, one per model.
2. **Convergence line chart**: functional accuracy vs. N at fixed K=15 for both models.
3. **Stability vs. accuracy scatter**: does stability ratio predict accuracy?
4. One clear conclusion statement on whether H1 or H2 is supported.

---

## 🔬 Experiment 2: Stress Test — *"Where Does RE:IN Break?"*

### Goal

Characterize RE:IN's computational limits as a function of network size and
connectivity. Produce a runtime curve. Also measure whether the degenerate
L-value problem scales with network size (i.e., does larger N mean more
unbalanced genes and therefore lower accuracy?).

### Model Selection

**20 models** from Biodivine, stratified by gene count:

| Bin | Gene count | Models | Notes                                   |
|-----|------------|--------|-----------------------------------------|
| A   | 5–10       | 5      | Include E1 models as anchors            |
| B   | 11–20      | 5      | Mix sparse and dense connectivity       |
| C   | 21–30      | 5      | Expect runtimes to rise sharply         |
| D   | 31–40      | 5      | Expect majority timeouts                |

Before running each model, record its structural properties:
- `n_genes`, `n_edges`
- `pct_unbalanced` — fraction of genes with only activators or only repressors

This lets you later test the correlation: *does higher unbalanced % predict
lower functional accuracy?*

### Fixed Parameters

| Parameter   | Value       |
|-------------|-------------|
| K           | 15          |
| N           | 10          |
| Timeout     | 300 seconds |
| Fresh tests | 200         |

### Per-model columns to collect

| Column               | Description                                       |
|----------------------|---------------------------------------------------|
| model_id             | Biodivine ID                                      |
| n_genes              | Total gene count                                  |
| n_edges              | Total edge count                                  |
| pct_unbalanced       | % genes with only activators or only repressors   |
| sat_status           | SAT / UNSAT / Timeout                             |
| runtime_s            | Wall-clock seconds                                |
| functional_accuracy  | % fresh trajectories matching ground truth        |
| gene_recovery        | % genes with correct truth table                  |
| degenerate_rate      | % genes with degenerate L-value                   |
| n_solutions          | Number of concrete networks found                 |

### Expected Outputs for Project Book

1. **Runtime curve**: genes (N) vs. runtime (seconds, log scale), points colored
   by SAT/UNSAT/Timeout.
2. **Degenerate rate vs. % unbalanced genes**: scatter plot showing the
   topology-to-failure relationship.
3. A threshold statement: *"RE:IN reliably synthesizes networks up to ~N genes
   within 300 seconds. Beyond this, state-space explosion dominates."*

---

## 📋 Metrics Reference

| Metric               | Formula                                           | Used in    |
|----------------------|---------------------------------------------------|------------|
| SAT Rate             | SAT models / total models run                     | E1, E2     |
| Functional Accuracy  | Matching fresh trajectories / 200                 | E1, E2, E3 |
| Gene Recovery        | Genes with correct truth table / total genes      | E1, E2, E3 |
| Stability Ratio      | Genes with variance=0 / total genes               | E1, E3     |
| Degenerate Rate      | Genes with degenerate L-value / total genes       | E1, E2, E3 |
| Runtime              | Wall-clock seconds from engine start to output    | E2         |
| Solution Count       | Concrete networks returned by SMT solver          | E1, E3     |

---

## 📁 Repository Structure

```
grn-rein-project/
├── README.md
├── GRN_EXPERIMENT_PLAN.md           ← this file
│
├── models/
│   ├── experiment_1_sanity/         ← 5 models (DONE)
│   ├── experiment_2_stress/         ← 20 models (pending)
│   └── experiment_3_observation/    ← 2 reference models × 15 runs (pending)
│
├── pipeline/
│   ├── aeon_to_rein.py
│   ├── rein_add_experiments.py
│   ├── ground_truth_extractor.py
│   ├── evaluate_results.py          ← implements all three metrics
│   └── analyzer.py
│
├── results/
│   ├── experiment_1/                ← DONE
│   ├── experiment_2/
│   └── experiment_3/
│
└── reports/
    ├── experiment_1_evaluation.csv  ← DONE
    ├── experiment_2_summary.csv
    └── experiment_3_heatmap.csv
```

---

_Last updated: March 2026 — v2.0 post E1 findings_
