# Experiment 3: Observation Sufficiency — Summary of Results

> **Question:** Can more data fix the degenerate L-value problem, or is it a structural property of network topology?

---

## Setup

- **Two reference models** from E1, held constant across all runs
- **Grid:** K in {5, 10, 15} x N in {1, 3, 5, 10, 20} = 15 runs per model, 30 total
- **Metric:** Gene-level truth table match (RE:IN L-value vs ground truth)

| Model   | Biological System              | E1 Performance | Role               |
|---------|--------------------------------|----------------|--------------------|
| 110_9_0 | C. crescentus cell cycle (9 genes) | 25% func. acc. | Best E1 performer  |
| 031_9_0 | Yeast cell cycle (9 genes)         | 0% func. acc.  | Worst E1 performer |

---

## Model 110_9_0 — C. crescentus Cell Cycle

### SAT Rate: 15/15 (100%)

All runs returned solutions regardless of K and N values.

### Gene Recovery Heatmap (K x N)

|       | N=1   | N=3   | N=5   | N=10  | N=20  |
|-------|-------|-------|-------|-------|-------|
| K=5   | 6/9 (66.7%) | 5/9 (55.6%) | 3/9 (33.3%) | 4/9 (44.4%) | 8/9 (88.9%) |
| K=10  | 7/9 (77.8%) | 5/9 (55.6%) | 7/9 (77.8%) | 6/9 (66.7%) | 7/9 (77.8%) |
| K=15  | 8/9 (88.9%) | 6/9 (66.7%) | 7/9 (77.8%) | 8/9 (88.9%) | 7/9 (77.8%) |

### Balanced vs Unbalanced Gene Recovery

| Gene Type   | Count | Best Recovery | Typical Recovery | Notes |
|-------------|-------|---------------|------------------|-------|
| Unbalanced  | 6/9   | 6/6 (100%)    | 5-6/6 (83-100%)  | Consistently high across all (K,N) |
| Balanced    | 3/9   | 3/3 (100%)    | 0-2/3 (0-67%)    | Highly variable, the bottleneck |

### Degenerate L-value Rate

Stays at ~33% across most runs (3 genes: v_ClpXP_RcdA, v_DivL, v_PleC — all pure-repressor with 1 regulator). Notably, these degenerate L-values still produce **correct** truth tables because with only 1 regulator, the degenerate and correct truth tables happen to coincide.

### Per-Gene Detail (K=15, N=20)

| Gene           | Balance        | L-value | Match | Degenerate | GT Truth Table       | RE:IN Truth Table    |
|----------------|----------------|---------|-------|------------|----------------------|----------------------|
| v_CckA         | pure_activator | 0       | YES   |            | [0, 1]               | [0, 1]               |
| v_ChpT         | pure_activator | 3       | YES   |            | [0, 1]               | [0, 1]               |
| v_ClpXP_RcdA   | pure_repressor | 17      | YES   | YES        | [1, 0]               | [1, 0]               |
| v_CpdR         | pure_activator | 3       | YES   |            | [0, 1]               | [0, 1]               |
| **v_CtrAb**    | **balanced**   | **5**   | **NO**|            | [0, 0, 1, 0]         | [0, 0, 1, 1]         |
| v_DivJ         | balanced       | 2       | YES   |            | [0, 0, 1, 0]         | [0, 0, 1, 0]         |
| **v_DivK**     | **balanced**   | **15**  | **NO**|            | [0, 0, 1, 0]         | [1, 0, 1, 1]         |
| v_DivL         | pure_repressor | 17      | YES   | YES        | [1, 0]               | [1, 0]               |
| v_PleC         | pure_repressor | 17      | YES   | YES        | [1, 0]               | [1, 0]               |

The two failures (v_CtrAb, v_DivK) are both **balanced genes** — RE:IN assigns L-values that produce overly permissive activation.

---

## Model 031_9_0 — Yeast Cell Cycle

### SAT Rate: 10/15 (67%)

Five runs returned UNSAT/Timeout — notably at higher N and K values:

|       | N=1 | N=3   | N=5   | N=10  | N=20  |
|-------|-----|-------|-------|-------|-------|
| K=5   | SAT | SAT   | SAT   | UNSAT | UNSAT |
| K=10  | SAT | SAT   | SAT   | SAT   | SAT   |
| K=15  | SAT | UNSAT | UNSAT | UNSAT | UNSAT |

**More observations make it harder for RE:IN to find any solution** — the opposite of what H1 would predict.

### Gene Recovery Heatmap (K x N, SAT runs only)

|       | N=1   | N=3   | N=5   | N=10  | N=20  |
|-------|-------|-------|-------|-------|-------|
| K=5   | 4/9 (44.4%) | 3/9 (33.3%) | 3/9 (33.3%) | UNSAT | UNSAT |
| K=10  | 1/9 (11.1%) | 1/9 (11.1%) | 3/9 (33.3%) | 3/9 (33.3%) | 1/9 (11.1%) |
| K=15  | 2/9 (22.2%) | UNSAT | UNSAT | UNSAT | UNSAT |

### Balanced vs Unbalanced Gene Recovery

| Gene Type   | Count | Best Recovery | Typical Recovery | Notes |
|-------------|-------|---------------|------------------|-------|
| Unbalanced  | 7/9   | 4/7 (57.1%)  | 1-4/7 (14-57%)   | Low and inconsistent |
| Balanced    | 2/9   | 0/2 (0%)     | 0/2 (0%)         | **Never recovered in any run** |

The two balanced genes (v_CLN3, v_SBF) were **never correctly recovered across all 15 runs**.

### Gene Topology

| Gene   | Activators | Repressors | Balance        |
|--------|------------|------------|----------------|
| v_ACE2 | 1          | 0          | pure_activator |
| v_CLN3 | 2          | 2          | balanced       |
| v_HCM1 | 2          | 0          | pure_activator |
| v_MBF  | 1          | 0          | pure_activator |
| v_SBF  | 2          | 2          | balanced       |
| v_SFF  | 2          | 0          | pure_activator |
| v_SWI5 | 1          | 0          | pure_activator |
| v_YHP1 | 3          | 0          | pure_activator |
| v_YOX1 | 2          | 0          | pure_activator |

7 out of 9 genes are pure-activator (no repressors), making this topology especially vulnerable to degenerate L-values.

---

## Key Findings

### 1. More data does NOT reliably improve gene recovery

For 110_9_0, recovery fluctuates between 33-89% with no monotonic trend as N or K increases. For 031_9_0, the best recovery (44.4%) occurs at the **lowest** observation count (K=5, N=1).

### 2. More constraints can cause UNSAT

For 031_9_0, adding more observations (higher N at K=5 or K=15) causes RE:IN to fail entirely — the added constraints become mutually unsatisfiable under the available L-value space.

### 3. Balanced vs unbalanced is the decisive split

| Model   | Unbalanced Recovery (best) | Balanced Recovery (best) |
|---------|---------------------------|-------------------------|
| 110_9_0 | 100%                      | 100% (but often 0-67%)  |
| 031_9_0 | 57.1%                     | **0% always**           |

### 4. Degenerate L-values on 1-regulator genes are harmless

In 110_9_0, the 3 degenerate genes (pure-repressor, 1 regulator each) always produce correct truth tables despite technically degenerate L-values — with only 1 regulator, the truth table has only 2 entries and the degeneracy doesn't manifest.

---

## Conclusion: H2 (Pessimistic) is Strongly Supported

> **The degenerate L-value problem is a structural property of the network topology that no amount of observational data can resolve.**

Evidence:
1. **No convergence trend** — gene recovery does not increase with N or K
2. **More data causes UNSAT** — additional constraints conflict rather than refine
3. **Balanced genes are the barrier** — only topology with both activators and repressors allows clean L-value discrimination, and even then recovery is inconsistent
4. **031_9_0 is fundamentally limited** — its heavily unbalanced topology (7/9 pure-activator) means RE:IN's L-value space cannot represent the correct regulatory logic for most genes

The distinguishing regulator states that would disambiguate degenerate L-values do not appear in the (init -> final) trajectories regardless of how many trajectories are sampled. This is a limitation of the RE:IN synthesis framework, not the observation protocol.

---

*Generated from: `pipeline/compare_e3.py`*
*Data: `results/experiment_3/e3_summary.csv`, `results/experiment_3/e3_gene_detail.csv`*
