#!/usr/bin/env python3
"""
Experiment 3 — Compare RE:IN L-value outputs against ground truth truth tables.

For each (model, K, N) run:
  - Load ground truth JSON (gene topology + truth tables)
  - Load L-values CSV from RE:IN
  - Compute RE:IN's implied truth table per gene
  - Compare to ground truth
  - Report: SAT status, gene recovery, degenerate L-value detection

Outputs:
  - Per-run summary CSV  (e3_summary.csv)
  - Per-gene detail CSV  (e3_gene_detail.csv)
  - Console report
"""

import json, csv, itertools, os, sys
from pathlib import Path


# ── RE:IN regulatory condition predicates (from compare_results.py) ──

def AllActivators(state, activators):
    return len(activators) > 0 and all(state[a] == 1 for a in activators)

def NoActivators(state, activators):
    return all(state[a] == 0 for a in activators)

def AllRepressors(state, repressors):
    return len(repressors) > 0 and all(state[r] == 1 for r in repressors)

def NoRepressors(state, repressors):
    return all(state[r] == 0 for r in repressors)


def evaluate_R(i, state, activators, repressors):
    AA = AllActivators(state, activators)
    NA = NoActivators(state, activators)
    AR = AllRepressors(state, repressors)
    NR = NoRepressors(state, repressors)
    rules = [
        AA and NR,                                      # R0
        (not NA) and NR,                                # R1
        AA and (not AR),                                # R2
        ((not NA) and NR) or (AA and (not AR)),         # R3
        AA,                                             # R4
        AA or (NR and (not NA)),                        # R5
        (not NA) and (not AR),                          # R6
        AA or ((not NA) and (not AR)),                  # R7
        not NA,                                         # R8
        NR,                                             # R9
        NR or (AA and (not AR)),                        # R10
        NR or ((not NA) and (not AR)),                  # R11
        not AR,                                         # R12
        NR or AA,                                       # R13
        (NR or AA) or ((not NA) and (not AR)),          # R14
        AA or (not AR),                                 # R15
        (not NA) or NR,                                 # R16
        (not NA) or (not AR),                           # R17
    ]
    if i <= 17:
        return rules[i]
    count_A = sum(1 for a in activators if state[a] == 1)
    count_R = sum(1 for r in repressors if state[r] == 1)
    cur = state.get('__self__', 0)
    if i == 18:
        return count_A > count_R or (count_A == count_R and cur == 1)
    if i == 19:
        return count_A > count_R
    return False


def compute_rein_truth_table(rein_L, gene_data):
    activators = sorted([r['name'] for r in gene_data['regulators']
                         if r['edge_signs'][0] == 'positive'])
    repressors = sorted([r['name'] for r in gene_data['regulators']
                         if r['edge_signs'][0] == 'negative'])
    reg_names = sorted([r['name'] for r in gene_data['regulators']])
    result = []
    for combo in itertools.product([0, 1], repeat=len(reg_names)):
        state = dict(zip(reg_names, combo))
        result.append(int(evaluate_R(rein_L, state, activators, repressors)))
    return result


def classify_gene_balance(gene_data):
    """Return 'balanced', 'pure_activator', or 'pure_repressor'."""
    activators = [r for r in gene_data['regulators'] if r['edge_signs'][0] == 'positive']
    repressors = [r for r in gene_data['regulators'] if r['edge_signs'][0] == 'negative']
    if activators and repressors:
        return 'balanced'
    elif activators and not repressors:
        return 'pure_activator'
    elif repressors and not activators:
        return 'pure_repressor'
    return 'none'


def is_degenerate(rein_L, gene_data):
    """Check if an L-value is degenerate for this gene's topology."""
    activators = [r['name'] for r in gene_data['regulators'] if r['edge_signs'][0] == 'positive']
    repressors = [r['name'] for r in gene_data['regulators'] if r['edge_signs'][0] == 'negative']

    # Degenerate conditions from the experiment plan:
    # Pure-repressor gene (no activators): NoActivators always True, ¬NoActivators always False
    # Pure-activator gene (no repressors): NoRepressors always True
    if not activators and not repressors:
        return False

    # R-conditions that contain NoRepressors (vacuously True when no repressors)
    nr_conditions = {9, 10, 11, 13, 16}
    # R-conditions that contain ¬NoActivators (vacuously False when no activators)
    not_na_conditions = {1, 3, 5, 6, 7, 8, 16, 17}

    if not repressors and rein_L in nr_conditions:
        return True
    if not activators and rein_L in not_na_conditions:
        return True
    return False


def main():
    base = Path(__file__).resolve().parent.parent / 'results' / 'experiment_3'

    gene_detail_rows = []
    summary_rows = []

    for model_id in ['110_9_0', '031_9_0']:
        model_dir = base / model_id
        if not model_dir.exists():
            print(f"WARNING: {model_dir} not found", file=sys.stderr)
            continue

        for kn_dir in sorted(model_dir.iterdir()):
            if not kn_dir.is_dir():
                continue
            variant = kn_dir.name  # e.g. K15_N20
            parts = variant.split('_')
            K = int(parts[0][1:])
            N = int(parts[1][1:])

            inner_dir = kn_dir / model_id
            gt_path = inner_dir / f'{model_id}_ground_truth.json'
            csv_path = inner_dir / f'{model_id}_{variant}_l_values.csv'

            if not gt_path.exists():
                print(f"WARNING: ground truth not found: {gt_path}", file=sys.stderr)
                continue
            if not csv_path.exists():
                print(f"WARNING: L-values CSV not found: {csv_path}", file=sys.stderr)
                continue

            gt = json.load(open(gt_path))

            # Parse L-values
            rein_lv = {}
            with open(csv_path) as f:
                reader = csv.DictReader(f)
                for row in reader:
                    rein_lv[row['gene']] = int(row['L_value'])

            # Determine SAT status
            sat_status = 'SAT' if rein_lv else 'UNSAT/Timeout'

            genes_total = 0
            genes_matched = 0
            genes_degenerate = 0
            balanced_total = 0
            balanced_matched = 0
            unbalanced_total = 0
            unbalanced_matched = 0

            for gene in gt['genes']:
                if gene['is_input']:
                    continue
                genes_total += 1
                name = gene['name']
                gt_tt = list(gene['truth_table'].values())
                balance = classify_gene_balance(gene)

                rein_L = rein_lv.get(name)
                if rein_L is None:
                    gene_detail_rows.append({
                        'model_id': model_id, 'K': K, 'N': N,
                        'gene': name, 'balance': balance,
                        'rein_L': '', 'match': False,
                        'degenerate': False,
                        'gt_tt': str(gt_tt), 'rein_tt': '',
                    })
                    if balance == 'balanced':
                        balanced_total += 1
                    else:
                        unbalanced_total += 1
                    continue

                rein_tt = compute_rein_truth_table(rein_L, gene)
                match = (gt_tt == rein_tt)
                degen = is_degenerate(rein_L, gene)

                if match:
                    genes_matched += 1
                if degen:
                    genes_degenerate += 1

                if balance == 'balanced':
                    balanced_total += 1
                    if match:
                        balanced_matched += 1
                else:
                    unbalanced_total += 1
                    if match:
                        unbalanced_matched += 1

                gene_detail_rows.append({
                    'model_id': model_id, 'K': K, 'N': N,
                    'gene': name, 'balance': balance,
                    'rein_L': rein_L, 'match': match,
                    'degenerate': degen,
                    'gt_tt': str(gt_tt), 'rein_tt': str(rein_tt),
                })

            gene_recovery = genes_matched / genes_total if genes_total > 0 else 0
            degenerate_rate = genes_degenerate / genes_total if genes_total > 0 else 0
            balanced_recovery = balanced_matched / balanced_total if balanced_total > 0 else 0
            unbalanced_recovery = unbalanced_matched / unbalanced_total if unbalanced_total > 0 else 0

            summary_rows.append({
                'model_id': model_id,
                'K': K,
                'N': N,
                'sat_status': sat_status,
                'gene_recovery': f'{gene_recovery:.1%}',
                'gene_recovery_raw': f'{genes_matched}/{genes_total}',
                'degenerate_rate': f'{degenerate_rate:.1%}',
                'balanced_recovery': f'{balanced_recovery:.1%}' if balanced_total > 0 else 'N/A',
                'unbalanced_recovery': f'{unbalanced_recovery:.1%}' if unbalanced_total > 0 else 'N/A',
                'balanced_detail': f'{balanced_matched}/{balanced_total}',
                'unbalanced_detail': f'{unbalanced_matched}/{unbalanced_total}',
            })

    # Sort by model, K, N
    summary_rows.sort(key=lambda r: (r['model_id'], r['K'], r['N']))
    gene_detail_rows.sort(key=lambda r: (r['model_id'], r['K'], r['N'], r['gene']))

    # ── Console output ──
    print("\n" + "=" * 100)
    print("EXPERIMENT 3 — Observation Sufficiency: Ground Truth vs RE:IN L-values")
    print("=" * 100)

    for model_id in ['110_9_0', '031_9_0']:
        model_rows = [r for r in summary_rows if r['model_id'] == model_id]
        print(f"\n{'─' * 80}")
        print(f"Model: {model_id}")
        print(f"{'─' * 80}")
        print(f"{'K':>4} {'N':>4}  {'SAT':<14} {'Gene Recovery':<16} {'Degen Rate':<12} {'Balanced':<12} {'Unbalanced':<12}")
        print(f"{'─' * 4} {'─' * 4}  {'─' * 14} {'─' * 16} {'─' * 12} {'─' * 12} {'─' * 12}")
        for r in model_rows:
            print(f"{r['K']:>4} {r['N']:>4}  {r['sat_status']:<14} "
                  f"{r['gene_recovery_raw']:>5} ({r['gene_recovery']:>5})  "
                  f"{r['degenerate_rate']:>10}  "
                  f"{r['balanced_detail']:>5} ({r['balanced_recovery']:>5})  "
                  f"{r['unbalanced_detail']:>5} ({r['unbalanced_recovery']:>5})")

    # ── Heatmap view ──
    print(f"\n{'=' * 60}")
    print("HEATMAP: Gene Recovery % (K rows × N columns)")
    print(f"{'=' * 60}")
    N_vals = [1, 3, 5, 10, 20]
    K_vals = [5, 10, 15]
    for model_id in ['110_9_0', '031_9_0']:
        print(f"\nModel: {model_id}")
        kn_header = 'K\\N'
        print(f"{kn_header:>6}", end='')
        for n in N_vals:
            print(f"{'N=' + str(n):>10}", end='')
        print()
        for k in K_vals:
            print(f"{'K=' + str(k):>6}", end='')
            for n in N_vals:
                match = [r for r in summary_rows
                         if r['model_id'] == model_id and r['K'] == k and r['N'] == n]
                if match:
                    val = match[0]['gene_recovery']
                    sat = match[0]['sat_status']
                    cell = val if sat == 'SAT' else 'UNSAT'
                else:
                    cell = '—'
                print(f"{cell:>10}", end='')
            print()

    # ── Per-gene detail for key runs ──
    print(f"\n{'=' * 100}")
    print("PER-GENE DETAIL (K=15, N=20 — maximum observation runs)")
    print(f"{'=' * 100}")
    for model_id in ['110_9_0', '031_9_0']:
        print(f"\nModel: {model_id}")
        detail = [r for r in gene_detail_rows
                  if r['model_id'] == model_id and r['K'] == 15 and r['N'] == 20]
        print(f"  {'Gene':<20} {'Balance':<16} {'L':>4} {'Match':>6} {'Degen':>6}  GT_TT  →  REIN_TT")
        print(f"  {'─' * 90}")
        for d in detail:
            m = 'YES' if d['match'] else 'NO'
            dg = 'YES' if d['degenerate'] else ''
            print(f"  {d['gene']:<20} {d['balance']:<16} {str(d['rein_L']):>4} {m:>6} {dg:>6}  "
                  f"{d['gt_tt']}  →  {d['rein_tt']}")

    # ── Save CSVs ──
    summary_path = base / 'e3_summary.csv'
    with open(summary_path, 'w', newline='') as f:
        w = csv.DictWriter(f, fieldnames=[
            'model_id', 'K', 'N', 'sat_status', 'gene_recovery',
            'gene_recovery_raw', 'degenerate_rate',
            'balanced_recovery', 'unbalanced_recovery',
            'balanced_detail', 'unbalanced_detail',
        ])
        w.writeheader()
        w.writerows(summary_rows)
    print(f"\nSaved: {summary_path}")

    detail_path = base / 'e3_gene_detail.csv'
    with open(detail_path, 'w', newline='') as f:
        w = csv.DictWriter(f, fieldnames=[
            'model_id', 'K', 'N', 'gene', 'balance',
            'rein_L', 'match', 'degenerate', 'gt_tt', 'rein_tt',
        ])
        w.writeheader()
        w.writerows(gene_detail_rows)
    print(f"Saved: {detail_path}")


if __name__ == '__main__':
    main()
