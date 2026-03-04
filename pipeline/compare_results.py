#!/usr/bin/env python3
"""
Compare RE:IN L-value output against ground truth truth tables.

For each gene, takes RE:IN's L-value, computes the truth table that L-value
produces using the gene's topology (activators/repressors), then compares
it to the pre-computed truth table in ground_truth.json.

Usage:
  python3 compare_results.py
"""

import json, csv, itertools, os, sys
from pathlib import Path


# ---------------------------------------------------------------------------
# RE:IN regulatory condition predicates
# ---------------------------------------------------------------------------

def AllActivators(state, activators):
    return len(activators) > 0 and all(state[a] == 1 for a in activators)

def NoActivators(state, activators):
    return all(state[a] == 0 for a in activators)

def AllRepressors(state, repressors):
    return len(repressors) > 0 and all(state[r] == 1 for r in repressors)

def NoRepressors(state, repressors):
    return all(state[r] == 0 for r in repressors)


def evaluate_R_raw(i, state, activators, repressors):
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
    # R18, R19 use counts
    count_A = sum(1 for a in activators if state[a] == 1)
    count_R = sum(1 for r in repressors if state[r] == 1)
    cur = state.get('__self__', 0)
    if i == 18:
        return count_A > count_R or (count_A == count_R and cur == 1)
    if i == 19:
        return count_A > count_R
    return False


def evaluate_R(i, state, activators, repressors):
    # No wrapper — just raw R conditions.
    # The wrapper exists in RE:IN's synthesis process to prevent degenerate
    # UNSAT results, but should not be applied when evaluating whether a
    # given L-value matches a truth table.
    return evaluate_R_raw(i, state, activators, repressors)


# ---------------------------------------------------------------------------
# Truth table computation from L-value + gene topology
# ---------------------------------------------------------------------------

def compute_rein_truth_table(rein_L, gene_data):
    """Given an L-value and gene topology, compute the truth table RE:IN implies."""
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


# ---------------------------------------------------------------------------
# Main comparison
# ---------------------------------------------------------------------------

def main():
    import argparse
    parser = argparse.ArgumentParser(description='Compare RE:IN L-values against ground truth')
    parser.add_argument('results_dir', nargs='?', default=None,
                        help='Results directory (default: results/experiment_1)')
    args = parser.parse_args()

    if args.results_dir:
        base = Path(args.results_dir)
    else:
        base = Path(__file__).resolve().parent.parent / 'results' / 'experiment_1'

    # Auto-discover model subdirectories
    models = sorted([d.name for d in base.iterdir() if d.is_dir()])

    rows = []
    for model_id in models:
        model_dir = base / model_id
        gt_path = model_dir / f'{model_id}_ground_truth.json'
        csv_path = model_dir / f'{model_id}_l_values_extracted.csv'

        if not gt_path.exists():
            print(f"WARNING: ground truth not found: {gt_path}", file=sys.stderr)
            continue
        if not csv_path.exists():
            print(f"WARNING: L-values CSV not found: {csv_path}", file=sys.stderr)
            continue

        gt = json.load(open(gt_path))

        # Parse RE:IN L-values CSV (4 columns: model_id, solution_index, gene, L_value)
        rein_lv = {}
        for line in open(csv_path).read().strip().split('\n')[1:]:
            if line.strip():
                parts = line.split(',')
                if len(parts) == 4:
                    rein_lv[parts[2]] = int(parts[3])

        for gene in gt['genes']:
            name = gene['name']
            if gene['is_input']:
                continue

            gt_tt = list(gene['truth_table'].values())
            rein_L = rein_lv.get(name)

            if rein_L is None:
                rows.append({
                    'model': model_id, 'gene': name,
                    'rein_L': 'MISSING', 'match': False,
                    'gt_tt': gt_tt, 'rein_tt': [],
                })
                continue

            rein_tt = compute_rein_truth_table(rein_L, gene)
            match = (gt_tt == rein_tt)
            rows.append({
                'model': model_id, 'gene': name,
                'rein_L': rein_L, 'match': match,
                'gt_tt': gt_tt, 'rein_tt': rein_tt,
            })

    # Print results
    print(f"\n{'Model':<12} {'Gene':<25} {'L':>4} {'Match':>6}  GT_TT  REIN_TT")
    print('-' * 90)
    for r in rows:
        match_str = 'YES' if r['match'] else 'NO'
        print(f"{r['model']:<12} {r['gene']:<25} {str(r['rein_L']):>4} {match_str:>6}  "
              f"{r['gt_tt']}  {r['rein_tt']}")

    matched = sum(1 for r in rows if r['match'])
    total = len(rows)
    print(f"\nRecovery ratio: {matched}/{total} = {matched/total:.1%}" if total > 0 else "\nNo genes compared.")

    # Save CSV
    out_path = base / f'{base.name}_comparison.csv'
    with open(out_path, 'w', newline='') as f:
        w = csv.DictWriter(f, fieldnames=['model', 'gene', 'rein_L', 'match', 'gt_tt', 'rein_tt'])
        w.writeheader()
        w.writerows(rows)
    print(f"Saved {out_path}")


if __name__ == '__main__':
    main()
