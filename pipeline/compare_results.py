#!/usr/bin/env python3
"""
Compare RE:IN L-value output against ground truth truth tables.

For each gene, takes RE:IN's L-value, computes the truth table that L-value
produces using the gene's topology (activators/repressors), then compares
it to the pre-computed truth table in ground_truth.json.

Usage:
  python3 compare_results.py
"""

import json, csv, itertools, os, sys, random, re
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
# Simulation helpers (for functional accuracy)
# ---------------------------------------------------------------------------

def _aeon_to_python(expr: str) -> str:
    py = expr.replace("!", "not ")
    py = py.replace("&", " and ")
    py = py.replace("|", " or ")
    return re.sub(r"\s+", " ", py).strip()


def simulate_gt(gt_json, initial_state, K):
    """Simulate ground-truth model for K steps, return final state."""
    update_functions = {}
    input_nodes = []
    for gene in gt_json['genes']:
        if gene['is_input']:
            input_nodes.append(gene['name'])
        else:
            update_functions[gene['name']] = gene['update_function']

    current = dict(initial_state)
    for _ in range(K):
        next_state = {}
        for g in input_nodes:
            next_state[g] = current[g]
        for g, expr in update_functions.items():
            py_expr = _aeon_to_python(expr)
            result = eval(py_expr, {"__builtins__": {}}, current)
            next_state[g] = int(bool(result))
        current = next_state
    return current


def simulate_rein(gt_json, rein_lv, initial_state, K):
    """Simulate RE:IN synthesized model for K steps using L-values."""
    input_nodes = []
    non_input_genes = []
    for gene in gt_json['genes']:
        if gene['is_input']:
            input_nodes.append(gene['name'])
        else:
            non_input_genes.append(gene)

    current = dict(initial_state)
    for _ in range(K):
        next_state = {}
        for g in input_nodes:
            next_state[g] = current[g]
        for gene in non_input_genes:
            name = gene['name']
            L = rein_lv.get(name)
            if L is None:
                # No solution for this gene — keep current value
                next_state[name] = current[name]
                continue
            activators = sorted([r['name'] for r in gene['regulators']
                                 if r['edge_signs'][0] == 'positive'])
            repressors = sorted([r['name'] for r in gene['regulators']
                                 if r['edge_signs'][0] == 'negative'])
            state_for_gene = {r['name']: current[r['name']] for r in gene['regulators']}
            next_state[name] = int(evaluate_R(L, state_for_gene, activators, repressors))
        current = next_state
    return current


def functional_accuracy(gt_json, rein_lv, K, n_samples=200, seed=42):
    """Simulate both models from fresh random initial states, compare final states."""
    rng = random.Random(seed)
    all_genes = sorted([g['name'] for g in gt_json['genes']])
    matches = 0
    for _ in range(n_samples):
        state = {g: rng.randint(0, 1) for g in all_genes}
        gt_final = simulate_gt(gt_json, state, K)
        rein_final = simulate_rein(gt_json, rein_lv, state, K)
        if gt_final == rein_final:
            matches += 1
    return matches / n_samples


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
    parser.add_argument('--K', type=int, default=15,
                        help='Trajectory length for functional accuracy (default: 15)')
    args = parser.parse_args()

    if args.results_dir:
        base = Path(args.results_dir)
    else:
        base = Path(__file__).resolve().parent.parent / 'results' / 'experiment_1'

    # Auto-discover model subdirectories
    models = sorted([d.name for d in base.iterdir() if d.is_dir()])

    rows = []
    model_func_acc = {}  # model_id -> functional accuracy

    for model_id in models:
        model_dir = base / model_id
        gt_path = model_dir / f'{model_id}_ground_truth.json'

        # Try multiple CSV naming conventions
        csv_path = model_dir / f'{model_id}_l_values_extracted.csv'
        if not csv_path.exists():
            # E3-style naming: model_id_K{k}_N{n}_l_values.csv
            candidates = list(model_dir.glob(f'{model_id}_*_l_values.csv'))
            if candidates:
                csv_path = candidates[0]

        if not gt_path.exists():
            print(f"WARNING: ground truth not found: {gt_path}", file=sys.stderr)
            continue
        if not csv_path.exists():
            print(f"WARNING: L-values CSV not found for {model_id}", file=sys.stderr)
            continue

        gt = json.load(open(gt_path))

        # Parse RE:IN L-values CSV — handle both formats
        rein_lv = {}
        with open(csv_path) as f:
            reader = csv.DictReader(f)
            for row in reader:
                gene_key = row.get('gene', row.get('Gene', ''))
                lv_key = row.get('L_value', row.get('l_value', ''))
                if not gene_key or not lv_key:
                    # Fallback: positional 4-column format
                    break
                rein_lv[gene_key] = int(lv_key)

        if not rein_lv:
            # Fallback: 4-column positional CSV (model_id, solution_index, gene, L_value)
            for line in open(csv_path).read().strip().split('\n')[1:]:
                if line.strip():
                    parts = line.split(',')
                    if len(parts) == 4:
                        rein_lv[parts[2]] = int(parts[3])

        # Compute functional accuracy for this model
        has_any_lvalues = bool(rein_lv)
        if has_any_lvalues:
            func_acc = functional_accuracy(gt, rein_lv, K=args.K)
            model_func_acc[model_id] = func_acc
            print(f"{model_id}: functional_accuracy = {func_acc:.1%} (K={args.K}, 200 samples)")
        else:
            model_func_acc[model_id] = None
            print(f"{model_id}: no L-values (UNSAT/Timeout)")

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
                    'functional_accuracy': model_func_acc[model_id],
                })
                continue

            rein_tt = compute_rein_truth_table(rein_L, gene)
            match = (gt_tt == rein_tt)
            rows.append({
                'model': model_id, 'gene': name,
                'rein_L': rein_L, 'match': match,
                'gt_tt': gt_tt, 'rein_tt': rein_tt,
                'functional_accuracy': model_func_acc[model_id],
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

    # Print functional accuracy summary
    print(f"\nFunctional Accuracy (K={args.K}, 200 fresh states):")
    for mid, acc in model_func_acc.items():
        if acc is not None:
            print(f"  {mid}: {acc:.1%}")
        else:
            print(f"  {mid}: N/A (no solution)")

    # Save CSV
    out_path = base / f'{base.name}_comparison.csv'
    with open(out_path, 'w', newline='') as f:
        w = csv.DictWriter(f, fieldnames=['model', 'gene', 'rein_L', 'match',
                                          'gt_tt', 'rein_tt', 'functional_accuracy'])
        w.writeheader()
        w.writerows(rows)
    print(f"Saved {out_path}")


if __name__ == '__main__':
    main()
