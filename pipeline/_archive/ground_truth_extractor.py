"""
AEON Boolean Network Simulator

Parses AEON format Boolean network files and simulates their synchronous
dynamics to produce trajectories. Part of the RE:IN inference testing pipeline.
"""

import json
import random
import re
from itertools import product


def parse_aeon(filepath):
    """Parse an AEON file and extract update functions and input nodes.

    Reads the regulatory graph edges and update function definitions from an
    AEON file. Genes that appear in edges but have no update function ($gene: ...)
    line are classified as input nodes.

    Args:
        filepath: Path to the .aeon file.

    Returns:
        A tuple (update_functions, input_nodes) where:
        - update_functions: dict mapping gene name to its Boolean expression string.
          Example: {"v_Dif": "!v_Cactus", "v_Targets": "v_Dif | v_Dorsal"}
        - input_nodes: sorted list of gene names that appear in edges but have
          no update function defined. These preserve their value during simulation.
    """
    update_functions = {}
    graph_genes = set()

    with open(filepath, "r") as f:
        for line in f:
            line = line.strip()
            if not line:
                continue

            # Update function line: $gene_name: expression
            if line.startswith("$"):
                match = re.match(r"^\$(\S+):\s*(.+)$", line)
                if match:
                    gene = match.group(1)
                    expr = match.group(2).strip()
                    update_functions[gene] = expr
                continue

            # Regulatory edge line: source -> target, source -| target, source -? target
            edge_match = re.match(r"^(\S+)\s+(-\||->|-\?)\s+(\S+)$", line)
            if edge_match:
                source = edge_match.group(1)
                target = edge_match.group(3)
                graph_genes.add(source)
                graph_genes.add(target)

    input_nodes = sorted(graph_genes - set(update_functions.keys()))
    return update_functions, input_nodes


def evaluate_state(update_functions, input_nodes, current_state):
    """Compute the next state via synchronous update.

    For input nodes, their value is copied unchanged. For all other genes,
    their Boolean expression is evaluated against the current state to produce
    the next value. All computations use the same current_state (no partial updates).

    Args:
        update_functions: dict mapping gene name to Boolean expression string.
        input_nodes: list of gene names that are inputs (no update function).
        current_state: dict mapping every gene name to 0 or 1.

    Returns:
        A new dict representing the state at t+1.
    """
    next_state = {}

    # Input nodes preserve their value
    for gene in input_nodes:
        next_state[gene] = current_state[gene]

    # Evaluate update functions for all non-input genes
    for gene, expr in update_functions.items():
        py_expr = _aeon_expr_to_python(expr)
        # Use current_state as the namespace so gene names resolve to 0/1 values
        result = eval(py_expr, {"__builtins__": {}}, current_state)
        next_state[gene] = int(bool(result))

    return next_state


def simulate(update_functions, input_nodes, initial_state, K=15):
    """Run synchronous simulation for K steps from an initial state.

    Args:
        update_functions: dict mapping gene name to Boolean expression string.
        input_nodes: list of input gene names.
        initial_state: dict mapping every gene name to 0 or 1.
        K: number of simulation steps (default 15).

    Returns:
        A list of K+1 state dicts: [state_t0, state_t1, ..., state_tK].
    """
    trajectory = [dict(initial_state)]
    current = dict(initial_state)

    for _ in range(K):
        current = evaluate_state(update_functions, input_nodes, current)
        trajectory.append(dict(current))

    return trajectory


def generate_trajectories(aeon_filepath, N=5, K=15):
    """Generate N trajectories from random initial states.

    Parses the AEON file once, then generates N trajectories each starting from
    a different random initial state where each gene is independently assigned
    0 or 1 with equal probability.

    Args:
        aeon_filepath: path to the .aeon file.
        N: number of trajectories to generate (default 5).
        K: number of simulation steps per trajectory (default 15).

    Returns:
        A tuple (trajectories, input_nodes) where:
        - trajectories: list of N trajectories, each a list of K+1 state dicts.
        - input_nodes: list of input gene names.
    """
    update_functions, input_nodes = parse_aeon(aeon_filepath)
    all_genes = sorted(set(list(update_functions.keys()) + input_nodes))

    trajectories = []
    for _ in range(N):
        initial_state = {gene: random.randint(0, 1) for gene in all_genes}
        traj = simulate(update_functions, input_nodes, initial_state, K)
        trajectories.append(traj)

    return trajectories, input_nodes


def expression_to_truth_table(expression, regulators):
    """Build a truth table for a Boolean expression over its regulators.

    Iterates over all 2^len(regulators) combinations of 0/1 values and evaluates
    the expression for each. The resulting truth table can be compared against
    RE:IN's recovered L-values.

    Args:
        expression: Boolean expression string in AEON syntax (e.g. "v_Dif | v_Dorsal").
        regulators: list of gene names that appear in the expression.

    Returns:
        A dict mapping each input combination (as a tuple of sorted (gene, value)
        pairs) to the output value (0 or 1).
    """
    py_expr = _aeon_expr_to_python(expression)
    truth_table = {}

    for values in product([0, 1], repeat=len(regulators)):
        assignment = dict(zip(regulators, values))
        result = eval(py_expr, {"__builtins__": {}}, assignment)
        key = tuple(sorted(assignment.items()))
        truth_table[key] = int(bool(result))

    return truth_table


def parse_aeon_edges(filepath):
    """Parse an AEON file and extract regulatory edges with their signs.

    Args:
        filepath: Path to the .aeon file.

    Returns:
        A list of dicts, each with keys: src, dst, sign.
        sign is one of: "positive" (->), "negative" (-|), "unknown" (-?).
    """
    edges = []
    with open(filepath, "r") as f:
        for line in f:
            line = line.strip()
            if not line or line.startswith("$"):
                continue
            edge_match = re.match(r"^(\S+)\s+(-\||->|-\?)\s+(\S+)$", line)
            if edge_match:
                source = edge_match.group(1)
                edge_type = edge_match.group(2)
                target = edge_match.group(3)
                sign = {"->": "positive", "-|": "negative", "-?": "unknown"}[edge_type]
                edges.append({"src": source, "dst": target, "sign": sign})
    return edges


def extract_regulators_from_expr(expression):
    """Extract gene names that appear in a Boolean expression.

    Finds all identifiers matching the v_NAME pattern used in AEON files.

    Args:
        expression: Boolean expression string.

    Returns:
        Sorted list of unique gene names found in the expression.
    """
    return sorted(set(re.findall(r"\bv_[A-Za-z0-9_]+\b", expression)))


def generate_ground_truth(aeon_filepath, N=5, K=1):
    """Generate complete ground truth data from an AEON model.

    Produces truth tables, regulators, edge signs, and simulated experiments
    purely from the AEON file — no RE:IN involved. The output can be compared
    against RE:IN's inferred results to measure inference accuracy.

    Args:
        aeon_filepath: Path to the .aeon file.
        N: Number of experiments (random trajectories) to generate.
        K: Number of simulation steps per experiment.

    Returns:
        A dict containing:
        - source_file: path to the AEON file
        - total_genes: total number of genes
        - total_inputs: number of input nodes
        - input_nodes: list of input gene names
        - genes: list of per-gene dicts with truth tables, regulators, edges
        - experiments: list of INIT/FIN experiment dicts with full trajectories
    """
    update_functions, input_nodes = parse_aeon(aeon_filepath)
    edges = parse_aeon_edges(aeon_filepath)
    all_genes = sorted(set(list(update_functions.keys()) + input_nodes))

    genes_info = []
    for gene in all_genes:
        is_input = gene in input_nodes

        if is_input:
            genes_info.append({
                "name": gene,
                "is_input": True,
                "update_function": None,
                "regulators": [],
                "truth_table": {},
                "edges_in": [],
            })
            continue

        expr = update_functions[gene]
        regs = extract_regulators_from_expr(expr)
        tt = expression_to_truth_table(expr, regs)

        # JSON-serializable truth table
        tt_serializable = {}
        for combo, output in sorted(tt.items()):
            key_str = ", ".join(f"{name}={val}" for name, val in combo)
            tt_serializable[key_str] = output

        # Edges pointing to this gene (from the AEON regulatory graph)
        gene_edges = [e for e in edges if e["dst"] == gene]

        # Build regulator list with sign info from edges
        reg_info = []
        for r in regs:
            edge_signs = [e["sign"] for e in gene_edges if e["src"] == r]
            reg_info.append({
                "name": r,
                "edge_signs": edge_signs if edge_signs else ["implicit"],
            })

        genes_info.append({
            "name": gene,
            "is_input": False,
            "update_function": expr,
            "regulators": reg_info,
            "num_regulators": len(regs),
            "truth_table": tt_serializable,
            "truth_table_size": len(tt_serializable),
            "edges_in": gene_edges,
        })

    # Generate experiments
    experiments = []
    for i in range(N):
        init = {gene: random.randint(0, 1) for gene in all_genes}
        traj = simulate(update_functions, input_nodes, init, K)
        experiments.append({
            "id": i,
            "K": K,
            "init_state": traj[0],
            "final_state": traj[-1],
            "full_trajectory": traj,
        })

    return {
        "source_file": str(aeon_filepath),
        "total_genes": len(all_genes),
        "total_update_functions": len(update_functions),
        "total_inputs": len(input_nodes),
        "input_nodes": input_nodes,
        "genes": genes_info,
        "experiment_params": {"N": N, "K": K},
        "experiments": experiments,
    }


def _aeon_expr_to_python(expr):
    """Convert an AEON Boolean expression to a Python-evaluable expression.

    Replaces AEON operators: ! -> not , & -> and, | -> or.
    Handles the case where ! may appear with or without a space before the operand.

    Args:
        expr: AEON Boolean expression string.

    Returns:
        Python-evaluable Boolean expression string.
    """
    # Replace ! with not (add space to handle cases like !v_X or !(expr))
    py_expr = expr.replace("!", "not ")
    py_expr = py_expr.replace("&", " and ")
    py_expr = py_expr.replace("|", " or ")
    # Clean up any double spaces
    py_expr = re.sub(r"\s+", " ", py_expr).strip()
    return py_expr


if __name__ == "__main__":
    import os
    import sys

    # --- Configuration ---
    example_path = os.path.join(os.path.dirname(__file__), "aeon", "007-05-0.aeon")
    output_dir = os.path.join(os.path.dirname(__file__), "ground_truth")
    N = 5   # number of experiments
    K = 1   # steps per experiment (K=1 matches rein_k1 benchmarks)

    # Allow command-line override: python3 aeon_simulator.py <aeon_file> [N] [K]
    if len(sys.argv) >= 2:
        example_path = sys.argv[1]
    if len(sys.argv) >= 3:
        N = int(sys.argv[2])
    if len(sys.argv) >= 4:
        K = int(sys.argv[3])

    if not os.path.exists(example_path):
        print(f"File not found: {example_path}")
        sys.exit(1)

    # --- Generate ground truth ---
    model_name = os.path.splitext(os.path.basename(example_path))[0]
    print(f"Generating ground truth for: {model_name}")
    print(f"  Source: {example_path}")
    print(f"  Params: N={N}, K={K}")

    gt = generate_ground_truth(example_path, N=N, K=K)

    # --- Save to JSON ---
    os.makedirs(output_dir, exist_ok=True)
    out_path = os.path.join(output_dir, f"{model_name}_ground_truth.json")
    with open(out_path, "w") as f:
        json.dump(gt, f, indent=2)
    print(f"  Saved: {out_path}")

    # --- Print summary ---
    print(f"\n{'=' * 60}")
    print(f"Ground Truth: {model_name}")
    print(f"{'=' * 60}")
    print(f"  Genes: {gt['total_genes']} ({gt['total_update_functions']} functions, {gt['total_inputs']} inputs)")
    print(f"  Input nodes: {gt['input_nodes']}")

    print(f"\n  Gene truth tables:")
    for g in gt["genes"]:
        if g["is_input"]:
            print(f"    {g['name']}: INPUT (self-loop)")
            continue
        regs_str = ", ".join(r["name"] for r in g["regulators"])
        print(f"    {g['name']}: {g['update_function']}")
        print(f"      regulators: [{regs_str}]  ({g['truth_table_size']} rows)")
        for combo, val in g["truth_table"].items():
            print(f"        {combo} => {val}")

    print(f"\n  Experiments (N={N}, K={K}):")
    all_genes = sorted(set(
        [g["name"] for g in gt["genes"]]
    ))
    for exp in gt["experiments"]:
        init_vals = "".join(str(exp["init_state"][g]) for g in all_genes)
        fin_vals = "".join(str(exp["final_state"][g]) for g in all_genes)
        print(f"    Exp {exp['id']}: {init_vals} -> {fin_vals}")
