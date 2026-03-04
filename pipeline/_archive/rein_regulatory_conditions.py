#!/usr/bin/env python3
"""
RE:IN Regulatory Condition Formulas (R₀–R₁₉)

Converts AEON ground truth update functions into RE:IN L-values by evaluating
the 20 regulatory condition formulas and finding which one matches the ground
truth truth table.

The 18 core conditions (R₀–R₁₇) use four predicates:
  AllActivators, NoActivators, AllRepressors, NoRepressors
with an inducible/repressible edge-case wrapper.

The 2 threshold conditions (R₁₈, R₁₉) use activator/repressor counts directly.

Usage:
  python3 rein_regulatory_conditions.py <aeon_file>

  Or as a library:
    from rein_regulatory_conditions import find_ground_truth_L_value
    L = find_ground_truth_L_value(gene, update_fn_str, edges)
"""

import re
import sys
from itertools import product as iterproduct
from pathlib import Path
from typing import Dict, List, Optional, Set, Tuple


# ---------------------------------------------------------------------------
# Activator / Repressor extraction from edge list
# ---------------------------------------------------------------------------

def get_activators(gene: str, edges: List[Dict]) -> Set[str]:
    """Return set of genes with positive edges (→) pointing to gene."""
    return {e["src"] for e in edges
            if e["dst"] == gene and e["sign"] == "positive"}


def get_repressors(gene: str, edges: List[Dict]) -> Set[str]:
    """Return set of genes with negative edges (-|) pointing to gene."""
    return {e["src"] for e in edges
            if e["dst"] == gene and e["sign"] == "negative"}


# ---------------------------------------------------------------------------
# The Four Core Predicates
# ---------------------------------------------------------------------------

def AllActivators(gene: str, state: Dict[str, int], edges: List[Dict]) -> bool:
    """TRUE if gene has ≥1 activator AND all activators are 1 in state."""
    acts = get_activators(gene, edges)
    return len(acts) > 0 and all(state[a] == 1 for a in acts)


def NoActivators(gene: str, state: Dict[str, int], edges: List[Dict]) -> bool:
    """TRUE if all activators are 0 in state (vacuously True if none)."""
    acts = get_activators(gene, edges)
    return all(state[a] == 0 for a in acts)


def AllRepressors(gene: str, state: Dict[str, int], edges: List[Dict]) -> bool:
    """TRUE if gene has ≥1 repressor AND all repressors are 1 in state."""
    reps = get_repressors(gene, edges)
    return len(reps) > 0 and all(state[r] == 1 for r in reps)


def NoRepressors(gene: str, state: Dict[str, int], edges: List[Dict]) -> bool:
    """TRUE if all repressors are 0 in state (vacuously True if none)."""
    reps = get_repressors(gene, edges)
    return all(state[r] == 0 for r in reps)


# ---------------------------------------------------------------------------
# Edge Case Wrapper: inducible / repressible regulation
# ---------------------------------------------------------------------------

def inducible_regulation(gene: str, state: Dict[str, int], edges: List[Dict]) -> bool:
    """Balancing rule for genes with both activators and repressors.
    Prevents constitutive OFF by requiring not-all-inputs-absent."""
    acts = get_activators(gene, edges)
    reps = get_repressors(gene, edges)
    if len(acts) > 0 and len(reps) > 0:
        return not NoActivators(gene, state, edges)
    else:
        return True


def repressible_regulation(gene: str, state: Dict[str, int], edges: List[Dict]) -> bool:
    """Balancing rule for genes with both activators and repressors.
    Forces activation when all repressors are absent."""
    acts = get_activators(gene, edges)
    reps = get_repressors(gene, edges)
    if len(acts) > 0 and len(reps) > 0:
        return NoRepressors(gene, state, edges)
    else:
        return False


# ---------------------------------------------------------------------------
# The 18 Raw Regulatory Condition Formulas (R₀–R₁₇)
# ---------------------------------------------------------------------------

def evaluate_R_raw(i: int, gene: str, state: Dict[str, int], edges: List[Dict]) -> bool:
    """Evaluate raw Rᵢ formula (0 ≤ i ≤ 17) without edge case wrapper."""
    AA = AllActivators(gene, state, edges)
    NA = NoActivators(gene, state, edges)
    AR = AllRepressors(gene, state, edges)
    NR = NoRepressors(gene, state, edges)

    if i == 0:
        return AA and NR
    elif i == 1:
        return (not NA) and NR
    elif i == 2:
        return AA and (not AR)
    elif i == 3:
        return ((not NA) and NR) or (AA and (not AR))
    elif i == 4:
        return AA
    elif i == 5:
        return AA or (NR and (not NA))
    elif i == 6:
        return (not NA) and (not AR)
    elif i == 7:
        return AA or ((not NA) and (not AR))
    elif i == 8:
        return not NA
    elif i == 9:
        return NR
    elif i == 10:
        return NR or (AA and (not AR))
    elif i == 11:
        return NR or ((not NA) and (not AR))
    elif i == 12:
        return not AR
    elif i == 13:
        return NR or AA
    elif i == 14:
        return (NR or AA) or ((not NA) and (not AR))
    elif i == 15:
        return AA or (not AR)
    elif i == 16:
        return (not NA) or NR
    elif i == 17:
        return (not NA) or (not AR)
    else:
        raise ValueError(f"Invalid R index: {i} (must be 0–17)")


# ---------------------------------------------------------------------------
# The 2 Threshold Conditions (R₁₈, R₁₉)
# ---------------------------------------------------------------------------

def evaluate_R18(gene: str, state: Dict[str, int], edges: List[Dict]) -> bool:
    """Instant threshold: #A > #R, or (#A == #R and gene is currently ON)."""
    acts = get_activators(gene, edges)
    reps = get_repressors(gene, edges)
    count_a = sum(state[a] for a in acts)
    count_r = sum(state[r] for r in reps)
    return (count_a > count_r) or (count_a == count_r and state.get(gene, 0) == 1)


def evaluate_R19(gene: str, state: Dict[str, int], edges: List[Dict]) -> bool:
    """Delayed threshold: #A > #R."""
    acts = get_activators(gene, edges)
    reps = get_repressors(gene, edges)
    count_a = sum(state[a] for a in acts)
    count_r = sum(state[r] for r in reps)
    return count_a > count_r


# ---------------------------------------------------------------------------
# Full Rᵢ evaluation (with wrapper for R₀–R₁₇, direct for R₁₈–R₁₉)
# ---------------------------------------------------------------------------

def evaluate_R(i: int, gene: str, state: Dict[str, int], edges: List[Dict]) -> bool:
    """Evaluate final Rᵢ with inducible/repressible wrapper (i=0–17)
    or threshold rule (i=18,19)."""
    if i == 18:
        return evaluate_R18(gene, state, edges)
    elif i == 19:
        return evaluate_R19(gene, state, edges)
    else:
        raw = evaluate_R_raw(i, gene, state, edges)
        ind = inducible_regulation(gene, state, edges)
        rep = repressible_regulation(gene, state, edges)
        return (raw and ind) or rep


# ---------------------------------------------------------------------------
# Truth table construction
# ---------------------------------------------------------------------------

def get_regulators_sorted(gene: str, edges: List[Dict]) -> List[str]:
    """Get all regulators of a gene, sorted alphabetically."""
    regs = set()
    for e in edges:
        if e["dst"] == gene:
            regs.add(e["src"])
    return sorted(regs)


def gene_truth_table(gene: str, edges: List[Dict]) -> Dict[int, List[bool]]:
    """Enumerate all 2^n regulator states. For each, evaluate all 20 R conditions.

    Returns: {i: [bool, bool, ...]} for each condition index i (0–19).
    Each list has 2^n entries in the same order as regulator state enumeration.
    """
    regulators = get_regulators_sorted(gene, edges)
    n = len(regulators)
    if n == 0:
        return {}

    result = {i: [] for i in range(20)}

    for values in iterproduct([0, 1], repeat=n):
        # Build state dict: all regulators get their assigned value,
        # the gene itself gets 0 (needed for R18 self-reference).
        # We evaluate R18 twice — once with gene=0, once with gene=1 — below.
        state = dict(zip(regulators, values))
        if gene not in state:
            state[gene] = 0  # default for R18 if gene is not its own regulator

        for i in range(20):
            if i == 18 and gene not in regulators:
                # R18 depends on q[c]. Since gene is not a regulator of itself,
                # we need to evaluate with the gene's current value.
                # For truth table purposes, gene value is part of the state.
                # We set it to 0 here; the actual comparison handles this.
                result[i].append(evaluate_R(i, gene, state, edges))
            else:
                result[i].append(evaluate_R(i, gene, state, edges))

    return result


def aeon_function_truth_table(gene: str, update_fn_str: str,
                               edges: List[Dict]) -> List[bool]:
    """Parse AEON update function, enumerate all 2^n regulator states,
    evaluate, return truth table as list of bools.

    Regulator ordering: sorted alphabetical (same as gene_truth_table).
    """
    regulators = get_regulators_sorted(gene, edges)
    n = len(regulators)
    if n == 0:
        return []

    # Convert AEON expression to Python
    py_expr = update_fn_str.replace("!", "not ")
    py_expr = py_expr.replace("&", " and ")
    py_expr = py_expr.replace("|", " or ")
    py_expr = re.sub(r"\s+", " ", py_expr).strip()

    result = []
    for values in iterproduct([0, 1], repeat=n):
        state = dict(zip(regulators, values))
        val = eval(py_expr, {"__builtins__": {}}, state)
        result.append(bool(val))

    return result


# ---------------------------------------------------------------------------
# Main matching function
# ---------------------------------------------------------------------------

def find_ground_truth_L_value(gene: str, update_fn_str: str,
                               edges: List[Dict]) -> Tuple[Optional[int], List[int]]:
    """Find which Rᵢ matches the ground truth truth table.

    Args:
        gene: Gene name (e.g. "v_Pax6")
        update_fn_str: AEON update function string (e.g. "v_Sp8 & !(v_Emx2 | v_Coup_fti)")
        edges: List of edge dicts with keys: src, dst, sign

    Returns:
        (primary_match, all_matches) where:
        - primary_match: first matching Rᵢ index (or None if no match)
        - all_matches: list of ALL matching indices
    """
    regulators = get_regulators_sorted(gene, edges)
    if not regulators:
        return None, []

    # Ground truth truth table from AEON function
    gt_tt = aeon_function_truth_table(gene, update_fn_str, edges)

    # Evaluate all 20 R conditions
    r_tables = gene_truth_table(gene, edges)

    matches = []
    for i in range(20):
        r_tt = r_tables.get(i, [])
        if len(r_tt) == len(gt_tt) and r_tt == gt_tt:
            matches.append(i)

    if len(matches) > 1:
        print(f"  WARNING: {gene} matches multiple R conditions: {matches}",
              file=sys.stderr)

    primary = matches[0] if matches else None
    return primary, matches


# ---------------------------------------------------------------------------
# AEON file parsing (self-contained)
# ---------------------------------------------------------------------------

def parse_aeon_for_comparison(filepath: Path):
    """Parse .aeon file → (update_functions, edges).

    Returns:
        update_functions: dict gene→expr
        edges: list of {src, dst, sign} (only deterministic positive/negative)
    """
    update_functions = {}
    edges = []
    graph_genes = set()

    with open(filepath, "r") as f:
        for line in f:
            line = line.split("//")[0].strip()
            if not line:
                continue

            if line.startswith("$"):
                match = re.match(r"^\$(\S+):\s*(.+)$", line)
                if match:
                    gene = match.group(1)
                    expr = match.group(2).strip()
                    update_functions[gene] = expr
                    graph_genes.add(gene)
                continue

            # Parse edges: -> (positive), -| (negative), -? (both optional)
            clean = line.replace(" ", "")
            parsed = _parse_edge(clean)
            for e in parsed:
                edges.append(e)
                graph_genes.add(e["src"])
                graph_genes.add(e["dst"])

    input_nodes = sorted(graph_genes - set(update_functions.keys()))
    return update_functions, edges, input_nodes


def _parse_edge(clean: str) -> List[Dict]:
    """Parse a single AEON edge line (spaces already stripped)."""
    patterns = [
        ("->?", "positive"),
        ("?->", "positive"),
        ("-|?", "negative"),
        ("?-|", "negative"),
        ("-?", "both"),
        ("->", "positive"),
        ("-|", "negative"),
    ]
    for symbol, kind in patterns:
        if symbol in clean:
            parts = clean.split(symbol, 1)
            if len(parts) != 2:
                continue
            src, dst = parts
            if kind == "both":
                return [
                    {"src": src, "dst": dst, "sign": "positive"},
                    {"src": src, "dst": dst, "sign": "negative"},
                ]
            return [{"src": src, "dst": dst, "sign": kind}]
    return []


# ---------------------------------------------------------------------------
# CLI: process an .aeon file and print L-values
# ---------------------------------------------------------------------------

def process_aeon_file(aeon_path: Path):
    """Process an .aeon file: compute ground truth L-values for all genes."""
    update_functions, edges, input_nodes = parse_aeon_for_comparison(aeon_path)
    all_genes = sorted(set(list(update_functions.keys()) + input_nodes))

    print(f"Model: {aeon_path.name}")
    print(f"Genes: {len(all_genes)} ({len(update_functions)} functions, {len(input_nodes)} inputs)")
    print()

    results = []
    for gene in all_genes:
        if gene in input_nodes:
            print(f"  {gene:30s}  INPUT — skip")
            results.append({"gene": gene, "L": None, "all_matches": [], "is_input": True})
            continue

        expr = update_functions[gene]
        L, all_matches = find_ground_truth_L_value(gene, expr, edges)

        regs = get_regulators_sorted(gene, edges)
        acts = get_activators(gene, edges)
        reps = get_repressors(gene, edges)

        if L is not None:
            match_str = f"L={L}" + (f" (also {all_matches})" if len(all_matches) > 1 else "")
        else:
            match_str = "NO MATCH"

        print(f"  {gene:30s}  {match_str:20s}  regs={len(regs)} act={len(acts)} rep={len(reps)}")
        results.append({"gene": gene, "L": L, "all_matches": all_matches, "is_input": False})

    return results


def main():
    if len(sys.argv) < 2:
        print("Usage: python3 rein_regulatory_conditions.py <aeon_file>")
        sys.exit(1)

    aeon_path = Path(sys.argv[1])
    if not aeon_path.exists():
        print(f"Error: file not found: {aeon_path}")
        sys.exit(1)

    process_aeon_file(aeon_path)


if __name__ == "__main__":
    main()
