import json
import re
from pathlib import Path
from typing import Dict, List, Set


# =========================
# AEON parsing
# =========================

def parse_aeon(file_path: Path):
    variables: Set[str] = set()
    update_functions: Dict[str, str] = {}
    regulations: List[Dict] = []

    with file_path.open() as f:
        for raw_line in f:
            line = raw_line.split("//")[0].strip()
            if not line:
                continue

            # Update function
            if line.startswith("$"):
                var, expr = line[1:].split(":", 1)
                var = var.strip()
                expr = expr.strip()
                update_functions[var] = expr
                variables.add(var)
                variables |= extract_vars_from_expr(expr)
                continue

            # Regulation parsing
            regs = parse_regulation(line)
            if regs:
                regulations.extend(regs)
                for r in regs:
                    variables.add(r["src"])
                    variables.add(r["dst"])

    return variables, update_functions, regulations


def extract_vars_from_expr(expr: str) -> Set[str]:
    return set(re.findall(r"[A-Za-z_][A-Za-z0-9_]*", expr))


def parse_regulation(line: str):
    """
    Parse AEON regulation line.
    Optional: ->? (positive), -|? (negative), or -? (optional, sign unknown = both).
    Deterministic: -> (positive), -| (negative).
    Maps to RE:IN 'optional' or plain.
    """
    line = line.replace(" ", "")

    # Order matters: check longer / optional symbols first; -?? and -? (unknown) before -> and -|
    patterns = [
        ("-??", "both", True),  # double-? variant seen in some AEON files (same semantics as -?)
        ("->?", "positive", True),
        ("?->", "positive", True),
        ("-|?", "negative", True),
        ("?-|", "negative", True),
        ("-?", "both", True),   # optional, unknown whether activating or inhibiting
        ("->", "positive", False),
        ("-|", "negative", False),
    ]

    for (symbol, kind, optional) in patterns:
        if symbol in line:
            parts = line.split(symbol, 1)
            if len(parts) != 2:
                continue
            src, dst = parts
            if kind == "both":
                return [
                    {"src": src, "dst": dst, "sign": "positive", "optional": optional},
                    {"src": src, "dst": dst, "sign": "negative", "optional": optional},
                ]
            return [{
                "src": src,
                "dst": dst,
                "sign": kind,
                "optional": optional,
            }]

    return []


def get_regulation_summary(aeon_file: Path) -> Dict:
    """
    Parse AEON file and return optional vs deterministic regulation counts and lists.
    Edge display: ->? / -|? (optional), -? (optional unknown), -> / -| (deterministic).
    """
    _, _, regulations = parse_aeon(aeon_file)
    optional_list: List[Dict] = []
    deterministic_list: List[Dict] = []
    # Pairs (src,dst) that came from "-?" (optional unknown) - show once as "A -? B"
    seen_unknown: set = set()
    for r in regulations:
        key = (r["src"], r["dst"])
        if r.get("optional"):
            if key in seen_unknown:
                continue  # already added as "A -? B"
            signs_for_key = [x.get("sign") for x in regulations if x.get("src") == r["src"] and x.get("dst") == r["dst"] and x.get("optional")]
            if "positive" in signs_for_key and "negative" in signs_for_key:
                seen_unknown.add(key)
                optional_list.append({"src": r["src"], "dst": r["dst"], "sign": "unknown", "edge": f"{r['src']} -? {r['dst']}"})
                continue
            edge = f"{r['src']} ->? {r['dst']}" if r.get("sign") == "positive" else f"{r['src']} -|? {r['dst']}"
        else:
            edge = f"{r['src']} -> {r['dst']}" if r.get("sign") == "positive" else f"{r['src']} -| {r['dst']}"
        entry = {"src": r["src"], "dst": r["dst"], "sign": r.get("sign", ""), "edge": edge}
        if r.get("optional"):
            optional_list.append(entry)
        else:
            deterministic_list.append(entry)
    return {
        "optional": optional_list,
        "deterministic": deterministic_list,
        "optional_count": len(optional_list),
        "deterministic_count": len(deterministic_list),
    }


# =========================
# Input detection
# =========================

def detect_inputs(variables, update_functions, regulations):
    incoming = {v: 0 for v in variables}
    for r in regulations:
        incoming[r["dst"]] += 1

    return [
        v for v in variables
        if v not in update_functions and incoming[v] == 0
    ]


def add_input_self_loops(regulations, inputs):
    for v in inputs:
        regulations.append({
            "src": v,
            "dst": v,
            "sign": "positive",
            "optional": True,
        })


# =========================
# Writing RE:IN
# =========================

def write_rein_base(rein_path: Path, variables, regulations, dynamics_mode: str = "sync"):
    """
    Write RE:IN base file. Regulations: optional (AEON with ?) -> 'sign optional;',
    mandatory (no ?) -> 'sign;' only.
    """
    dynamics_mode = (dynamics_mode or "sync").lower().strip()
    if dynamics_mode not in ("sync", "async"):
        dynamics_mode = "sync"
    mode_comment = "Synchronous dynamics" if dynamics_mode == "sync" else "Asynchronous dynamics"

    with rein_path.open("w") as f:
        f.write(f"// {mode_comment}\n")
        f.write(f"directive updates {dynamics_mode};\n\n")
        f.write("// Default regulation conditions\n")
        f.write("directive regulation legacy;\n\n")

        parts = [f"{v}[-+] (0..17)" for v in sorted(variables)]
        f.write("; ".join(parts) + ";\n\n")

        for r in regulations:
            # AEON optional (->?, -|?, -?) -> RE:IN "optional"; deterministic (->, -|) -> no keyword
            if r.get("optional", False):
                f.write(f"{r['src']} {r['dst']} {r['sign']} optional;\n")
            else:
                f.write(f"{r['src']} {r['dst']} {r['sign']};\n")



# =========================
# Writing JSON
# =========================

def write_json(json_path: Path, variables, update_functions, inputs, regulations):
    data = {
        "variables": sorted(variables),
        "update_functions": update_functions,
        "inputs": sorted(inputs),
        "regulations": regulations
    }

    with json_path.open("w") as f:
        json.dump(data, f, indent=2)


# =========================
# Public pipeline function
# =========================

def aeon_to_rein_pipeline(aeon_file: Path, rein_file: Path, json_file: Path, dynamics_mode: str = "sync"):
    """
    Convert AEON to RE:IN. dynamics_mode: 'sync' or 'async'.
    Optional regulations (AEON with ?) -> RE:IN optional; mandatory -> no optional keyword.
    """
    variables, updates, regulations = parse_aeon(aeon_file)
    inputs = detect_inputs(variables, updates, regulations)
    add_input_self_loops(regulations, inputs)

    write_rein_base(rein_file, variables, regulations, dynamics_mode)
    write_json(json_file, variables, updates, inputs, regulations)
