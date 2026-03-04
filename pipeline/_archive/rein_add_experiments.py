import json
import random
import re
from pathlib import Path
from typing import Dict

from aeon_to_rein import aeon_to_rein_pipeline


# =========================
# Configuration (ONLY HERE)
# =========================

BASE_DIR = Path(__file__).resolve().parent
MODEL_NAME = "example"

AEON_FILE = BASE_DIR / "aeon" / f"{MODEL_NAME}.aeon"
REIN_FILE = BASE_DIR / "rein" / f"{MODEL_NAME}.rein"
JSON_FILE = BASE_DIR / "rein" / f"{MODEL_NAME}.json"

# Trajectory length (K) and number of experiments (N) - overridden by analyze_aeon.py
K = 1
N = 5
DYNAMICS_MODE = "sync"


# =========================
# JSON loading
# =========================

def load_json(json_path: Path):
    with json_path.open() as f:
        return json.load(f)


# =========================
# Simulation
# =========================

def random_initial_state(variables):
    return {v: random.randint(0, 1) for v in variables}


def eval_update_function(expr: str, state: Dict[str, int]) -> int:
    python_expr = expr.replace("!", " not ").replace("&", " and ").replace("|", " or ")

    for var, value in state.items():
        # Escape var so names with regex metacharacters (e.g. ?) do not break the pattern
        escaped = re.escape(var)
        python_expr = re.sub(rf"\b{escaped}\b", str(bool(value)), python_expr)

    return int(eval(python_expr))


def simulate(initial_state, update_functions, K):
    state = initial_state.copy()

    for _ in range(K):
        next_state = state.copy()
        for var, expr in update_functions.items():
            next_state[var] = eval_update_function(expr, state)
        state = next_state

    return state


# =========================
# REIN experiment blocks
# =========================

def format_state_block(name, state):
    lines = [f"{var} = {val}" for var, val in sorted(state.items())]
    return f"${name} := {{\n  " + " and\n  ".join(lines) + "\n};\n"


def generate_experiment_block(i: int, init_state: Dict[str, int], final_state: Dict[str, int], trajectory_k: int) -> str:
    """Generate #ExperimentN block with step 0..K alignment. K = trajectory length."""
    block = ""
    block += f"// Experiment {i}\n"
    block += f"#Experiment{i}[0] |= $INIT{i} \"initial state\";\n"
    block += f"#Experiment{i}[0] |= $NoKnockDowns \"no knockdowns\";\n"
    block += f"#Experiment{i}[0] |= $NoOverExpression \"no overexpression\";\n"
    block += f"#Experiment{i}[{trajectory_k}] |= $FIN{i} \"final state\";\n\n"
    block += format_state_block(f"INIT{i}", init_state) + "\n"
    block += format_state_block(f"FIN{i}", final_state) + "\n"
    return block


def generate_global_blocks(variables):
    block = "$NoKnockDowns := {\n  "
    block += " and\n  ".join(f"KO({v}) = 0" for v in variables)
    block += "\n};\n\n"

    block += "$NoOverExpression := {\n  "
    block += " and\n  ".join(f"FE({v}) = 0" for v in variables)
    block += "\n};\n\n"

    return block


# =========================
# Main
# =========================

def main():
    # Stage 1: AEON → REIN base + JSON (with dynamics mode)
    aeon_to_rein_pipeline(AEON_FILE, REIN_FILE, JSON_FILE, DYNAMICS_MODE)

    # Stage 2: Dynamic experiments - use K (trajectory length) and N (count) from module
    data = load_json(JSON_FILE)
    trajectory_k = max(1, K)  # Ensure at least 1 step
    num_experiments = max(1, N)

    # Generate exactly N experiment blocks (allow repeated initial states so N is always met)
    with REIN_FILE.open("a") as f:
        f.write("\n\n// =====================\n")
        f.write("// Experiments\n")
        f.write("// =====================\n\n")

        for i in range(num_experiments):
            init_state = random_initial_state(data["variables"])
            final_state = simulate(init_state, data["update_functions"], trajectory_k)
            f.write(generate_experiment_block(i, init_state, final_state, trajectory_k))

        f.write(generate_global_blocks(data["variables"]))

    print("✔ Full pipeline completed successfully")


if __name__ == "__main__":
    main()
