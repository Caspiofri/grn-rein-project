"""
VM Jupyter Notebook — L-Value Extraction Code (F#)

This file stores the F# code added to the RE:IN Jupyter notebook on the VM
that extracts L-values from the RE:IN engine solutions. Having it here in the
git repo ensures it is versioned and not lost.

The code below produces a CSV with columns:
  solution_index, gene, L_value

How it works:
  - solutions = ReinAPI.Enumerate <limit> <model>
  - Each solution has a .solution property (Option<Solution>)
  - solution.vars is a Map<string, string> containing all SMT solver variables
  - Variables named "L_v_<gene>" hold the L-value for that gene
  - The L-value is RE:IN's internal encoding of the Boolean update function

To inspect a solution's structure:
  solutions.[0].species   -- shows reg_conds (L-value), lVar name per gene
  solutions.[0].solution  -- shows the full SMT solution with vars map
"""

# =============================================================================
# F# CODE — ADD THESE CELLS TO THE BOTTOM OF THE RE:IN JUPYTER NOTEBOOK
# =============================================================================

FSHARP_CELL_1_LOAD_AND_SOLVE = """
// Cell 1: Load model and enumerate solutions
let model = ReinAPI.LoadFile "./path/to/model.rein"
ReinAPI.CheckAndPrint model

let solutions = ReinAPI.Enumerate 65 model
printfn "Solution count: %d" (Seq.length solutions)
"""

FSHARP_CELL_2_EXTRACT_L_VALUES = """
// Cell 2: Extract L-values to CSV
open System.IO

let outputDir = "./output"
Directory.CreateDirectory(outputDir) |> ignore

let lines = System.Collections.Generic.List<string>()
lines.Add("solution_index,gene,L_value")

// Loop over ALL solutions
for i in 0 .. solutions.Length - 1 do
    match solutions.[i].solution with
    | Some s ->
        s.vars |> Map.iter (fun key value ->
            if key.StartsWith("L_") then
                lines.Add(sprintf "%d,%s,%s" i (key.Replace("L_","")) value))
    | None -> ()

File.WriteAllLines(Path.Combine(outputDir, "l_values_extracted.csv"), lines)
printfn "Saved %d rows across %d solutions" (lines.Count - 1) solutions.Length
"""

FSHARP_CELL_3_INSPECT_SPECIES = """
// Cell 3 (optional): Inspect species details — shows reg_conds and lVar
printfn "%A" solutions.[0].species
"""

FSHARP_CELL_4_INSPECT_SOLUTION_VARS = """
// Cell 4 (optional): Print all solution variables
match solutions.[0].solution with
| Some s -> s.vars |> Map.iter (fun k v -> printfn "%s = %s" k v)
| None -> printfn "No solution"
"""

# =============================================================================
# END OF F# CODE
# =============================================================================
