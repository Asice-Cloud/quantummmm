#!/usr/bin/env python3
"""
Simple Anyon F/R helper and scaffold

Features:
- Provides built-in F and R data for common models (Ising, Fibonacci).
- Parses simple fusion rules.
- Scaffolding functions to build pentagon/hexagon equations (TODO: general solver).

This is a starting point — solving general pentagon/hexagon equations
requires a nonlinear solver and careful gauge fixing. The provided
built-in models show the expected F and R matrices.
"""
from __future__ import annotations
import cmath
import json
import argparse
from typing import Dict, Tuple, List, Any
try:
    import sympy as sp
except Exception:
    sp = None


def format_matrix(mat):
    """Format a matrix (nested lists/tuples or sympy.Matrix) into aligned rows/columns."""
    try:
        from sympy import Matrix
    except Exception:
        Matrix = None

    if Matrix is not None and isinstance(mat, Matrix):
        rows = [[mat[i, j] for j in range(mat.cols)] for i in range(mat.rows)]
    else:
        rows = [list(r) for r in mat]

    # convert entries to strings with modest precision
    str_rows = []
    for r in rows:
        srow = []
        for x in r:
            if isinstance(x, complex):
                s = f"{x.real:.6g}+{x.imag:.6g}j"
            elif isinstance(x, float):
                s = f"{x:.6g}"
            else:
                s = str(x)
            srow.append(s)
        str_rows.append(srow)

    # compute column widths
    if not str_rows:
        return "[]"
    cols = len(str_rows[0])
    widths = [0] * cols
    for r in str_rows:
        for j, v in enumerate(r):
            widths[j] = max(widths[j], len(v))

    # build lines
    lines = []
    for r in str_rows:
        pieces = [v.rjust(widths[j]) for j, v in enumerate(r)]
        lines.append("[ " + "  ".join(pieces) + " ]")
    return "\n".join(lines)

Anyon = str
FusionRules = Dict[Tuple[Anyon, Anyon], List[Anyon]]

# Examples (Ising/Fibonacci) are kept in a separate module for clarity.
import os
import sys

# Ensure project root is on sys.path so `python3 tools/anyon_solver.py` can
# import the `tools` package when executed from the project root.
_this_dir = os.path.dirname(__file__)
_project_root = os.path.abspath(os.path.join(_this_dir, '..'))
if _project_root not in sys.path:
    sys.path.insert(0, _project_root)

try:
    from tools.examples import ising_model, fibonacci_model
except Exception:
    # If the module is executed as a script (python tools/anyon_solver.py),
    # the package context may not be set; try a relative import fallback.
    try:
        from .examples import ising_model, fibonacci_model
    except Exception:
        # final fallback: clear error helpers
        def ising_model():
            raise RuntimeError("tools.examples not available; restore ising_model or adjust PYTHONPATH")

        def fibonacci_model():
            raise RuntimeError("tools.examples not available; restore fibonacci_model or adjust PYTHONPATH")


def get_standard_model(name: str):
    name = name.lower()
    if name in ("ising", "majorana"):
        return ising_model()
    if name in ("fibonacci", "golden"):
        return fibonacci_model()
    raise ValueError("Unknown model")


def print_model(name: str):
    rules, F, R = get_standard_model(name)
    print("Model:", name)
    print("Fusion rules:")
    for (a, b), clist in rules.items():
        print(f"  {a} x {b} -> {clist}")
    print("\nF matrices (sample):")
    for k, mat in F.items():
        print(f"  {k} =")
        try:
            print(format_matrix(mat))
        except Exception:
            print(mat)
    print("\nR phases:")
    for k, v in R.items():
        print(f"  {k} = {v}")


def dump_pentagon_strings(fusion_rules: FusionRules, out=None):
    """
    For each 5-tuple of labels (a,b,c,d,e) and total s, produce two explicit
    index-sum strings representing the two pentagon paths. This does not
    attempt to solve equations, but gives a concrete algebraic expansion
    (symbol names) suitable for feeding into a symbolic solver later.
    """
    # produce structured JSON-like entries for each relevant F block
    Fs = gen_F_symbols(fusion_rules)
    entries = []
    for key, info in Fs.items():
        d, a, b, c = key
        E = info['E']
        Flist = info['Flist']
        mat = info['mat']
        # flatten symbol names
        sym_names = [str(s) for s in mat]
        lhs_expr = ' + '.join(sym_names) if sym_names else '0'
        entry = {
            'labels': [a, b, c],
            'd': d,
            'E': E,
            'Flist': Flist,
            'symbols': sym_names,
            'lhs_expr': lhs_expr,
            'rhs_expr': 'TODO',
        }
        entries.append(entry)

    # output formats
    if out is None:
        print(json.dumps(entries, indent=2))
    else:
        if out.endswith('.json'):
            with open(out, 'w') as f:
                json.dump(entries, f, indent=2)
        elif out.endswith('.py'):
            # write a python file that constructs sympy symbols and expressions
            if sp is None:
                raise RuntimeError('sympy required to dump .py expressions')
            with open(out, 'w') as f:
                f.write('import sympy as sp\n\n')
                # collect unique symbol names
                all_syms = sorted({s for e in entries for s in e['symbols']})
                if all_syms:
                    f.write('symbols = sp.symbols("' + ' '.join(all_syms) + '")\n')
                    f.write('# symbols is a tuple matching the names above\n')
                f.write('\nexpressions = []\n')
                for i, e in enumerate(entries):
                    if not e['symbols']:
                        f.write(f'expressions.append(0)\n')
                        continue
                    # build sum by indexing into symbols tuple
                    terms = []
                    for s in e['symbols']:
                        idx = all_syms.index(s)
                        terms.append(f'symbols[{idx}]')
                    f.write('expressions.append(' + ' + '.join(terms) + ')\n')
        else:
            with open(out, 'w') as f:
                f.write(json.dumps(entries, indent=2))



def build_pentagon_equations(fusion_rules: FusionRules):
    """
    Scaffolding: construct symbolic pentagon equations from fusion_rules.
    Full implementation is nontrivial; this function is a placeholder that
    documents where to build index sums and equation sets.
    """
    if sp is None:
        raise RuntimeError("sympy is required to build pentagon equations; please install sympy")
    # Placeholder: this function will be implemented to return a list of
    # sympy equations representing the pentagon constraints for the provided
    # fusion rules. For now, provide a helpful message.
    raise NotImplementedError("Pentagon equation builder not yet implemented")


def build_unitarity_equations(fusion_rules: FusionRules, out_py: str = None):
    """
    Build unitarity equations F^
    and optionally dump a SymPy .py script that defines the symbols and the
    equations Eq(F.H*F - I, 0) for each F block.
    This helps to verify and to provide inputs to a numerical solver.
    """
    if sp is None:
        raise RuntimeError("sympy required to build unitarity equations")
    Fs = gen_F_symbols(fusion_rules)
    eqs = []
    # collect unique symbols and mapping
    for key, info in Fs.items():
        mat: sp.Matrix = info['mat']
        # build unitary condition: mat.H * mat - I = 0
        # note: mat may be rectangular; require square for unitarity
        if mat.shape[0] != mat.shape[1]:
            continue
        I = sp.eye(mat.shape[0])
        eq = (mat.H * mat - I)
        eqs.append({'key': key, 'mat': mat, 'eq': eq})

    if out_py is None:
        # print summary
        for item in eqs:
            d,a,b,c = item['key']
            print(f"Unitary eq for F^{d}_{{{a}{b}{c}}}: matrix shape={item['mat'].shape}")
    else:
        # write a .py file that constructs these symbols and equations
        with open(out_py, 'w') as f:
            f.write('import sympy as sp\n')
            f.write('\n')
            # collect all unique symbol names
            all_syms = []
            for item in eqs:
                for s in item['mat']:
                    name = str(s)
                    if name not in all_syms:
                        all_syms.append(name)
            if all_syms:
                f.write('symbols = sp.symbols("' + ' '.join(all_syms) + '")\n')
                f.write('\n')
            # helper to map name->symbol
            f.write('sym_map = {')
            for i, name in enumerate(all_syms):
                f.write(f'"{name}": symbols[{i}],')
            f.write('}\n\n')
            f.write('eqs = []\n')
            for item in eqs:
                d,a,b,c = item['key']
                mat = item['mat']
                rows, cols = mat.shape
                f.write(f'# F^{d}_{{{a}{b}{c}}} ({rows}x{cols})\n')
                f.write('M = sp.Matrix([')
                row_strs = []
                for r in range(rows):
                    elems = []
                    for cidx in range(cols):
                        name = str(mat[r, cidx])
                        elems.append(f'sym_map["{name}"]')
                    row_strs.append('[' + ','.join(elems) + ']')
                f.write(','.join(row_strs) + '])\n')
                f.write('eqs.append(sp.Eq(M.H*M, sp.eye(%d)))\n' % rows)
            f.write('\n# eqs is a list of SymPy Eq objects representing unitarity constraints\n')
    return eqs


def build_hexagon_equations(fusion_rules: FusionRules):
    if sp is None:
        raise RuntimeError("sympy is required to build hexagon equations; please install sympy")
    # Placeholder for hexagon construction
    raise NotImplementedError("Hexagon equation builder not yet implemented")


def gen_F_symbols(fusion_rules: FusionRules):
    """
    Generate sympy symbols for all needed F^d_{a b c} matrices (for N_{ab}^c in {0,1}).
    Returns dict keyed by (d,a,b,c) mapping to a structure with keys:
      - 'E': list of e values (rows)
      - 'Flist': list of f values (cols)
      - 'mat': sympy.Matrix of symbols of shape (len(E), len(Flist))
    The sets are defined by the standard fusion constraints:
      E = { e | N_ab^e = 1 and N_ec^d = 1 }
      Flist = { f | N_bc^f = 1 and N_af^d = 1 }
    """
    if sp is None:
        raise RuntimeError("sympy required for symbol generation")
    Fs = {}
    # collect labels present in fusion_rules
    labels = set()
    for (x, y), outs in fusion_rules.items():
        labels.add(x); labels.add(y)
        for o in outs: labels.add(o)
    labels = sorted(labels)
    for a in labels:
        for b in labels:
            for c in labels:
                for d in labels:
                    E = [e for e in fusion_rules.get((a, b), []) if d in fusion_rules.get((e, c), [])]
                    Flist = [f for f in fusion_rules.get((b, c), []) if d in fusion_rules.get((a, f), [])]
                    if not E or not Flist:
                        continue
                    # create matrix of fresh symbols
                    mat = sp.Matrix([[sp.symbols(f'F_{d}_{a}_{b}_{c}_{e}_{f}') for f in Flist] for e in E])
                    Fs[(d, a, b, c)] = {'E': E, 'Flist': Flist, 'mat': mat}
    return Fs


def pentagon_skeleton(fusion_rules: FusionRules, labels: Tuple[str, str, str, str, str], total: str):
    """
    Return a human-readable skeleton of the pentagon equation for five anyons
    labels = (a,b,c,d,e) and total is the total charge 'total'. This function
    does not expand the full algebraic sum but shows the summed index structure
    for the two canonical paths; it is meant to help formulate the symbolic
    equations using the output of `gen_F_symbols`.
    """
    a, b, c, d, e = labels
    # Two canonical paths connecting the same start/end trees produce products
    # of three F matrices vs two F matrices; here we return skeleton strings.
    path1 = (
        "P1: Sum_{g,h} F^{?}_{%s %s %s}[g, _] * F^{?}_{%s g %s}[h, _] * F^{?}_{%s %s %s}[_, _]"
        % (a, b, c, a, d, b, c, d)
    )
    path2 = "P2: Sum_{x,y} F... (other path product)"
    return {'path1': path1, 'path2': path2, 'note': 'Use gen_F_symbols to obtain index sets and build full sums.'}


def main():
    p = argparse.ArgumentParser()
    p.add_argument("model", nargs="?", default="ising", help="standard model: ising or fibonacci")
    # `--gen` enabled by default; use `--no-gen` to disable.
    p.add_argument("--gen", dest="gen", action="store_true", help="Generate and pretty-print F symbols from fusion rules")
    p.add_argument("--no-gen", dest="gen", action="store_false", help="Don't generate F symbols; just print model")
    p.set_defaults(gen=True)
    args = p.parse_args()
    if args.gen:
        # Try to generate F symbols; if sympy is not available, fall back
        # to printing the model to avoid crashing the CLI.
        if sp is None:
            print("SymPy not available; falling back to printing model.")
            print_model(args.model)
            return
        rules, _, _ = get_standard_model(args.model)
        Fs = gen_F_symbols(rules)
        if not Fs:
            print("No generated F symbols for model", args.model)
            return
        print(f"Generated F symbols for model: {args.model}")
        for key, info in Fs.items():
            d, a, b, c = key
            print(f"\nF^{d}_{{{a}{b}{c}}} (rows={len(info['E'])}, cols={len(info['Flist'])}):")
            try:
                print(format_matrix(info['mat']))
            except Exception:
                print(info['mat'])
    else:
        print_model(args.model)


if __name__ == "__main__":
    main()
