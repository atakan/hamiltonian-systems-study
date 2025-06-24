import sympy as sp
import numpy as np
from sympy.printing.pycode import pycode

def define_symbols(n):
    q_syms = sp.Matrix(sp.symbols(f'q0:{n}'))
    p_syms = sp.Matrix(sp.symbols(f'p0:{n}'))
    return q_syms, p_syms

def build_hamiltonian(T_func, U_func, q_syms, p_syms):
    q_array = np.array(q_syms).astype(object)
    p_array = np.array(p_syms).astype(object)
    return T_func(q_array, p_array) + U_func(q_array, p_array)

def compute_derivatives(H, q_syms, p_syms):
    qdot = [sp.diff(H, pi) for pi in p_syms]
    pdot = [-sp.diff(H, qi) for qi in q_syms]
    return qdot, pdot

def generate_numpy_function_code(name, exprs, q_syms, p_syms):
    """Generate a NumPy-based Python function."""
    func_lines = [
        "import numpy as np",
        f"def {name}(q, p):"
    ]
    for i, expr in enumerate(exprs):
        code = pycode(expr)
        for j, sym in enumerate(q_syms):
            code = code.replace(pycode(sym), f"q[{j}]")
        for j, sym in enumerate(p_syms):
            code = code.replace(pycode(sym), f"p[{j}]")
        func_lines.append(f"    out_{i} = {code}")
    func_lines.append(f"    return np.array([{', '.join(f'out_{i}' for i in range(len(exprs)))}])")
    return "\n".join(func_lines)

def hamiltonian_to_numpy_code(T_func, U_func, n):
    q_syms, p_syms = define_symbols(n)
    H = build_hamiltonian(T_func, U_func, q_syms, p_syms)
    qdot_syms, pdot_syms = compute_derivatives(H, q_syms, p_syms)

    qdot_code = generate_numpy_function_code("qdot", qdot_syms, q_syms, p_syms)
    pdot_code = generate_numpy_function_code("pdot", pdot_syms, q_syms, p_syms)

    print("# === Python function for qdot(q, p) ===")
    print(qdot_code)
    print("\n# === Python function for pdot(q, p) ===")
    print(pdot_code)

# === Example symbolic T and U ===
def T_example(q, p):
    return 0.5 * np.dot(p, p)

def U_example(q, p):
    return 0.5 * np.dot(q, q)

# === Generate functions for dimension n ===
hamiltonian_to_numpy_code(T_example, U_example, n=2)
