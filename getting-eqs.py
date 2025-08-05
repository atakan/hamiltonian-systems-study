import sys
import argparse
import sympy as sp
import numpy as np
from sympy.parsing.sympy_parser import parse_expr

def main():
    parser = argparse.ArgumentParser(description='Symbolic differentiation of Hamiltonian')
    parser.add_argument('-i', '--input', required=True, help='Input file containing Hamiltonian')
    parser.add_argument('-o', '--output', help='Output file for generated functions')
    parser.add_argument('-v', '--verbose', action='store_true', help='Print progress messages')
    args = parser.parse_args()

    if args.verbose:
        print("Reading input file...", file=sys.stderr)

    try:
        with open(args.input, 'r') as f:
            lines = f.readlines()
            N = int(lines[0].strip())
            H_expr = lines[1].strip()
    except Exception as e:
        print(f"Error reading input file: {e}", file=sys.stderr)
        sys.exit(1)

    if args.verbose:
        print(f"Parsing Hamiltonian with N={N}...", file=sys.stderr)

    # Create symbols
    q = sp.symbols([f'q{i+1}' for i in range(N)])
    p = sp.symbols([f'p{i+1}' for i in range(N)])
    t = sp.symbols('t')

    # Parse the Hamiltonian expression
    try:
        H = parse_expr(H_expr)
    except Exception as e:
        print(f"Error parsing Hamiltonian expression: {e}", file=sys.stderr)
        sys.exit(1)

    if args.verbose:
        print("Computing derivatives...", file=sys.stderr)

    # Compute dqdt (∂H/∂p)
    dqdt_exprs = [sp.diff(H, p_i) for p_i in p]
    
    # Compute dpdt (-∂H/∂q)
    dpdt_exprs = [-sp.diff(H, q_i) for q_i in q]

    # Generate the Python functions
    def generate_function(name, exprs):
        lines = []
        lines.append(f"def {name}(q, p, t, N):")
        lines.append(f"    return np.asarray([")
        for i, expr in enumerate(exprs):
            # Convert SymPy expr to Python/Numpy code
            code = sp.pycode(expr).replace('math.', 'np.')
            lines.append(f"        {code},")
        lines.append("    ])")
        return '\n'.join(lines)

    dqdt_func = generate_function('dqdt', dqdt_exprs)
    dpdt_func = generate_function('dpdt', dpdt_exprs)

    output_content = f"{dqdt_func}\n\n{dpdt_func}\n"

    if args.output:
        if args.verbose:
            print(f"Writing output to {args.output}...", file=sys.stderr)
        try:
            with open(args.output, 'w') as f:
                f.write(output_content)
        except Exception as e:
            print(f"Error writing output file: {e}", file=sys.stderr)
            sys.exit(1)
    else:
        print(output_content)

    if args.verbose:
        print("Successfully generated dqdt() and dpdt().", file=sys.stderr)

if __name__ == '__main__':
    main()
