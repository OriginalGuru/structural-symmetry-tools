#!/usr/bin/env python3
"""
run_findsym.py

Generate a FINDSYM input file from a POSCAR and run findsym, writing
a human-readable log and a CIF file.

This script calls poscar_to_findsym.py logic internally, then runs
findsym as a subprocess. findsym must be on PATH and the ISODATA
environment variable must be set.

Required files in the working directory
----------------------------------------
    POSCAR      VASP input structure

Required environment
--------------------
    findsym     On PATH (e.g. $WORK/software/isotropy/findsym)
    ISODATA     Set to the isotropy data directory with trailing slash
                (e.g. /work/01253/guru/ls6/software/isotropy/)

Usage
-----
    python run_findsym.py [options]

Options
-------
    --input FILENAME        Use an existing findsym input file instead of
                            generating one from POSCAR
    --output-stem STEM      Stem for output filenames (default: findsym)
                            Produces: <stem>.log and <stem>.cif
    --lattol FLOAT          Lattice tolerance in Angstrom (default: 0.00001)
    --postol FLOAT          Atomic position tolerance in Angstrom (default: 0.001)
    --keep-input            Keep the generated findsym_input.txt (default: deleted)

Output
------
    findsym.log     Human-readable space group identification
    findsym.cif     CIF file of the identified structure

Example
-------
    python scripts/run_findsym.py
    python scripts/run_findsym.py --postol 0.01 --output-stem distorted
"""

import os
import sys
import argparse
import subprocess
import shutil

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from vasp_phonon_tools import support
from scripts.poscar_to_findsym import write_findsym_input


def check_findsym_available():
    """Exit with a clear message if findsym is not on PATH or ISODATA is unset."""
    if shutil.which('findsym') is None:
        print('Error: findsym not found on PATH.')
        print('       Add the isotropy directory to PATH, e.g.:')
        print('       export PATH=$PATH:$WORK/software/isotropy')
        sys.exit(1)

    if not os.environ.get('ISODATA'):
        print('Error: ISODATA environment variable is not set.')
        print('       Set it to the isotropy data directory with a trailing slash, e.g.:')
        print('       export ISODATA=$WORK/software/isotropy/')
        sys.exit(1)


def split_log_and_cif(stdout_text):
    """
    Split findsym stdout into the human-readable log section and the CIF section.

    findsym writes both to stdout. The CIF block begins with a line
    starting with '# CIF file'.

    Returns
    -------
    log_text : str
    cif_text : str or None
    """
    lines     = stdout_text.splitlines(keepends=True)
    cif_start = None

    for i, line in enumerate(lines):
        if line.startswith('# CIF file'):
            cif_start = i
            break

    if cif_start is None:
        return stdout_text, None

    log_text = ''.join(lines[:cif_start])
    cif_text = ''.join(lines[cif_start:])
    return log_text, cif_text


def extract_space_group(log_text):
    """
    Extract the space group line from findsym log text for summary printing.
    Returns the line as a string, or None if not found.
    """
    for line in log_text.splitlines():
        if line.startswith('Space Group:'):
            return line.strip()
    return None


def main():
    parser = argparse.ArgumentParser(
        description='Run FINDSYM on a POSCAR and write log and CIF output.'
    )
    parser.add_argument('--input', default=None,
                        help='Use an existing findsym input file instead of '
                             'generating one from POSCAR')
    parser.add_argument('--output-stem', default='findsym',
                        help='Stem for output filenames (default: findsym)')
    parser.add_argument('--lattol', type=float, default=0.00001,
                        help='Lattice tolerance in Angstrom (default: 0.00001)')
    parser.add_argument('--postol', type=float, default=0.001,
                        help='Atomic position tolerance in Angstrom (default: 0.001)')
    parser.add_argument('--keep-input', action='store_true',
                        help='Keep the generated findsym_input.txt after running')
    args = parser.parse_args()

    check_findsym_available()

    path = os.getcwd()

    # --- Generate or use existing input file ---
    generated_input = False
    if args.input:
        input_path = os.path.join(path, args.input)
        if not os.path.isfile(input_path):
            print(f'Error: specified input file not found: {input_path}')
            sys.exit(1)
        print(f'Using existing input file: {args.input}')
    else:
        poscar_path = os.path.join(path, 'POSCAR')
        if not os.path.isfile(poscar_path):
            print(f'Error: POSCAR not found in {path}')
            sys.exit(1)
        input_path = os.path.join(path, 'findsym_input.txt')
        poscar_lines = support.read_file_lines(poscar_path)
        write_findsym_input(
            poscar_lines, input_path,
            lattol=args.lattol,
            postol=args.postol
        )
        generated_input = True
        print(f'Generated: findsym_input.txt')

    # --- Run findsym ---
    print(f'Running findsym...')
    try:
        result = subprocess.run(
            ['findsym', input_path],
            capture_output=True,
            text=True
        )
    except Exception as e:
        print(f'Error running findsym: {e}')
        sys.exit(1)

    if result.returncode != 0:
        print(f'Error: findsym exited with code {result.returncode}')
        if result.stderr:
            print(result.stderr)
        sys.exit(1)

    # --- Split and write output ---
    log_text, cif_text = split_log_and_cif(result.stdout)

    log_path = os.path.join(path, f'{args.output_stem}.log')
    with open(log_path, 'w') as f:
        f.write(log_text)
    print(f'Written: {args.output_stem}.log')

    if cif_text:
        cif_path = os.path.join(path, f'{args.output_stem}.cif')
        with open(cif_path, 'w') as f:
            f.write(cif_text)
        print(f'Written: {args.output_stem}.cif')
    else:
        print('Warning: no CIF block found in findsym output.')

    # --- Print space group summary ---
    sg_line = extract_space_group(log_text)
    if sg_line:
        print()
        print(f'Result: {sg_line}')

    # --- Clean up generated input file unless --keep-input ---
    if generated_input and not args.keep_input:
        os.remove(input_path)


main()
