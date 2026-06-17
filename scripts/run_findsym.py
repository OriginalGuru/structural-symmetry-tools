#!/usr/bin/env python3
"""
run_findsym.py

Generate a FINDSYM input file from a POSCAR and run findsym, writing
a human-readable log and a CIF file.

findsym must be on PATH and ISODATA must be set. The script changes
into the directory containing the POSCAR before calling findsym, so
that the input path passed to the Fortran executable stays short.

Usage
-----
    python run_findsym.py [POSCAR] [options]

Arguments
---------
    POSCAR                  Path to POSCAR file (default: ./POSCAR)

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

    Both files are written to the same directory as the POSCAR.

Examples
--------
    # Run on POSCAR in current directory
    python scripts/run_findsym.py

    # Run on a displaced POSCAR
    python scripts/run_findsym.py POSCAR_phonon_30_-0.086THz_0.1Ang

    # Looser tolerance, custom output stem
    python scripts/run_findsym.py POSCAR_sam_3_0.1Ang --postol 0.01 --output-stem distorted_sam3
"""

import os
import sys
import argparse
import subprocess
import shutil

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from vasp_phonon_tools import support

# Import write_findsym_input directly — avoids module path issues
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from poscar_to_findsym import write_findsym_input


def check_findsym_available():
    if shutil.which('findsym') is None:
        print('Error: findsym not found on PATH.')
        print('       export PATH=$PATH:$WORK/software/isotropy')
        sys.exit(1)
    if not os.environ.get('ISODATA'):
        print('Error: ISODATA environment variable is not set.')
        print('       export ISODATA=$WORK/software/isotropy/')
        sys.exit(1)


def split_log_and_cif(stdout_text):
    """Split findsym stdout into log section and CIF section."""
    lines     = stdout_text.splitlines(keepends=True)
    cif_start = None
    for i, line in enumerate(lines):
        if line.startswith('# CIF file'):
            cif_start = i
            break
    if cif_start is None:
        return stdout_text, None
    return ''.join(lines[:cif_start]), ''.join(lines[cif_start:])


def extract_space_group(log_text):
    for line in log_text.splitlines():
        if line.startswith('Space Group:'):
            return line.strip()
    return None


def main():
    parser = argparse.ArgumentParser(
        description='Run FINDSYM on a POSCAR and write log and CIF output.'
    )
    parser.add_argument('poscar', nargs='?', default='POSCAR',
                        help='Path to POSCAR file (default: ./POSCAR)')
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

    poscar_path = os.path.abspath(args.poscar)
    if not os.path.isfile(poscar_path):
        print(f'Error: POSCAR file not found: {poscar_path}')
        sys.exit(1)

    # All output goes alongside the POSCAR
    work_dir = os.path.dirname(poscar_path)

    # --- Generate or locate input file ---
    generated_input = False
    if args.input:
        input_filename = args.input
        input_path     = os.path.join(work_dir, input_filename)
        if not os.path.isfile(input_path):
            print(f'Error: specified input file not found: {input_path}')
            sys.exit(1)
        print(f'Using existing input file: {input_filename}')
    else:
        input_filename = 'findsym_input.txt'
        input_path     = os.path.join(work_dir, input_filename)
        poscar_lines   = support.read_file_lines(poscar_path)
        write_findsym_input(poscar_lines, input_path,
                            lattol=args.lattol, postol=args.postol)
        generated_input = True
        print(f'Generated: {input_filename}')

    # --- Run findsym from the work directory so the path stays short ---
    print('Running findsym...')
    try:
        result = subprocess.run(
            ['findsym', input_filename],
            capture_output=True,
            text=True,
            cwd=work_dir          # <-- key fix: short filename, correct cwd
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

    log_path = os.path.join(work_dir, f'{args.output_stem}.log')
    with open(log_path, 'w') as f:
        f.write(log_text)
    print(f'Written: {log_path}')

    if cif_text:
        cif_path = os.path.join(work_dir, f'{args.output_stem}.cif')
        with open(cif_path, 'w') as f:
            f.write(cif_text)
        print(f'Written: {cif_path}')
    else:
        print('Warning: no CIF block found in findsym output.')

    sg_line = extract_space_group(log_text)
    if sg_line:
        print()
        print(f'Result: {sg_line}')

    if generated_input and not args.keep_input:
        os.remove(input_path)


main()
