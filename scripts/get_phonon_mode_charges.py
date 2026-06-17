#!/usr/bin/env python3
"""
get_phonon_mode_charges.py

Compute phonon frequencies, reduced masses, mode effective charges,
and symmetry-adapted mode decomposition from a VASP DFPT calculation.

Required files in the working directory
----------------------------------------
    vasprun.xml     VASP output: Hessian, Born effective charges (optional),
                    species masses and ZVAL
    POSCAR          VASP input: structure used in the DFPT run
    symmetry_basis  ISODISTORT output: symmetry-adapted mode basis vectors
    symmetry_list   ISODISTORT output: irreducible representation labels

Usage
-----
    python get_phonon_mode_charges.py [-s] [-d N] [-n]

Options
-------
    -s, --summary-only   Print only the summary table; omit the per-mode
                         SAM decomposition blocks.
    -d, --decimals N     Number of digits after the decimal point in the
                         numeric columns (default: 12).
    -n, --no-symmetry    Omit the Symmetry / SAM Indices columns from the
                         summary table.

Output
------
    Printed to stdout. Redirect to a file with:
        python get_phonon_mode_charges.py > output.txt
"""

import os
import sys
import math
import argparse
import numpy as np

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from vasp_phonon_tools import parsers, support


def check_required_files(path):
    """Exit with a clear message if any required input file is missing."""
    required = {
        "vasprun.xml":    os.path.isfile(os.path.join(path, "vasprun.xml")),
        "POSCAR":         os.path.isfile(os.path.join(path, "POSCAR")),
        "symmetry_basis": os.path.isfile(os.path.join(path, "symmetry_basis")),
        "symmetry_list":  os.path.isfile(os.path.join(path, "symmetry_list")),
    }
    missing = [k for k, present in required.items() if not present]
    if missing:
        print("\nError: missing required input files:")
        for name in missing:
            print(f"  {name}")
        print("\nSee README.md for file format details.")
        sys.exit(1)


def strip_irrep(label):
    """Reduce a full SAM label to its irrep token.

    'P-3m1[1/2,0,0]M2-(a;b;c)[Re1:b:dsp]A2u(c)' -> 'M2-'
    'P-3m1[0,0,0]GM1+(a)[S1:d:dsp]A1(a)'        -> 'GM1+'

    The irrep is the token after the k-point bracket and before its
    order-parameter direction '(...)'.
    """
    return label.split("]", 1)[1].split("(", 1)[0]


def mode_symmetry(row):
    """Per-mode symmetry assignment from a phonon's SAM decomposition.

    Uses the already-filtered SAM indices (row[8]) and their labels
    (row[9]; row[9][k] corresponds to row[8][k]). Indices are grouped
    by irrep, preserving first-appearance order.

    Returns
    -------
    irreps  : 'M2-'              or  'GM1+; M2-'
    indices : '27, 30, 33, 36'   or  '1; 27, 30, 33, 36'
    """
    groups, pos = [], {}
    for k, sam_idx in enumerate(row[8]):
        irrep = strip_irrep(row[9][k][1])
        if irrep not in pos:
            pos[irrep] = len(groups)
            groups.append([irrep, []])
        groups[pos[irrep]][1].append(sam_idx)
    irreps  = "; ".join(g[0] for g in groups)
    indices = "; ".join(", ".join(str(i) for i in g[1]) for g in groups)
    return irreps, indices


def parse_args():
    p = argparse.ArgumentParser(
        description="Phonon frequencies, mode effective charges, and "
                    "symmetry-adapted mode decomposition from a VASP DFPT run."
    )
    p.add_argument(
        "-s", "--summary-only", action="store_true",
        help="Print only the summary table; omit the per-mode "
             "SAM decomposition blocks.",
    )
    p.add_argument(
        "-d", "--decimals", type=int, default=12, metavar="N",
        help="Digits after the decimal point in the numeric columns "
             "(default: 12).",
    )
    p.add_argument(
        "-n", "--no-symmetry", action="store_true",
        help="Omit the Symmetry / SAM Indices columns from the summary table.",
    )
    return p.parse_args()


def main(args):
    path = os.getcwd()
    vasprun_path = os.path.join(path, "vasprun.xml")
    poscar_path  = os.path.join(path, "POSCAR")

    check_required_files(path)

    # --- Read structure ---
    poscar_lines = support.read_file_lines(poscar_path)
    atom_counts  = support.get_atom_counts(poscar_lines)

    # --- Read from vasprun.xml ---
    species_masses = parsers.parse_species_masses(vasprun_path)
    species_zvals  = parsers.parse_species_zvals(vasprun_path)
    hessian        = parsers.get_dynamical_matrix(vasprun_path)

    # --- Born effective charges: use calculated values or fall back to ZVAL ---
    bec = parsers.get_born_effective_charges(vasprun_path)
    bec_approximate = False
    if bec is None:
        print("\nWarning: born_charges not found in vasprun.xml.")
        print("         Constructing diagonal BEC from ZVAL (rigid-ion approximation).")
        print("         Phonon frequencies and eigenvectors are unaffected.")
        print("         Mode effective charges (Z*) will be approximate.\n")
        bec = support.build_bec_from_zvals(species_zvals, atom_counts)
        bec_approximate = True

    # --- Phonon eigenmodes ---
    frequencies, eigenvectors = support.solve_dynamical_eigenmodes(hessian)
    reduced_masses, displacements = support.get_real_space_displacements(
        atom_counts, species_masses, eigenvectors
    )

    # --- Mode effective charges ---
    mec = support.get_mode_effective_charges(atom_counts, bec, displacements)

    # --- Symmetry-adapted mode decomposition ---
    sam_basis, sam_labels = support.read_symmetry_adapted_modes(path)
    phonon_decomp = support.decompose_phonons_into_sam(sam_basis, eigenvectors)
    output_list = support.assign_phonon_symmetry(
        phonon_decomp, sam_labels, frequencies, reduced_masses, mec
    )

    # -------------------------------------------------------------------
    # Output
    # -------------------------------------------------------------------

    COL_DECIMALS = max(0, args.decimals)
    SEP   = "  "
    SYM_W = 14         # width of the irrep column in the summary

    NUM_HEADERS = [
        ("Frequency",      "(THz)"),
        ("Force Constant", "(eV/A2)"),
        ("Reduced Mass",   "(AMU)"),
        ("Z*_x",           "(e)"),
        ("Z*_y",           "(e)"),
        ("Z*_z",           "(e)"),
        ("|Z*|",           "(e)"),
    ]

    # Numeric field formatting. Width adapts to the requested decimal
    # precision and to the largest integer part actually present, and is
    # never narrower than the widest column header.
    int_digits = max(
        (len(str(int(abs(v)))) for row in output_list for v in row[1:8]),
        default=1,
    )
    content = 1 + int_digits + (1 + COL_DECIMALS if COL_DECIMALS > 0 else 0)
    NUM_W = max(content + 1, max(len(name) for name, _ in NUM_HEADERS))
    COL = "%+{}.{}f".format(NUM_W, COL_DECIMALS)

    # Pre-compute digit width for the index column.
    predecimal = max(len(str(int(row[0]))) for row in output_list)
    IDX_W = max(len("Phonon"), predecimal + 1)

    def fmt_idx(row):
        return ("+" + str(row[0]).rjust(predecimal, "0")).ljust(IDX_W)

    def print_header(with_symmetry=False):
        line1 = "Phonon".ljust(IDX_W)
        line2 = "Mode".ljust(IDX_W)
        for name, unit in NUM_HEADERS:
            line1 += SEP + name.ljust(NUM_W)
            line2 += SEP + unit.ljust(NUM_W)
        if with_symmetry:
            line1 += SEP + "Symmetry".ljust(SYM_W) + SEP + "SAM Indices"
            line2 += SEP + "(irrep)".ljust(SYM_W)  + SEP + "(index list)"
        print(line1.rstrip())
        print(line2.rstrip())

    def print_phonon_row(row, with_symmetry=False):
        # Force constant carries the sign of the frequency, mirroring the
        # imaginary-frequency-as-negative convention used elsewhere.
        fc = math.copysign(row[2], row[1])
        values = [row[1], fc, row[3], row[4], row[5], row[6], row[7]]
        line = fmt_idx(row) + "".join(SEP + (COL % v) for v in values)
        if with_symmetry:
            irreps, indices = mode_symmetry(row)
            line += SEP + irreps.ljust(SYM_W) + SEP + indices
        print(line)

    # Summary table (with per-mode symmetry assignment columns)
    show_symmetry = not args.no_symmetry
    if bec_approximate:
        print("[Z* approximate: rigid-ion]")
    print_header(with_symmetry=show_symmetry)
    for row in output_list:
        print_phonon_row(row, with_symmetry=show_symmetry)

    print()
    print()

    # Detailed table with SAM decomposition per mode
    if not args.summary_only:
        for i, row in enumerate(output_list):
            print_header()
            print_phonon_row(row)
            print()
            print()
            print("\t    Symmetry Adapted Mode Decomposition (amplitudes below 10^-2 suppressed)")
            print("\t    Mode Index" + "".rjust(predecimal + 4) +
                  "Amplitude         Irreducible Representation")
            for sam_idx in row[8]:
                amp = phonon_decomp[i][sam_idx - 1]
                sign = " " if amp >= 0 else ""
                print(f"\t    {str(sam_idx).rjust(predecimal, '0')}"
                      f"            {sign}{amp:.12f}"
                      f"   {row[9][row[8].index(sam_idx)][1] if sam_idx in row[8] else ''}")
            print()
            print()


main(parse_args())
