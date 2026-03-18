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
    python get_phonon_mode_charges.py

Output
------
    Printed to stdout. Redirect to a file with:
        python get_phonon_mode_charges.py > output.txt
"""

import os
import sys
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


def main():
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

    COL = "%+16.12f"

    def print_header():
        bec_note = "  [Z* approximate: rigid-ion]" if bec_approximate else ""
        print(
            "Phonon  "
            "  Frequency       "
            "  Force Constant  "
            "  Reduced Mass    "
            "  Z*_x            "
            "  Z*_y            "
            "  Z*_z            "
            "  |Z*|"
            + bec_note
        )
        print(
            "Mode    "
            "  (THz)           "
            "  (eV/A2)         "
            "  (AMU)           "
            "  (e)             "
            "  (e)             "
            "  (e)             "
            "  (e)"
        )

    # Pre-compute digit widths for index column alignment
    predecimal = max(len(str(int(row[0]))) for row in output_list)

    def print_phonon_row(row):
        idx = "+" + str(row[0]).rjust(predecimal, "0")
        values = "  ".join(COL % row[k] for k in range(1, 8))
        print(idx + "  " + values)

    # Summary table
    print_header()
    for row in output_list:
        print_phonon_row(row)

    print()
    print()

    # Detailed table with SAM decomposition per mode
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


main()
