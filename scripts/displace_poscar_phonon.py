#!/usr/bin/env python3
"""
displace_poscar_phonon.py

Apply a phonon eigenvector displacement to a POSCAR and write the
distorted structure to disk.

Required files in the working directory
----------------------------------------
    vasprun.xml     VASP output: Hessian and species masses
    POSCAR          VASP input: structure used in the DFPT run

Usage
-----
    python displace_poscar_phonon.py <mode_index> <amplitude>

Arguments
---------
    mode_index  1-based phonon mode index (highest frequency = 1)
    amplitude   Displacement amplitude in Angstrom (may be negative)

Output
------
    POSCAR_phonon_<mode_index>_<amplitude>Ang
    Written to the current working directory.

Example
-------
    python displace_poscar_phonon.py 5 0.1
    -> POSCAR_phonon_5_0.1Ang
"""

import os
import sys
import numpy as np

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from vasp_phonon_tools import parsers, support


def main():
    if len(sys.argv) != 3:
        print("Usage: python displace_poscar_phonon.py <mode_index> <amplitude>")
        print("  mode_index : 1-based phonon mode index (highest frequency = 1)")
        print("  amplitude  : displacement in Angstrom")
        sys.exit(1)

    try:
        mode_index = int(sys.argv[1])
        amplitude  = float(sys.argv[2])
    except ValueError:
        print("Error: mode_index must be an integer and amplitude must be a float.")
        sys.exit(1)

    path         = os.getcwd()
    vasprun_path = os.path.join(path, "vasprun.xml")
    poscar_path  = os.path.join(path, "POSCAR")

    for f in (vasprun_path, poscar_path):
        if not os.path.isfile(f):
            print(f"Error: required file not found: {f}")
            sys.exit(1)

    poscar_lines   = support.read_file_lines(poscar_path)
    atom_counts    = support.get_atom_counts(poscar_lines)
    species_masses = parsers.parse_species_masses(vasprun_path)
    hessian        = parsers.get_dynamical_matrix(vasprun_path)

    frequencies, eigenvectors = support.solve_dynamical_eigenmodes(hessian)
    n_modes = len(frequencies)

    if mode_index < 1 or mode_index > n_modes:
        print(f"Error: mode_index {mode_index} out of range. "
              f"Valid range: 1 to {n_modes}.")
        sys.exit(1)

    _, displacements = support.get_real_space_displacements(
        atom_counts, species_masses, eigenvectors
    )

    # mode_index is 1-based; displacements rows are ordered highest-frequency first
    disp_vector = displacements[mode_index - 1]

    freq = frequencies[mode_index - 1]
    label = f"phonon mode {mode_index} ({freq:+.6f} THz) | amplitude: {amplitude} Ang"
    distorted = support.displace_poscar_atoms(
        poscar_lines, disp_vector, amplitude, label
    )

    out_name = f"POSCAR_phonon_{mode_index}_{freq:+.3f}THz_{amplitude}Ang"
    support.write_poscar(distorted, out_name)
    print(f"Written: {out_name}")
    print(f"  Mode {mode_index}: {freq:+.6f} THz")
    print(f"  Amplitude: {amplitude} Ang")


main()
