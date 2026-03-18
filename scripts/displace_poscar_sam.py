#!/usr/bin/env python3
"""
displace_poscar_sam.py

Apply a symmetry-adapted mode (SAM) displacement to a POSCAR and write
the distorted structure to disk. Supports single SAMs or linear
combinations of SAMs.

Required files in the working directory
----------------------------------------
    POSCAR          VASP input structure
    symmetry_basis  ISODISTORT output: SAM basis vectors
    symmetry_list   ISODISTORT output: irreducible representation labels

COORDINATE NOTE: This script assumes SAM vectors in symmetry_basis are
in CARTESIAN coordinates. If your ISODISTORT export used Direct
(fractional) coordinates, set --direct on the command line and the
POSCAR lattice vectors will be used to convert automatically.

Usage
-----
Single SAM:
    python displace_poscar_sam.py <sam_index> <amplitude> [--direct]

Linear combination (comma-separated amplitude:index pairs):
    python displace_poscar_sam.py "1.0:3,1.0:5" <amplitude> [--direct]

Arguments
---------
    sam_index   1-based SAM index, or a quoted combination string
    amplitude   Displacement amplitude in Angstrom
    --direct    Flag: treat symmetry_basis vectors as Direct coordinates

Output
------
    POSCAR_sam_<label>_<amplitude>Ang
    Written to the current working directory.

Examples
--------
    python displace_poscar_sam.py 3 0.2
    python displace_poscar_sam.py "1.0:3,1.0:5" 0.2
    python displace_poscar_sam.py 3 0.2 --direct
"""

import os
import sys
import numpy as np

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from vasp_phonon_tools import support


def parse_sam_selection(selection_str, n_sam):
    """
    Parse a SAM selection string into a weight vector.

    Single index:   "3"         -> weight 1.0 on SAM 3
    Combination:    "1.0:3,0.5:5" -> weights on SAMs 3 and 5

    Returns
    -------
    numpy.ndarray, shape (n_sam,)
    """
    weights = np.zeros(n_sam)
    if ":" in selection_str:
        for part in selection_str.split(","):
            part = part.strip()
            amp_str, idx_str = part.split(":")
            amp = float(amp_str)
            idx = int(idx_str)
            if idx < 1 or idx > n_sam:
                print(f"Error: SAM index {idx} out of range (1 to {n_sam}).")
                sys.exit(1)
            weights[idx - 1] = amp
    else:
        idx = int(selection_str)
        if idx < 1 or idx > n_sam:
            print(f"Error: SAM index {idx} out of range (1 to {n_sam}).")
            sys.exit(1)
        weights[idx - 1] = 1.0
    return weights


def main():
    if len(sys.argv) < 3:
        print("Usage: python displace_poscar_sam.py <sam_index|combination> "
              "<amplitude> [--direct]")
        sys.exit(1)

    selection_str = sys.argv[1]
    try:
        amplitude = float(sys.argv[2])
    except ValueError:
        print("Error: amplitude must be a float.")
        sys.exit(1)

    use_direct = "--direct" in sys.argv

    path        = os.getcwd()
    poscar_path = os.path.join(path, "POSCAR")

    if not os.path.isfile(poscar_path):
        print(f"Error: POSCAR not found in {path}")
        sys.exit(1)

    poscar_lines = support.read_file_lines(poscar_path)
    sam_basis, sam_labels = support.read_symmetry_adapted_modes(path)
    n_sam = len(sam_basis)

    weights = parse_sam_selection(selection_str, n_sam)

    # Build the combined displacement vector
    sam_array = np.array(sam_basis)
    disp_vector = weights @ sam_array  # shape (3*N,)

    if np.linalg.norm(disp_vector) == 0:
        print("Error: resulting displacement vector is zero.")
        sys.exit(1)

    # Convert from Direct to Cartesian if needed
    if use_direct:
        lattice = support.get_lattice_vectors(poscar_lines)
        disp_vector = support.direct_to_cartesian_displacement(disp_vector, lattice)
        coord_note = "_direct-converted"
    else:
        coord_note = ""

    # Build filename label — index only, shell-safe
    if ":" in selection_str:
        indices = [part.split(':')[1].strip() for part in selection_str.split(',')]
        index_label = '+'.join(indices)
    else:
        index_label = selection_str.strip()

    direct_note = "_direct" if use_direct else ""
    filename_label = f"{index_label}{direct_note}"

    # Build POSCAR title label — full human-readable detail
    if ":" in selection_str:
        parts = selection_str.split(',')
        combo_str = ' + '.join(
            f"SAM {p.split(':')[1].strip()} (w={float(p.split(':')[0]):.3g})"
            for p in parts
        )
        title_label = f"{combo_str} | amplitude: {amplitude} Ang"
    else:
        idx = int(selection_str.strip())
        irrep = sam_labels[idx - 1][1] if len(sam_labels[idx - 1]) > 1 \
                else sam_labels[idx - 1][0]
        title_label = f"SAM {idx}: {irrep} | amplitude: {amplitude} Ang"
    if use_direct:
        title_label += " | Direct->Cartesian converted"

    distorted = support.displace_poscar_atoms(
        poscar_lines, disp_vector, amplitude, title_label
    )

    out_name = f"POSCAR_sam_{filename_label}_{amplitude}Ang"
    support.write_poscar(distorted, out_name)
    print(f"Written: {out_name}")
    if ":" not in selection_str:
        idx = int(selection_str)
        print(f"  SAM {idx}: {sam_labels[idx - 1]}")
    print(f"  Amplitude: {amplitude} Ang")
    if use_direct:
        print("  Coordinate conversion: Direct -> Cartesian applied.")


main()
