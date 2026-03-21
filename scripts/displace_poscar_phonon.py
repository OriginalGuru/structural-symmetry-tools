#!/usr/bin/env python3
"""
displace_poscar.py

Apply a phonon-based displacement to a POSCAR and write the distorted
structure to disk. Three modes of operation:

  Mode A  Single phonon eigenvector.
  Mode B  Explicit linear combination: --op-direction gives unnormalized
          weights over the listed phonon eigenvectors.
  Mode C  SAM-projected direction: --sam-indices triggers this mode.
          --op-direction gives an unnormalized target in the named SAM
          subspace; the script solves for the linear combination of the
          listed phonon eigenvectors that best matches it.

In all modes --amplitude sets the total displacement magnitude in Angstrom.

Required files in the working directory
----------------------------------------
    vasprun.xml         VASP DFPT output: Hessian and species masses
    POSCAR              Structure used in the DFPT run
    symmetry_basis      ISODISTORT SAM basis vectors        (Mode C only)
    symmetry_list       ISODISTORT SAM irrep labels         (Mode C only)
    (or symmetry_file combining both)

All files may be softlinks.

Usage
-----
    # Mode A: single phonon
    python displace_poscar.py --modes 31 --amplitude 0.1

    # Mode B: explicit direction over named phonons (unnormalized weights)
    python displace_poscar.py --modes 31,32,33 --op-direction 1,0,-2 --amplitude 0.1

    # Mode C: SAM-projected direction (--sam-indices triggers Mode C)
    python displace_poscar.py \\
        --modes 31,32,33 \\
        --sam-indices 19,20,21 \\
        --op-direction 1,1,0 \\
        --amplitude 0.1

Options
-------
    --modes STR         Comma-separated 1-based phonon mode indices.
    --amplitude FLOAT   Total displacement amplitude in Angstrom.
    --op-direction STR  Unnormalized direction. Mode B: one weight per
                        --modes entry. Mode C: one component per
                        --sam-indices entry.
    --sam-indices STR   Comma-separated 1-based SAM indices defining the
                        target subspace. Triggers Mode C. All indices must
                        share one irrep label.
    --sym-path PATH     Directory containing symmetry files (default: cwd).
    --vasprun PATH      Path to vasprun.xml (default: ./vasprun.xml).
    --poscar PATH       Path to POSCAR (default: ./POSCAR).

Output
------
    Mode A:  POSCAR_phonon_<N>_<freq>THz_<amp>Ang
    Mode B:  POSCAR_phonons_<N1>-<N2>-..._<amp>Ang
    Mode C:  POSCAR_<irrep>_op<dir>_<amp>Ang

    The POSCAR title line encodes all physical quantities for traceability.
    Mode C also prints the per-mode coefficients so the result can be
    reproduced exactly via Mode B.

Coordinate convention
---------------------
SAM basis vectors are assumed to be in Cartesian coordinates, consistent
with get_phonon_mode_charges.py and the ISODISTORT export convention used
in this project. No coordinate conversion is applied to SAM vectors.
"""

import os
import sys
import argparse
import numpy as np

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from vasp_phonon_tools import parsers, support


# ---------------------------------------------------------------------------
# SAM label parsing
# ---------------------------------------------------------------------------

def extract_irrep_label(label_string):
    """
    Extract the irrep token from a symmetry_list label string.

    For example, from 'P-3m1[1/2,0,0]M1-(a;b;c)[Re1:b:dsp]Eu(a)'
    returns 'M1-'.

    Parameters
    ----------
    label_string : str

    Returns
    -------
    str
    """
    kpt_end = label_string.index(']')
    after   = label_string[kpt_end + 1:]
    paren   = after.index('(')
    return after[:paren]


def get_irrep_for_sam(sam_labels, sam_index_1based):
    """
    Return the irrep label for a single SAM (1-based index).

    Parameters
    ----------
    sam_labels : list of list of str
        From support.read_symmetry_adapted_modes().
    sam_index_1based : int

    Returns
    -------
    str
    """
    tokens = sam_labels[sam_index_1based - 1]
    return extract_irrep_label(tokens[1])


def find_all_sams_for_irrep(sam_labels, irrep_label):
    """
    Return all 1-based SAM indices whose irrep label matches irrep_label.

    Parameters
    ----------
    sam_labels : list of list of str
    irrep_label : str

    Returns
    -------
    list of int
    """
    result = []
    for i, tokens in enumerate(sam_labels):
        try:
            if extract_irrep_label(tokens[1]) == irrep_label:
                result.append(i + 1)
        except (IndexError, ValueError):
            continue
    return result


# ---------------------------------------------------------------------------
# Displacement construction
# ---------------------------------------------------------------------------

def build_displacement_mode_a(displacements, mode_index):
    """
    Return the unit displacement vector for a single phonon mode.

    Parameters
    ----------
    displacements : numpy.ndarray, shape (n_modes, 3*N)
    mode_index : int
        1-based.

    Returns
    -------
    numpy.ndarray, shape (3*N,)
    """
    return displacements[mode_index - 1].copy()


def build_displacement_mode_b(displacements, mode_indices, op_direction):
    """
    Build a displacement vector as a weighted sum of phonon eigenvectors.

    op_direction gives unnormalized weights over the listed modes. The
    combined vector is returned un-normalized; displace_poscar_atoms
    handles normalization and amplitude scaling.

    Parameters
    ----------
    displacements : numpy.ndarray, shape (n_modes, 3*N)
    mode_indices : list of int
        1-based phonon mode indices.
    op_direction : list of float
        One weight per mode_index entry (unnormalized).

    Returns
    -------
    numpy.ndarray, shape (3*N,)
    """
    combined = np.zeros(displacements.shape[1])
    for idx, weight in zip(mode_indices, op_direction):
        combined += weight * displacements[idx - 1]
    return combined


def build_displacement_mode_c(displacements, mode_indices,
                               sam_basis, sam_labels,
                               sam_subspace_indices, op_direction):
    """
    Find the linear combination of degenerate phonon eigenvectors whose
    SAM projection best matches a target direction in a user-specified
    SAM subspace.

    Algorithm
    ---------
    1. Verify all sam_subspace_indices share one irrep label.
    2. Find all SAMs sharing that label (full irrep subspace).
    3. Build projection matrix P of shape (n_modes, n_full_sam) where
       P[i, j] = dot(displacement of mode i, SAM basis vector j).
    4. Construct target vector t in the full SAM subspace: op_direction
       components fill the sam_subspace_indices slots, zeros elsewhere.
    5. Solve P^T @ c = t (least-squares) for coefficients c over modes.
    6. Return the combined (un-normalized) displacement and metadata.

    SAM vectors are used as-is (Cartesian coordinates assumed).

    Parameters
    ----------
    displacements : numpy.ndarray, shape (n_modes, 3*N)
        Real-space unit displacement vectors, Cartesian.
    mode_indices : list of int
        1-based indices of the degenerate phonon modes to mix.
    sam_basis : list of list of float
        Shape (n_sam, 3*N). Cartesian coordinates.
    sam_labels : list of list of str
    sam_subspace_indices : list of int
        1-based SAM indices defining the target direction subspace.
    op_direction : list of float
        Unnormalized target direction, one component per sam_subspace_indices.

    Returns
    -------
    combined : numpy.ndarray, shape (3*N,)
        Combined displacement vector (not yet normalized or amplitude-scaled).
    coefficients : list of (int, float)
        (mode_index, raw_coefficient) pairs for reporting and reproducibility.
    irrep_label : str
    full_sam_indices : list of int
        All SAM indices belonging to the identified irrep.
    """
    # --- Verify all sam_subspace_indices share one irrep label ---
    irrep_labels_found = [get_irrep_for_sam(sam_labels, s)
                          for s in sam_subspace_indices]

    if len(set(irrep_labels_found)) != 1:
        print("Error: --sam-indices do not all share the same irrep label.")
        for s, lbl in zip(sam_subspace_indices, irrep_labels_found):
            print(f"  SAM {s}: {lbl}")
        sys.exit(1)

    irrep_label = irrep_labels_found[0]

    # --- Find all SAMs with this irrep label ---
    full_sam_indices = find_all_sams_for_irrep(sam_labels, irrep_label)
    if not full_sam_indices:
        print(f"Error: no SAMs found with irrep label '{irrep_label}'.")
        sys.exit(1)

    # --- Build SAM matrix for the full irrep subspace ---
    sam_array        = np.array(sam_basis)
    full_sam_matrix  = sam_array[[s - 1 for s in full_sam_indices]].copy()

    for i in range(len(full_sam_matrix)):
        norm = np.linalg.norm(full_sam_matrix[i])
        if norm < 1e-10:
            print(f"Error: SAM {full_sam_indices[i]} has zero norm.")
            sys.exit(1)
        full_sam_matrix[i] /= norm

    # --- Build projection matrix P: shape (n_modes_subset, n_full) ---
    n_full         = len(full_sam_indices)
    n_modes_subset = len(mode_indices)
    P = np.zeros((n_modes_subset, n_full))
    for i, midx in enumerate(mode_indices):
        u = displacements[midx - 1]
        for j in range(n_full):
            P[i, j] = np.dot(u, full_sam_matrix[j])

    # --- Build target vector in full SAM subspace ---
    full_idx_map = {s: pos for pos, s in enumerate(full_sam_indices)}
    target = np.zeros(n_full)
    for comp, s in zip(op_direction, sam_subspace_indices):
        if s not in full_idx_map:
            print(f"Error: SAM {s} is in --sam-indices but not in the full "
                  f"irrep subspace for '{irrep_label}'.")
            sys.exit(1)
        target[full_idx_map[s]] = comp

    if np.linalg.norm(target) < 1e-10:
        print("Error: --op-direction maps to a zero vector in the SAM subspace.")
        sys.exit(1)

    # --- Solve P^T @ c = target ---
    c, _, _, _ = np.linalg.lstsq(P.T, target, rcond=None)

    # --- Build combined displacement vector ---
    combined = np.zeros(sam_array.shape[1])
    for i, midx in enumerate(mode_indices):
        combined += c[i] * displacements[midx - 1]

    coefficients = list(zip(mode_indices, c.tolist()))
    return combined, coefficients, irrep_label, full_sam_indices


# ---------------------------------------------------------------------------
# Output filename and title line construction
# ---------------------------------------------------------------------------

def format_amplitude(amp):
    """Compact amplitude string without trailing zeros."""
    return f"{amp:.6f}".rstrip('0').rstrip('.')


def fmt_float(x):
    """Compact float: integer-valued floats drop the decimal point."""
    return str(int(x)) if x == int(x) else str(x)


def make_output_name_mode_a(mode_index, freq, amplitude):
    return (f"POSCAR_phonon_{mode_index}_{freq:+.3f}THz"
            f"_{format_amplitude(amplitude)}Ang")


def make_output_name_mode_b(mode_indices, amplitude):
    indices_str = "-".join(str(m) for m in mode_indices)
    return f"POSCAR_phonons_{indices_str}_{format_amplitude(amplitude)}Ang"


def make_output_name_mode_c(irrep_label, op_direction, amplitude):
    dir_str    = "".join(fmt_float(x) for x in op_direction)
    irrep_safe = irrep_label.replace('-', 'm').replace('+', 'p')
    return f"POSCAR_{irrep_safe}_op{dir_str}_{format_amplitude(amplitude)}Ang"


def make_title_mode_a(original_title, mode_index, freq, amplitude):
    return (f"{original_title.strip()} | "
            f"phonon mode: {mode_index} ({freq:+.6f} THz) | "
            f"amplitude: {format_amplitude(amplitude)} Ang")


def make_title_mode_b(original_title, mode_indices, op_direction,
                      frequencies, amplitude):
    parts = [f"{m}({frequencies[m-1]:+.3f}THz):{fmt_float(w)}"
             for m, w in zip(mode_indices, op_direction)]
    return (f"{original_title.strip()} | "
            f"phonon modes: {', '.join(parts)} | "
            f"amplitude: {format_amplitude(amplitude)} Ang")


def make_title_mode_c(original_title, irrep_label, mode_indices, frequencies,
                      coefficients, full_sam_indices, sam_subspace_indices,
                      op_direction, amplitude):
    dir_str     = ",".join(fmt_float(x) for x in op_direction)
    sams_str    = ",".join(str(s) for s in full_sam_indices)
    sub_str     = ",".join(str(s) for s in sam_subspace_indices)
    coeff_parts = [f"{midx}({frequencies[midx-1]:+.3f}THz):{coeff:+.6f}"
                   for midx, coeff in coefficients]
    return (f"{original_title.strip()} | "
            f"irrep: {irrep_label} | "
            f"phonon modes: {','.join(str(m) for m in mode_indices)} | "
            f"SAMs: {sams_str} | "
            f"op-direction: {dir_str} (in SAMs {sub_str}) | "
            f"amplitude: {format_amplitude(amplitude)} Ang | "
            f"coefficients: {', '.join(coeff_parts)}")


# ---------------------------------------------------------------------------
# Argument parsing and mode determination
# ---------------------------------------------------------------------------

def parse_args():
    parser = argparse.ArgumentParser(
        description='Displace a POSCAR along a phonon eigenvector or '
                    'order-parameter direction.'
    )
    parser.add_argument('--modes', required=True,
                        help='Comma-separated 1-based phonon mode indices.')
    parser.add_argument('--amplitude', type=float, required=True,
                        help='Total displacement amplitude in Angstrom.')
    parser.add_argument('--op-direction', default=None,
                        help='Unnormalized direction. Mode B: one weight per '
                             '--modes entry. Mode C: one component per '
                             '--sam-indices entry.')
    parser.add_argument('--sam-indices', default=None,
                        help='Comma-separated 1-based SAM indices. '
                             'Triggers Mode C.')
    parser.add_argument('--sym-path', default=None,
                        help='Directory containing symmetry files (default: cwd).')
    parser.add_argument('--vasprun', default=None,
                        help='Path to vasprun.xml (default: ./vasprun.xml).')
    parser.add_argument('--poscar', default=None,
                        help='Path to POSCAR (default: ./POSCAR).')
    return parser.parse_args()


def determine_mode(args):
    """
    Return mode string ('A', 'B', or 'C') and parsed lists.
    Exits with a clear message on invalid combinations.

    Returns
    -------
    mode : str
    mode_indices : list of int
    sam_indices : list of int or None
    op_direction : list of float or None
    """
    mode_indices = [int(x) for x in args.modes.split(',')]
    has_op_dir   = args.op_direction is not None
    has_sam_idx  = args.sam_indices  is not None

    # --- Mode C ---
    if has_sam_idx:
        if not has_op_dir:
            print("Error: --sam-indices requires --op-direction.")
            sys.exit(1)
        if len(mode_indices) < 2:
            print("Error: Mode C (--sam-indices) requires at least two --modes.")
            sys.exit(1)
        sam_indices  = [int(x)   for x in args.sam_indices.split(',')]
        op_direction = [float(x) for x in args.op_direction.split(',')]
        if len(op_direction) != len(sam_indices):
            print(f"Error: --op-direction has {len(op_direction)} components "
                  f"but --sam-indices has {len(sam_indices)} entries.")
            sys.exit(1)
        return 'C', mode_indices, sam_indices, op_direction

    # --- Mode B ---
    if has_op_dir:
        if len(mode_indices) < 2:
            print("Error: --op-direction with multiple modes requires at least "
                  "two entries in --modes.")
            sys.exit(1)
        op_direction = [float(x) for x in args.op_direction.split(',')]
        if len(op_direction) != len(mode_indices):
            print(f"Error: --op-direction has {len(op_direction)} components "
                  f"but --modes has {len(mode_indices)} entries.")
            sys.exit(1)
        return 'B', mode_indices, None, op_direction

    # --- Mode A ---
    if len(mode_indices) != 1:
        print("Error: multiple --modes supplied without --op-direction. "
              "Provide --op-direction weights or add --sam-indices.")
        sys.exit(1)
    return 'A', mode_indices, None, None


def resolve_paths(args):
    cwd          = os.getcwd()
    vasprun_path = args.vasprun  or os.path.join(cwd, 'vasprun.xml')
    poscar_path  = args.poscar   or os.path.join(cwd, 'POSCAR')
    sym_path     = args.sym_path or cwd
    for p in (vasprun_path, poscar_path):
        if not os.path.isfile(p):
            print(f"Error: required file not found: {p}")
            sys.exit(1)
    return vasprun_path, poscar_path, sym_path


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

def main():
    args = parse_args()
    mode, mode_indices, sam_indices, op_direction = determine_mode(args)
    vasprun_path, poscar_path, sym_path = resolve_paths(args)

    # --- Load structure and solve dynamical matrix ---
    poscar_lines   = support.read_file_lines(poscar_path)
    atom_counts    = support.get_atom_counts(poscar_lines)
    species_masses = parsers.parse_species_masses(vasprun_path)
    hessian        = parsers.get_dynamical_matrix(vasprun_path)

    frequencies, eigenvectors = support.solve_dynamical_eigenmodes(hessian)
    n_modes = len(frequencies)

    for m in mode_indices:
        if m < 1 or m > n_modes:
            print(f"Error: mode index {m} out of range (1–{n_modes}).")
            sys.exit(1)

    _, displacements = support.get_real_space_displacements(
        atom_counts, species_masses, eigenvectors
    )

    original_title = poscar_lines[0].rstrip('\n')
    amplitude      = args.amplitude

    # -----------------------------------------------------------------------
    # Mode A: single phonon
    # -----------------------------------------------------------------------
    if mode == 'A':
        m    = mode_indices[0]
        freq = frequencies[m - 1]
        disp = build_displacement_mode_a(displacements, m)

        label    = make_title_mode_a(original_title, m, freq, amplitude)
        out_name = make_output_name_mode_a(m, freq, amplitude)

        distorted = support.displace_poscar_atoms(
            poscar_lines, disp, amplitude, label
        )
        support.write_poscar(distorted, out_name)

        print(f"Written: {out_name}")
        print(f"  Phonon mode {m}: {freq:+.6f} THz")
        print(f"  Amplitude: {amplitude} Ang")

    # -----------------------------------------------------------------------
    # Mode B: explicit weighted combination
    # -----------------------------------------------------------------------
    elif mode == 'B':
        disp = build_displacement_mode_b(displacements, mode_indices, op_direction)
        if np.linalg.norm(disp) < 1e-10:
            print("Error: combined displacement vector has zero norm. "
                  "Check --op-direction weights.")
            sys.exit(1)

        label    = make_title_mode_b(original_title, mode_indices,
                                     op_direction, frequencies, amplitude)
        out_name = make_output_name_mode_b(mode_indices, amplitude)

        distorted = support.displace_poscar_atoms(
            poscar_lines, disp, amplitude, label
        )
        support.write_poscar(distorted, out_name)

        print(f"Written: {out_name}")
        print(f"  Amplitude: {amplitude} Ang")
        for m, w in zip(mode_indices, op_direction):
            print(f"  Mode {m} ({frequencies[m-1]:+.6f} THz): weight {fmt_float(w)}")

    # -----------------------------------------------------------------------
    # Mode C: SAM-projected OP direction
    # -----------------------------------------------------------------------
    elif mode == 'C':
        sam_basis, sam_labels = support.read_symmetry_adapted_modes(sym_path)

        for s in sam_indices:
            if s < 1 or s > len(sam_labels):
                print(f"Error: SAM index {s} out of range (1–{len(sam_labels)}).")
                sys.exit(1)

        disp, coefficients, irrep_label, full_sam_indices = \
            build_displacement_mode_c(
                displacements, mode_indices,
                sam_basis, sam_labels,
                sam_indices, op_direction
            )

        if np.linalg.norm(disp) < 1e-10:
            print("Error: resulting displacement vector has zero norm. "
                  "Check that --modes span the requested OP direction.")
            sys.exit(1)

        label    = make_title_mode_c(
            original_title, irrep_label, mode_indices, frequencies,
            coefficients, full_sam_indices, sam_indices, op_direction, amplitude
        )
        out_name = make_output_name_mode_c(irrep_label, op_direction, amplitude)

        distorted = support.displace_poscar_atoms(
            poscar_lines, disp, amplitude, label
        )
        support.write_poscar(distorted, out_name)

        print(f"Written: {out_name}")
        print(f"  Irrep: {irrep_label}")
        print(f"  Full SAM subspace: {full_sam_indices}")
        print(f"  OP direction {op_direction} in SAMs {sam_indices}")
        print(f"  Amplitude: {amplitude} Ang")
        print(f"  Per-mode coefficients (reproducible via Mode B --op-direction):")
        for midx, coeff in coefficients:
            print(f"    mode {midx} ({frequencies[midx-1]:+.6f} THz): {coeff:+.10f}")


main()
