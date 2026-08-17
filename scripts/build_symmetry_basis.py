"""
build_symmetry_basis.py
=======================
Extract symmetry-adapted mode (SAM) basis vectors from an ISODISTORT output
file and a VASP POSCAR, producing two output files:

  symmetry_basis  -- N_modes x N_dof matrix (one SAM per row), unit-normalized,
                     in Cartesian coordinates, in POSCAR atom order.
  symmetry_list   -- ordered list of SAM labels (one per line: index label).

Usage
-----
  python build_symmetry_basis.py <iso_output.txt> <POSCAR> [options]

Options
-------
  --basis  PATH   Output path for symmetry_basis  (default: symmetry_basis)
  --list   PATH   Output path for symmetry_list   (default: symmetry_list)
  --tol    FLOAT  Fractional-coordinate tolerance for atom matching (default: 0.02)
  --skip-strain   Exclude strain modes (default: True)

Output format
-------------
  symmetry_basis : plain text, one row per SAM, space-separated floats.
                   Each row is a unit-normalized Cartesian displacement vector
                   in POSCAR DOF order (atom0_x, atom0_y, atom0_z, atom1_x, ...).

  symmetry_list  : plain text, one line per SAM:
                   "<1-based index>  <full ISODISTORT label>"

Normalization convention
------------------------
  Each SAM vector is normalized to unit length. With unit SAMs and VASP
  eigenvectors that are also unit-normalized, the decomposition amplitudes
  are direct projections and their squared sum equals 1 for a complete basis.

Coordinate convention
---------------------
  ISODISTORT gives displacement vectors in fractional coordinates of the
  superstructure cell. This script converts them to Cartesian by multiplying
  each atom's (dx, dy, dz) by the lattice matrix A (rows = lattice vectors).
"""

import argparse
import sys
import numpy as np


# ---------------------------------------------------------------------------
# Parsing
# ---------------------------------------------------------------------------

def parse_poscar(path):
    """
    Read a VASP POSCAR file.

    Returns
    -------
    A : (3, 3) array
        Lattice matrix, rows are lattice vectors (Angstrom).
    positions : (N, 3) array
        Fractional coordinates of all atoms.
    """
    with open(path) as f:
        lines = f.readlines()

    scale = float(lines[1])
    a1 = [float(x) for x in lines[2].split()]
    a2 = [float(x) for x in lines[3].split()]
    a3 = [float(x) for x in lines[4].split()]
    A = np.array([a1, a2, a3]) * scale

    # Handle both old (no species line) and new POSCAR formats
    try:
        counts = [int(x) for x in lines[6].split()]
        coord_line = 7
    except ValueError:
        counts = [int(x) for x in lines[7].split()]
        coord_line = 8

    n_atoms = sum(counts)

    # Skip "Direct" / "Cartesian" header line
    coord_line += 1

    positions = []
    for i in range(n_atoms):
        xyz = [float(x) for x in lines[coord_line + i].split()[:3]]
        positions.append(xyz)

    return A, np.array(positions)


def parse_isotropy(path, skip_strain=True):
    """
    Parse an ISODISTORT complete-mode-details output file.

    Returns
    -------
    lattice_params : dict with keys a, b, c, alpha, beta, gamma
    superstructure_positions : (N, 3) array of fractional coords (undistorted)
    modes : list of dicts, each with keys:
        'label'   : full ISODISTORT label string
        'normfactor' : float
        'atoms'   : list of (atom_name, x, y, z, dx, dy, dz)
    """
    with open(path) as f:
        lines = f.readlines()

    # -- Lattice parameters from superstructure block -------------------------
    lattice_params = None
    superstructure_positions = []
    in_superstructure = False

    for i, line in enumerate(lines):
        if 'Undistorted superstructure' in line:
            in_superstructure = True
            # Next line has lattice params
            lp_line = lines[i + 1]
            lp = {}
            for token in lp_line.split(','):
                k, v = token.strip().split('=')
                lp[k.strip()] = float(v)
            lattice_params = lp
            continue

        if in_superstructure:
            if line.strip().startswith('atom'):
                continue  # header
            if line.strip() == '' or 'Distorted' in line:
                in_superstructure = False
                continue
            parts = line.split()
            if len(parts) >= 5:
                try:
                    superstructure_positions.append(
                        (parts[0],
                         float(parts[2]), float(parts[3]), float(parts[4]))
                    )
                except ValueError:
                    pass

    # -- Mode blocks ----------------------------------------------------------
    modes = []
    current_label = None
    current_atoms = []
    current_nf = None

    for line in lines:
        if 'normfactor' in line:
            # Save previous mode
            if current_label is not None:
                modes.append({
                    'label': current_label,
                    'normfactor': current_nf,
                    'atoms': current_atoms,
                })
            label_part, nf_part = line.split('normfactor')
            current_label = label_part.strip()
            current_nf = float(nf_part.strip().split()[1])
            current_atoms = []
        elif current_label is not None:
            parts = line.split()
            if len(parts) == 7:
                try:
                    current_atoms.append((
                        parts[0],
                        float(parts[1]), float(parts[2]), float(parts[3]),
                        float(parts[4]), float(parts[5]), float(parts[6]),
                    ))
                except ValueError:
                    pass

    # Save last mode
    if current_label is not None:
        modes.append({
            'label': current_label,
            'normfactor': current_nf,
            'atoms': current_atoms,
        })

    if skip_strain:
        modes = [m for m in modes if 'strain' not in m['label']]

    return lattice_params, superstructure_positions, modes


# ---------------------------------------------------------------------------
# Atom matching
# ---------------------------------------------------------------------------

def match_atoms(iso_positions, poscar_positions, tol=0.02):
    """
    For each ISOTROPY superstructure atom, find the index of the matching
    POSCAR atom by minimum fractional-coordinate distance (periodic).

    Returns
    -------
    iso_to_poscar : dict  {iso_index: poscar_index}
    """
    poscar_pos = np.array(poscar_positions)
    iso_to_poscar = {}
    unmatched = []

    for iso_idx, (name, x, y, z) in enumerate(iso_positions):
        iso_pos = np.array([x, y, z]) % 1.0
        delta = (poscar_pos - iso_pos + 0.5) % 1.0 - 0.5
        dists = np.linalg.norm(delta, axis=1)
        best = int(np.argmin(dists))
        if dists[best] < tol:
            iso_to_poscar[iso_idx] = best
        else:
            unmatched.append((iso_idx, name, (x, y, z), dists[best]))

    if unmatched:
        print(f"WARNING: {len(unmatched)} ISOTROPY atoms could not be matched "
              f"to POSCAR atoms (tolerance={tol}):")
        for iso_idx, name, pos, d in unmatched:
            print(f"  iso[{iso_idx}] {name} {pos} -> min_dist={d:.4f}")

    return iso_to_poscar


# ---------------------------------------------------------------------------
# Basis construction
# ---------------------------------------------------------------------------

def build_basis(modes, iso_positions, poscar_positions, A, tol=0.02):
    """
    Build the SAM basis matrix.

    Each row corresponds to one SAM. Displacements are placed into the
    POSCAR DOF slots (3*poscar_idx + component), converted to Cartesian,
    then unit-normalized.

    Returns
    -------
    basis : (n_modes, 3*n_atoms) array
    labels : list of str
    warnings : list of str
    """
    n_atoms = len(poscar_positions)
    n_dof = 3 * n_atoms
    poscar_pos = np.array(poscar_positions)

    # Pre-build position lookup for speed
    basis = []
    labels = []
    warnings = []

    for mode in modes:
        label = mode['label']
        row = np.zeros(n_dof)
        unmatched_atoms = []

        for atom_name, ax, ay, az, dx, dy, dz in mode['atoms']:
            # Skip zero-displacement entries
            if abs(dx) + abs(dy) + abs(dz) < 1e-10:
                continue

            # Match to POSCAR atom
            iso_pos = np.array([ax, ay, az]) % 1.0
            delta = (poscar_pos - iso_pos + 0.5) % 1.0 - 0.5
            dists = np.linalg.norm(delta, axis=1)
            best = int(np.argmin(dists))

            if dists[best] > tol:
                unmatched_atoms.append(
                    f"{atom_name}({ax:.4f},{ay:.4f},{az:.4f}) dist={dists[best]:.4f}"
                )
                continue

            row[3*best + 0] = dx
            row[3*best + 1] = dy
            row[3*best + 2] = dz

        if unmatched_atoms:
            warnings.append(
                f"Mode '{label}': unmatched atoms: {', '.join(unmatched_atoms)}"
            )

        # Fractional -> Cartesian: multiply each atom's (dx,dy,dz) by A
        cart = row.reshape(n_atoms, 3) @ A

        # Unit normalize
        vec = cart.flatten()
        norm = np.linalg.norm(vec)
        if norm < 1e-10:
            warnings.append(f"Mode '{label}': zero-norm vector, skipping normalization.")
        else:
            vec = vec / norm

        basis.append(vec)
        labels.append(label)

    return np.array(basis), labels


# ---------------------------------------------------------------------------
# Validation
# ---------------------------------------------------------------------------

def validate_basis(basis, labels):
    """
    Check cross-irrep block orthogonality and within-block norms.
    Prints a summary.
    """
    import re

    def gm_label(label):
        m = re.search(r'(GM\d+[+-])', label)
        return m.group(1) if m else 'unknown'

    gm = [gm_label(l) for l in labels]
    unique_gm = list(dict.fromkeys(gm))

    print("\nValidation:")
    print(f"  Basis shape: {basis.shape}")
    print(f"  Unique irreps: {unique_gm}")

    # Row norms
    norms = np.linalg.norm(basis, axis=1)
    if np.allclose(norms, 1.0, atol=1e-6):
        print("  Row norms: all 1.0 â")
    else:
        bad = np.where(np.abs(norms - 1.0) > 1e-6)[0]
        print(f"  Row norms: {len(bad)} rows not unit-normalized: {bad}")

    # Cross-block orthogonality
    max_cross = 0.0
    bad_pairs = []
    for i, u in enumerate(unique_gm):
        for j, v in enumerate(unique_gm):
            if j <= i:
                continue
            idx_u = [k for k, x in enumerate(gm) if x == u]
            idx_v = [k for k, x in enumerate(gm) if x == v]
            block = basis[idx_u] @ basis[idx_v].T
            mx = float(np.abs(block).max())
            if mx > max_cross:
                max_cross = mx
            if mx > 1e-6:
                bad_pairs.append(f"{u}x{v}: {mx:.2e}")

    if bad_pairs:
        print(f"  Cross-block orthogonality FAILED:")
        for bp in bad_pairs:
            print(f"    {bp}")
    else:
        print(f"  Cross-block orthogonality: max={max_cross:.2e} â")


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

def main():
    parser = argparse.ArgumentParser(description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument('isotropy', help='ISODISTORT complete-mode-details output file')
    parser.add_argument('poscar',   help='VASP POSCAR file')
    parser.add_argument('--basis',  default='symmetry_basis',
                        help='Output path for symmetry_basis (default: symmetry_basis)')
    parser.add_argument('--list',   default='symmetry_list',
                        help='Output path for symmetry_list (default: symmetry_list)')
    parser.add_argument('--tol',    type=float, default=0.02,
                        help='Fractional-coord tolerance for atom matching (default: 0.02)')
    parser.add_argument('--skip-strain', action='store_true', default=True,
                        help='Exclude strain modes (default: True)')
    args = parser.parse_args()

    # -- Parse inputs ---------------------------------------------------------
    print(f"Reading POSCAR: {args.poscar}")
    A, poscar_positions = parse_poscar(args.poscar)
    n_atoms = len(poscar_positions)
    print(f"  {n_atoms} atoms, lattice:\n{A}")

    print(f"\nReading ISOTROPY file: {args.isotropy}")
    lattice_params, iso_positions, modes = parse_isotropy(
        args.isotropy, skip_strain=args.skip_strain)
    print(f"  {len(iso_positions)} superstructure atoms")
    print(f"  {len(modes)} displacement modes")
    if len(iso_positions) != n_atoms:
        print(f"  WARNING: atom count mismatch â ISOTROPY has {len(iso_positions)}, "
              f"POSCAR has {n_atoms}")

    # -- Build basis ----------------------------------------------------------
    print(f"\nBuilding basis (atom match tolerance={args.tol})...")
    basis, labels = build_basis(
        modes, iso_positions, poscar_positions, A, tol=args.tol)

    # -- Validate -------------------------------------------------------------
    validate_basis(basis, labels)

    # -- Write outputs --------------------------------------------------------
    print(f"\nWriting {args.basis}...")
    with open(args.basis, 'w') as f:
        for row in basis:
            f.write('  '.join(f'{v:.10f}' for v in row) + '\n')

    print(f"Writing {args.list}...")
    with open(args.list, 'w') as f:  # Note: args.list shadows builtin; fine here
        for i, label in enumerate(labels, 1):
            f.write(f'{i:4d}  {label}\n')

    print(f"\nDone. {len(labels)} SAMs written.")


if __name__ == '__main__':
    main()

