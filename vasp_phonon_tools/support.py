"""
support.py

Core physics and POSCAR manipulation functions for vasp-phonon-tools.

Naming conventions
------------------
- All functions are lowercase with underscores.
- Functions that read files are prefixed read_ or validate_.
- Functions that compute physical quantities are prefixed get_ or solve_.
- Functions that modify a POSCAR are prefixed poscar_ or displace_.
- Functions that write files are prefixed write_.

Coordinate conventions
----------------------
All displacement vectors (phonon eigenvectors, symmetry-adapted mode vectors)
are in CARTESIAN coordinates with components ordered as
[x1, y1, z1, x2, y2, z2, ...] matching POSCAR atom order.

ISODISTORT can export symmetry-adapted modes in either Cartesian or
Direct (fractional) coordinates. Before passing SAM vectors to any
function here, the caller must convert them to Cartesian using
direct_to_cartesian_displacement() if needed.

POSCAR files may be in Direct or Cartesian format. Functions that
operate on atomic positions handle both via poscar_to_cartesian() and
poscar_to_direct().
"""

import numpy as np
import os
import copy
import sys

np.set_printoptions(precision=20, threshold=np.inf, suppress=True, linewidth=10000)


# ---------------------------------------------------------------------------
# File reading
# ---------------------------------------------------------------------------

def read_file_lines(path):
    """
    Read a plain-text file and return its lines as a list of strings.

    Parameters
    ----------
    path : str

    Returns
    -------
    list of str
    """
    try:
        with open(path, "r") as f:
            return f.readlines()
    except FileNotFoundError:
        print(f"Error: file not found: {path}")
        sys.exit(1)


# ---------------------------------------------------------------------------
# POSCAR validation and basic extraction
# ---------------------------------------------------------------------------

def validate_poscar(lines):
    """
    Check a POSCAR (supplied as a list of lines from read_file_lines) for
    common formatting errors. Strips trailing velocity lines if present.

    Parameters
    ----------
    lines : list of str

    Returns
    -------
    list of str
        Validated and trimmed POSCAR lines.
    """
    if len(lines) < 10:
        print("Error: POSCAR has fewer than 10 lines; too short to be valid.")
        sys.exit(1)

    try:
        scale = float(lines[1])
    except ValueError:
        print("Error: POSCAR line 2 (scale factor) is not a float.")
        sys.exit(1)

    if scale != 1.0:
        print("Warning: POSCAR scale factor is not 1. All code assumes scale=1.")
        sys.exit(1)

    for i in (2, 3, 4):
        if len(lines[i].split()) != 3:
            print(f"Error: POSCAR line {i+1} (lattice vector) does not have 3 columns.")
            sys.exit(1)
        try:
            [float(x) for x in lines[i].split()]
        except ValueError:
            print(f"Error: POSCAR line {i+1} lattice vector contains non-float values.")
            sys.exit(1)

    if len(lines[5].split()) != len(lines[6].split()):
        print("Error: POSCAR species names (line 6) and counts (line 7) have different lengths.")
        sys.exit(1)

    try:
        atom_counts = [int(x) for x in lines[6].split()]
    except ValueError:
        print("Error: POSCAR line 7 atom counts are not integers.")
        sys.exit(1)

    coord_type = lines[7].split()[0].lower()
    if coord_type not in ("direct", "cartesian"):
        print(f"Error: POSCAR line 8 coordinate type '{lines[7].split()[0]}' "
              "not recognised. Expected 'Direct' or 'Cartesian'.")
        sys.exit(1)

    n_atoms = sum(atom_counts)
    for i in range(n_atoms):
        idx = 8 + i
        try:
            cols = lines[idx].split()
        except IndexError:
            print(f"Error: POSCAR ended before all {n_atoms} atomic positions were read.")
            sys.exit(1)
        if len(cols) != 3:
            print(f"Error: POSCAR atom position line {idx+1} does not have 3 columns.")
            sys.exit(1)
        try:
            [float(x) for x in cols]
        except ValueError:
            print(f"Error: POSCAR atom position line {idx+1} contains non-float values.")
            sys.exit(1)

    # Trim velocity lines if present
    end = 8 + n_atoms
    if len(lines) > end and len(lines[end].split()) != 0:
        print("Error: unexpected content after atomic positions in POSCAR.")
        sys.exit(1)

    return lines[:end]


def get_atom_counts(poscar_lines):
    """
    Return the number of atoms per species from a validated POSCAR.

    Parameters
    ----------
    poscar_lines : list of str

    Returns
    -------
    list of int
    """
    poscar_lines = validate_poscar(poscar_lines)
    return [int(x) for x in poscar_lines[6].split()]


def get_cell_volume(poscar_lines):
    """
    Compute the unit-cell volume from the POSCAR lattice vectors.

    Parameters
    ----------
    poscar_lines : list of str

    Returns
    -------
    float
        Volume in Angstrom^3. Exits if volume is non-positive.
    """
    poscar_lines = validate_poscar(poscar_lines)
    lattice = np.array([[float(x) for x in poscar_lines[i].split()] for i in (2, 3, 4)])
    volume = np.linalg.det(lattice)
    if volume <= 0:
        print("Error: unit-cell volume is non-positive. Use a right-handed coordinate system.")
        sys.exit(1)
    return volume


def get_lattice_vectors(poscar_lines):
    """
    Return the 3x3 lattice vector matrix from a validated POSCAR.

    Parameters
    ----------
    poscar_lines : list of str

    Returns
    -------
    numpy.ndarray, shape (3, 3)
        Rows are lattice vectors a, b, c.
    """
    poscar_lines = validate_poscar(poscar_lines)
    return np.array([[float(x) for x in poscar_lines[i].split()] for i in (2, 3, 4)])


# ---------------------------------------------------------------------------
# Masses and BEC (species-level operations)
# ---------------------------------------------------------------------------

def expand_species_to_atoms(species_values, atom_counts):
    """
    Expand a per-species list to a per-atom list using POSCAR atom counts.

    Parameters
    ----------
    species_values : list
        One value per species, in POSCAR species order.
    atom_counts : list of int
        Number of atoms per species, from get_atom_counts().

    Returns
    -------
    list
        One value per atom, repeated according to atom_counts.

    Example
    -------
    expand_species_to_atoms([87.62, 47.88, 16.0], [2, 2, 6])
    -> [87.62, 87.62, 47.88, 47.88, 16.0, 16.0, 16.0, 16.0, 16.0, 16.0]
    """
    if len(species_values) != len(atom_counts):
        print("Error: length of species_values does not match number of species "
              "in atom_counts.")
        sys.exit(1)
    result = []
    for value, count in zip(species_values, atom_counts):
        result.extend([value] * count)
    return result


def build_bec_from_zvals(species_zvals, atom_counts):
    """
    Construct diagonal Born effective charge tensors from ZVAL.

    This is the rigid-ion approximation: electronic screening is ignored.
    Phonon frequencies and eigenvectors are unaffected. Mode effective
    charges computed from these tensors will be approximate.

    Parameters
    ----------
    species_zvals : list of float
        ZVAL per species from parse_species_zvals().
    atom_counts : list of int
        From get_atom_counts().

    Returns
    -------
    list of list of float, shape [3*N_atoms][3]
        Diagonal BEC tensors. Each atom contributes three rows:
        [Z, 0, 0], [0, Z, 0], [0, 0, Z]
    """
    bec = []
    for z, count in zip(species_zvals, atom_counts):
        for _ in range(count):
            bec.extend([[z, 0.0, 0.0],
                        [0.0, z, 0.0],
                        [0.0, 0.0, z]])
    return bec


# ---------------------------------------------------------------------------
# Phonon eigenmodes
# ---------------------------------------------------------------------------

def solve_eigenmodes(matrix):
    """
    Solve the symmetric eigenvalue problem for a real square matrix.

    Parameters
    ----------
    matrix : numpy.ndarray, shape (n, n)

    Returns
    -------
    eigenvalues : numpy.ndarray, shape (n,)
    eigenvectors : numpy.ndarray, shape (n, n)
        Column k corresponds to eigenvalue k.
    """
    if matrix.shape[0] != matrix.shape[1]:
        print(f"Error: matrix is not square (shape {matrix.shape}).")
        sys.exit(1)
    return np.linalg.eigh(matrix)


def solve_dynamical_eigenmodes(hessian):
    """
    Symmetrise the Hessian, diagonalise it, and return phonon frequencies
    and mass-weighted eigenvectors.

    The Hessian from vasprun.xml is force-symmetrised as
    H_sym = -(H + H^T) / 2  (the minus sign follows VASP convention for
    the dynamical matrix sign). Eigenvalues are converted to frequencies
    in THz. Imaginary frequencies (negative eigenvalues) are returned as
    negative THz values.

    Parameters
    ----------
    hessian : numpy.ndarray, shape (3*N, 3*N)
        From parsers.get_dynamical_matrix().

    Returns
    -------
    frequencies : numpy.ndarray, shape (3*N,)
        Signed frequencies in THz. Ordered highest to lowest.
    eigenvectors : numpy.ndarray, shape (3*N, 3*N)
        Column k is the mass-weighted eigenvector for frequencies[k].
        Cartesian components ordered [x1,y1,z1, x2,y2,z2, ...].
    """
    sym = (hessian + hessian.T) * -0.5
    eigenvalues, eigenvectors = solve_eigenmodes(sym)

    # Reverse to highest-first order
    eigenvalues = eigenvalues[::-1]
    eigenvectors = np.fliplr(eigenvectors)

    # Track negative eigenvalues before sqrt
    negative = eigenvalues < 0
    frequencies = np.sqrt(np.abs(eigenvalues))
    frequencies[negative] *= -1

    return frequencies, eigenvectors


def get_real_space_displacements(atom_counts, species_masses, eigenvectors):
    """
    Convert mass-weighted eigenvectors to real-space displacements in Angstrom,
    and compute the reduced mass of each mode.

    The mass-weighted eigenvectors e_{k,alpha} have units of 1/sqrt(AMU).
    Real-space displacements are u_{k,alpha} = e_{k,alpha} / sqrt(m_k),
    giving units of Angstrom. The normalisation condition becomes
    sum_k  m_k |u_k|^2 = 1, which defines the reduced mass as its inverse.

    Parameters
    ----------
    atom_counts : list of int
        From get_atom_counts().
    species_masses : list of float
        Per-species masses from parsers.parse_species_masses().
    eigenvectors : numpy.ndarray, shape (3*N, 3*N)
        From solve_dynamical_eigenmodes(). Each column is one mode.

    Returns
    -------
    reduced_masses : list of float
        Reduced mass per mode in AMU.
    displacements : numpy.ndarray, shape (3*N, 3*N)
        Row k is the real-space displacement vector for mode k,
        in Cartesian coordinates (Angstrom), normalised to unit length.
    """
    if len(atom_counts) == 0 or len(species_masses) == 0 or len(eigenvectors) == 0:
        print("Error: empty input to get_real_space_displacements.")
        sys.exit(1)
    if sum(atom_counts) * 3 != len(eigenvectors):
        print("Error: atom count does not match eigenvector dimension.")
        sys.exit(1)

    # Build diagonal inverse-mass matrix (per degree of freedom)
    inv_sqrt_mass = []
    for mass, count in zip(species_masses, atom_counts):
        inv_sqrt_mass.extend([1.0 / np.sqrt(mass)] * (3 * count))
    inv_sqrt_mass = np.diag(inv_sqrt_mass)

    # u = M^{-1/2} e  (columns = modes)
    disp_cols = np.dot(inv_sqrt_mass, eigenvectors)

    # Reduced mass and normalisation: work column-by-column
    reduced_masses = []
    for k in range(disp_cols.shape[1]):
        col = disp_cols[:, k]
        norm_sq = np.dot(col, col)
        reduced_masses.append(1.0 / norm_sq)
        disp_cols[:, k] = col / np.linalg.norm(col)

    # Return as rows (one displacement vector per row)
    return reduced_masses, disp_cols.T


# ---------------------------------------------------------------------------
# Mode effective charges
# ---------------------------------------------------------------------------

def get_mode_effective_charges(atom_counts, bec, displacements):
    """
    Compute mode effective charges for each phonon mode.

    Based on equation (53) of Gonze & Lee, PRB 55, 10355 (1997).

    Parameters
    ----------
    atom_counts : list of int
    bec : list of list of float, shape [3*N_atoms][3]
        Born effective charges from parsers.get_born_effective_charges()
        or build_bec_from_zvals().
    displacements : numpy.ndarray, shape (3*N, 3*N)
        Real-space displacement vectors, one per row, from
        get_real_space_displacements(). Units: Angstrom.

    Returns
    -------
    list of list of float
        Shape [n_modes][3]. Mode effective charge vector (x, y, z) in
        units of electron charge e.
    """
    n_atoms = sum(atom_counts)
    if len(bec) != 3 * n_atoms or len(displacements) != 3 * n_atoms:
        print("Error: atom count does not match BEC or displacement dimensions.")
        sys.exit(1)

    mec_list = []
    for n in range(len(displacements)):
        u = displacements[n]
        norm = np.dot(u, u)
        mec = []
        for alpha in range(3):
            z_star = sum(
                u[3 * k + beta] * bec[3 * k + alpha][beta]
                for k in range(n_atoms)
                for beta in range(3)
            )
            if abs(z_star / norm**0.5) < 1e-10:
                z_star = 0.0
            mec.append(z_star / norm**0.5)
        mec_list.append(mec)
    return mec_list


# ---------------------------------------------------------------------------
# Symmetry-adapted modes
# ---------------------------------------------------------------------------

def read_symmetry_adapted_modes(path):
    """
    Read symmetry-adapted mode (SAM) basis vectors and their irrep labels.

    Looks first for a combined 'symmetry_file', then falls back to the
    separate 'symmetry_basis' and 'symmetry_list' files generated by
    ISODISTORT (https://iso.byu.edu).

    COORDINATE WARNING: ISODISTORT can export SAM vectors in either
    Cartesian or Direct (fractional) coordinates. This function reads
    whatever is in the file. If the vectors are in Direct coordinates,
    convert them with direct_to_cartesian_displacement() before use.

    Parameters
    ----------
    path : str
        Directory containing the symmetry files.

    Returns
    -------
    basis : list of list of float
        Each inner list is one normalised SAM basis vector of length 3*N.
    labels : list of list of str
        Each inner list is the label tokens for the corresponding SAM.
    """
    combined = os.path.join(path, "symmetry_file")
    basis_path = os.path.join(path, "symmetry_basis")
    list_path = os.path.join(path, "symmetry_list")

    basis = []
    labels = []

    if os.path.isfile(combined):
        lines = read_file_lines(combined)
        # First block: labels (until blank line), second block: vectors
        lines = lines[1:]  # skip header
        i = 0
        while i < len(lines) and lines[i].split():
            labels.append(lines[i].split())
            i += 1
        i += 1  # skip blank line
        while i < len(lines):
            if lines[i].split():
                basis.append(lines[i].split())
            i += 1
    else:
        if not os.path.isfile(basis_path):
            print(f"Error: symmetry_basis not found in {path}")
            sys.exit(1)
        if not os.path.isfile(list_path):
            print(f"Error: symmetry_list not found in {path}")
            sys.exit(1)
        for line in read_file_lines(basis_path):
            if line.split():
                basis.append(line.split())
        for line in read_file_lines(list_path):
            if line.split():
                labels.append(line.split())

    if not basis or not labels:
        print("Error: symmetry basis or label file is empty.")
        sys.exit(1)

    # Convert to float and normalise
    basis = [[float(x) for x in row] for row in basis]
    labels = [[str(x) for x in row] for row in labels]

    return basis, labels


def direct_to_cartesian_displacement(displacement_vector, lattice_vectors):
    """
    Convert a displacement vector from Direct (fractional) to Cartesian
    coordinates.

    SAM vectors from ISODISTORT may be in Direct coordinates, while all
    physics functions here expect Cartesian. Use this function to convert
    before decomposing phonons or displacing a POSCAR.

    The conversion is:
        v_cart[k] = L @ v_direct[k]   for each atom k
    where L is the 3x3 matrix of lattice vectors (rows = a, b, c).

    Parameters
    ----------
    displacement_vector : array-like, shape (3*N,)
        Displacement vector in Direct coordinates, ordered
        [x1,y1,z1, x2,y2,z2, ...].
    lattice_vectors : numpy.ndarray, shape (3, 3)
        From get_lattice_vectors(). Rows are lattice vectors a, b, c.

    Returns
    -------
    numpy.ndarray, shape (3*N,)
        Displacement vector in Cartesian coordinates (Angstrom).
    """
    v = np.array(displacement_vector)
    n_atoms = len(v) // 3
    if len(v) % 3 != 0:
        print("Error: displacement vector length is not a multiple of 3.")
        sys.exit(1)
    v_matrix = v.reshape(n_atoms, 3)
    cart_matrix = v_matrix @ lattice_vectors
    return cart_matrix.reshape(3 * n_atoms)


def assign_phonon_symmetry(phonon_decomp, labels, frequencies,
                           reduced_masses, mec=None, tolerance=1e-2):
    """
    Match each phonon to its dominant symmetry-adapted modes and assemble
    the output table.

    Parameters
    ----------
    phonon_decomp : numpy.ndarray, shape (n_modes, n_sam)
        SAM amplitudes per phonon, from decompose_phonons_into_sam().
    labels : list of list of str
        SAM labels from read_symmetry_adapted_modes().
    frequencies : numpy.ndarray
        Phonon frequencies in THz.
    reduced_masses : list of float
        Per-mode reduced masses in AMU.
    mec : list of list of float or None
        Mode effective charges. If None, zeros are used.
    tolerance : float
        SAM amplitudes below this threshold are suppressed (default 1e-2).

    Returns
    -------
    list
        One entry per phonon. Each entry is:
        [index, frequency, force_constant, reduced_mass,
         mec_x, mec_y, mec_z, mec_magnitude,
         sam_indices (list), sam_labels (list)]
    """
    if mec is None:
        mec = [[0.0, 0.0, 0.0]] * len(frequencies)

    output = []
    for i in range(len(phonon_decomp)):
        freq = frequencies[i]
        mu = reduced_masses[i]
        # Force constant: k = mu * (2*pi*nu)^2 converted to eV/A^2
        force_constant = mu * freq * freq * 0.00409164963432946585666737210006

        mec_x, mec_y, mec_z = mec[i]
        mec_mag = (mec_x**2 + mec_y**2 + mec_z**2)**0.5

        # SAM indices (1-based) with amplitude above tolerance
        sam_indices = [j + 1 for j, amp in enumerate(phonon_decomp[i])
                       if abs(amp) > tolerance]

        # Unique irrep labels for those SAMs (preserve order, no duplicates)
        seen = []
        sam_labels = []
        for j in sam_indices:
            lbl = labels[j - 1]
            if lbl not in seen:
                seen.append(lbl)
                sam_labels.append(lbl)

        output.append([
            str(i + 1), freq, force_constant, mu,
            mec_x, mec_y, mec_z, mec_mag,
            sam_indices, sam_labels
        ])

    return output


def decompose_phonons_into_sam(sam_basis, eigenvectors):
    """
    Decompose phonon eigenvectors into the symmetry-adapted mode basis
    using least-squares projection.

    Parameters
    ----------
    sam_basis : list of list of float
        From read_symmetry_adapted_modes(). Shape (n_sam, 3*N).
    eigenvectors : numpy.ndarray, shape (3*N, 3*N)
        Mass-weighted eigenvectors from solve_dynamical_eigenmodes().
        Each column is one mode.

    Returns
    -------
    numpy.ndarray, shape (n_modes, n_sam)
        Row k gives the SAM amplitudes for phonon mode k.
    """
    sam = np.array(sam_basis)
    # Solve: sam.T @ coeffs = eigenvectors  (one solution per mode column)
    coeffs, _, _, _ = np.linalg.lstsq(sam.T, eigenvectors, rcond=None)
    # Transpose so rows = modes
    return coeffs.T


# ---------------------------------------------------------------------------
# POSCAR coordinate conversion
# ---------------------------------------------------------------------------

def poscar_to_cartesian(poscar_lines):
    """
    Convert a POSCAR from Direct to Cartesian coordinates.

    Parameters
    ----------
    poscar_lines : list of str

    Returns
    -------
    list of str
        POSCAR lines with Cartesian atomic positions.
    """
    poscar_lines = validate_poscar(poscar_lines)
    if poscar_lines[7].split()[0].lower() == "cartesian":
        print("Error: POSCAR is already in Cartesian coordinates.")
        sys.exit(1)

    lattice = get_lattice_vectors(poscar_lines)
    n_atoms = sum(get_atom_counts(poscar_lines))
    direct = np.array([[float(x) for x in poscar_lines[8 + i].split()]
                       for i in range(n_atoms)])
    cartesian = direct @ lattice

    result = list(poscar_lines)
    result[7] = "Cartesian\n"
    for i in range(n_atoms):
        result[8 + i] = " ".join(str(x) for x in cartesian[i])
    return result


def poscar_to_direct(poscar_lines):
    """
    Convert a POSCAR from Cartesian to Direct coordinates.

    Parameters
    ----------
    poscar_lines : list of str

    Returns
    -------
    list of str
        POSCAR lines with Direct atomic positions.
    """
    poscar_lines = validate_poscar(poscar_lines)
    if poscar_lines[7].split()[0].lower() == "direct":
        print("Error: POSCAR is already in Direct coordinates.")
        sys.exit(1)

    lattice = get_lattice_vectors(poscar_lines)
    n_atoms = sum(get_atom_counts(poscar_lines))
    cart = np.array([[float(x) for x in poscar_lines[8 + i].split()]
                     for i in range(n_atoms)])
    direct = cart @ np.linalg.inv(lattice)

    result = list(poscar_lines)
    result[7] = "Direct\n"
    for i in range(n_atoms):
        result[8 + i] = " ".join(str(x) for x in direct[i])
    return result


# ---------------------------------------------------------------------------
# POSCAR displacement
# ---------------------------------------------------------------------------

def displace_poscar_atoms(poscar_lines, displacement_vector, amplitude, label):
    """
    Apply a displacement to a POSCAR and return the distorted structure.

    The displacement vector must be in Cartesian coordinates. If your
    vector is in Direct coordinates (e.g. from ISODISTORT), convert it
    first with direct_to_cartesian_displacement().

    The vector is normalised to unit length before scaling by amplitude,
    so the total displacement magnitude equals abs(amplitude).

    Parameters
    ----------
    poscar_lines : list of str
        Input POSCAR.
    displacement_vector : array-like, shape (3*N,)
        Displacement pattern in Cartesian coordinates.
    amplitude : float
        Displacement magnitude in Angstrom.
    label : str
        Appended to the POSCAR title line to identify the distortion.

    Returns
    -------
    list of str
        Distorted POSCAR in Direct coordinates, with 16-decimal precision.
    """
    poscar_lines = validate_poscar(poscar_lines)
    n_atoms = sum(get_atom_counts(poscar_lines))

    if len(displacement_vector) != 3 * n_atoms:
        print(f"Error: displacement vector has length {len(displacement_vector)}, "
              f"expected {3 * n_atoms} (3 * {n_atoms} atoms).")
        sys.exit(1)

    # Work in Cartesian
    working = list(poscar_lines)
    was_direct = working[7].split()[0].lower() == "direct"
    if was_direct:
        working = poscar_to_cartesian(working)

    positions = np.array([[float(x) for x in working[8 + i].split()]
                          for i in range(n_atoms)])

    disp = np.array(displacement_vector, dtype=float)
    disp = disp / np.linalg.norm(disp) * amplitude

    disp_matrix = disp.reshape(n_atoms, 3)
    new_positions = positions + disp_matrix

    # Amplitude sanity check
    actual = np.linalg.norm(new_positions - positions)
    if abs(actual - abs(amplitude)) > 1e-8:
        print("Warning: displacement amplitude check failed. Proceed with caution.")

    for i in range(n_atoms):
        working[8 + i] = " ".join(str(x) for x in new_positions[i])

    result = poscar_to_direct(working)
    result[0] = (poscar_lines[0].rstrip()
                 + f" | displaced: {label}, amplitude: {amplitude} Ang\n")

    # Format positions to 16 decimal places
    for i in range(n_atoms):
        coords = [float(x) for x in result[8 + i].split()]
        result[8 + i] = (" {: .16f}  {: .16f}  {: .16f}\n".format(*coords))

    return result


def write_poscar(poscar_lines, filename):
    """
    Write a POSCAR to disk with consistent column-aligned formatting.

    Parameters
    ----------
    poscar_lines : list of str
    filename : str
        Output filename. Written to the current working directory.
    """
    if not filename.strip():
        print("Error: output filename is empty.")
        sys.exit(1)

    output_path = os.path.join(os.getcwd(), filename)
    with open(output_path, "w") as f:
        for i, line in enumerate(poscar_lines):
            if i >= 8:
                parts = line.split()
                formatted = []
                for p in parts:
                    val = float(p)
                    s = f"{val: .16f}"
                    formatted.append(s)
                f.write("  " + "   ".join(formatted) + "\n")
            else:
                f.write(line)


def decompose_poscar_displacement(parent_poscar, child_poscar,
                                  displacement_basis, cutoff=1e-3):
    """
    Find the amplitudes of each basis vector in the displacement between
    two POSCARs (parent and child/distorted structure).

    Useful for identifying which symmetry-adapted modes or phonons are
    present in a known structural distortion.

    Parameters
    ----------
    parent_poscar : list of str
    child_poscar : list of str
    displacement_basis : list of array-like
        Basis vectors in Cartesian coordinates, each of length 3*N.
    cutoff : float
        Amplitudes below this value are not printed.
    """
    parent = validate_poscar(list(parent_poscar))
    child = validate_poscar(list(child_poscar))
    basis = np.array(displacement_basis)

    # Normalise basis
    for i in range(len(basis)):
        basis[i] /= np.linalg.norm(basis[i])

    lattice_p = get_lattice_vectors(parent)
    lattice_c = get_lattice_vectors(child)
    n_atoms = sum(get_atom_counts(parent))

    if sum(get_atom_counts(child)) != n_atoms:
        print("Error: parent and child POSCARs have different atom counts.")
        sys.exit(1)

    # Convert both to Cartesian
    def to_cart(poscar, lattice):
        if poscar[7].split()[0].lower() == "cartesian":
            return np.array([[float(x) for x in poscar[8 + i].split()]
                             for i in range(n_atoms)])
        direct = np.array([[float(x) for x in poscar[8 + i].split()]
                           for i in range(n_atoms)])
        return direct @ lattice

    cart_p = to_cart(parent, lattice_p)
    cart_c = to_cart(child, lattice_c)
    diff = (cart_c - cart_p).reshape(n_atoms * 3)

    amplitudes, _, _, _ = np.linalg.lstsq(basis.T, diff, rcond=None)
    for i, amp in enumerate(amplitudes):
        if abs(amp) > cutoff:
            print(f"  {i + 1:4d}   {amp:+.12f}")
