"""
parsers.py

All reading from vasprun.xml lives here. Each public function is independent:
it opens the file, finds its block, parses it, validates the result, and
returns a plain Python/numpy object. No function depends on another.

The shared helper _find_xml_block() handles the streaming parse so that
error messages and file-handling are consistent across all callers.

Blocks targeted:
    atomtypes   -> species masses and ZVAL (always present)
    born_charges -> per-atom BEC tensors (present only after LEPSILON/LCALCEPS)
    hessian     -> force-constant / dynamical matrix (present only after IBRION=8)
"""

import xml.etree.ElementTree as ET
import numpy as np
import sys


# ---------------------------------------------------------------------------
# Private helper
# ---------------------------------------------------------------------------

def _find_xml_block(vasprun_path, tag, name_attr, required=True):
    """
    Stream vasprun.xml with iterparse and return the first Element whose
    tag matches `tag` and whose 'name' attribute matches `name_attr`.

    Parameters
    ----------
    vasprun_path : str
        Path to vasprun.xml.
    tag : str
        XML element tag to match, e.g. 'array' or 'varray'.
    name_attr : str
        Value of the element's 'name' attribute, e.g. 'atomtypes'.
    required : bool
        If True, exit with an error message when the block is absent.
        If False, return None when the block is absent.

    Returns
    -------
    xml.etree.ElementTree.Element or None
    """
    try:
        context = ET.iterparse(vasprun_path, events=("end",))
    except FileNotFoundError:
        print(f"Error: vasprun.xml not found at {vasprun_path}")
        sys.exit(1)
    except ET.ParseError as e:
        print(f"Error: vasprun.xml could not be parsed: {e}")
        sys.exit(1)

    for _event, elem in context:
        if elem.tag == tag and elem.attrib.get("name") == name_attr:
            return elem

    if required:
        print(f"Error: <{tag} name=\"{name_attr}\"> block not found in {vasprun_path}")
        print(f"       Check that the VASP run completed and produced this quantity.")
        sys.exit(1)

    return None


# ---------------------------------------------------------------------------
# Public parsers
# ---------------------------------------------------------------------------

def parse_species_masses(vasprun_path):
    """
    Read atomic masses from the atomtypes block of vasprun.xml.

    Returns one mass per species (not per atom), in POSCAR species order.
    Expand to per-atom values downstream using atom counts from the POSCAR.

    Parameters
    ----------
    vasprun_path : str

    Returns
    -------
    list of float
        Atomic masses in AMU, one entry per species.

    Example
    -------
    For SrTiO3 (species Sr, Ti, O) returns [87.62, 47.88, 16.0]
    """
    elem = _find_xml_block(vasprun_path, "array", "atomtypes", required=True)

    masses = []
    for rc in elem.iter("rc"):
        cells = rc.findall("c")
        # Row layout: atomspertype | element | mass | valence | pseudopotential
        if len(cells) < 3:
            print("Error: atomtypes row has fewer columns than expected.")
            sys.exit(1)
        masses.append(float(cells[2].text))

    if not masses:
        print("Error: no species rows found in atomtypes block.")
        sys.exit(1)

    return masses


def parse_species_zvals(vasprun_path):
    """
    Read pseudopotential valence charges (ZVAL) from the atomtypes block.

    Returns one ZVAL per species (not per atom), in POSCAR species order.
    Used as a fallback when Born effective charges are absent.

    Parameters
    ----------
    vasprun_path : str

    Returns
    -------
    list of float
        ZVAL per species.

    Example
    -------
    For SrTiO3 (Sr_sv, Ti_sv, O) returns [10.0, 12.0, 6.0]
    """
    elem = _find_xml_block(vasprun_path, "array", "atomtypes", required=True)

    zvals = []
    for rc in elem.iter("rc"):
        cells = rc.findall("c")
        # Row layout: atomspertype | element | mass | valence | pseudopotential
        if len(cells) < 4:
            print("Error: atomtypes row has fewer columns than expected.")
            sys.exit(1)
        zvals.append(float(cells[3].text))

    if not zvals:
        print("Error: no species rows found in atomtypes block.")
        sys.exit(1)

    return zvals


def get_born_effective_charges(vasprun_path):
    """
    Read Born effective charge tensors from vasprun.xml.

    Returns None (without error) when the born_charges block is absent,
    which is the normal case for runs without LEPSILON=.TRUE. or LCALCEPS=.TRUE.
    The caller is responsible for falling back to ZVAL-based diagonal tensors.

    Parameters
    ----------
    vasprun_path : str

    Returns
    -------
    list of list of float, shape [3*N_atoms][3], or None
        BEC tensors flattened row-wise: for each atom, three rows of its
        3x3 tensor are appended in order. Matches the format expected by
        get_mode_effective_charges().

    Notes
    -----
    Each atom's tensor is stored as a <set> containing three <v> elements
    (one per Cartesian row of the 3x3 tensor).
    """
    elem = _find_xml_block(vasprun_path, "array", "born_charges", required=False)

    if elem is None:
        return None

    bec = []
    for atom_set in elem.iter("set"):
        for v in atom_set.findall("v"):
            row = [float(x) for x in v.text.split()]
            if len(row) != 3:
                print("Error: Born effective charge row does not have 3 components.")
                sys.exit(1)
            bec.append(row)

    if not bec:
        print("Error: born_charges block found but contains no data.")
        sys.exit(1)

    if len(bec) % 3 != 0:
        print(f"Error: born_charges has {len(bec)} rows, expected a multiple of 3.")
        sys.exit(1)

    return bec


def get_dynamical_matrix(vasprun_path):
    """
    Read the Hessian matrix from vasprun.xml and return it as a numpy array.

    The Hessian is stored in VASP's vasprun.xml in mass-weighted Cartesian
    coordinates (units THz^2, i.e. eV / (Angstrom^2 * AMU) after VASP's
    internal normalisation). Components are ordered as
    [x1, y1, z1, x2, y2, z2, ...] matching POSCAR atom order.

    This function does NOT diagonalise the matrix. Diagonalisation,
    symmetrisation, and frequency extraction are handled by
    solve_dynamical_eigenmodes() in support.py.

    Parameters
    ----------
    vasprun_path : str

    Returns
    -------
    numpy.ndarray, shape [3*N_atoms, 3*N_atoms]
        Raw Hessian matrix.
    """
    elem = _find_xml_block(vasprun_path, "varray", "hessian", required=True)

    rows = []
    for v in elem.findall("v"):
        row = [float(x) for x in v.text.split()]
        rows.append(row)

    if not rows:
        print("Error: hessian block found but contains no data.")
        sys.exit(1)

    matrix = np.array(rows)

    if matrix.shape[0] != matrix.shape[1]:
        print(f"Error: Hessian is not square: shape is {matrix.shape}.")
        sys.exit(1)

    return matrix
