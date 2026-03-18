# vasp-phonon-tools

Analysis of phonon modes and symmetry-adapted displacements from VASP DFPT calculations.

## What it does

- Computes phonon frequencies, reduced masses, and force constants from the VASP Hessian
- Calculates mode effective charges (Z\*) using Born effective charges from VASP or a rigid-ion fallback
- Decomposes phonon eigenvectors into symmetry-adapted modes (SAMs) from ISODISTORT
- Generates distorted POSCAR files along phonon modes or SAMs for use in further calculations

## Requirements

- Python >= 3.9
- NumPy >= 1.21

No other dependencies. The standard library `xml.etree.ElementTree` handles all XML parsing.

## Installation

```bash
# Clone the repository
git clone https://github.com/yourname/vasp-phonon-tools.git
cd vasp-phonon-tools

# Create and activate a virtual environment (recommended)
python3 -m venv venv
source venv/bin/activate

# Install dependencies
pip install -r requirements.txt
```

See `docs/lonestar6.md` for TACC Lonestar6 specific instructions.

## Required input files

| File | Source | Contents |
|---|---|---|
| `vasprun.xml` | VASP output | Hessian, Born effective charges (optional), masses, ZVAL |
| `POSCAR` | VASP input | Structure used in the DFPT run |
| `symmetry_basis` | ISODISTORT | Symmetry-adapted mode basis vectors |
| `symmetry_list` | ISODISTORT | Irreducible representation labels |

All four files must be in the working directory. `vasprun.xml` and `OUTCAR` are produced
by a DFPT run with `IBRION=8` and (optionally) `LEPSILON=.TRUE.` in the INCAR.

Born effective charges are optional. If absent from `vasprun.xml`, diagonal tensors are
constructed from the pseudopotential valence charge (ZVAL). Frequencies and eigenvectors
are unaffected; mode effective charges will be approximate.

## Usage

### Phonon mode charges and SAM decomposition

Run from the directory containing the input files:

```bash
python scripts/get_phonon_mode_charges.py > output.txt
```

### Displace a POSCAR along a phonon mode

```bash
python scripts/displace_poscar_phonon.py <mode_index> <amplitude>
```

- `mode_index`: 1-based, ordered highest frequency first
- `amplitude`: displacement in Ångström (may be negative)

```bash
# Example: displace along mode 5 by 0.1 Å
python scripts/displace_poscar_phonon.py 5 0.1
```

### Validate output against a reference
```bash
python scripts/validate_output.py output.txt examples/SrTiO3/expected_output.txt
```

Checks frequency, force constant, reduced mass, and |Z*| for each mode.
Sign-dependent quantities (Z*_x/y/z, SAM amplitudes) are intentionally
skipped — these are not invariant across platforms or numpy versions.
See `docs/coordinate_conventions.md` for explanation.

### Displace a POSCAR along a symmetry-adapted mode

Single SAM:
```bash
python scripts/displace_poscar_sam.py <sam_index> <amplitude>
```

Linear combination:
```bash
python scripts/displace_poscar_sam.py "1.0:3,1.0:5" <amplitude>
```

If your `symmetry_basis` vectors are in Direct (fractional) coordinates rather than
Cartesian, add `--direct`:
```bash
python scripts/displace_poscar_sam.py 3 0.2 --direct
```

See `docs/coordinate_conventions.md` for details on this distinction.

## Example

The `examples/SrTiO3/` directory contains input files for a 10-atom SrTiO₃ supercell
(tetragonal, 2×1×1). See `examples/SrTiO3/README.md` for expected output.

## Code structure

```
vasp_phonon_tools/
    parsers.py    # All vasprun.xml reading (independent functions, iterparse)
    support.py    # Physics and POSCAR operations
scripts/
    get_phonon_mode_charges.py    # Main analysis script
    displace_poscar_phonon.py     # Phonon -> POSCAR
    displace_poscar_sam.py        # SAM -> POSCAR
docs/
    coordinate_conventions.md
    lonestar6.md
    usage.md
```

## Coordinate conventions

All displacement vectors (phonon eigenvectors, SAM vectors passed to POSCAR
displacement functions) must be in **Cartesian coordinates**.

ISODISTORT can export SAM vectors in either Cartesian or Direct coordinates.
Use `--direct` with `displace_poscar_sam.py`, or call
`support.direct_to_cartesian_displacement()` directly, to convert before use.

See `docs/coordinate_conventions.md` for full details.
