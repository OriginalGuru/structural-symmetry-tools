# structural-symmetry-tools

Analysis of phonon modes, symmetry-adapted displacements, and space group
identification for crystal structures from VASP DFPT calculations.

## What it does

- Computes phonon frequencies, reduced masses, and force constants from the VASP Hessian
- Calculates mode effective charges (Z\*) using Born effective charges from VASP or a rigid-ion fallback
- Decomposes phonon eigenvectors into symmetry-adapted modes (SAMs) from ISODISTORT
- Generates distorted POSCAR files along phonon modes or SAMs for use in further calculations
- Identifies space groups using FINDSYM from the ISOTROPY Software Suite

## Requirements

- Python >= 3.9
- NumPy >= 1.21

No other Python dependencies. The standard library `xml.etree.ElementTree` handles all XML parsing.

For space group identification, `findsym` from the ISOTROPY Software Suite must be
installed separately and available on PATH. See `docs/lonestar6.md`.

## Installation

```bash
# Clone the repository
git clone https://github.com/OriginalGuru/structural-symmetry-tools.git
cd structural-symmetry-tools

# Create and activate a virtual environment (recommended)
python3 -m venv venv
source venv/bin/activate

# Install dependencies
pip install -r requirements.txt
```

See `docs/lonestar6.md` for TACC Lonestar6 specific instructions.

## Required input files

### Phonon analysis

| File | Source | Contents |
|---|---|---|
| `vasprun.xml` | VASP output | Hessian, Born effective charges (optional), masses, ZVAL |
| `POSCAR` | VASP input | Structure used in the DFPT run |
| `symmetry_basis` | ISODISTORT | Symmetry-adapted mode basis vectors |
| `symmetry_list` | ISODISTORT | Irreducible representation labels |

`vasprun.xml` requires a DFPT run with `IBRION=8` in the INCAR. Born effective
charges are optional (`LEPSILON=.TRUE.`). If absent, diagonal tensors are
constructed from ZVAL. Frequencies and eigenvectors are unaffected; mode
effective charges will be approximate.

### Space group identification

| File | Source |
|---|---|
| `POSCAR` | Any VASP structure file — reference or displaced |

## Usage

### Phonon mode charges and SAM decomposition

Run from the directory containing the input files:

```bash
python scripts/get_phonon_mode_charges.py > output.txt
```

### Validate output against a reference

```bash
python scripts/validate_output.py output.txt examples/SrTiO3/expected_output.txt
```

Checks frequency, force constant, reduced mass, and |Z\*| for each mode.
Sign-dependent quantities (Z\*_x/y/z, SAM amplitudes) are intentionally
skipped — these are not invariant across platforms or numpy versions.
See `docs/coordinate_conventions.md` for explanation.

### Displace a POSCAR along a phonon mode

```bash
python scripts/displace_poscar_phonon.py <mode_index> <amplitude>
```

- `mode_index`: 1-based, ordered highest frequency first
- `amplitude`: displacement in Angstrom (may be negative)

```bash
# Displace along mode 30 by 0.1 Angstrom
python scripts/displace_poscar_phonon.py 30 0.1
# Output: POSCAR_phonon_30_-0.086THz_0.1Ang
```

The frequency is embedded in the output filename. The POSCAR title line records
the full mode description and amplitude.

### Displace a POSCAR along a symmetry-adapted mode

Single SAM:
```bash
python scripts/displace_poscar_sam.py <sam_index> <amplitude>
```

Linear combination:
```bash
python scripts/displace_poscar_sam.py "1.0:3,1.0:5" <amplitude>
```

```bash
# Examples
python scripts/displace_poscar_sam.py 3 0.1
# Output: POSCAR_sam_3_0.1Ang

python scripts/displace_poscar_sam.py "1.0:9,1.0:12" 0.1
# Output: POSCAR_sam_9+12_0.1Ang
```

If your `symmetry_basis` vectors are in Direct (fractional) coordinates rather than
Cartesian, add `--direct`:
```bash
python scripts/displace_poscar_sam.py 3 0.2 --direct
```

See `docs/coordinate_conventions.md` for details.

### Identify space group with FINDSYM

Generate a FINDSYM input file and run in one step:
```bash
python scripts/run_findsym.py
```

Run on any specific POSCAR file, including displaced structures:
```bash
python scripts/run_findsym.py POSCAR_phonon_30_-0.086THz_0.1Ang
python scripts/run_findsym.py POSCAR_sam_3_0.1Ang --postol 0.01 --output-stem distorted_sam3
```

Generate the input file only (to inspect before running):
```bash
python scripts/poscar_to_findsym.py
python scripts/poscar_to_findsym.py POSCAR_sam_3_0.1Ang --postol 0.01
findsym findsym_input.txt
```

Output: `findsym.log` and `findsym.cif` written alongside the POSCAR.
Requires `findsym` on PATH and `ISODATA` set. See `docs/lonestar6.md`.

## Examples

| Directory | Structure | Atoms | Description |
|---|---|---|---|
| `examples/SrTiO3/` | SrTiO3 | 10 | Tetragonal supercell, phonon + SAM analysis |
| `examples/ReS2_164_slab/` | ReS2 slab | 12 | Hexagonal supercell, space group identification |

## Code structure

```
vasp_phonon_tools/
    parsers.py    # All vasprun.xml reading (independent functions, iterparse)
    support.py    # Physics and POSCAR operations
scripts/
    get_phonon_mode_charges.py    # Phonon frequencies, Z*, SAM decomposition
    validate_output.py            # Validate output against a reference
    displace_poscar_phonon.py     # Apply phonon displacement to POSCAR
    displace_poscar_sam.py        # Apply SAM displacement to POSCAR
    poscar_to_findsym.py          # Convert POSCAR to FINDSYM input
    run_findsym.py                # Run FINDSYM and write log + CIF
examples/
    SrTiO3/                       # SrTiO3 10-atom supercell
    ReS2_164_slab/                # ReS2 12-atom slab
docs/
    coordinate_conventions.md
    lonestar6.md
```

## Coordinate conventions

All displacement vectors (phonon eigenvectors, SAM vectors passed to POSCAR
displacement functions) must be in **Cartesian coordinates**.

ISODISTORT can export SAM vectors in either Cartesian or Direct coordinates.
Use `--direct` with `displace_poscar_sam.py`, or call
`support.direct_to_cartesian_displacement()` directly, to convert before use.

See `docs/coordinate_conventions.md` for full details.
