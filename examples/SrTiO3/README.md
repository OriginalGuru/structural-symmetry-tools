# Example: SrTiO3 (10-atom supercell)

Tetragonal SrTiO3, 2x1x1 supercell (2 Sr, 2 Ti, 6 O, 10 atoms, 30 degrees of freedom).
Space group Pm-3m (No. 221). Lattice parameters: a = b = 3.9 Ang, c = 7.8 Ang.

## Files

| File | Description |
|---|---|
| `POSCAR` | 10-atom structure in Direct coordinates |
| `symmetry_basis` | 30 SAM basis vectors from ISODISTORT (Gamma and X point modes) |
| `symmetry_list` | Irreducible representation labels for each SAM |
| `expected_output.txt` | Reference output from `get_phonon_mode_charges.py` |
| `findsym_input.txt` | Reference FINDSYM input generated from POSCAR |

`vasprun.xml` is not included (too large for the repository). To reproduce the
output, run a VASP DFPT calculation with `IBRION=8` and place `vasprun.xml` in
this directory alongside the other files.

## Running the phonon analysis

```bash
cd examples/SrTiO3
# place your vasprun.xml here
python ../../scripts/get_phonon_mode_charges.py > output.txt
python ../../scripts/validate_output.py output.txt expected_output.txt
```

A passing validation prints:
```
PASSED -- all 30 modes within tolerance 1.00e-06
```

## Generating displaced structures

```bash
# Displace along the soft mode (mode 30, imaginary)
python ../../scripts/displace_poscar_phonon.py 30 0.1

# Displace along a single SAM
python ../../scripts/displace_poscar_sam.py 3 0.1

# Displace along a linear combination of SAMs
python ../../scripts/displace_poscar_sam.py "1.0:9,1.0:12" 0.1
```

## Space group identification

```bash
# Identify space group of the reference structure
python ../../scripts/run_findsym.py
# Expected: Space Group: 221  Oh-1  Pm-3m

# Check symmetry of a displaced structure (looser tolerance)
python ../../scripts/run_findsym.py POSCAR_phonon_30_-0.086THz_0.1Ang --postol 0.01 --output-stem distorted
```

## Notes on this calculation

- Born effective charges are **not present** in the reference `vasprun.xml`
  (the run used `IBRION=8` without `LEPSILON=.TRUE.`). Mode effective charges
  are therefore computed using the rigid-ion approximation (diagonal BEC from ZVAL).
  The warning in `expected_output.txt` is expected.

- Three near-zero frequencies (modes 25-27) and three negative frequencies
  (modes 28-30) are present. The negative modes indicate structural instabilities
  at this geometry, consistent with the known soft-mode behaviour of SrTiO3.

- The SAM basis covers both Gamma-point (GM) and X-point modes. Phonons with
  significant X-point SAM character reflect zone-boundary instabilities.

## Coordinate convention for symmetry_basis

The vectors in `symmetry_basis` were exported from ISODISTORT in **Cartesian**
coordinates. No `--direct` flag is needed when using `displace_poscar_sam.py`.
