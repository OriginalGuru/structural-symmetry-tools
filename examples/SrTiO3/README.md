# Example: SrTiO₃ (10-atom supercell)

Tetragonal SrTiO₃, 2×1×1 supercell (2 Sr, 2 Ti, 6 O, 10 atoms, 30 degrees of freedom).
Space group I4/mcm. Lattice parameters: a = b = 3.9 Å, c = 7.8 Å.

## Files

| File | Description |
|---|---|
| `POSCAR` | 10-atom structure in Direct coordinates |
| `symmetry_basis` | 30 SAM basis vectors from ISODISTORT (Γ and X point modes) |
| `symmetry_list` | Irreducible representation labels for each SAM |
| `expected_output.txt` | Reference output from `get_phonon_mode_charges.py` |

`vasprun.xml` is not included (too large for the repository). To reproduce the
output, run a VASP DFPT calculation with `IBRION=8` and place `vasprun.xml` in
this directory alongside the other files.

## Running the example

```bash
cd examples/SrTiO3
# place your vasprun.xml here
python ../../scripts/get_phonon_mode_charges.py > output.txt
diff output.txt expected_output.txt
```

## Notes on this calculation

- Born effective charges are **not present** in the reference `vasprun.xml`
  (the run used `IBRION=8` without `LEPSILON=.TRUE.`). Mode effective charges
  are therefore computed using the rigid-ion approximation (diagonal BEC from ZVAL).
  The warning in `expected_output.txt` is expected.

- Three near-zero frequencies (modes 25–27) and three negative frequencies
  (modes 28–30) are present. The negative modes indicate structural instabilities
  at this geometry, consistent with the known soft-mode behaviour of SrTiO₃.

- The SAM basis covers both Γ-point (GM) and X-point modes. Phonons with
  significant X-point SAM character reflect zone-boundary instabilities.

## Coordinate convention for symmetry_basis

The vectors in `symmetry_basis` were exported from ISODISTORT in **Cartesian**
coordinates. No `--direct` flag is needed when using `displace_poscar_sam.py`.
