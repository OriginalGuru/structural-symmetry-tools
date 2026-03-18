# Coordinate conventions

This document describes the coordinate systems used throughout vasp-phonon-tools
and where conversions are needed.

## Phonon eigenvectors

The Hessian in `vasprun.xml` is written by VASP in **Cartesian coordinates**,
in units of THz² (eV / Å² / AMU after mass-weighting). Components are ordered as

    [x₁, y₁, z₁, x₂, y₂, z₂, ...]

matching the POSCAR atom order. Diagonalising the Hessian yields eigenvectors in
the same Cartesian, mass-weighted basis.

`get_real_space_displacements()` removes the mass-weighting (multiplies by M⁻¹/²)
to give real-space displacements in **Angstrom**, still in Cartesian coordinates.

## Symmetry-adapted mode vectors

ISODISTORT can export symmetry-adapted mode (SAM) vectors in either:

- **Cartesian coordinates** — directly compatible with the phonon eigenvectors
- **Direct (fractional) coordinates** — components along the lattice vectors a, b, c

The files `symmetry_basis` and `symmetry_list` carry no metadata indicating which
convention was used. **You must know which export option you selected in ISODISTORT.**

### How to check

If your structure is cubic or orthogonal (lattice vectors along x, y, z), Cartesian
and Direct vectors are proportional and the distinction may not matter numerically.
For non-orthogonal lattices, the vectors will differ and using the wrong one will
give incorrect SAM decomposition amplitudes and incorrect POSCAR displacements.

### Converting Direct to Cartesian

The conversion for a displacement vector is:

    v_cart[k] = L @ v_direct[k]   for each atom k

where L is the 3×3 matrix of lattice vectors (rows = a, b, c).

In code:

```python
lattice = support.get_lattice_vectors(poscar_lines)
v_cart  = support.direct_to_cartesian_displacement(v_direct, lattice)
```

When using `displace_poscar_sam.py`, pass `--direct` to apply this conversion
automatically.

## POSCAR atomic positions

POSCAR files may store atomic positions in either Direct or Cartesian format
(indicated by line 8). Functions that operate on atomic positions handle both:

- `poscar_to_cartesian()` converts Direct -> Cartesian
- `poscar_to_direct()` converts Cartesian -> Direct
- `displace_poscar_atoms()` always works internally in Cartesian and converts
  back to Direct for output

## Summary table

| Quantity | Coordinate system | Conversion needed? |
|---|---|---|
| Hessian (vasprun.xml) | Cartesian | No |
| Phonon eigenvectors | Cartesian (mass-weighted) | No |
| Real-space displacements | Cartesian (Å) | No |
| SAM vectors (symmetry_basis) | **Unknown — depends on ISODISTORT export** | Maybe |
| POSCAR atomic positions | Direct or Cartesian (line 8 of POSCAR) | Handled internally |
| Input to `displace_poscar_atoms()` | Must be Cartesian | Convert if needed |
