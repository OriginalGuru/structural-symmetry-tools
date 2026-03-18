#!/usr/bin/env python3
"""
poscar_to_findsym.py

Convert a VASP POSCAR file to a FINDSYM keyword-format input file.

The output uses !latticeBasisVectors (Cartesian, Angstrom) and
!atomPosition (Direct/fractional coordinates). Cartesian POSCARs are
converted to Direct before writing. The centering is always set to P
since POSCAR positions are already fully explicit.

Usage
-----
    python poscar_to_findsym.py [POSCAR] [options]

Arguments
---------
    POSCAR                  Path to POSCAR file (default: ./POSCAR)

Options
-------
    --output FILENAME       Output filename (default: findsym_input.txt)
                            Written to the same directory as the POSCAR
    --lattol FLOAT          Lattice tolerance in Angstrom (default: 0.00001)
    --postol FLOAT          Atomic position tolerance in Angstrom (default: 0.001)
    --title STRING          Title line (default: taken from POSCAR line 1)

Examples
--------
    # Use POSCAR in current directory
    python scripts/poscar_to_findsym.py

    # Use a specific POSCAR file
    python scripts/poscar_to_findsym.py POSCAR_phonon_30_-0.086THz_0.1Ang

    # Looser tolerance for a distorted structure
    python scripts/poscar_to_findsym.py POSCAR_sam_3_0.1Ang --postol 0.01
"""

import os
import sys
import argparse

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from vasp_phonon_tools import support


def write_findsym_input(poscar_lines, output_path,
                        lattol=0.00001, postol=0.001, title=None):
    """
    Convert validated POSCAR lines to a FINDSYM keyword-format input file.

    Parameters
    ----------
    poscar_lines : list of str
    output_path : str
        Full path to write the findsym input file.
    lattol : float
    postol : float
    title : str or None
    """
    poscar_lines = support.validate_poscar(poscar_lines)

    title_line = title if title else poscar_lines[0].strip()
    lattice    = support.get_lattice_vectors(poscar_lines)
    counts     = support.get_atom_counts(poscar_lines)
    species    = poscar_lines[5].split()
    n_atoms    = sum(counts)

    atom_types = []
    for symbol, count in zip(species, counts):
        atom_types.extend([symbol] * count)

    coord_type = poscar_lines[7].split()[0].lower()
    if coord_type == 'cartesian':
        working = support.poscar_to_direct(list(poscar_lines))
    else:
        working = list(poscar_lines)

    positions = []
    for i in range(n_atoms):
        coords = [float(x) for x in working[8 + i].split()]
        positions.append(coords)

    with open(output_path, 'w') as f:
        f.write('!useKeyWords\n\n')
        f.write('!title\n')
        f.write(f'{title_line}\n\n')
        f.write('!latticeTolerance\n')
        f.write(f'{lattol}\n\n')
        f.write('!atomicPositionTolerance\n')
        f.write(f'{postol}\n\n')
        f.write('!latticeBasisVectors\n')
        for vec in lattice:
            f.write(f'  {vec[0]:20.16f}  {vec[1]:20.16f}  {vec[2]:20.16f}\n')
        f.write('\n')
        f.write('!unitCellCentering\nP\n\n')
        f.write('!atomCount\n')
        f.write(f'{n_atoms}\n\n')
        f.write('!atomType\n')
        f.write('  ' + '  '.join(atom_types) + '\n\n')
        f.write('!atomPosition\n')
        for coords in positions:
            f.write(f'  {coords[0]:20.16f}  {coords[1]:20.16f}  {coords[2]:20.16f}\n')
        f.write('\n')


def main():
    parser = argparse.ArgumentParser(
        description='Convert POSCAR to FINDSYM keyword-format input file.'
    )
    parser.add_argument('poscar', nargs='?', default='POSCAR',
                        help='Path to POSCAR file (default: ./POSCAR)')
    parser.add_argument('--output', default='findsym_input.txt',
                        help='Output filename (default: findsym_input.txt)')
    parser.add_argument('--lattol', type=float, default=0.00001,
                        help='Lattice tolerance in Angstrom (default: 0.00001)')
    parser.add_argument('--postol', type=float, default=0.001,
                        help='Atomic position tolerance in Angstrom (default: 0.001)')
    parser.add_argument('--title', default=None,
                        help='Title line (default: from POSCAR line 1)')
    args = parser.parse_args()

    poscar_path = os.path.abspath(args.poscar)
    if not os.path.isfile(poscar_path):
        print(f'Error: POSCAR file not found: {poscar_path}')
        sys.exit(1)

    # Write output alongside the POSCAR
    output_path = os.path.join(os.path.dirname(poscar_path), args.output)

    poscar_lines = support.read_file_lines(poscar_path)
    write_findsym_input(poscar_lines, output_path,
                        lattol=args.lattol, postol=args.postol, title=args.title)

    print(f'Written: {output_path}')
    print(f'  Lattice tolerance        : {args.lattol} Ang')
    print(f'  Atomic position tolerance: {args.postol} Ang')
    print(f'  Run with: findsym {args.output}')


main()
