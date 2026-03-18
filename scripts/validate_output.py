#!/usr/bin/env python3
"""
validate_output.py

Validate the output of get_phonon_mode_charges.py against a reference file
by comparing only the physically invariant quantities:

    - Phonon frequency (signed, THz)
    - Force constant (eV/Å²)
    - Reduced mass (AMU)
    - |Z*| (mode effective charge magnitude, e)

Quantities that are NOT checked because they are not invariant across
platforms or numpy/LAPACK versions:

    - Z*_x, Z*_y, Z*_z  (sign depends on eigenvector sign convention)
    - SAM amplitudes     (sign and mixing within degenerate subspaces vary)
    - SAM indices        (ordering within degenerate subspaces varies)

See docs/coordinate_conventions.md for a full explanation.

Usage
-----
    python validate_output.py <output_file> <reference_file> [--tol 1e-6]

Arguments
---------
    output_file     Output from get_phonon_mode_charges.py to validate
    reference_file  Reference output to compare against
    --tol           Absolute tolerance for floating point comparison
                    (default: 1e-6)

Exit codes
----------
    0   All checks passed
    1   One or more checks failed
    2   Could not parse one or both files

Example
-------
    python ../../scripts/validate_output.py output.txt expected_output.txt
"""

import sys
import argparse


# ---------------------------------------------------------------------------
# Parsing
# ---------------------------------------------------------------------------

def parse_summary_table(lines):
    """
    Extract the first summary table from output lines.

    The summary table is the first block of mode rows (lines starting with
    '+' followed by digits) before the detailed per-mode section begins.
    Returns a list of dicts, one per mode, with keys:
        index, frequency, force_constant, reduced_mass, zstar_mag
    """
    modes = []
    in_summary = False

    for line in lines:
        stripped = line.strip()

        # Summary table rows start with +NN
        if stripped.startswith('+') and stripped[1:3].isdigit():
            parts = stripped.split()
            if len(parts) < 8:
                continue
            try:
                modes.append({
                    'index':          int(parts[0].lstrip('+')),
                    'frequency':      float(parts[1]),
                    'force_constant': float(parts[2]),
                    'reduced_mass':   float(parts[3]),
                    'zstar_mag':      float(parts[7]),
                })
                in_summary = True
            except (ValueError, IndexError):
                continue

        elif in_summary and stripped == '':
            # First blank line after summary rows signals end of summary table
            break

    return modes


def parse_file(path):
    """
    Read an output file and return parsed mode data.

    Returns
    -------
    list of dict, or None if parsing failed
    """
    try:
        with open(path, 'r') as f:
            lines = f.readlines()
    except FileNotFoundError:
        print(f"Error: file not found: {path}")
        return None

    modes = parse_summary_table(lines)

    if not modes:
        print(f"Error: no mode data found in {path}")
        print("       Is this a valid get_phonon_mode_charges.py output file?")
        return None

    return modes


# ---------------------------------------------------------------------------
# Comparison
# ---------------------------------------------------------------------------

def compare(output_modes, reference_modes, tol):
    """
    Compare two lists of parsed mode dicts.

    Returns
    -------
    passed : bool
    messages : list of str
    """
    messages = []
    passed = True

    if len(output_modes) != len(reference_modes):
        messages.append(
            f"FAIL: mode count mismatch — "
            f"output has {len(output_modes)}, reference has {len(reference_modes)}"
        )
        return False, messages

    fields = [
        ('frequency',      'THz'),
        ('force_constant', 'eV/Å²'),
        ('reduced_mass',   'AMU'),
        ('zstar_mag',      'e'),
    ]

    for i, (out, ref) in enumerate(zip(output_modes, reference_modes)):
        mode_num = ref['index']
        for field, unit in fields:
            diff = abs(out[field] - ref[field])
            if diff > tol:
                messages.append(
                    f"FAIL: mode {mode_num:2d}  {field:<16s}  "
                    f"output={out[field]:+.12f}  "
                    f"reference={ref[field]:+.12f}  "
                    f"diff={diff:.3e}  [{unit}]"
                )
                passed = False

    return passed, messages


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

def main():
    parser = argparse.ArgumentParser(
        description="Validate get_phonon_mode_charges.py output against a reference."
    )
    parser.add_argument('output_file',    help="Output file to validate")
    parser.add_argument('reference_file', help="Reference file to compare against")
    parser.add_argument('--tol', type=float, default=1e-6,
                        help="Absolute tolerance (default: 1e-6)")
    args = parser.parse_args()

    print(f"Validating : {args.output_file}")
    print(f"Reference  : {args.reference_file}")
    print(f"Tolerance  : {args.tol:.2e}")
    print()

    output_modes    = parse_file(args.output_file)
    reference_modes = parse_file(args.reference_file)

    if output_modes is None or reference_modes is None:
        sys.exit(2)

    print(f"Modes found in output    : {len(output_modes)}")
    print(f"Modes found in reference : {len(reference_modes)}")
    print()
    print("Checking: frequency, force_constant, reduced_mass, |Z*|")
    print("Skipping: Z*_x/y/z, SAM amplitudes (not invariant across platforms)")
    print()

    passed, messages = compare(output_modes, reference_modes, args.tol)

    if passed:
        print(f"PASSED — all {len(output_modes)} modes within tolerance {args.tol:.2e}")
        sys.exit(0)
    else:
        for msg in messages:
            print(msg)
        print()
        print(f"FAILED — {len(messages)} check(s) did not pass")
        sys.exit(1)


main()
