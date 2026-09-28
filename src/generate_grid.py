#!/usr/bin/env python3
"""Generate the ideal periodic 72x72 Ace-Gly-NMe torsional grid."""

from __future__ import annotations

import argparse
from pathlib import Path

import numpy as np
from rdkit import Chem
from rdkit.Chem import rdMolTransforms


def cyclic_error(measured: float, target: float) -> float:
    """Signed angular error in [-180, 180) degrees."""
    return (measured - target + 180.0) % 360.0 - 180.0


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser()
    parser.add_argument("--reference", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--step", type=float, default=5.0)
    parser.add_argument("--phi", type=int, nargs=4, default=(5, 7, 9, 11))
    parser.add_argument("--psi", type=int, nargs=4, default=(7, 9, 11, 13))
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    if 360.0 % args.step != 0.0:
        raise ValueError("The angular step must divide 360 degrees exactly.")

    supplier = Chem.SDMolSupplier(str(args.reference), removeHs=False)
    reference = next((mol for mol in supplier if mol is not None), None)
    if reference is None or reference.GetNumConformers() != 1:
        raise ValueError("A valid SDF reference geometry with one conformer is required.")

    phi_atoms = tuple(i - 1 for i in args.phi)
    psi_atoms = tuple(i - 1 for i in args.psi)
    angles = np.arange(0.0, 360.0, args.step)
    args.output.parent.mkdir(parents=True, exist_ok=True)

    writer = Chem.SDWriter(str(args.output))
    max_phi_error = 0.0
    max_psi_error = 0.0
    count = 0

    for phi in angles:
        for psi in angles:
            mol = Chem.Mol(reference)
            conf = mol.GetConformer()
            rdMolTransforms.SetDihedralDeg(conf, *phi_atoms, float(phi))
            rdMolTransforms.SetDihedralDeg(conf, *psi_atoms, float(psi))

            phi_measured = rdMolTransforms.GetDihedralDeg(conf, *phi_atoms) % 360.0
            psi_measured = rdMolTransforms.GetDihedralDeg(conf, *psi_atoms) % 360.0
            phi_error = abs(cyclic_error(phi_measured, phi))
            psi_error = abs(cyclic_error(psi_measured, psi))
            max_phi_error = max(max_phi_error, phi_error)
            max_psi_error = max(max_psi_error, psi_error)

            mol.SetProp("grid_index", str(count))
            mol.SetProp("phi_target_deg", f"{phi:.8f}")
            mol.SetProp("psi_target_deg", f"{psi:.8f}")
            mol.SetProp("phi_measured_deg", f"{phi_measured:.8f}")
            mol.SetProp("psi_measured_deg", f"{psi_measured:.8f}")
            writer.write(mol)
            count += 1

    writer.close()
    print(f"Generated {count} configurations.")
    print(f"Maximum cyclic error in phi: {max_phi_error:.10f} degrees")
    print(f"Maximum cyclic error in psi: {max_psi_error:.10f} degrees")


if __name__ == "__main__":
    main()
