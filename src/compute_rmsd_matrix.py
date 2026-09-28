#!/usr/bin/env python3
"""Compute the all-atom pairwise RMSD matrix after proper Kabsch alignment."""

from __future__ import annotations

import argparse
from pathlib import Path

import numpy as np
from numpy.lib.format import open_memmap
from rdkit import Chem


def load_coordinates(path: Path) -> np.ndarray:
    supplier = Chem.SDMolSupplier(str(path), removeHs=False)
    coordinates = []
    atom_count = None
    for mol in supplier:
        if mol is None:
            raise ValueError("Invalid molecule encountered in the SDF file.")
        if atom_count is None:
            atom_count = mol.GetNumAtoms()
        if mol.GetNumAtoms() != atom_count:
            raise ValueError("The atom count is not constant across configurations.")
        coordinates.append(np.asarray(mol.GetConformer().GetPositions(), dtype=np.float64))
    if not coordinates:
        raise ValueError("No configurations were read.")
    return np.stack(coordinates)


def kabsch_rmsd(x: np.ndarray, y: np.ndarray) -> float:
    p = x - x.mean(axis=0)
    q = y - y.mean(axis=0)
    covariance = p.T @ q
    u, _, vt = np.linalg.svd(covariance)
    rotation = vt.T @ u.T
    if np.linalg.det(rotation) < 0.0:
        vt[-1, :] *= -1.0
        rotation = vt.T @ u.T
    aligned = (rotation @ p.T).T
    return float(np.sqrt(np.mean(np.sum((aligned - q) ** 2, axis=1))))


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser()
    parser.add_argument("--sdf", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    xyz = load_coordinates(args.sdf)
    n = xyz.shape[0]
    args.output.parent.mkdir(parents=True, exist_ok=True)
    distances = open_memmap(args.output, mode="w+", dtype=np.float64, shape=(n, n))

    for i in range(n):
        distances[i, i] = 0.0
        for j in range(i + 1, n):
            value = kabsch_rmsd(xyz[i], xyz[j])
            distances[i, j] = value
            distances[j, i] = value
        if (i + 1) % 50 == 0 or i + 1 == n:
            print(f"Completed row {i + 1}/{n}", flush=True)

    distances.flush()
    print(f"Maximum RMSD: {float(np.max(distances)):.8f} angstrom")


if __name__ == "__main__":
    main()
