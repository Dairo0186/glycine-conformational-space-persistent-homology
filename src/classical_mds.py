#!/usr/bin/env python3
"""Three-dimensional classical MDS from a precomputed RMSD matrix."""

from __future__ import annotations

import argparse
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.sparse.linalg import LinearOperator, eigsh


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser()
    parser.add_argument("--matrix", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--dimensions", type=int, default=3)
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    distances = np.load(args.matrix, mmap_mode="r")
    if distances.ndim != 2 or distances.shape[0] != distances.shape[1]:
        raise ValueError("A square distance matrix is required.")

    squared = np.asarray(distances**2)
    n = squared.shape[0]

    def gram_product(vector: np.ndarray) -> np.ndarray:
        centered = vector - vector.mean()
        product = squared @ centered
        product -= product.mean()
        return -0.5 * product

    gram = LinearOperator((n, n), matvec=gram_product, dtype=np.float64)
    eigenvalues, eigenvectors = eigsh(gram, k=args.dimensions, which="LA")
    order = np.argsort(eigenvalues)[::-1]
    eigenvalues = eigenvalues[order]
    eigenvectors = eigenvectors[:, order]
    if np.any(eigenvalues <= 0):
        raise ValueError("The requested leading MDS eigenvalues are not all positive.")

    coordinates = eigenvectors * np.sqrt(eigenvalues)
    frame = pd.DataFrame(coordinates, columns=[f"MDS{i + 1}" for i in range(args.dimensions)])
    args.output.parent.mkdir(parents=True, exist_ok=True)
    frame.to_csv(args.output, index=False)
    print("Leading positive eigenvalues:", eigenvalues)


if __name__ == "__main__":
    main()
