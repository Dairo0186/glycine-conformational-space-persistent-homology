#!/usr/bin/env python3
"""Persistent homology of deterministic subsamples of a periodic torsional grid."""

from __future__ import annotations

import argparse
from pathlib import Path

import numpy as np
import pandas as pd
from ripser import ripser


def one_dimensional_indices(full_size: int, sample_size: int) -> np.ndarray:
    return np.floor(full_size * np.arange(sample_size) / sample_size).astype(int)


def grid_indices(full_size: int, sample_size: int) -> np.ndarray:
    selected = one_dimensional_indices(full_size, sample_size)
    return np.asarray([i * full_size + j for i in selected for j in selected], dtype=int)


def intervals_frame(diagrams: list[np.ndarray]) -> pd.DataFrame:
    rows = []
    for dimension, diagram in enumerate(diagrams):
        for birth, death in diagram:
            persistence = death - birth if np.isfinite(death) else np.inf
            rows.append(
                {
                    "dimension": dimension,
                    "birth_angstrom": birth,
                    "death_angstrom": death,
                    "persistence_angstrom": persistence,
                }
            )
    return pd.DataFrame(rows)


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser()
    parser.add_argument("--matrix", type=Path, required=True)
    parser.add_argument("--grid-size", type=int, default=72)
    parser.add_argument("--resolutions", type=int, nargs="+", default=(25, 36))
    parser.add_argument("--output-dir", type=Path, required=True)
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    distances = np.load(args.matrix, mmap_mode="r")
    expected = args.grid_size**2
    if distances.shape != (expected, expected):
        raise ValueError(f"Expected a {expected} x {expected} distance matrix.")

    args.output_dir.mkdir(parents=True, exist_ok=True)
    for resolution in args.resolutions:
        indices = grid_indices(args.grid_size, resolution)
        reduced = np.asarray(distances[np.ix_(indices, indices)])
        result = ripser(
            reduced,
            distance_matrix=True,
            coeff=2,
            maxdim=2,
        )
        frame = intervals_frame(result["dgms"])
        stem = f"{resolution}x{resolution}"
        np.save(args.output_dir / f"{stem}_indices.npy", indices)
        np.save(args.output_dir / f"{stem}_rmsd.npy", reduced)
        for dimension, diagram in enumerate(result["dgms"]):
            np.save(args.output_dir / f"{stem}_H{dimension}.npy", diagram)
        frame.to_csv(args.output_dir / f"{stem}_intervals.csv", index=False)
        print(f"Completed persistent homology for {stem}.")


if __name__ == "__main__":
    main()
