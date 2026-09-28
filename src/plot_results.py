#!/usr/bin/env python3
"""Create the persistence diagram and three-dimensional MDS representation."""

from __future__ import annotations

import argparse
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser()
    parser.add_argument("--persistence", type=Path, required=True)
    parser.add_argument("--mds", type=Path, required=True)
    parser.add_argument("--output-dir", type=Path, required=True)
    return parser.parse_args()


def plot_persistence(frame: pd.DataFrame, output: Path) -> None:
    finite = frame[np.isfinite(frame["death_angstrom"])].copy()
    maximum = float(max(finite["birth_angstrom"].max(), finite["death_angstrom"].max()))
    colors = {0: "#1f77b4", 1: "#ff7f0e", 2: "#2ca02c"}

    fig, ax = plt.subplots(figsize=(6.2, 5.4))
    ax.plot([0, maximum], [0, maximum], "--", color="0.55", linewidth=1)
    for dimension, group in finite.groupby("dimension"):
        ax.scatter(
            group["birth_angstrom"],
            group["death_angstrom"],
            s=28,
            alpha=0.8,
            color=colors.get(int(dimension), "black"),
            label=rf"$H_{int(dimension)}$",
        )
    ax.set_xlabel("Birth (Å)")
    ax.set_ylabel("Death (Å)")
    ax.legend(frameon=False)
    ax.set_aspect("equal", adjustable="box")
    fig.tight_layout()
    fig.savefig(output, dpi=300)
    plt.close(fig)


def plot_mds(frame: pd.DataFrame, output: Path) -> None:
    required = ["MDS1", "MDS2", "MDS3"]
    if not set(required).issubset(frame.columns):
        raise ValueError("The MDS file must contain MDS1, MDS2, and MDS3.")

    fig = plt.figure(figsize=(7, 6))
    ax = fig.add_subplot(111, projection="3d")
    ax.scatter(frame["MDS1"], frame["MDS2"], frame["MDS3"], s=3, alpha=0.55)
    ax.set_xlabel("MDS1")
    ax.set_ylabel("MDS2")
    ax.set_zlabel("MDS3")
    fig.tight_layout()
    fig.savefig(output, dpi=300)
    plt.close(fig)


def main() -> None:
    args = parse_args()
    args.output_dir.mkdir(parents=True, exist_ok=True)
    plot_persistence(pd.read_csv(args.persistence), args.output_dir / "persistence_36x36.png")
    plot_mds(pd.read_csv(args.mds), args.output_dir / "mds_3d.png")


if __name__ == "__main__":
    main()
