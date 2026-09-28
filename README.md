# Persistent Homology of the Ace–Gly–NMe Conformational Space

Reproducible code and selected results accompanying the study **“Topological Characterization of the Conformational Space of the Ace–Gly–NMe Glycine Model Using Persistent Homology.”**

## Overview

The ideal all-atom conformational space of Ace–Gly–NMe is sampled on a periodic `72 × 72` grid over the backbone dihedral angles φ and ψ. The 5184 configurations (19 atoms, embedded in R^57) are compared using all-atom RMSD after optimal Kabsch alignment. Vietoris–Rips persistent homology is evaluated on deterministic `25 × 25` and `36 × 36` subsamples.

The dominant signature recovered at both resolutions is

```text
(β0, β1, β2) = (1, 2, 1),
```

which is compatible with the topology of the two-dimensional torus `S¹ × S¹`.

## Repository structure

```text
.
├── README.md
├── CITATION.cff
├── LICENSE
├── requirements.txt
├── src/
│   ├── generate_grid.py
│   ├── compute_rmsd_matrix.py
│   ├── persistent_homology.py
│   ├── classical_mds.py
│   └── plot_results.py
└── results/
    ├── README.md
    └── topological_consistency.csv
```

## Reproducible workflow

Use Python 3.12.2 and install the dependencies:

```bash
python -m pip install -r requirements.txt
```

Place a correctly atom-ordered Ace–Gly–NMe reference geometry at
`data/raw/reference_ace_gly_nme.sdf`. The default torsions use the one-based atom numbering reported in the paper:

- φ = D(5, 7, 9, 11)
- ψ = D(7, 9, 11, 13)

Run:

```bash
python src/generate_grid.py \
  --reference data/raw/reference_ace_gly_nme.sdf \
  --output data/processed/glycine_grid_72x72.sdf

python src/compute_rmsd_matrix.py \
  --sdf data/processed/glycine_grid_72x72.sdf \
  --output data/processed/glycine_all_atom_rmsd.npy

python src/persistent_homology.py \
  --matrix data/processed/glycine_all_atom_rmsd.npy \
  --grid-size 72 \
  --resolutions 25 36 \
  --output-dir results/persistence

python src/classical_mds.py \
  --matrix data/processed/glycine_all_atom_rmsd.npy \
  --output results/mds_coordinates.csv

python src/plot_results.py \
  --persistence results/persistence/36x36_intervals.csv \
  --mds results/mds_coordinates.csv \
  --output-dir figures
```

## Large-data policy

The full `5184 × 5184` RMSD matrix is approximately 268 MB as CSV and is intentionally not tracked by GitHub. It can be regenerated with the scripts above. Selected numerical results needed to verify the main topological conclusions are stored in `results/`.

## Main result

| Resolution | Configurations | Dominant H1 classes | Dominant H2 classes | Stable Betti window (Å) |
|---|---:|---:|---:|---|
| 25 × 25 | 625 | 2 | 1 | [0.250913, 1.134388) |
| 36 × 36 | 1296 | 2 | 1 | [0.167302, 1.137197) |

## Software

- Python 3.12.2
- RDKit 2024.09.6
- NumPy
- pandas
- SciPy
- Ripser.py 0.6.12
- Matplotlib

## License

Code is released under the MIT License. Numerical results retain their scholarly attribution and should be cited using `CITATION.cff`.
