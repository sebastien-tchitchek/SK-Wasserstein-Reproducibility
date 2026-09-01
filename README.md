# SK Wasserstein Reproducibility Repo

This repository contains the complete software environment and the end-to-end experimental workflow used for the SK Wasserstein evaluation on Ubuntu 24.04.

- `SK-Wasserstein-Reproducibility/` downloads the pinned TTK-ParaView 5.13.0 source from GitHub and builds it together with the included custom TTK 1.3.0 source containing the SK Wasserstein implementation.
- `SK_WASSERSTEIN_REPRODUCTION/` downloads the benchmark data, computes the normalized persistence diagrams and distance matrices, and reproduces the figures, and tables of the "Experimental evaluation" section.

The custom TTK source is included. The installer automatically downloads `topology-tool-kit/ttk-paraview` at tag `v5.13.0` and verifies commit `8b383eeb53821e71f08593ac431a1b3b1855ccac` before building.

## 1. Install ParaView and TTK

From the root of this repository:

```bash
cd SK-Wasserstein-Reproducibility
chmod +x install.sh
. ./install.sh
```

The installer downloads, builds, and installs:

- TTK-ParaView 5.13.0 from tag `v5.13.0` at commit `8b383eeb53821e71f08593ac431a1b3b1855ccac`;
- custom TTK 1.3.0;
- `TTK PersistenceDiagramDistanceMatrix` with the SK Wasserstein implementation.

ParaView can then be launched with:

```bash
paraview
```

## 2. Reproduce the experiments

Return to the repository root and run:

```bash
cd ../SK_WASSERSTEIN_REPRODUCTION
chmod +x install.sh
bash install.sh
```

The workflow performs the following operations:

1. downloads the scalar-field benchmark;
2. computes classical TTK persistence diagrams in all dimensions with the Discrete Morse Sandwich backend;
3. applies one common normalization per collection in `[0,1]^2`;
4. computes the SK Wasserstein, planar-surrogate, and numerical 2-Wasserstein distance matrices;
5. reproduces the figures and tables of the "Experimental evaluation" section.

The generated outputs are available under:

```text
SK_WASSERSTEIN_REPRODUCTION/results/
├── normalized_diagrams
├── matrices
├── analysis
└── computation_time.txt
```

If the computation is interrupted, run the same command again from `SK_WASSERSTEIN_REPRODUCTION/`:

```bash
bash install.sh
```