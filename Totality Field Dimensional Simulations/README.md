# Totality Field Dimensional Simulations

This repository contains numerical simulations of nonlinear wave dynamics across 1D to 5D spatial domains, developed within the framework of the **Totality Field Theory (TFT)**. These simulations aim to explore how stable structures — analogous to particles — can spontaneously emerge from interference patterns in a fundamental existence field, denoted as $\Psi(X, \Phi)$.

---

## 🔬 Purpose

The primary objective is to investigate **dimensional dependence of structure formation** in nonlinear existence fields governed by a generalized wave equation:

\[
\frac{\partial^2 \Psi}{\partial \Phi^2} = c^2 \nabla^2 \Psi - \lambda \Psi^3
\]

This system mimics core principles of TFT, where **space, time, matter, and even consciousness** are emergent from the evolving interference topology of the field $\Psi$.

---

## 🧠 Files Overview

Each script corresponds to a simulation in different spatial dimensions. They share a common structure:

| File | Description |
|------|-------------|
| `1D.py` | Simulates and analyzes interference structures in 1 spatial dimension. Includes fractal analysis of localized peaks. |
| `2D.py` | Extends the model to 2D grids, generating 2D structure maps and energy histograms. |
| `3D.py` | Simulates volumetric (voxel-based) structure formation. Outputs fractal estimates and energy distributions. |
| `4D.py` | Simulates 4D field evolution with projection analysis. Computationally intensive. |
| `5D.py` | Explores the breakdown of coherent structure formation in 5D — useful for testing dimensional stability in TFT. |

---

## 📈 Features

- **Spontaneous Structure Emergence**: From random Gaussian initial fields, stable localized peaks emerge without external forcing.
- **Peak Tracking**: Counts and records number of coherent peaks per timestep.
- **Structure Analysis**:
  - Volume (cell count)
  - Energy (integrated amplitude squared)
  - Fractal dimension estimate \( D_f = \log(V) / \log(\text{BoxSize}) \)
- **Visualization**: Automatically generates histograms of structure distributions.

---

## ⚙️ How to Run

Each script is standalone and can be run with Python 3.8+:

```bash
python 3D.py

