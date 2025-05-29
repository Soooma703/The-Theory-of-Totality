# Detailed Structure Analysis in Totality Field Theory (TFT)

This repository provides a comprehensive simulation and analysis pipeline for the nonlinear evolution of the Totality Field $\Psi(X, \Phi)$ in 3D. It focuses on identifying and analyzing particle-like interference structures that emerge at the critical nonlinearity coefficient **λ = 0.57**, as predicted by Totality Field Theory.

---

## 📌 Purpose

This script allows you to:

- Simulate the nonlinear wave evolution of the existence field Ψ over time
- Detect localized high-amplitude interference regions
- Quantify and analyze:
  - Structure volumes
  - Peak amplitudes
  - Energy distributions
  - Spatial separations
  - Centroid trajectories of dominant structures

---

## 🧪 Governing Equation

The field evolves according to the 3D nonlinear wave equation:
∂²Ψ/∂Φ² = c² ∇²Ψ − λ Ψ³

where `Ψ = Ψ(x, y, z, Φ)` and λ = 0.57 is the critical value empirically identified for stable structure emergence in 3D.

---

## 🔍 Features

- **3D wave evolution** over 2000 timesteps
- **Interference peak detection** using smoothed thresholding
- **Volume and amplitude histograms** of all localized structures
- **Energy computation** per structure: `E ∝ A² × Volume`
- **Nearest-neighbor distance analysis**
- **Trajectory tracking** of the largest structure's centroid

---

## 📂 Output Files

| Filename | Description |
|----------|-------------|
| `peak_count_evolution.png` | Number of localized peaks vs time |
| `visualize_slice.png` | 2D slice visualization of final field |
| `volume_distribution.png` | Histogram of structure volumes |
| `amplitude_distribution.png` | Histogram of structure amplitudes |
| `energy_distribution.png` | Histogram of structure energies |
| `structure_distance_distribution.png` | Histogram of nearest-neighbor distances |
| `largest_structure_trajectory.png` | Motion of the largest structure over time |

---

## 🚀 How to Run

Install dependencies:

```bash
pip install numpy matplotlib scipy tqdm numba
python detailed_analysis.py
```
The simulation will take a few minutes (depending on machine power). Outputs will be saved in the current directory.

## ⚙️ Parameter Options

You can modify these key parameters at the top of `detailed_analysis.py` to control simulation behavior:

- `lambda_c = 0.57`  
  → Nonlinearity coefficient. This controls the self-interaction strength of the field Ψ.

- `N = 256`  
  → Number of grid points along each spatial axis. Increasing this improves resolution but increases computation time.

- `Nt = 2000`  
  → Number of simulation time steps (Φ iterations). More steps reveal long-term structure dynamics.

- `threshold_ratio = 0.3`  
  → Relative amplitude threshold for identifying localized structures. Expressed as a fraction of the maximum amplitude.

- `smooth_sigma = 1.0`  
  → Standard deviation for Gaussian smoothing applied before structure detection. Helps reduce noise and isolate meaningful peaks.

---

## 🧠 Scientific Context

According to the Totality Field Theory (TFT), the field Ψ(X, Φ) is the fundamental substrate from which all physical phenomena emerge. In this framework:

- Localized interference nodes in Ψ are interpreted as **proto-particles**, arising without external potentials.
- These stable structures form when nonlinear self-interaction balances with dispersive spreading.
- The dimensionality of the simulation plays a critical role: in particular, **3D supports stable structure formation**, while higher dimensions tend to destabilize coherence.

This simulation explores how and why such structures form, evolve, and persist over time — providing a computational probe into the foundations of matter and geometry.

---

## 📝 License

This project is released under the [MIT License](https://en.wikipedia.org/wiki/MIT_License).

You are free to use, modify, and distribute this code for research, education, or derivative work, provided proper attribution is given to the author and original framework (Totality Field Theory).


