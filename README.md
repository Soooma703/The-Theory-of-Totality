# Totality Field Theory (TFT) - Computational Validation Package

This repository contains the full computational framework to reproduce the key numerical validation results of the **Totality Field Theory (TFT)** as presented in:

> Ohno (2025). The Totality Field Theory: Unifying Spacetime and Observation without Time.  

The code is fully open source and aims to maximize transparency, reproducibility, and testability of the theoretical predictions.

---

## 🧩 Project Structure

### 📂 1. Black Hole X-ray FFT Analysis

Located in: [`TFT_Xray_Analysis/`](./TFT_Xray_Analysis/)

This module performs the **definitive spectral validation** of the TFT by:

- Simulating the nonlinear wave equation to derive ε and a
- Applying TFT correction to observed black hole X-ray spectra
- Comparing FFT power spectra between model and observation
- Performing bootstrap confidence and phase-shift sensitivity tests

Data is provided for:
- Cygnus X-1
- MAXI J1820+070
- XTE J1550-564

This fully reproduces Figures 2-5 of the paper.

---

### 📂 2. Lambda Stability Scan

Located in: [`TFT_Xray_LambdaScan/`](./TFT_Xray_LambdaScan/)

This module reproduces **Appendix A** of the paper by:

- Scanning candidate values of λ
- Running long-time field simulations to calculate total energy fluctuation
- Determining the unique energy-stable self-interaction coefficient λ = 0.0860

This provides the only first-principles prediction of λ without external fitting.

---

## 🔧 Requirements

Each module contains its own requirements.txt.

You will need:
- Python ≥3.8
- numpy, matplotlib, scipy, astropy, tqdm (depending on module)

---

## 📜 Notes

This repository was designed for maximum reproducibility for research and peer review purposes.  
If you use this code, please cite the original paper:

> Ohno (2025). The Totality Field Theory: Unifying Spacetime and Observation without Time.

We welcome questions, feedback, or collaborative proposals.

