# Totality Field Theory: λ Determination by Energy Stability

This repository contains code to reproduce **Appendix A: Determination of the Self-Interaction Coefficient λ** of the Totality Field Theory (TFT) paper.

We scan a wide range of candidate values of λ to find the unique value that minimizes long-term energy fluctuation of the nonlinear existential field Ψ(X,Φ), following:

    δE = max(E[Ψ(Φ)]) - min(E[Ψ(Φ)])

---

## 🧪 What it Does

1. Simulates the one-dimensional nonlinear wave equation:

       ∂²Ψ/∂Φ² - c² ∂²Ψ/∂X² + λ Ψ³ = 0

2. Calculates total energy at each time step:

       E[Ψ] = ∫ [ (1/2)(∂Ψ/∂Φ)² + (1/2)c²(∂Ψ/∂X)² + (1/4)λΨ⁴ ] dX

3. Selects λ minimizing energy fluctuation.

This reproduces the theoretical procedure of the TFT paper (Appendix A).

---

## 📁 Repository Structure

```plaintext
lambda_determination/
├── lambda_scan.py   ← Coarse scan over wide range (λ ∈ [0.01, 1.0], Δλ = 0.0005)
├── lambda_scan_highres.py  ← Precision scan around theoretical value (λ ∈ [0.06, 0.1], Δλ = 0.0001)
├── README.md               ← This file

---

## ▶️ How to Run

Install dependencies:

```bash
pip install numpy matplotlib tqdm
python lambda_scan.py
