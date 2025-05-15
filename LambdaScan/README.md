# Totality Field Theory: λ Determination by Energy Stability

This repository contains code to reproduce **Appendix A: Determination of the Self-Interaction Coefficient λ** of the Totality Field Theory (TFT) paper.

We scan a wide range of candidate values of λ to find the unique value that minimizes long-term energy fluctuation of the nonlinear existential field Ψ(X,Φ), following:

\[
\delta E = \max_\Phi E[\Psi(\Phi)] - \min_\Phi E[\Psi(\Phi)]
\]

---

## 🧪 What it Does

1. Simulates the one-dimensional nonlinear wave equation:
\[
\frac{\partial^2 \Psi}{\partial \Phi^2} - c^2 \frac{\partial^2 \Psi}{\partial X^2} + \lambda \Psi^3 = 0
\]
2. Calculates total energy at each time step:
\[
E[\Psi] = \int \left[ \frac{1}{2} \left( \frac{\partial \Psi}{\partial \Phi} \right)^2 + \frac{1}{2} c^2 \left( \frac{\partial \Psi}{\partial X} \right)^2 + \frac{1}{4} \lambda \Psi^4 \right] dX
\]
3. Selects λ minimizing energy fluctuation.

This reproduces the theoretical procedure of the TFT paper (Appendix A).

---

## ▶️ How to Run

Install dependencies:
```bash
pip install numpy matplotlib tqdm
python lambda_scan.py


