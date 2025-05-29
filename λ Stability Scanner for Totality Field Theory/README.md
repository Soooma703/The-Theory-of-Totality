# λ Stability Scanner for Totality Field Theory

This repository contains a simulation tool for estimating the critical nonlinear coefficient λ (lambda) at which **stable interference structures** spontaneously emerge in the Totality Field Theory (TFT). This scanning approach was used to identify **λ ≈ 0.57** as a threshold value for coherent structure formation in 3D.

---

## 🔬 Purpose

In Totality Field Theory, the emergence of matter-like structures is driven by nonlinear interference dynamics of the existence field Ψ(X, Φ). This script scans a range of λ values to find where **long-lived, spatially coherent formations** begin to appear. It does this by evaluating fluctuations in the field’s internal energy.

---

## 📈 Methodology

The script simulates the 3D nonlinear wave equation:

∂²Ψ/∂Φ² - c² ∇²Ψ + λ Ψ³ = 0

For each candidate λ:

1. Initialize Ψ with random Gaussian peaks
2. Run time evolution over Nt = 500 steps
3. Compute gradient-based internal energy E(Φ)
4. Record δE = max(E) − min(E) during the evolution
5. Analyze δE and its derivative to find:
   - Local minima (most stable region)
   - Local maxima (instability onset)

The point where **energy fluctuations become minimal** is interpreted as a **critical λ_c** where structure becomes statistically stable.

---

## 📂 Files

| File | Description |
|------|-------------|
| `lambda_calcurate.py` | Main simulation + analysis script. Generates plots and prints out estimated λ_c. |
| `deltaE_vs_lambda.png` | δE vs λ plot showing global energy fluctuation across the scan. |
| `d_deltaE_d_lambda.png` | First derivative plot showing rate of change of energy fluctuation. |
| `d_deltaE_d_lambda_peaks.png` | Peaks in the derivative (local instabilities). |

---

## ⚙️ How to Run

Make sure you have Python 3.8+ and the following libraries:

```bash
pip install numpy matplotlib scipy tqdm numba
python lambda_calcurate.py
```
---

## ⚙️ Parameter Configuration

You can modify the range of λ values scanned by adjusting the following parameter in the script:

```python
lambda_range = (0, 3.5)
```
