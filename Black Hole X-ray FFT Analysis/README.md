# Totality Field Theory: Simulated ε Applied to X-ray FFT Spectrum

This repository contains code to validate the **Totality Field Theory** by:

1. Simulating a nonlinear wave equation to derive the amplitude fluctuation parameter ε
2. Applying this ε to model the fine FFT power spectrum of the black hole **Cygnus X-1**
3. Comparing the resulting modulation to actual observed **X-ray spectral data** (FITS)

---

## 🧪 What it Does

- Numerically integrates:　∂²Ψ/∂Φ² − c² ∇²Ψ + λ Ψ³ = 0
- Estimates the time-averaged amplitude fluctuation as:　ε = ϵ=mean[σ(Ψ(t))]
- Uses this ε to construct a correction function:　α(E) = 1 + ε · sin(aE)
- Optimizes the parameter `a` to **maximize the Pearson correlation** between the observed and modeled FFT power spectra.

---

## 📚 Data Source

The X-ray spectrum of the black hole binary system **Cygnus X-1** used in this project is publicly available from the **NASA HEASARC (High Energy Astrophysics Science Archive Research Center)** archive.

- [Main Archive and Examples (XSPEC)](https://heasarc.gsfc.nasa.gov/xamin/)
- [Direct file link (FITS format)](https://heasarc.gsfc.nasa.gov/FTP/heao1/data//a2/spectra/a2_h3l_cygx1a_p.pha.Z)

File: `a2_h3l_cygx1a_p.pha`

This file contains energy-resolved photon count data, typically used for spectral fitting in high-energy astrophysics. It is based on RXTE (Rossi X-ray Timing Explorer) PCA observations of Cygnus X-1.

For publication use, please cite the original data source accordingly (e.g., RXTE mission, PCA instrument).

---


## ▶️ How to Run

Install required packages:

```bash
pip install -r requirements.txt
python nonlinear_fit_and_xray_fft.py

