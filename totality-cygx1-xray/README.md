# Totality Field Theory: Simulated ε Applied to X-ray FFT Spectrum

This repository contains code to validate the Totality Field Theory by:

1. Simulating a nonlinear wave equation to derive the amplitude fluctuation ε
2. Applying this ε to model the fine FFT power spectrum of the black hole Cygnus X-1
3. Comparing the resulting modulation with actual observed X-ray spectrum data

## 🧪 What it Does

- Numerically integrates:  
  $$ \frac{\partial^2 \Psi}{\partial \Phi^2} - c^2 \nabla^2 \Psi + \lambda \Psi^3 = 0 $$
- Estimates:
  $$ \epsilon = \langle \sigma[\Psi(t)] \rangle $$
- Fits:
  $$ \alpha(E) = 1 + \epsilon \sin(aE) $$
- Optimizes $a$ to maximize the FFT power spectrum correlation

## ▶️ How to Run

```bash
pip install -r requirements.txt
python nonlinear_fit_and_xray_fft.py
