# Totality Field Theory: Simulated ε Applied to Black Hole X-ray FFT Spectra

This repository contains code to validate the **Totality Field Theory (TFT)** by:

1. Simulating a nonlinear wave equation to derive the amplitude fluctuation parameter **ε**
2. Applying this ε to model the fine FFT power spectrum of black hole X-ray emission spectra
3. Comparing the resulting modulation to actual observed **X-ray spectral data** for three independent black hole systems

---

## 🧪 What it Does

- Numerically integrates:　∂²Ψ/∂Φ² − c² ∇²Ψ + λ Ψ³ = 0
- Estimates the time-averaged amplitude fluctuation as:　ε = mean[σ(Ψ(Φ))]
- Uses this ε to construct a correction function:　α(E) = 1 + ε · sin(aE)
- Optimizes the parameter `a` to **maximize the Pearson correlation** between the observed and modeled FFT power spectra.

---

## 📚 Data Sources

This repository uses X-ray spectral data of three well-studied black hole binary systems, publicly available from the **NASA HEASARC (High Energy Astrophysics Science Archive Research Center)**:

1. **Cygnus X-1**  
   - Format: `.pha` (Pulse Height Analyzer spectrum)  
   - [Data file](https://heasarc.gsfc.nasa.gov/FTP/heao1/data/a2/spectra/a2_h3l_cygx1a_p.pha.Z): `a2_h3l_cygx1a_p.pha`

2. **MAXI J1820+070**  
   - Format: `.fits` (event data or spectral data)  
   - [Data file](https://heasarc.gsfc.nasa.gov/FTP/chandra/data/byobsid/7/20207/primary/acisf20207N002_evt2.fits.gz): `acisf20207N002_evt2.fits`

3. **XTE J1550-564**  
   - Format: `.pi` (Pulse Invariant spectrum, equivalent to `.pha`)  
   - [Data file](https://heasarc.gsfc.nasa.gov/FTP/asca/data/rev2/15606000/spectra/ad15606000s010102_1.pi.gz): `ad15606000s010102_1.pi`

These files contain energy-resolved photon count data used for spectral analysis in high-energy astrophysics.

**Please cite the original missions and instruments (e.g., RXTE, Chandra, ASCA) when publishing results.**

---

## 🛠️ File Conversion and Analysis

### File Conversion Tools (`file conversion codes` folder)

You can convert the original data files into analysis-ready CSV format using the following scripts:

- `fits_to_csv.py`  
   Converts `.fits` files (e.g., `acisf20207N002_evt2.fits`) to CSV  
   
- `pi_to_csv.py`  
   Converts `.pi` or `.pha` files (e.g., `a2_h3l_cygx1a_p.pha`, `ad15606000s010102_1.pi`) to CSV  

This step prepares the data for use in the main TFT analysis.

### Analysis Tools (`Analysis` folder)

Once the data is converted, you can run the Totality Field Theory spectral validation using:

- `csv_analysis.py`  
   Performs full TFT FFT modeling on pre-converted CSV data (recommended for MAXI J1820+070)

- `pha_analysis.py`  
   Directly analyzes `.pha` or `.pi` spectral files using the TFT framework (recommended for Cygnus X-1, XTE J1550-564)

These scripts will:
1. Derive ε and a from the nonlinear wave simulation
2. Apply the α(E) correction to the spectrum
3. Compute FFT power spectra
4. Perform statistical validation including Pearson correlation and bootstrap confidence intervals
5. Generate ready-to-publish plots for scientific papers

---

## ▶️ How to Run

Install required packages:

```bash
pip install -r requirements.txt
