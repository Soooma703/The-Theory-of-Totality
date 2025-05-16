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

You can convert the original data files into analysis-ready CSV or pha format using the following scripts:

- `fits_to_csv.py`  
   Converts `.fits` files (e.g., `acisf20207N002_evt2.fits`) to CSV  
   
- `pi_to_pha.py`  
   Converts `.pi` or `.pha` files (e.g., `a2_h3l_cygx1a_p.pha`, `ad15606000s010102_1.pi`) to CSV  

This step prepares the data for use in the main TFT analysis.

### Analysis Tools (`Analysis` folder)

Once the data is converted, you can run the Totality Field Theory spectral validation using:

- `csv_analysis.py`  
   Performs full TFT FFT modeling on pre-converted CSV data (recommended for MAXI J1820+070)

- `pha_analysis.py`  
   Directly analyzes `.pha` or `.pi` spectral files using the TFT framework (recommended for Cygnus X-1, XTE J1550-564)
  
- `alpha_frequency_sensitivity.py`  
    Performs the Energy Mapping and Scale Fixing test (reproduces Figure 1 of the paper).

These scripts will:

1. Simulate a nonlinear wave equation to estimate the fluctuation amplitude **ε** and dominant frequency parameter **a** from the theoretical dynamics of the Totality Field Ψ(X,Φ).
2. Load real black hole X-ray spectral data from either converted CSV files or original `.pha` / `.pi` formats.
3. Uniformly interpolate the energy channels and compute the observed FFT power spectrum.
4. Apply the Totality Field correction model α(E) = 1 + ε · sin(aE) to the spectral data.
5. Recompute the FFT power spectrum of the corrected spectrum and statistically compare it to the original using:
   - **Pearson correlation coefficient**
   - **Empirical p-value estimation using random noise trials**
   - **Bootstrap resampling to calculate 95% confidence intervals**
6. Test the sensitivity of the correlation to phase shifts in α(E) by scanning over a wide range of phase values and generating a diagnostic plot.
7. Output publication-quality plots including:
   - The comparison of observed vs. modeled FFT spectra
   - The phase shift sensitivity curve confirming the deterministic structure predicted by the theory.


---

## ▶️ How to Run

Install required packages:

```bash
pip install numpy matplotlib pandas scipy astropy
