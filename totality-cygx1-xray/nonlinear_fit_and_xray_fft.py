import numpy as np
import matplotlib.pyplot as plt
from scipy.interpolate import interp1d
from scipy.stats import pearsonr
import pandas as pd
from astropy.io import fits

# -------------------------------
# 1. Nonlinear Wave Simulation → ε推定
# -------------------------------
print("🔄 Running nonlinear wave simulation to estimate ε...")

L = 100
Nx = 500
dx = L / Nx
x = np.linspace(0, L, Nx)
c = 1.0
dt = 0.05
Nt = 200
lambda_nl = 0.1

u = np.exp(-((x - L / 2) ** 2) / 4)
u_new = u.copy()
u_old = u.copy()
std_amplitudes = []

for t in range(Nt):
    nonlinear_term = lambda_nl * u**3
    u_new[1:-1] = (
        2 * u[1:-1]
        - u_old[1:-1]
        + (c * dt / dx) ** 2 * (u[2:] - 2 * u[1:-1] + u[:-2])
        - dt**2 * nonlinear_term[1:-1]
    )
    u_old = u.copy()
    u = u_new.copy()
    std_amplitudes.append(np.std(u))

epsilon_estimate = np.mean(std_amplitudes)
print(f"✅ Estimated ε from nonlinear wave simulation: {epsilon_estimate:.5f}")

# -------------------------------
# 2. Load FITS X-ray Spectrum
# -------------------------------
print("📥 Loading X-ray spectrum FITS file...")

fits_path = "a2_h3l_cygx1a_p.pha"  # 適宜変更
hdul = fits.open(fits_path)
data = hdul[1].data

channel = data["CHANNEL"]
counts = data["COUNTS"]
df = pd.DataFrame({"channel": channel, "counts": counts})
df.to_csv("cygx1_spectrum.csv", index=False)
print("✅ Saved to cygx1_spectrum.csv")

# -------------------------------
# 3. FFT Preparation (uniform spacing)
# -------------------------------
chan_uniform = np.linspace(df["channel"].min(), df["channel"].max(), 512)
interp = interp1d(df["channel"], df["counts"], kind="linear", fill_value="extrapolate")
counts_uniform = interp(chan_uniform)

fft_obs = np.fft.fft(counts_uniform - np.mean(counts_uniform))
freq = np.fft.fftfreq(len(chan_uniform), d=(chan_uniform[1] - chan_uniform[0]))
power_obs = np.abs(fft_obs) ** 2
mask = freq > 0

# -------------------------------
# 4. Fit a using estimated ε
# -------------------------------
print("🔎 Searching for best-fit a using ε from simulation...")


def alpha_sin(E, epsilon, a):
    return 1 + epsilon * np.sin(a * E)


a_values = np.linspace(0.01, 1.0, 500)
best_corr = -1
best_a = None

for a in a_values:
    alpha = alpha_sin(chan_uniform, epsilon_estimate, a)
    model_counts = counts_uniform / alpha
    fft_model = np.fft.fft(model_counts - np.mean(model_counts))
    power_model = np.abs(fft_model) ** 2
    corr, _ = pearsonr(power_obs[mask], power_model[mask])
    if corr > best_corr:
        best_corr = corr
        best_a = a

print(f"🔍 Best-fit a = {best_a:.5f}")
print(f"📊 Max Pearson correlation (ε = {epsilon_estimate:.5f}): {best_corr:.6f}")

# -------------------------------
# 5. Visualization
# -------------------------------
alpha_best = alpha_sin(chan_uniform, epsilon_estimate, best_a)
model_counts_best = counts_uniform / alpha_best
fft_model_best = np.fft.fft(model_counts_best - np.mean(model_counts_best))
power_model_best = np.abs(fft_model_best) ** 2

plt.figure(figsize=(10, 6))
plt.plot(freq[mask], power_obs[mask], label="Observed FFT Power")
plt.plot(
    freq[mask],
    power_model_best[mask],
    "--",
    label=f"Model: ε={epsilon_estimate:.5f}, a={best_a:.3f}",
)
plt.xlabel("Frequency [1/channel]")
plt.ylabel("Power")
plt.title("FFT Spectrum of Cygnus X-1 with Simulated ε Model")
plt.grid(True)
plt.legend()
plt.tight_layout()
plt.savefig("FFT_Xray_CygX1_alpha_simulated_eps.png")
plt.show()
