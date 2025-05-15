import numpy as np
import matplotlib.pyplot as plt
from scipy.stats import pearsonr
from scipy.fft import fft, fftfreq
from scipy.interpolate import interp1d
from astropy.io import fits

# === 固定された理論パラメータ ===
epsilon = 0.106
a_theory = 0.100

# === FITSファイル読み込み（Cygnus X-1など）===
fits_path = "ad15606000s010102_1.pha"  # 実際のファイル名に変更してください
hdul = fits.open(fits_path)
data = hdul[1].data
channel = data["CHANNEL"]
counts = data["COUNTS"]

# === チャネル補間とFFT ===
chan_uniform = np.linspace(channel.min(), channel.max(), 512)
interp = interp1d(channel, counts, kind="linear", fill_value="extrapolate")
counts_uniform = interp(chan_uniform)
fft_obs = fft(counts_uniform - np.mean(counts_uniform))
freq = fftfreq(len(chan_uniform), d=(chan_uniform[1] - chan_uniform[0]))
power_obs = np.abs(fft_obs) ** 2
mask = freq > 0

# === k_RXTE の範囲を変えながら相関係数を計算 ===
k_nominal = 3.02  # keV·s
k_range = np.linspace(k_nominal - 0.02, k_nominal + 0.02, 100)
correlations = []

for k in k_range:
    a_eff = a_theory * k
    alpha = 1 + epsilon * np.sin(a_eff * chan_uniform)
    model_counts = counts_uniform / alpha
    fft_model = fft(model_counts - np.mean(model_counts))
    power_model = np.abs(fft_model) ** 2
    r, _ = pearsonr(power_obs[mask], power_model[mask])
    correlations.append(r)

# === 図示 ===
plt.figure(figsize=(8, 5))
plt.plot(k_range, correlations, label="Pearson correlation")
plt.axvline(
    k_nominal, color="gray", linestyle="--", label="Nominal $k_{\\mathrm{RXTE}}$"
)
plt.axvspan(
    k_nominal - 0.01, k_nominal + 0.01, color="gray", alpha=0.2, label="$\\pm 1\\sigma$"
)
plt.xlabel("Energy mapping coefficient $k_{\\mathrm{RXTE}}$ [keV·s]")
plt.ylabel("Pearson correlation $r$")
plt.ylim(0.98, 1.001)
plt.title("Sensitivity of Spectral Match to Energy Scaling")
plt.grid(True)
plt.legend()
plt.tight_layout()
plt.savefig("alpha_frequency_sensitivity.png")
plt.show()
