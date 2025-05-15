import numpy as np
import matplotlib.pyplot as plt
from scipy.interpolate import interp1d
from scipy.stats import pearsonr
from scipy.fft import fft, fftfreq
import pandas as pd

# === 1. Nonlinear Wave Simulation → ε と a を理論から導出 ===
L = 100
Nx = 500
dx = L / Nx
x = np.linspace(0, L, Nx)
c = 1.0
dt = 0.05
Nt = 200
lambda_nl = 0.0860

u = np.exp(-((x - L / 2) ** 2) / 4)
u_new = u.copy()
u_old = u.copy()

std_amplitudes = []
psi_time_series = []
probe_index = Nx // 2  # X固定点

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
    psi_time_series.append(u[probe_index])

epsilon_estimate = np.mean(std_amplitudes)
fft_vals = fft(psi_time_series - np.mean(psi_time_series))
freqs = fftfreq(Nt, d=dt)
power = np.abs(fft_vals) ** 2
mask = freqs > 0
a_estimate = freqs[mask][np.argmax(power[mask])]

print(f"Estimated ε = {epsilon_estimate:.5f}")
print(f"Estimated a = {a_estimate:.5f}")

# === 2. CSVから X線スペクトルを読み込み ===
csv_path = "maxi1820_spectrum.csv"  # ← CSVパスに変更
df = pd.read_csv(csv_path)
channel = df["CHANNEL"].values
counts = df["COUNTS"].values

# === 3. チャネル補間 → FFT ===
chan_uniform = np.linspace(channel.min(), channel.max(), 512)
interp = interp1d(channel, counts, kind="linear", fill_value="extrapolate")
counts_uniform = interp(chan_uniform)

fft_obs = fft(counts_uniform - np.mean(counts_uniform))
freq = fftfreq(len(chan_uniform), d=(chan_uniform[1] - chan_uniform[0]))
power_obs = np.abs(fft_obs) ** 2
mask_freq = freq > 0


# === 4. Totalityモデル α(E) 適用 ===
def alpha_sin(E, epsilon, a):
    return 1 + epsilon * np.sin(a * E)


alpha = alpha_sin(chan_uniform, epsilon_estimate, a_estimate)
model_counts = counts_uniform / alpha
fft_model = fft(model_counts - np.mean(model_counts))
power_model = np.abs(fft_model) ** 2

# === 5. 相関とp値 ===
true_corr, _ = pearsonr(power_obs[mask_freq], power_model[mask_freq])
print(f"Pearson correlation (theoretical): {true_corr:.6f}")

n_trials = 1000
random_corrs = []
for _ in range(n_trials):
    noise = np.random.normal(0, 1, size=counts_uniform.shape)
    fft_noise = fft(noise - np.mean(noise))
    power_noise = np.abs(fft_noise) ** 2
    corr, _ = pearsonr(power_obs[mask_freq], power_noise[mask_freq])
    random_corrs.append(corr)

random_corrs = np.array(random_corrs)
p_empirical = np.sum(random_corrs >= true_corr) / n_trials
print(f"Empirical p-value (from {n_trials} trials): p = {p_empirical:.5f}")

# === 6. 図示 ===
plt.figure(figsize=(10, 6))
plt.plot(freq[mask_freq], power_obs[mask_freq], label="Observed FFT Power")
plt.plot(freq[mask_freq], power_model[mask_freq], "--", label="Model Prediction")
plt.xlabel("Frequency [1/channel]")
plt.ylabel("Power")
plt.title("FFT Spectrum with Totality Model (from CSV)")
plt.legend()
plt.grid(True)
plt.tight_layout()
plt.savefig("FFT_Xray_MAXI_J1820+070_theoretical_alpha.png")
plt.show()
# -------------------------------
# 6.1 ブートストラップによる信頼区間の推定
# -------------------------------
bootstrap_trials = 1000
bootstrap_corrs = []
rng = np.random.default_rng()

# ブートストラップ：mask_freqで選ばれたインデックスからリサンプリング
indices = np.arange(np.sum(mask_freq))
for _ in range(bootstrap_trials):
    resample_idx = rng.choice(indices, size=len(indices), replace=True)
    corr_bs, _ = pearsonr(
        power_obs[mask_freq][resample_idx], power_model[mask_freq][resample_idx]
    )
    bootstrap_corrs.append(corr_bs)

# 平均と信頼区間（95% CI）
bootstrap_corrs = np.array(bootstrap_corrs)
mean_corr = np.mean(bootstrap_corrs)
ci_lower = np.percentile(bootstrap_corrs, 2.5)
ci_upper = np.percentile(bootstrap_corrs, 97.5)

print(f"Bootstrap mean correlation: {mean_corr:.6f}")
print(f"95% confidence interval: [{ci_lower:.6f}, {ci_upper:.6f}]")
# -------------------------------
# 6.2 全位相ズレテスト（phi ∈ [0, 200]）
# -------------------------------
phase_range = np.linspace(0, 100, 10000)  # 0〜200ラジアンを1000分割
abs_correlations = []

for phi in phase_range:
    alpha_phase = 1 + epsilon_estimate * np.sin(a_estimate * chan_uniform + phi)
    model_counts_phase = counts_uniform / alpha_phase
    fft_model_phase = fft(model_counts_phase - np.mean(model_counts_phase))
    power_model_phase = np.abs(fft_model_phase) ** 2
    r, _ = pearsonr(power_obs[mask_freq], power_model_phase[mask_freq])
    abs_correlations.append(np.abs(r))  # 相関の絶対値を記録

# 図示
plt.figure(figsize=(10, 6))
plt.plot(phase_range, abs_correlations, color="blue")
plt.xlabel("Phase shift $\\phi$ [rad]")
plt.ylabel("Absolute Pearson correlation $|r|$")
plt.ylim(0.999, 1.001)  # 精密さを強調
plt.title(
    "Correlation vs Phase Shift in $\\alpha(E) = 1 + \\epsilon \\sin(aE + \\phi)$"
)
plt.grid(True)
plt.tight_layout()
plt.savefig("alpha_phase_shift_sensitivity.png")
plt.show()
