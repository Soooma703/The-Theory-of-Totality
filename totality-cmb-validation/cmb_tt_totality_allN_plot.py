import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
from scipy.optimize import minimize
from scipy.stats import pearsonr, chi2
import re

# -------------------------------
# 1. データ読み込み
# -------------------------------
file_path = "COM_PowerSpect_CMB-TT-full_R3.01.txt"

with open(file_path, "r") as f:
    lines = f.readlines()
data_lines = [
    line.strip() for line in lines if re.match(r"^\s*\d+(\.\d+)?[eE]?[+-]?\d*\s", line)
]
data = [list(map(float, line.split())) for line in data_lines]
df = pd.DataFrame(data, columns=["ell", "Dl", "dDl_low", "dDl_high"])
df["dDl"] = (df["dDl_low"] + df["dDl_high"]) / 2

# 対象データ：ℓ > 2
df = df[df["ell"] > 2].reset_index(drop=True)
ell = df["ell"].values
Dl_obs = df["Dl"].values
dDl = df["dDl"].values
n = len(ell)

# -------------------------------
# 2. 共分散行列（対角近似）
# -------------------------------
cov = np.diag(dDl**2)
cov_inv = np.diag(1.0 / dDl**2)

# -------------------------------
# 3. ΛCDM: 多項式フィット
# -------------------------------
coeffs = np.polyfit(ell, Dl_obs, deg=6)
Dl_smooth = np.polyval(coeffs, ell)
residual_lcdm = Dl_obs - Dl_smooth
chi2_lcdm = residual_lcdm @ cov_inv @ residual_lcdm
k_lcdm = 7
aic_lcdm = chi2_lcdm + 2 * k_lcdm
bic_lcdm = chi2_lcdm + k_lcdm * np.log(n)
p_lcdm = chi2.sf(chi2_lcdm, df=n - k_lcdm)
r_lcdm, _ = pearsonr(Dl_smooth, Dl_obs)

print(f"[ΛCDM]")
print(f"  χ²      = {chi2_lcdm:.2f}")
print(f"  DoF     = {n - k_lcdm}")
print(f"  p-value = {p_lcdm:.6f}")
print(f"  AIC     = {aic_lcdm:.2f}")
print(f"  BIC     = {bic_lcdm:.2f}")
print(f"  r       = {r_lcdm:.6f}\n")


# -------------------------------
# 4. Totality モデル
# -------------------------------
def alpha_multi(ell, epsilons, as_, phis):
    total = 1.0
    for eps, a, phi in zip(epsilons, as_, phis):
        total += eps * np.sin(a * ell + phi)
    return total


def make_model(N):
    def model(params):
        epsilons = params[:N]
        as_ = params[N : 2 * N]
        phis = params[2 * N :]
        return Dl_smooth * alpha_multi(ell, epsilons, as_, phis)

    def chi2_fn(params):
        model_vals = model(params)
        delta = Dl_obs - model_vals
        return delta @ cov_inv @ delta

    # 初期値と範囲
    x0 = [0.05] * N + [0.05] * N + [0.0] * N
    bounds = [(-0.2, 0.2)] * N + [(0.005, 0.1)] * N + [(-np.pi, np.pi)] * N

    result = minimize(chi2_fn, x0=x0, bounds=bounds)
    model_vals = model(result.x)
    chi2_val = chi2_fn(result.x)
    k = 3 * N
    dof = n - k
    aic = chi2_val + 2 * k
    bic = chi2_val + k * np.log(n)
    p = chi2.sf(chi2_val, df=dof)
    r, _ = pearsonr(model_vals, Dl_obs)

    return model_vals, result.x, chi2_val, aic, bic, p, dof, r


# -------------------------------
# 5. プロットと統計出力
# -------------------------------
plt.figure(figsize=(12, 6))
plt.errorbar(
    ell, Dl_obs, yerr=dDl, fmt=".", label="Planck Data", color="black", alpha=0.4
)
plt.plot(ell, Dl_smooth, label="ΛCDM (smooth fit)", color="red", linewidth=2)

colors = ["blue", "green", "purple", "orange", "brown"]
for i, N in enumerate(range(1, 6)):
    model_vals, params, chi2_val, aic_val, bic_val, p_val, dof_val, r_val = make_model(
        N
    )
    eps_str = ", ".join([f"{e:.3f}" for e in params[:N]])
    a_str = ", ".join([f"{a:.3f}" for a in params[N : 2 * N]])
    phi_str = ", ".join([f"{p:.2f}" for p in params[2 * N :]])

    print(f"[Totality N={N}]")
    print(f"  χ²      = {chi2_val:.2f}")
    print(f"  DoF     = {dof_val}")
    print(f"  p-value = {p_val:.6f}")
    print(f"  AIC     = {aic_val:.2f}")
    print(f"  BIC     = {bic_val:.2f}")
    print(f"  r       = {r_val:.6f}")
    print(f"  ε       = [{eps_str}]")
    print(f"  a       = [{a_str}]")
    print(f"  φ       = [{phi_str}]\n")

    label = f"Totality N={N}"
    plt.plot(ell, model_vals, linestyle="--", color=colors[i], label=label)

plt.xlabel(r"Multipole moment $\ell$")
plt.ylabel(r"$D_\ell\ [\mu K^2]$")
plt.title("CMB TT Power Spectrum: ΛCDM vs Totality Models (N=1~5)")
plt.legend(fontsize=9)
plt.grid(True)
plt.tight_layout()
plt.savefig("cmb_tt_totality_allN_diagcov.png")
plt.show()
