from astropy.io import fits
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd

# === 1. Load evt2.fits ===
evt_path = "acisf20207N002_evt2.fits"
hdul = fits.open(evt_path)
data = hdul["EVENTS"].data

# === 2. ENERGYカラム or PIカラムからヒストグラム（PHA）を作成 ===
# 通常は PI (Pulse Invariant) を使用
energies = data["PI"]  # または data["ENERGY"]（単位 keV の場合もあり）

# === 3. バイニング設定とヒストグラム作成 ===
nbins = 1024
min_pi = np.min(energies)
max_pi = np.max(energies)

hist, bin_edges = np.histogram(energies, bins=nbins, range=(min_pi, max_pi))
bin_centers = 0.5 * (bin_edges[:-1] + bin_edges[1:])

# === 4. CSVとして保存 ===
df = pd.DataFrame({"CHANNEL": bin_centers, "COUNTS": hist})
df.to_csv("maxi1820_spectrum.csv", index=False)
print("✅ Saved: maxi1820_spectrum.csv")

# === 5. Optional: Plot the spectrum ===
plt.figure(figsize=(10, 5))
plt.plot(bin_centers, hist, drawstyle="steps-mid")
plt.xlabel("PI Channel")
plt.ylabel("Counts")
plt.title("Extracted Spectrum from evt2.fits")
plt.grid(True)
plt.tight_layout()
plt.savefig("maxi1820_spectrum_plot.png")
plt.show()
