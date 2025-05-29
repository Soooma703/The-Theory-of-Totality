import numpy as np
import matplotlib.pyplot as plt
from scipy.ndimage import gaussian_filter, label, find_objects
from tqdm import tqdm


# === 初期化 ===
def initialize_field(D=5, N=16, L=10.0, seed=42):
    np.random.seed(seed)
    dx = L / N
    coords = [np.linspace(0, L, N, endpoint=False) for _ in range(D)]
    mesh = np.meshgrid(*coords, indexing="ij")
    u = np.zeros([N] * D)
    for _ in range(10):
        center = np.random.uniform(0, L, D)
        w = np.random.uniform(0.1 * L, 0.3 * L)
        A = np.random.uniform(0.5, 1.0)
        gauss = sum((m - c) ** 2 for m, c in zip(mesh, center)) / (2 * w**2)
        u += A * np.exp(-gauss)
    return u, dx


# === ラプラシアン（±1近傍）===
def laplacian_5d(u, dx):
    lap = np.zeros_like(u)
    for axis in range(5):
        lap += np.roll(u, 1, axis=axis)
        lap += np.roll(u, -1, axis=axis)
    lap -= 10 * u
    return lap / dx**2


# === 時間更新 ===
def update_field(u, u_old, dx, dt, lambda_nl=0.57, c=1.0):
    lap = laplacian_5d(u, dx)
    u_new = 2 * u - u_old + dt**2 * (c**2 * lap - lambda_nl * u**3)
    return u_new


# === ピークカウント ===
def count_peaks(u, threshold_ratio=0.3, sigma=0.5):
    u_smooth = gaussian_filter(u, sigma=sigma)
    threshold = threshold_ratio * np.max(u_smooth)
    binary = u_smooth > threshold
    labeled, num = label(binary)
    return num


# === シミュレーション ===
def simulate_D5(steps=500, dt=0.01, lambda_nl=0.57):
    u, dx = initialize_field()
    u_old = u.copy()
    peak_history = []
    snapshots = []

    for _ in tqdm(range(steps), desc="Simulating D=5"):
        u_new = update_field(u, u_old, dx, dt, lambda_nl)
        u_old, u = u, u_new
        peak_count = count_peaks(u)
        peak_history.append(peak_count)
        snapshots.append(u.copy())

    return peak_history, snapshots


# === 構造解析 ===
def analyze_d5_structures(u, threshold_ratio=0.3, sigma=0.5):
    u_smooth = gaussian_filter(u, sigma=sigma)
    threshold = threshold_ratio * np.max(u_smooth)
    binary = u_smooth > threshold
    labeled, num_features = label(binary)
    regions = find_objects(labeled)

    volumes, energies, box_sizes, fractal_dims = [], [], [], []

    for i, slc in enumerate(regions, start=1):
        if slc is None:
            continue
        region = labeled[slc] == i
        volume = np.sum(region)
        peak = np.max(u_smooth[slc][region])
        energy = (peak**2) * volume
        box_size = np.prod([s.stop - s.start for s in slc])

        if volume > 0 and box_size > 0:
            log_box = np.log(box_size)
            if log_box > 1e-8:
                fractal = np.log(volume) / log_box
                fractal_dims.append(fractal)

        volumes.append(volume)
        energies.append(energy)
        box_sizes.append(box_size)

    return (
        np.array(volumes),
        np.array(energies),
        np.array(box_sizes),
        np.array(fractal_dims),
    )


# === 実行部 ===
if __name__ == "__main__":
    print("Simulating D=5...")
    peak_series, snapshots = simulate_D5()

    print("Analyzing final structure...")
    volumes, energies, box_sizes, fractals = analyze_d5_structures(snapshots[-1])

    # === 可視化 ===
    plt.figure(figsize=(14, 4))

    plt.subplot(1, 3, 1)
    plt.hist(volumes, bins=20, color="skyblue", edgecolor="black")
    plt.title("D=5 Structure Volumes")
    plt.xlabel("Volume (voxel count)")
    plt.ylabel("Count")

    plt.subplot(1, 3, 2)
    plt.hist(energies, bins=20, color="gold", edgecolor="black")
    plt.title("D=5 Structure Energies")
    plt.xlabel("Energy (amplitude² × volume)")

    plt.subplot(1, 3, 3)
    valid_fractals = fractals[~np.isnan(fractals)]
    if len(valid_fractals) > 0:
        plt.hist(valid_fractals, bins=20, color="salmon", edgecolor="black")
        plt.title("D=5 Fractal Dimension Estimates")
        plt.xlabel("Estimated D_f = log(V)/log(BoxSize)")
    else:
        plt.title("No valid fractal dimension")

    plt.tight_layout()
    plt.savefig("D5_structure_analysis.png")
    plt.show()

    # === 数値出力 ===
    print(f"\n構造数: {len(volumes)}")
    print(
        f"体積 平均={np.mean(volumes):.2f}, 最大={np.max(volumes)}, 最小={np.min(volumes)}"
    )
    print(f"エネルギー 平均={np.mean(energies):.2f}, 最大={np.max(energies):.2f}")
    if len(valid_fractals) > 0:
        print(
            f"フラクタル次元 平均={np.mean(valid_fractals):.2f}, 範囲=({np.min(valid_fractals):.2f} ～ {np.max(valid_fractals):.2f})"
        )
    else:
        print("フラクタル次元を評価できる構造がありません。")
