import numpy as np
import matplotlib.pyplot as plt
from scipy.ndimage import gaussian_filter, label, find_objects
from tqdm import tqdm


# === 初期化 ===
def initialize_field(N=128, L=10.0, seed=42):
    np.random.seed(seed)
    dx = L / N
    x = np.linspace(0, L, N, endpoint=False)
    y = np.linspace(0, L, N, endpoint=False)
    X, Y = np.meshgrid(x, y, indexing="ij")
    u = np.zeros_like(X)
    for _ in range(10):
        cx, cy = np.random.uniform(0, L, 2)
        w = np.random.uniform(0.1 * L, 0.3 * L)
        A = np.random.uniform(0.5, 1.0)
        u += A * np.exp(-((X - cx) ** 2 + (Y - cy) ** 2) / (2 * w**2))
    return u, dx


# === 時間更新（2Dラプラシアン）===
def update_field(u, u_old, dx, dt, lambda_nl=0.57, c=1.0):
    u_new = np.zeros_like(u)
    u_new[1:-1, 1:-1] = (
        2 * u[1:-1, 1:-1]
        - u_old[1:-1, 1:-1]
        + dt**2
        * (
            c**2
            * (
                u[2:, 1:-1]
                + u[:-2, 1:-1]
                + u[1:-1, 2:]
                + u[1:-1, :-2]
                - 4 * u[1:-1, 1:-1]
            )
            / dx**2
            - lambda_nl * u[1:-1, 1:-1] ** 3
        )
    )
    return u_new


# === ピークカウント ===
def count_peaks(u, threshold_ratio=0.3, sigma=1.0):
    u_smooth = gaussian_filter(u, sigma=sigma)
    threshold = threshold_ratio * np.max(u_smooth)
    binary = u_smooth > threshold
    labeled, num = label(binary)
    return num


# === シミュレーション本体 ===
def simulate_D2(steps=500, dt=0.01, lambda_nl=0.57):
    u, dx = initialize_field()
    u_old = u.copy()
    peak_history = []
    snapshots = []

    for _ in tqdm(range(steps), desc="Simulating D=2"):
        u_new = update_field(u, u_old, dx, dt, lambda_nl)
        u_old, u = u, u_new
        peak_count = count_peaks(u)
        peak_history.append(peak_count)
        snapshots.append(u.copy())

    return peak_history, snapshots


# === 構造解析 ===
def analyze_d2_structures(u, threshold_ratio=0.3, sigma=1.0):
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
    print("Simulating D=2...")
    peak_series, snapshots = simulate_D2()

    print("Analyzing final structure...")
    volumes, energies, box_sizes, fractals = analyze_d2_structures(snapshots[-1])

    # === 可視化 ===
    plt.figure(figsize=(14, 4))

    plt.subplot(1, 3, 1)
    plt.hist(volumes, bins=20, color="skyblue", edgecolor="black")
    plt.title("D=2 Structure Volumes")
    plt.xlabel("Volume (pixel count)")
    plt.ylabel("Count")

    plt.subplot(1, 3, 2)
    plt.hist(energies, bins=20, color="gold", edgecolor="black")
    plt.title("D=2 Structure Energies")
    plt.xlabel("Energy (amplitude² × volume)")

    plt.subplot(1, 3, 3)
    valid_fractals = fractals[~np.isnan(fractals)]
    if len(valid_fractals) > 0:
        plt.hist(valid_fractals, bins=20, color="salmon", edgecolor="black")
        plt.title("D=2 Fractal Dimension Estimates")
        plt.xlabel("Estimated D_f = log(V)/log(BoxSize)")
    else:
        plt.title("No valid fractal dimension")

    plt.tight_layout()
    plt.savefig("D2_structure_analysis.png")
    plt.show()

    # === 統計出力 ===
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
