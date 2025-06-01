import numpy as np
from numba import njit, prange
from tqdm import tqdm
import matplotlib.pyplot as plt
from scipy.ndimage import gaussian_filter, label, find_objects
from scipy.spatial.distance import pdist, squareform


# === Initialization ===
def initialize_field(N, L, seed=42):
    np.random.seed(seed)
    dx = L / N
    x = np.linspace(0, L, N, endpoint=False)
    y = np.linspace(0, L, N, endpoint=False)
    z = np.linspace(0, L, N, endpoint=False)
    X, Y, Z = np.meshgrid(x, y, z, indexing="ij")
    field = np.zeros((N, N, N))
    for _ in range(10):
        cx, cy, cz = np.random.uniform(0, L, 3)
        w = np.random.uniform(0.1 * L, 0.3 * L)
        A = np.random.uniform(0.5, 1.0)
        field += A * np.exp(
            -((X - cx) ** 2 + (Y - cy) ** 2 + (Z - cz) ** 2) / (2 * w**2)
        )
    return field, dx


# === Time update ===
@njit(parallel=True)
def update(u, u_old, dx, dt, lambda_nl, c):
    Nx, Ny, Nz = u.shape
    u_new = np.empty_like(u)
    for i in prange(1, Nx - 1):
        for j in range(1, Ny - 1):
            for k in range(1, Nz - 1):
                lap = (
                    u[i + 1, j, k]
                    + u[i - 1, j, k]
                    + u[i, j + 1, k]
                    + u[i, j - 1, k]
                    + u[i, j, k + 1]
                    + u[i, j, k - 1]
                    - 6 * u[i, j, k]
                ) / dx**2
                nonlin = lambda_nl * u[i, j, k] ** 3
                u_new[i, j, k] = (
                    2 * u[i, j, k] - u_old[i, j, k] + dt**2 * (c**2 * lap - nonlin)
                )
    return u_new


# === Slice visualization ===
def visualize_slice(u, slice_axis="z", index=None, cmap="viridis"):
    if slice_axis == "z":
        if index is None:
            index = u.shape[2] // 2
        plt.imshow(u[:, :, index], cmap=cmap, origin="lower")
    elif slice_axis == "y":
        if index is None:
            index = u.shape[1] // 2
        plt.imshow(u[:, index, :], cmap=cmap, origin="lower")
    elif slice_axis == "x":
        if index is None:
            index = u.shape[0] // 2
        plt.imshow(u[index, :, :], cmap=cmap, origin="lower")
    else:
        raise ValueError("Invalid axis")
    plt.colorbar()
    plt.title(f"Slice at {slice_axis}={index}")
    plt.tight_layout()
    plt.savefig("visualize_slice")
    plt.show()


# === Peak detection ===
def count_localized_peaks(u, threshold_ratio=0.3, smooth_sigma=1.0):
    u_smooth = gaussian_filter(u, sigma=smooth_sigma)
    threshold = threshold_ratio * np.max(u_smooth)
    binary = u_smooth > threshold
    labeled, num_features = label(binary)
    return num_features


# === Time evolution and full analysis ===
def run_simulation(lambda_nl, N=256, L=10.0, dt=0.01, Nt=2000):
    c = 1.0
    u, dx = initialize_field(N, L)
    u_old = u.copy()

    snapshots = []
    peak_counts = []

    for t in tqdm(range(Nt), desc=f"Simulating λ={lambda_nl:.5f}"):
        u_new = update(u, u_old, dx, dt, lambda_nl, c)
        u_old, u = u, u_new

        if t % 10 == 0:
            snapshots.append(u.copy())
            num_peaks = count_localized_peaks(u)
            peak_counts.append(num_peaks)

    return snapshots, peak_counts


# === Execution ===
if __name__ == "__main__":
    lambda_c = 0.57
    snapshots, peak_counts = run_simulation(lambda_c)

    # Visualize a slice at the final time step
    visualize_slice(snapshots[-1], slice_axis="z")

    # Plot evolution of peak count
    time_points = [i * 10 for i in range(len(peak_counts))]
    plt.figure(figsize=(8, 5))
    plt.plot(time_points, peak_counts, marker="o")
    plt.xlabel("Time step")
    plt.ylabel("Number of localized peaks")
    plt.title(f"Temporal evolution of particle-like structures (λ={lambda_c})")
    plt.grid(True)
    plt.tight_layout()
    plt.savefig("peak_count_evolution.png")
    plt.show()
    print(f"最終時刻におけるピーク数: {peak_counts[-1]}")


# === Analyze structure volume and amplitude ===
def analyze_structures(u, threshold_ratio=0.3, smooth_sigma=1.0):
    u_smooth = gaussian_filter(u, sigma=smooth_sigma)
    threshold = threshold_ratio * np.max(u_smooth)
    binary = u_smooth > threshold
    labeled, num_features = label(binary)
    regions = find_objects(labeled)

    volumes = []
    amplitudes = []

    for i, slc in enumerate(regions, start=1):
        if slc is None:
            continue
        region = labeled[slc] == i
        volume = np.sum(region)
        peak = np.max(u_smooth[slc][region])
        volumes.append(volume)
        amplitudes.append(peak)

    return volumes, amplitudes


# === Statistical analysis ===
volumes, amplitudes = analyze_structures(snapshots[-1], threshold_ratio=0.3)

# === Histogram of volume distribution ===
plt.figure(figsize=(7, 4))
plt.hist(volumes, bins=20, color="skyblue", edgecolor="black")
plt.xlabel("Structure volume (voxel count)")
plt.ylabel("Count")
plt.title("Distribution of structure volumes")
plt.tight_layout()
plt.savefig("volume_distribution.png")
plt.show()

# === Histogram of amplitude distribution ===
plt.figure(figsize=(7, 4))
plt.hist(amplitudes, bins=20, color="salmon", edgecolor="black")
plt.xlabel("Peak amplitude")
plt.ylabel("Count")
plt.title("Distribution of structure peak amplitudes")
plt.tight_layout()
plt.savefig("amplitude_distribution.png")
plt.show()

# === Numerical statistics output ===
print(f"Number of structures: {len(volumes)}")
print(
    f"Volume: mean = {np.mean(volumes):.2f}, max = {np.max(volumes)}, min = {np.min(volumes)}"
)
print(
    f"Amplitude: mean = {np.mean(amplitudes):.3f}, max = {np.max(amplitudes):.3f}, min = {np.min(amplitudes):.3f}"
)


# === Compute centroid, energy, and position of structures ===
def analyze_structures_detailed(u, threshold_ratio=0.3, smooth_sigma=1.0):
    u_smooth = gaussian_filter(u, sigma=smooth_sigma)
    threshold = threshold_ratio * np.max(u_smooth)
    binary = u_smooth > threshold
    labeled, num_features = label(binary)
    regions = find_objects(labeled)

    volumes, amplitudes, energies, centers = [], [], [], []

    for i, slc in enumerate(regions, start=1):
        if slc is None:
            continue
        region = labeled[slc] == i
        volume = np.sum(region)
        peak = np.max(u_smooth[slc][region])
        energy = (peak**2) * volume

        # centroid calculation
        coords = np.argwhere(region)
        abs_coords = coords + [s.start for s in slc]
        center = np.mean(abs_coords, axis=0)

        volumes.append(volume)
        amplitudes.append(peak)
        energies.append(energy)
        centers.append(center)

    return volumes, amplitudes, energies, np.array(centers)


# === Analysis and visualization ===
volumes, amplitudes, energies, centers = analyze_structures_detailed(
    snapshots[-1], threshold_ratio=0.3
)

# --- Histogram of energy distribution
plt.figure(figsize=(7, 4))
plt.hist(energies, bins=20, color="gold", edgecolor="black")
plt.xlabel("Structure energy (amplitude² × volume)")
plt.ylabel("Count")
plt.title("Distribution of structure energies")
plt.tight_layout()
plt.savefig("energy_distribution.png")
plt.show()

# --- Distribution of inter-structure distances (nearest neighbor)
if len(centers) >= 2:
    dist_matrix = squareform(pdist(centers))
    np.fill_diagonal(dist_matrix, np.inf)
    nearest_distances = np.min(dist_matrix, axis=1)

    plt.figure(figsize=(7, 4))
    plt.hist(nearest_distances, bins=20, color="gray", edgecolor="black")
    plt.xlabel("Nearest neighbor distance")
    plt.ylabel("Count")
    plt.title("Distribution of structure separations")
    plt.tight_layout()
    plt.savefig("structure_distance_distribution.png")
    plt.show()
else:
    print("Not enough structures to compute inter-structure distances.")


# --- Track the centroid of the largest structure during time evolution
def track_largest_structure_centroid(snapshots, threshold_ratio=0.3, smooth_sigma=1.0):
    trajectory = []

    for u in snapshots:
        u_smooth = gaussian_filter(u, sigma=smooth_sigma)
        threshold = threshold_ratio * np.max(u_smooth)
        binary = u_smooth > threshold
        labeled, num_features = label(binary)
        regions = find_objects(labeled)

        max_vol = 0
        max_center = None

        for i, slc in enumerate(regions, start=1):
            if slc is None:
                continue
            region = labeled[slc] == i
            volume = np.sum(region)
            if volume > max_vol:
                coords = np.argwhere(region)
                abs_coords = coords + [s.start for s in slc]
                max_center = np.mean(abs_coords, axis=0)
                max_vol = volume

        if max_center is not None:
            trajectory.append(max_center)

    return np.array(trajectory)


trajectory = track_largest_structure_centroid(snapshots, threshold_ratio=0.3)

if len(trajectory) > 0:
    plt.figure(figsize=(7, 5))
    plt.plot(trajectory[:, 0], label="x")
    plt.plot(trajectory[:, 1], label="y")
    plt.plot(trajectory[:, 2], label="z")
    plt.xlabel("Snapshot index")
    plt.ylabel("Centroid coordinate")
    plt.title("Trajectory of the largest structure centroid")
    plt.legend()
    plt.tight_layout()
    plt.savefig("largest_structure_trajectory.png")
    plt.show()
else:
    print("No valid trajectory for largest structure.")
