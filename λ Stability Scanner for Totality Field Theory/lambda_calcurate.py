import numpy as np
from numba import njit, prange
from tqdm import tqdm
import matplotlib.pyplot as plt
from scipy.signal import argrelextrema


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


# === Time Evolution ===
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


# === Energy (gradient term only) ===
def compute_energy(u, dx):
    grad = np.gradient(u, dx)
    grad_squared = sum(g**2 for g in grad)
    return np.mean(grad_squared)


# === Evaluate δE across a λ-scan ===
def find_lambda_for_stable_structure(
    N=64, L=10.0, dt=0.01, Nt=500, lambda_range=(0, 3.5), num_lambda=100
):
    c = 1.0
    lambda_values = np.linspace(lambda_range[0], lambda_range[1], num_lambda)
    delta_E_list = []

    for lambda_nl in tqdm(lambda_values, desc="Scanning λ"):
        u, dx = initialize_field(N, L)
        u_old = u.copy()
        energies = []

        for t in range(Nt):
            u_new = update(u, u_old, dx, dt, lambda_nl, c)
            u_old, u = u, u_new

            if t % 10 == 0:
                E = compute_energy(u, dx)
                energies.append(E)

        delta_E = np.max(energies) - np.min(energies)
        delta_E_list.append((lambda_nl, delta_E))

    return lambda_values, delta_E_list


# === Plot and analyze ===
if __name__ == "__main__":
    lambdas, delta_Es = find_lambda_for_stable_structure()
    lambda_vals = [x[0] for x in delta_Es]
    delta_vals = [x[1] for x in delta_Es]

    plt.figure(figsize=(8, 5))
    plt.plot(lambda_vals, delta_vals, marker="o")
    plt.xlabel("λ (nonlinear coefficient)")
    plt.ylabel("δE (energy fluctuation)")
    plt.title("Stability Threshold in 3D: δE vs λ")
    plt.grid(True)
    plt.tight_layout()
    plt.savefig("deltaE_vs_lambda.png")
    plt.show()

    # === Compute the numerical derivative (fluctuation rate) ===
    delta_vals = np.array(delta_vals)
    lambda_vals = np.array(lambda_vals)
    derivative = np.gradient(delta_vals, lambda_vals)

    # === Plot fluctuation rate ===
    plt.figure(figsize=(8, 5))
    plt.plot(lambda_vals, derivative, marker="x", color="darkred")
    plt.xlabel("λ (nonlinear coefficient)")
    plt.ylabel("d(δE)/dλ (energy fluctuation rate)")
    plt.title("Energy Fluctuation Rate vs λ")
    plt.grid(True)
    plt.tight_layout()
    plt.savefig("d_deltaE_d_lambda.png")
    plt.show()

    # Identify λ with the minimum fluctuation rate
    min_idx = np.argmin(derivative)
    lambda_c = lambda_vals[min_idx]
    print(f"Estimated critical λ_c for structural stability ≈ {lambda_c:.5f}")

    # === Find peaks in the derivative (local maxima of d(δE)/dλ) ===
    peak_indices = argrelextrema(derivative, np.greater)[0]

    print("=== Local maxima in the first derivative (peaks in d(δE)/dλ) ===")
    for idx in peak_indices:
        print(
            f"λ ≈ {lambda_vals[idx]:.5f}, d(δE)/dλ ≈ {derivative[idx]:.5f}, δE ≈ {delta_vals[idx]:.5f}"
        )

    # Plot the derivative and mark the peaks
    plt.figure(figsize=(8, 5))
    plt.plot(lambda_vals, derivative, label="d(δE)/dλ", color="darkred")
    plt.scatter(
        lambda_vals[peak_indices],
        derivative[peak_indices],
        color="orange",
        label="Local max (peak)",
        zorder=5,
    )
    plt.xlabel("λ (nonlinear coefficient)")
    plt.ylabel("d(δE)/dλ")
    plt.title("Detection of Peaks in the Fluctuation Rate")
    plt.legend()
    plt.grid(True)
    plt.tight_layout()
    plt.savefig("d_deltaE_d_lambda_peaks.png")
    plt.show()
