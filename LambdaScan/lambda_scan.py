import numpy as np
import matplotlib.pyplot as plt
from multiprocessing import Pool, cpu_count
from tqdm import tqdm

# 定数定義（グローバル変数として共有される）
L = 100
Nx = 500
dx = L / Nx
x = np.linspace(0, L, Nx)
c = 1.0
dt = 0.01
Nt = 50000


def simulate_energy(lambda_nl):
    u = np.exp(-((x - L / 2) ** 2) / 4)
    u_new = u.copy()
    u_old = u.copy()
    energies = []

    for _ in range(Nt):
        nonlinear_term = lambda_nl * u**3
        u_new[1:-1] = (
            2 * u[1:-1]
            - u_old[1:-1]
            + (c * dt / dx) ** 2 * (u[2:] - 2 * u[1:-1] + u[:-2])
            - dt**2 * nonlinear_term[1:-1]
        )

        du_dphi = (u - u_old) / dt
        kinetic = 0.5 * du_dphi**2
        gradient = 0.5 * c**2 * ((np.roll(u, -1) - np.roll(u, 1)) / (2 * dx)) ** 2
        potential = 0.25 * lambda_nl * u**4
        energy_density = kinetic + gradient + potential
        total_energy = np.sum(energy_density) * dx
        energies.append(total_energy)

        u_old = u.copy()
        u = u_new.copy()

    variation = np.max(energies) - np.min(energies)
    return lambda_nl, variation


if __name__ == "__main__":
    # 探索範囲
    lambda_candidates = np.round(np.arange(0.01, 1.001, 0.0005), 4)

    # 並列処理プール
    with Pool(processes=cpu_count()) as pool:
        results = list(
            tqdm(
                pool.imap(simulate_energy, lambda_candidates),
                total=len(lambda_candidates),
            )
        )

    # 結果整理
    lambda_vals, variations = zip(*results)
    best_lambda = lambda_vals[np.argmin(variations)]

    print(f"\nBest λ (min energy variation): {best_lambda:.4f}")

    # オプション：可視化
    plt.plot(lambda_vals, variations)
    plt.xlabel("λ")
    plt.ylabel("Energy Variation")
    plt.title("Energy Stability vs λ")
    plt.grid(True)
    plt.tight_layout()
    plt.show()
