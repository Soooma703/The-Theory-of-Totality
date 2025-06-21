import numpy as np
import pandas as pd
import importlib.util
import os
from tqdm import tqdm

# === 定数設定 ===
dim_labels = ["1D", "2D", "3D", "4D", "5D"]
file_names = ["1D.py", "2D.py", "3D.py", "4D.py", "5D.py"]
module_names = ["mod1d", "mod2d", "mod3d", "mod4d", "mod5d"]
lambda_values = np.linspace(0.01, 2, 50)
steps = 500  # 個別実行と一致させる
global_seed = 42


# === 実行関数（逐次処理） ===
def run_single_simulation(dim, mod_name, file, lam):
    try:
        # 各条件に対してシード固定（dim と λ に基づいて変化）
        np.random.seed(global_seed + int(lam * 1000) + dim * 10000)

        # モジュール読み込み
        spec = importlib.util.spec_from_file_location(mod_name, file)
        mod = importlib.util.module_from_spec(spec)
        spec.loader.exec_module(mod)

        simulate_func = getattr(mod, f"simulate_D{dim}")
        analyze_func = getattr(mod, f"analyze_d{dim}_structures")

        _, snapshots = simulate_func(lambda_nl=lam, steps=steps)
        u = snapshots[-1]
        volumes, energies, _, fractals = analyze_func(u)

        return {
            "dimension": dim,
            "lambda": lam,
            "structure_count": len(volumes),
            "mean_volume": np.mean(volumes) if len(volumes) > 0 else 0,
            "mean_energy": np.mean(energies) if len(energies) > 0 else 0,
            "mean_fractal_dim": np.mean(fractals) if len(fractals) > 0 else 0,
        }

    except Exception as e:
        return {
            "dimension": dim,
            "lambda": lam,
            "structure_count": 0,
            "mean_volume": 0,
            "mean_energy": 0,
            "mean_fractal_dim": 0,
            "error": str(e),
        }


# === メイン処理 ===
if __name__ == "__main__":
    # フルパスで .py ファイル指定
    base_dir = os.path.dirname(os.path.abspath(__file__))
    full_paths = [os.path.join(base_dir, fname) for fname in file_names]

    results = []
    for dim, (mod_name, file) in enumerate(zip(module_names, full_paths), start=1):
        for lam in tqdm(lambda_values, desc=f"Simulating D={dim}"):
            result = run_single_simulation(dim, mod_name, file, lam)
            results.append(result)

    df = pd.DataFrame(results)
    df.to_csv("tft_dimension_scan_results.csv", index=False)
    print("完了：'tft_dimension_scan_results.csv' に保存されました。")
