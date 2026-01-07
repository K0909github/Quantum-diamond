"""vacancy_list（OVITO出力）から vacancy の深さ分布をヒストグラム表示する。

前提:

深さは `depth = SURFACE_Z - z` (Å) と定義する。
"""

from __future__ import annotations

import math
import shlex
from pathlib import Path


SURFACE_Z = 125.0  # Å  (diamond_substrate_250Å の上面)
BIN_WIDTH = 1  # Å


def read_positions(path: Path) -> list[tuple[float, float, float]]:
    x_idx: int | None = None
    y_idx: int | None = None
    z_idx: int | None = None

    positions: list[tuple[float, float, float]] = []
    for raw_line in path.read_text(encoding="utf-8", errors="ignore").splitlines():
        line = raw_line.strip()
        if not line:
            continue

        if line.startswith("#"):
            if "position.x" in line.lower():
                header = shlex.split(line.lstrip("#").strip())
                header_l = [h.lower() for h in header]
                try:
                    x_idx = header_l.index("position.x")
                    y_idx = header_l.index("position.y")
                    z_idx = header_l.index("position.z")
                except ValueError:
                    x_idx = y_idx = z_idx = None
            continue

        tokens = line.split()

        if x_idx is not None and y_idx is not None and z_idx is not None:
            if len(tokens) <= max(x_idx, y_idx, z_idx):
                continue
            try:
                x = float(tokens[x_idx])
                y = float(tokens[y_idx])
                z = float(tokens[z_idx])
            except ValueError:
                continue
            positions.append((x, y, z))
            continue

        if len(tokens) < 3:
            continue
        try:
            x = float(tokens[0])
            y = float(tokens[1])
            z = float(tokens[2])
        except ValueError:
            continue
        positions.append((x, y, z))

    return positions


def compute_depths(points: list[tuple[float, float, float]]) -> list[float]:
    return [SURFACE_Z - z for (_, _, z) in points]


def mean_pair_distance(points: list[tuple[float, float, float]]) -> float:
    n = len(points)
    if n < 2:
        return math.nan

    total = 0.0
    pairs = 0
    for i in range(n - 1):
        xi, yi, zi = points[i]
        for j in range(i + 1, n):
            xj, yj, zj = points[j]
            dx = xi - xj
            dy = yi - yj
            dz = zi - zj
            total += math.sqrt(dx * dx + dy * dy + dz * dz)
            pairs += 1
    return total / pairs


def main() -> None:
    data_dir = Path(__file__).parent
    path = data_dir / "vacancy_list"
    if not path.exists():
        path = data_dir / "vacancy_list.txt"
    if not path.exists():
        print("vacancy_list（または vacancy_list.txt）が見つかりません")
        return

    pts = read_positions(path)
    if not pts:
        print("vacancy座標が読み取れませんでした")
        return

    depths = [d for d in compute_depths(pts) if d >= 0.0]
    if not depths:
        print("深さ(depth)が計算できませんでした")
        return

    mean_distance = mean_pair_distance(pts)
    print(f"{path.name}: vacancies={len(depths)} mean_distance={mean_distance:.3f} Å")

    try:
        import matplotlib

        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
    except ModuleNotFoundError:
        print("matplotlib が必要です: pip install matplotlib")
        return

    max_depth = max(depths)
    xmax = int(math.ceil(max_depth / BIN_WIDTH) * BIN_WIDTH)
    edges = list(range(0, xmax + BIN_WIDTH, BIN_WIDTH))

    plt.figure(figsize=(7, 4))
    plt.hist(depths, bins=edges, range=(0.0, float(xmax)))
    plt.xlim(0.0, float(xmax))
    plt.xlabel("vacancy depth (Å)")
    plt.ylabel("count")
    plt.title("Vacancy depth histogram")
    plt.grid(True, alpha=0.3)
    plt.margins(x=0)
    plt.tight_layout()

    out_path = data_dir / "vacancy_depths.png"
    plt.savefig(out_path, dpi=200)
    print(f"Saved: {out_path}")


if __name__ == "__main__":
    main()

