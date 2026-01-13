"""N_list から N 原子の注入深さ分布（ヒストグラム）を作る。

対応フォーマット:
- LAMMPS dump 形式: 先頭に `ITEM:` があり、`ITEM: ATOMS ...` の後に座標が続く
- OVITO のテキスト出力風: 先頭に `#` のヘッダ行があり、列に `Position.X Position.Y Position.Z` が含まれる

深さは `depth = SURFACE_Z - z` (Å) と定義する。
グラフは深さヒストグラムのみを保存し、N-N 平均距離は数値だけ出力する。
"""

from __future__ import annotations

import math
from pathlib import Path
import shlex


SURFACE_Z = 125.0  # Å  (diamond_substrate_250Å の上面)
BIN_WIDTH = 5  # Å


def _read_positions_from_lammps_dump(path: Path) -> list[tuple[float, float, float]]:
    """LAMMPS dump 形式から (x, y, z) を読む。

    複数スナップショットが入っている場合は、最後の ITEM: ATOMS ブロックを採用する。
    """

    positions: list[tuple[float, float, float]] = []
    current: list[tuple[float, float, float]] = []
    reading = False
    x_idx = y_idx = z_idx = None

    for raw_line in path.read_text(encoding="utf-8", errors="ignore").splitlines():
        line = raw_line.strip()
        if not line:
            continue

        if line.startswith("ITEM:"):
            if line.startswith("ITEM: ATOMS"):
                cols = line.split()[2:]
                col_l = [c.lower() for c in cols]
                try:
                    x_idx = col_l.index("x")
                    y_idx = col_l.index("y")
                    z_idx = col_l.index("z")
                except ValueError:
                    x_idx = y_idx = z_idx = None
                current = []
                reading = True
            else:
                if reading:
                    positions = current
                reading = False
            continue

        if not reading or x_idx is None or y_idx is None or z_idx is None:
            continue

        tokens = line.split()
        if len(tokens) <= max(x_idx, y_idx, z_idx):
            continue
        try:
            x = float(tokens[x_idx])
            y = float(tokens[y_idx])
            z = float(tokens[z_idx])
        except ValueError:
            continue
        current.append((x, y, z))

    if reading and current:
        positions = current

    return positions


def _read_positions_from_ovito_table(path: Path) -> list[tuple[float, float, float]]:
    """`#` ヘッダ + 表形式（OVITOのエクスポート風）から (x, y, z) を読む。"""

    lines = path.read_text(encoding="utf-8", errors="ignore").splitlines()

    header_cols: list[str] | None = None
    for raw_line in lines:
        line = raw_line.strip()
        if not line:
            continue
        if not line.startswith("#"):
            # ヘッダが来る前にデータが始まっているケースもあり得るが、
            # この形式ではヘッダが必須とみなす。
            break

        # `# Position.X ...` のような列名行を探す
        candidate = line.lstrip("#").strip()
        if not candidate:
            continue
        try:
            cols = shlex.split(candidate)
        except ValueError:
            cols = candidate.split()
        lower_cols = [c.strip().strip('"').lower() for c in cols]
        if "position.x" in lower_cols and "position.y" in lower_cols and "position.z" in lower_cols:
            header_cols = lower_cols
            break

    if header_cols is None:
        return []

    x_idx = header_cols.index("position.x")
    y_idx = header_cols.index("position.y")
    z_idx = header_cols.index("position.z")

    positions: list[tuple[float, float, float]] = []
    for raw_line in lines:
        line = raw_line.strip()
        if not line or line.startswith("#"):
            continue
        tokens = line.split()
        if len(tokens) <= max(x_idx, y_idx, z_idx):
            continue
        try:
            x = float(tokens[x_idx])
            y = float(tokens[y_idx])
            z = float(tokens[z_idx])
        except ValueError:
            continue
        positions.append((x, y, z))

    return positions


def read_positions(path: Path) -> list[tuple[float, float, float]]:
    """`N_list`の入力形式を自動判別して (x, y, z) を返す。"""

    # 先頭の有効行だけ見て形式を判定
    first_meaningful: str | None = None
    for raw_line in path.read_text(encoding="utf-8", errors="ignore").splitlines():
        line = raw_line.strip()
        if not line:
            continue
        first_meaningful = line
        break

    if first_meaningful is None:
        return []

    if first_meaningful.startswith("ITEM:"):
        return _read_positions_from_lammps_dump(path)
    if first_meaningful.startswith("#"):
        return _read_positions_from_ovito_table(path)

    # どちらでもない場合は、従来互換のためdumpとして一応試す
    pts = _read_positions_from_lammps_dump(path)
    if pts:
        return pts
    return _read_positions_from_ovito_table(path)


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
    path = data_dir / "N_list"
    if not path.exists():
        print("N_list が見つかりません")
        return

    pts = read_positions(path)
    if not pts:
        print("N 座標が読み取れませんでした")
        return

    depths = [d for d in compute_depths(pts) if d >= 0.0]
    if not depths:
        print("深さ(depth)が計算できませんでした")
        return

    mean_distance = mean_pair_distance(pts)
    print(f"{path.name}: N={len(depths)} mean_distance={mean_distance:.3f} Å")

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
    plt.xlabel("nitrogen depth (Å)")
    plt.ylabel("count")
    plt.title("Nitrogen depth histogram")
    plt.grid(True, alpha=0.3)
    plt.margins(x=0)
    plt.tight_layout()

    out_path = data_dir / "nitrogen_depths.png"
    plt.savefig(out_path, dpi=200)
    print(f"Saved: {out_path}")


if __name__ == "__main__":
    main()

