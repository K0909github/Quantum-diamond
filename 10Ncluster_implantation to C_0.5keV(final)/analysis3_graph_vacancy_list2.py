"""2つの vacancy_list（OVITO出力）から vacancy の深さ分布を線グラフで重ね描きする。

前提:

深さは `depth = SURFACE_Z - z` (Å) と定義する。
"""

from __future__ import annotations

import argparse
import glob
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

        # 1) 単純な3列: "x y z"
        try:
            x = float(tokens[0])
            y = float(tokens[1])
            z = float(tokens[2])
            positions.append((x, y, z))
            continue
        except ValueError:
            
            pass

        # 2) XYZ/拡張XYZの典型: "C x y z" (以降に追加列があってもOK)
        if len(tokens) >= 4:
            try:
                x = float(tokens[1])
                y = float(tokens[2])
                z = float(tokens[3])
            except ValueError:
                continue
            positions.append((x, y, z))
            continue

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


def _candidate_files_in_dir(data_dir: Path) -> list[Path]:
    return [
        data_dir / "vacancy_list.xyz",
        data_dir / "vacancy_list",
        data_dir / "vacancy_list.txt",
    ]


def _resolve_inputs(raw_inputs: list[str] | None) -> list[Path]:
    if not raw_inputs:
        # 従来互換: スクリプトと同じフォルダから探す
        data_dir = Path(__file__).parent
        for cand in _candidate_files_in_dir(data_dir):
            if cand.exists():
                return [cand]
        return []

    resolved: list[Path] = []
    for raw in raw_inputs:
        # run_*/... のようにディレクトリ側にワイルドカードがある場合も扱えるように、先にglob展開する。
        if any(ch in raw for ch in ["*", "?", "[", "]"]):
            matches = sorted(glob.glob(raw, recursive=True))
            for m in matches:
                mp = Path(m)
                if mp.is_dir():
                    hit = next((c for c in _candidate_files_in_dir(mp) if c.exists()), None)
                    if hit is not None:
                        resolved.append(hit)
                elif mp.exists():
                    resolved.append(mp)
            continue

        p = Path(raw)

        if p.is_dir():
            hit = next((c for c in _candidate_files_in_dir(p) if c.exists()), None)
            if hit is not None:
                resolved.append(hit)
            continue

        if p.exists():
            resolved.append(p)
            continue

        # ワイルドカード等: 親から glob
        parent = p.parent if str(p.parent) else Path(".")
        pattern = p.name
        for m in sorted(parent.glob(pattern)):
            if m.is_dir():
                hit = next((c for c in _candidate_files_in_dir(m) if c.exists()), None)
                if hit is not None:
                    resolved.append(hit)
            elif m.exists():
                resolved.append(m)

    # 重複排除（順序維持）
    uniq: list[Path] = []
    seen: set[Path] = set()
    for p in resolved:
        rp = p.resolve()
        if rp in seen:
            continue
        seen.add(rp)
        uniq.append(p)
    return uniq


def _depths_from_inputs(raw_inputs: list[str]) -> list[float]:
    paths = _resolve_inputs(raw_inputs)
    if not paths:
        return []

    all_depths: list[float] = []
    for path in paths:
        pts = read_positions(path)
        if not pts:
            print(f"skip: {path} (vacancy座標が読み取れません)")
            continue

        depths = [d for d in compute_depths(pts) if d >= 0.0]
        if not depths:
            print(f"skip: {path} (深さ(depth)が計算できません)")
            continue

        mean_distance = mean_pair_distance(pts)
        print(f"{path}: vacancies={len(depths)} mean_distance={mean_distance:.3f} Å")
        all_depths.extend(depths)

    return all_depths


def _hist_counts(depths: list[float], *, bin_width: float, xmax: float) -> tuple[list[float], list[int]]:
    if not depths:
        return [], []
    if bin_width <= 0:
        raise ValueError("bin_width must be > 0")
    if xmax <= 0:
        raise ValueError("xmax must be > 0")

    n_bins = int(math.ceil(xmax / bin_width))
    edges = [i * bin_width for i in range(n_bins + 1)]
    centers = [(edges[i] + edges[i + 1]) / 2.0 for i in range(n_bins)]
    counts = [0 for _ in range(n_bins)]

    for d in depths:
        if d < 0:
            continue
        idx = int(d / bin_width)
        if 0 <= idx < n_bins:
            counts[idx] += 1

    return centers, counts


def main() -> None:
    global SURFACE_Z, BIN_WIDTH

    parser = argparse.ArgumentParser(
        description=(
            "2つの vacancy_list 群から vacancy 深さ分布を線グラフで重ね描きする。"
            "（例: total50keVの10クラスター と 70keV を同時表示）"
        )
    )
    parser.add_argument(
        "--a",
        nargs="+",
        required=True,
        help=(
            "データセットAの入力（ファイル/フォルダ/ワイルドカード、複数可）。"
            "フォルダ指定の場合は vacancy_list.xyz / vacancy_list / vacancy_list.txt を自動探索。"
        ),
    )
    parser.add_argument(
        "--b",
        nargs="+",
        required=True,
        help=(
            "データセットBの入力（ファイル/フォルダ/ワイルドカード、複数可）。"
            "フォルダ指定の場合は vacancy_list.xyz / vacancy_list / vacancy_list.txt を自動探索。"
        ),
    )
    parser.add_argument("--label-a", type=str, default=None, help="凡例に表示するデータセットAの名前")
    parser.add_argument("--label-b", type=str, default=None, help="凡例に表示するデータセットBの名前")
    parser.add_argument("--surface-z", type=float, default=SURFACE_Z, help="表面 z 座標 (Å)")
    parser.add_argument("--bin-width", type=float, default=float(BIN_WIDTH), help="ビン幅 (Å)")
    parser.add_argument(
        "--xmax",
        type=float,
        default=50.0,
        help="深さの表示上限 (Å)。0以下を指定するとデータから自動推定",
    )
    parser.add_argument("--out", type=str, default="vacancy_depths_overlay.png", help="出力画像名")
    args = parser.parse_args()
    SURFACE_Z = float(args.surface_z)
    BIN_WIDTH = float(args.bin_width)

    data_dir = Path(__file__).parent
    depths_a = _depths_from_inputs(args.a)
    depths_b = _depths_from_inputs(args.b)

    if not depths_a:
        print("データセットA: 有効な vacancy データがありませんでした")
    if not depths_b:
        print("データセットB: 有効な vacancy データがありませんでした")
    if not depths_a and not depths_b:
        return

    try:
        import matplotlib

        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
    except ModuleNotFoundError:
        print("matplotlib が必要です: pip install matplotlib")
        return

    max_depth = 0.0
    if depths_a:
        max_depth = max(max_depth, max(depths_a))
    if depths_b:
        max_depth = max(max_depth, max(depths_b))
    xmax_data = float(math.ceil(max_depth / BIN_WIDTH) * BIN_WIDTH) if max_depth > 0 else float(BIN_WIDTH)
    xmax = float(args.xmax) if float(args.xmax) > 0 else xmax_data

    label_a = args.label_a
    if label_a is None:
        label_a = str(Path(args.a[0]).name)
    label_b = args.label_b
    if label_b is None:
        label_b = str(Path(args.b[0]).name)

    plt.figure(figsize=(7, 4))
    if depths_a:
        x_a, y_a = _hist_counts(depths_a, bin_width=BIN_WIDTH, xmax=xmax)
        plt.plot(x_a, y_a, label=f"{label_a} (n={len(depths_a)})")
    if depths_b:
        x_b, y_b = _hist_counts(depths_b, bin_width=BIN_WIDTH, xmax=xmax)
        plt.plot(x_b, y_b, label=f"{label_b} (n={len(depths_b)})")

    plt.xlim(0.0, float(xmax))
    plt.xlabel("vacancy depth (Å)")
    plt.ylabel("count")
    plt.title("Vacancy depth distribution (overlay)")
    plt.grid(True, alpha=0.3)
    plt.legend()
    plt.margins(x=0)
    plt.tight_layout()

    out_arg = Path(args.out)
    out_path = out_arg if out_arg.is_absolute() else (data_dir / out_arg)
    out_path.parent.mkdir(parents=True, exist_ok=True)
    plt.savefig(out_path, dpi=200)
    print(f"Saved: {out_path}")


if __name__ == "__main__":
    main()

