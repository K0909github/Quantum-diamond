"""N_list から N 原子の注入深さ分布（ヒストグラム）を作る。

対応フォーマット:
- LAMMPS dump 形式: 先頭に `ITEM:` があり、`ITEM: ATOMS ...` の後に座標が続く
- OVITO のテキスト出力風: 先頭に `#` のヘッダ行があり、列に `Position.X Position.Y Position.Z` が含まれる

メモ:
- LAMMPS dump を直接渡す場合、`--atom-type 3` のように type フィルタを使うと
    「N 原子だけ」を抽出できます（dump に type 列がある前提）。

深さは `depth = SURFACE_Z - z` (Å) と定義する。
グラフは深さヒストグラムのみを保存し、N-N 平均距離は数値だけ出力する。
"""

from __future__ import annotations

import argparse
import math
from pathlib import Path
import shlex


SURFACE_Z = 125.0  # Å  (diamond_substrate_250Å の上面)
BIN_WIDTH = 5  # Å


def _candidate_n_files_in_dir(data_dir: Path) -> list[Path]:
    return [
        data_dir / "N_list",
        data_dir / "N_list.txt",
        data_dir / "N_list.xyz",
    ]


def _resolve_inputs(raw_inputs: list[str] | None) -> list[Path]:
    if not raw_inputs:
        # 従来互換: スクリプトと同じフォルダの N_list
        data_dir = Path(__file__).parent
        for cand in _candidate_n_files_in_dir(data_dir):
            if cand.exists():
                return [cand]
        return []

    resolved: list[Path] = []
    for raw in raw_inputs:
        p = Path(raw)

        if p.is_dir():
            hit = next((c for c in _candidate_n_files_in_dir(p) if c.exists()), None)
            if hit is not None:
                resolved.append(hit)
            continue

        if p.exists():
            resolved.append(p)
            continue

        parent = p.parent if str(p.parent) else Path(".")
        pattern = p.name
        for m in sorted(parent.glob(pattern)):
            if m.is_dir():
                hit = next((c for c in _candidate_n_files_in_dir(m) if c.exists()), None)
                if hit is not None:
                    resolved.append(hit)
            elif m.exists():
                resolved.append(m)

    uniq: list[Path] = []
    seen: set[Path] = set()
    for p in resolved:
        rp = p.resolve()
        if rp in seen:
            continue
        seen.add(rp)
        uniq.append(p)
    return uniq


def _read_positions_from_lammps_dump(path: Path, atom_type: int | None) -> list[tuple[float, float, float]]:
    """LAMMPS dump 形式から (x, y, z) を読む。

    複数スナップショットが入っている場合は、最後の ITEM: ATOMS ブロックを採用する。
    """

    positions: list[tuple[float, float, float]] = []
    current: list[tuple[float, float, float]] = []
    reading = False
    x_idx = y_idx = z_idx = None
    type_idx: int | None = None

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
                type_idx = col_l.index("type") if "type" in col_l else None
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

        if atom_type is not None and type_idx is not None:
            if len(tokens) <= type_idx:
                continue
            try:
                t = int(float(tokens[type_idx]))
            except ValueError:
                continue
            if t != atom_type:
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


def read_positions(path: Path, atom_type: int | None = None) -> list[tuple[float, float, float]]:
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
        return _read_positions_from_lammps_dump(path, atom_type=atom_type)
    if first_meaningful.startswith("#"):
        return _read_positions_from_ovito_table(path)

    # どちらでもない場合は、従来互換のためdumpとして一応試す
    pts = _read_positions_from_lammps_dump(path, atom_type=atom_type)
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
    global SURFACE_Z, BIN_WIDTH

    parser = argparse.ArgumentParser(description="N_list から N 原子の注入深さ分布をヒストグラム表示する")
    parser.add_argument(
        "inputs",
        nargs="*",
        help=(
            "入力ファイル/フォルダ/ワイルドカード（複数指定可）。"
            "フォルダ指定の場合は N_list / N_list.txt / N_list.xyz を自動探索。"
            "未指定なら、このスクリプトと同じフォルダから探します。"
        ),
    )
    parser.add_argument("--surface-z", type=float, default=SURFACE_Z, help="表面 z 座標 (Å)")
    parser.add_argument("--bin-width", type=float, default=float(BIN_WIDTH), help="ビン幅 (Å)")
    parser.add_argument("--atom-type", type=int, default=None, help="LAMMPS dump の type フィルタ（例: N=3）")
    parser.add_argument("--out", type=str, default="nitrogen_depths.png", help="出力画像名")
    args = parser.parse_args()
    SURFACE_Z = float(args.surface_z)
    BIN_WIDTH = float(args.bin_width)

    data_dir = Path(__file__).parent
    paths = _resolve_inputs(args.inputs)
    if not paths:
        print("N_list 入力が見つかりません（ファイル/フォルダ/パターンを指定してください）")
        return

    all_depths: list[float] = []
    for path in paths:
        pts = read_positions(path, atom_type=args.atom_type)
        if not pts:
            print(f"skip: {path} (N 座標が読み取れません)")
            continue

        depths = [d for d in compute_depths(pts) if d >= 0.0]
        if not depths:
            print(f"skip: {path} (深さ(depth)が計算できません)")
            continue

        mean_distance = mean_pair_distance(pts)
        print(f"{path}: N={len(depths)} mean_distance={mean_distance:.3f} Å")
        all_depths.extend(depths)

    if not all_depths:
        print("有効な N データがありませんでした")
        return

    try:
        import matplotlib

        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
    except ModuleNotFoundError:
        print("matplotlib が必要です: pip install matplotlib")
        return

    max_depth = max(all_depths)
    xmax = int(math.ceil(max_depth / BIN_WIDTH) * BIN_WIDTH)
    edges = list(range(0, xmax + BIN_WIDTH, BIN_WIDTH))

    plt.figure(figsize=(7, 4))
    plt.hist(all_depths, bins=edges, range=(0.0, float(xmax)))
    plt.xlim(0.0, float(xmax))
    plt.xlabel("nitrogen depth (Å)")
    plt.ylabel("count")
    plt.title(f"Nitrogen depth histogram (n={len(all_depths)})")
    plt.grid(True, alpha=0.3)
    plt.margins(x=0)
    plt.tight_layout()

    out_path = data_dir / args.out
    plt.savefig(out_path, dpi=200)
    print(f"Saved: {out_path}")


if __name__ == "__main__":
    main()

