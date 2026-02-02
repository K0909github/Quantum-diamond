"""N_list から N 原子の注入深さ分布（ヒストグラム）を作る。

対応フォーマット:
- LAMMPS dump 形式: 先頭に `ITEM:` があり、`ITEM: ATOMS ...` の後に座標が続く
    - 座標列は `x y z` だけでなく `xu yu zu` や `xs ys zs` でもOK
    - `xs ys zs`（scaled）の場合は `ITEM: BOX BOUNDS` を用いて実座標へ変換
- OVITO のテキスト出力風: 先頭に `#` のヘッダ行があり、列に `Position.X Position.Y Position.Z` が含まれる

メモ:
- LAMMPS dump を直接渡す場合、`--atom-type 3` のように type フィルタを使うと
    「N 原子だけ」を抽出できます（dump に type 列がある前提）。

深さは `depth = SURFACE_Z - z` (Å) と定義する。
グラフは深さヒストグラムのみを保存し、N-N 平均距離は数値だけ出力する。
"""

from __future__ import annotations

import argparse
import glob
import math
from pathlib import Path
import shlex


SURFACE_Z = 125.0  # Å  (diamond_substrate_250Å の上面)
BIN_WIDTH = 5  # Å
DEFAULT_MAX_DEPTH = 250.0  # Å  (基板厚)


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

    # 実行場所(cwd)が違っても動くように、相対パスは「スクリプトのあるフォルダ基準」でも探索する。
    base_dir = Path(__file__).parent

    resolved: list[Path] = []
    for raw in raw_inputs:
        # run_*/... のようにディレクトリ側にワイルドカードがある場合は pathlib だけだと扱いにくいので、
        # 文字列globで先に展開する。
        if any(ch in raw for ch in ["*", "?", "[", "]"]):
            matches = sorted(glob.glob(raw, recursive=True))
            if not matches:
                p_raw = Path(raw)
                if not p_raw.is_absolute():
                    matches = sorted(glob.glob(str(base_dir / raw), recursive=True))
            for m in matches:
                mp = Path(m)
                if mp.is_dir():
                    hit = next((c for c in _candidate_n_files_in_dir(mp) if c.exists()), None)
                    if hit is not None:
                        resolved.append(hit)
                elif mp.exists():
                    resolved.append(mp)
            continue

        p = Path(raw)
        if not p.is_absolute() and not p.exists() and not p.is_dir():
            # cwd基準で見つからない相対パスは、スクリプト基準でも探す
            p2 = base_dir / p
            if p2.exists() or p2.is_dir():
                p = p2

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

        # 親がcwd基準だとズレることがあるので、スクリプト基準でも同じパターンを試す
        if not resolved and not p.is_absolute():
            parent2 = (base_dir / parent).resolve()
            try:
                for m in sorted(parent2.glob(pattern)):
                    if m.is_dir():
                        hit = next((c for c in _candidate_n_files_in_dir(m) if c.exists()), None)
                        if hit is not None:
                            resolved.append(hit)
                    elif m.exists():
                        resolved.append(m)
            except FileNotFoundError:
                pass

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
    `x y z` / `xu yu zu` / `xs ys zs` をサポートする。
    """

    positions: list[tuple[float, float, float]] = []
    current: list[tuple[float, float, float]] = []
    reading_atoms = False
    x_idx = y_idx = z_idx = None
    type_idx: int | None = None
    coord_mode: str | None = None  # 'abs' | 'scaled'
    box_bounds: tuple[float, float, float, float, float, float] | None = None

    def _pick_idx(cols_lower: list[str], candidates: list[str]) -> int | None:
        for name in candidates:
            if name in cols_lower:
                return cols_lower.index(name)
        return None

    def _maybe_scaled_to_abs(value: float, lo: float, hi: float) -> float:
        return lo + value * (hi - lo)

    lines = path.read_text(encoding="utf-8", errors="ignore").splitlines()
    i = 0
    while i < len(lines):
        raw_line = lines[i]
        i += 1

        line = raw_line.strip()
        if not line:
            continue

        if line.startswith("ITEM:"):
            if line.startswith("ITEM: BOX BOUNDS"):
                # 次の3行が bounds
                if i + 2 >= len(lines):
                    continue

                def _parse_bounds_row(s: str) -> tuple[float, float] | None:
                    parts = s.split()
                    if len(parts) < 2:
                        return None
                    try:
                        lo = float(parts[0])
                        hi = float(parts[1])
                    except ValueError:
                        return None
                    return lo, hi

                bx = _parse_bounds_row(lines[i].strip())
                by = _parse_bounds_row(lines[i + 1].strip())
                bz = _parse_bounds_row(lines[i + 2].strip())
                i += 3
                if bx is None or by is None or bz is None:
                    box_bounds = None
                else:
                    box_bounds = (bx[0], bx[1], by[0], by[1], bz[0], bz[1])
                continue

            if line.startswith("ITEM: ATOMS"):
                cols = line.split()[2:]
                col_l = [c.lower() for c in cols]

                x_idx = _pick_idx(col_l, ["x", "xu", "xsu", "xus", "xs"])
                y_idx = _pick_idx(col_l, ["y", "yu", "ysu", "yus", "ys"])
                z_idx = _pick_idx(col_l, ["z", "zu", "zsu", "zus", "zs"])
                type_idx = col_l.index("type") if "type" in col_l else None

                if x_idx is not None and y_idx is not None and z_idx is not None:
                    x_name = col_l[x_idx]
                    y_name = col_l[y_idx]
                    z_name = col_l[z_idx]
                    if x_name.endswith("s") and y_name.endswith("s") and z_name.endswith("s"):
                        coord_mode = "scaled"
                    else:
                        coord_mode = "abs"
                else:
                    coord_mode = None

                current = []
                reading_atoms = True
                continue

            # 他のITEMに入ったら、直前のATOMSブロックを採用
            if reading_atoms:
                positions = current
            reading_atoms = False
            continue

        if not reading_atoms or x_idx is None or y_idx is None or z_idx is None or coord_mode is None:
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
            x_raw = float(tokens[x_idx])
            y_raw = float(tokens[y_idx])
            z_raw = float(tokens[z_idx])
        except ValueError:
            continue

        if coord_mode == "scaled":
            if box_bounds is None:
                continue
            xlo, xhi, ylo, yhi, zlo, zhi = box_bounds
            x = _maybe_scaled_to_abs(x_raw, xlo, xhi)
            y = _maybe_scaled_to_abs(y_raw, ylo, yhi)
            z = _maybe_scaled_to_abs(z_raw, zlo, zhi)
        else:
            x, y, z = x_raw, y_raw, z_raw

        current.append((x, y, z))

    if reading_atoms and current:
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


def mean(values: list[float]) -> float:
    if not values:
        return math.nan
    return sum(values) / len(values)


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
    parser.add_argument(
        "--out",
        type=str,
        default="nitrogen_depths.png",
        help=(
            "出力画像（パス可）。ファイル名のみならスクリプトと同じフォルダに保存。"
            "ディレクトリを含む相対パスならカレントディレクトリ基準で保存。"
        ),
    )
    parser.add_argument(
        "--max-depth",
        type=float,
        default=DEFAULT_MAX_DEPTH,
        help=(
            "深さの最大値 (Å)。この上限は描画だけでなく平均計算/集計にも適用します。"
            "（既定: 250）"
        ),
    )
    parser.add_argument(
        "--no-max-depth",
        action="store_true",
        help="深さ上限を無効化（描画/平均ともに全データを使用）",
    )
    args = parser.parse_args()
    SURFACE_Z = float(args.surface_z)
    BIN_WIDTH = float(args.bin_width)

    max_depth_limit: float | None
    if args.no_max_depth:
        max_depth_limit = None
    else:
        max_depth_limit = float(args.max_depth)
        if max_depth_limit <= 0:
            raise ValueError("--max-depth must be > 0")

    data_dir = Path(__file__).parent
    paths = _resolve_inputs(args.inputs)
    if not paths:
        print("N_list 入力が見つかりません（ファイル/フォルダ/パターンを指定してください）")
        print(f"cwd: {Path.cwd()}")
        print(f"script_dir: {data_dir}")
        print(f"inputs: {args.inputs!r}")
        if args.inputs:
            for raw in args.inputs:
                if any(ch in raw for ch in ["*", "?", "[", "]"]):
                    matches = sorted(glob.glob(raw, recursive=True))
                    if not matches:
                        p_raw = Path(raw)
                        if not p_raw.is_absolute():
                            matches = sorted(glob.glob(str(data_dir / raw), recursive=True))
                    print(f"- glob: {raw!r} -> {len(matches)} matches")
                else:
                    p = Path(raw)
                    p2 = (data_dir / p) if (not p.is_absolute()) else p
                    print(
                        f"- path: {raw!r} -> exists(cwd)={p.exists()} exists(script_dir)={p2.exists()} "
                        f"resolved(cwd)={p.resolve()} resolved(script_dir)={p2.resolve()}"
                    )
        return

    all_depths: list[float] = []
    for path in paths:
        pts = read_positions(path, atom_type=args.atom_type)
        if not pts:
            print(f"skip: {path} (N 座標が読み取れません)")
            continue

        depths = [d for d in compute_depths(pts) if d >= 0.0]
        if max_depth_limit is not None:
            depths = [d for d in depths if d <= max_depth_limit]
        if not depths:
            print(f"skip: {path} (深さ(depth)が計算できません)")
            continue

        mean_distance = mean_pair_distance(pts)
        mean_depth = mean(depths)
        print(
            f"{path}: N={len(depths)} mean_depth={mean_depth:.3f} Å mean_distance={mean_distance:.3f} Å"
        )
        all_depths.extend(depths)

    if not all_depths:
        print("有効な N データがありませんでした")
        return

    if max_depth_limit is None:
        print(f"ALL: N={len(all_depths)} mean_depth={mean(all_depths):.3f} Å")
    else:
        print(
            f"ALL (<= {max_depth_limit:g} Å): N={len(all_depths)} mean_depth={mean(all_depths):.3f} Å"
        )

    depths_for_plot = all_depths

    try:
        import matplotlib

        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
    except ModuleNotFoundError:
        print("matplotlib が必要です: pip install matplotlib")
        return

    max_depth = max(depths_for_plot)
    bin_width = float(BIN_WIDTH)
    if bin_width <= 0:
        raise ValueError("--bin-width must be > 0")

    xmax = math.ceil(max_depth / bin_width) * bin_width
    # 浮動小数の丸め誤差対策（例: 25.0000000004）
    xmax = round(float(xmax), 10)
    n_bins = int(round(xmax / bin_width))
    edges = [i * bin_width for i in range(n_bins + 1)]
    edges[-1] = xmax

    plt.figure(figsize=(7, 4))
    plt.hist(depths_for_plot, bins=edges, range=(0.0, float(xmax)))
    plt.xlim(0.0, float(xmax))
    plt.xlabel("nitrogen depth (Å)")
    plt.ylabel("count")
    title_n = len(depths_for_plot)
    if max_depth_limit is None:
        plt.title(f"Nitrogen depth histogram (n={title_n})")
    else:
        plt.title(f"Nitrogen depth histogram (n={title_n}, <= {max_depth_limit:g} Å)")
    plt.grid(True, alpha=0.3)
    plt.margins(x=0)
    plt.tight_layout()

    out_arg = Path(args.out)
    if out_arg.is_absolute():
        out_path = out_arg
    elif out_arg.parent != Path("."):
        out_path = Path.cwd() / out_arg
    else:
        # 従来互換: ファイル名だけ指定された場合はスクリプトのあるフォルダに出す
        out_path = data_dir / out_arg
    out_path.parent.mkdir(parents=True, exist_ok=True)
    plt.savefig(out_path, dpi=200)
    print(f"Saved: {out_path}")


if __name__ == "__main__":
    main()

