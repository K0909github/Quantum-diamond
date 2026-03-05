"""Compute nitrogen implantation depths and plot them.

This script reads LAMMPS data files written by `write_data final_data_run_*.data`,
extracts nitrogen atoms (type=3), computes implantation depth as:

    depth = surface_z - z

and outputs a per-run depth graph.
"""

from __future__ import annotations

import argparse
import math
import re
from pathlib import Path
from statistics import mean

DEFAULT_SURFACE_Z = 125.0  # Å  (diamond_substrate_250Å の上面)
DEFAULT_BIN_WIDTH = 10  # Å  (ヒストグラム表示を 0,10,20,... に揃える)


def read_atoms(path: Path) -> list[tuple[int, float]]:
    atoms: list[tuple[int, float]] = []
    reading = False
    atom_style: str | None = None

    with path.open("r", encoding="utf-8", errors="ignore") as handle:
        for raw_line in handle:
            line = raw_line.strip()
            if not line:
                if reading and atoms:
                    break
                continue

            lower = line.lower()
            if lower.startswith("atoms"):
                reading = True
                atom_style = None
                if "#" in line:
                    atom_style = line.split("#", 1)[1].strip().split()[0].lower()
                continue

            if reading:
                if line[0].isalpha():
                    break
                tokens = line.split()
                atom = _parse_atom(tokens, atom_style)
                if atom is None:
                    continue
                atoms.append(atom)

    return atoms


def _parse_atom(tokens: list[str], atom_style: str | None) -> tuple[int, float] | None:
    count = len(tokens)
    if count < 5:
        return None

    style = (atom_style or "").lower()
    has_images = count >= 8 and all(_looks_like_int(tok) for tok in tokens[-3:])

    try:
        if style.startswith("atomic") or not style:
            atom_type = int(tokens[1])
            z_coord = float(tokens[4])
        elif style.startswith("charge") or style.startswith("electron"):
            atom_type = int(tokens[1])
            z_coord = float(tokens[5])
        elif style.startswith("molecular"):
            atom_type = int(tokens[2])
            z_coord = float(tokens[5])
        elif style.startswith("full"):
            atom_type = int(tokens[2])
            z_coord = float(tokens[6])
        else:
            atom_type = int(tokens[1])
            z_coord = float(tokens[-4] if has_images else tokens[-1])
    except ValueError:
        return None

    return atom_type, z_coord


def _looks_like_int(token: str) -> bool:
    try:
        int(token)
        return True
    except ValueError:
        return False


def compute_depths(data_path: Path) -> list[float]:
    atoms = read_atoms(data_path)
    nitrogen_z = [z for atom_type, z in atoms if atom_type == 3]
    if not nitrogen_z:
        return []
    return [DEFAULT_SURFACE_Z - z for z in nitrogen_z]


def compute_depths_with_surface(
    data_path: Path,
    surface_z: float,
    include_above_surface: bool,
) -> list[float]:
    atoms = read_atoms(data_path)
    nitrogen_z = [z for atom_type, z in atoms if atom_type == 3]
    if not nitrogen_z:
        return []

    depths = [surface_z - z for z in nitrogen_z]
    if include_above_surface:
        return depths
    return [d for d in depths if d >= 0.0]


def _extract_run_index(path: Path) -> int | None:
    match = re.search(r"final_data_run_(\d+)\.data$", path.name)
    if not match:
        return None
    return int(match.group(1))


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--data-dir",
        type=Path,
        default=Path(__file__).parent,
        help="Directory containing final_data_run_*.data files.",
    )
    parser.add_argument(
        "--surface-z",
        type=float,
        default=DEFAULT_SURFACE_Z,
        help="Surface z position in Å used as depth=surface_z-z (default: 125.0).",
    )
    parser.add_argument(
        "--include-above-surface",
        action="store_true",
        help="Include nitrogen atoms above the surface (negative depths).",
    )
    parser.add_argument(
        "--out",
        type=Path,
        default=None,
        help="Output image path (default: <data-dir>/implantation_depths.png).",
    )
    parser.add_argument(
        "--show",
        action="store_true",
        help="Show the plot window after saving (requires GUI backend).",
    )
    parser.add_argument(
        "--bins",
        type=str,
        default="auto",
        help="Histogram bins passed to matplotlib (e.g. 'auto', 'fd', or an integer like '20').",
    )
    args = parser.parse_args()

    data_files = sorted(args.data_dir.glob("final_data_run_*.data"))
    if not data_files:
        print("No final_data_run_*.data files found.")
        return

    all_depths: list[float] = []

    for i, data_file in enumerate(data_files, start=1):
        depths = compute_depths_with_surface(
            data_file,
            surface_z=args.surface_z,
            include_above_surface=args.include_above_surface,
        )
        if not depths:
            print(f"Skipping {data_file.name}: no nitrogen atoms found.")
            continue

        shot_depth = mean(depths)
        all_depths.extend(depths)
        print(f"{data_file.name}: mean depth = {shot_depth:.3f} Å (n={len(depths)})")

    if not all_depths:
        print("No nitrogen depths computed.")
        return

    print("-" * 40)
    print(f"Overall average depth: {mean(all_depths):.3f} Å across {len(all_depths)} atoms")

    try:
        import matplotlib.pyplot as plt
    except ModuleNotFoundError:
        print("matplotlib が見つかりません。グラフ出力には matplotlib が必要です。")
        print("例: pip install matplotlib")
        return

    bins: str | int
    try:
        bins = int(args.bins)
    except ValueError:
        bins = args.bins

    min_depth = min(all_depths)
    max_depth = max(all_depths)
    if max_depth <= 0.0:
        max_depth = 1.0

    plt.figure(figsize=(7, 4))
    if min_depth >= 0.0:
        if isinstance(bins, str) and bins.lower() == "auto":
            # auto はビン境界が 20, 30... などのキリ番に揃わず、棒がズレて見えることがある。
            # デフォルト時は 10Å 刻みの境界に固定して見た目を安定させる。
            xmax = int(math.ceil(max_depth / DEFAULT_BIN_WIDTH) * DEFAULT_BIN_WIDTH)
            edges = list(range(0, xmax + DEFAULT_BIN_WIDTH, DEFAULT_BIN_WIDTH))
            plt.hist(all_depths, bins=edges, range=(0.0, float(xmax)))
            plt.xlim(0.0, float(xmax))
            plt.xticks(edges)
        else:
            xmax = math.ceil(max_depth)
            plt.hist(all_depths, bins=bins, range=(0.0, float(xmax)))
            plt.xlim(0.0, float(xmax))
    else:
        plt.hist(all_depths, bins=bins)
    plt.xlabel("implantation depth (Å)")
    plt.ylabel("count")
    plt.title("Nitrogen implantation depth histogram")
    plt.grid(True, alpha=0.3)
    plt.margins(x=0)

    out_path = args.out or (args.data_dir / "implantation_depths.png")
    out_path.parent.mkdir(parents=True, exist_ok=True)
    plt.tight_layout()
    plt.savefig(out_path, dpi=200)
    print(f"Saved plot: {out_path}")
    if args.show:
        plt.show()


if __name__ == "__main__":
    main()
