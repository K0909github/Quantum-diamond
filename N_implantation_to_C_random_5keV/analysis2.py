"""Compute nitrogen implantation depth assuming a fixed surface height."""

from __future__ import annotations

import argparse
from pathlib import Path
from statistics import mean

# diamond_substrate_250Å/in.diamond で作った基板(厚さ250Å)の上面は z=125Å
# 必要なら CLI オプションで上書きできます
DEFAULT_SURFACE_Z = 125.0  # Å


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


def compute_depths(data_path: Path, surface_z: float, include_above_surface: bool) -> list[float]:
    atoms = read_atoms(data_path)
    nitrogen_z = [z for atom_type, z in atoms if atom_type == 3]
    if not nitrogen_z:
        return []

    depths = [surface_z - z for z in nitrogen_z]
    if include_above_surface:
        return depths
    return [d for d in depths if d >= 0.0]


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
    args = parser.parse_args()

    data_files = sorted(args.data_dir.glob("final_data_run_*.data"))
    if not data_files:
        print("No final_data_run_*.data files found.")
        return

    all_depths: list[float] = []
    for data_file in data_files:
        depths = compute_depths(
            data_file,
            surface_z=args.surface_z,
            include_above_surface=args.include_above_surface,
        )
        if not depths:
            print(f"Skipping {data_file.name}: no nitrogen atoms found.")
            continue
        shot_mean = mean(depths)
        all_depths.extend(depths)
        print(f"{data_file.name}: mean depth = {shot_mean:.3f} Å (n={len(depths)})")

    if not all_depths:
        print("No nitrogen depths computed.")
        return

    print("-" * 40)
    print(f"Overall average depth: {mean(all_depths):.3f} Å across {len(all_depths)} atoms")


if __name__ == "__main__":
    main()
