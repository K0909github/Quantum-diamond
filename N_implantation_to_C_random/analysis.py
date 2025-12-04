"""Compute average implantation depth of nitrogen from final_data_run_*.data files."""

from __future__ import annotations

import argparse
from pathlib import Path
from statistics import mean
from typing import Iterable


def parse_atoms_section(lines: Iterable[str]) -> list[tuple[int, float]]:
	"""Return list of (type, z) tuples parsed from a LAMMPS data file."""

	atoms: list[tuple[int, float]] = []
	reading = False
	seen_data = False
	atom_style = None

	for line in lines:
		stripped = line.strip()
		if not stripped:
			if reading and seen_data:
				break  # blank line after data terminates section
			continue

		lower = stripped.lower()
		if lower.startswith("atoms"):
			reading = True
			seen_data = False
			atom_style = None
			if "#" in stripped:
				atom_style = stripped.split("#", 1)[1].strip().split()[0].lower()
			continue

		if reading:
			if stripped[0].isalpha():  # reached next section header
				break
			parts = stripped.split()
			atom_info = _extract_type_and_z(parts, atom_style)
			if atom_info is None:
				continue
			atoms.append(atom_info)
			seen_data = True

	return atoms


def _extract_type_and_z(parts: list[str], atom_style: str | None) -> tuple[int, float] | None:
	"""Handle common atom_style layouts (atomic, charge, molecular, full)."""

	cols = len(parts)
	if cols < 5:
		return None

	style = (atom_style or "").lower()

	try:
		if style.startswith("atomic"):
			# id type x y z (ix iy iz optional)
			return int(parts[1]), float(parts[4])
		if style.startswith("charge") or style.startswith("electron"):
			# id type q x y z
			if cols < 6:
				return None
			return int(parts[1]), float(parts[5])
		if style.startswith("molecular"):
			# id mol type x y z
			if cols < 6:
				return None
			return int(parts[2]), float(parts[5])
		if style.startswith("full"):
			# id mol type q x y z
			if cols < 7:
				return None
			return int(parts[2]), float(parts[6])

		# fallback heuristic: assume second entry is type, z just before image flags
		atom_type = int(parts[1])
		if cols >= 8 and all(_looks_like_int(tok) for tok in parts[-3:]):
			z_index = cols - 4
		else:
			z_index = cols - 1
		return atom_type, float(parts[z_index])
	except ValueError:
		return None


def _looks_like_int(token: str) -> bool:
	try:
		int(token)
		return True
	except ValueError:
		return False


def compute_depths(data_path: Path) -> list[float]:
	"""Compute injection depths for all nitrogen atoms in a data file."""

	with data_path.open("r", encoding="utf-8", errors="ignore") as fh:
		atoms = parse_atoms_section(fh)
	if not atoms:
		return []

	substrate_z = [z for atom_type, z in atoms if atom_type in {1, 2}]
	nitrogen_z = [z for atom_type, z in atoms if atom_type == 3]
	if not substrate_z or not nitrogen_z:
		return []

	surface_z = max(substrate_z)
	return [surface_z - z for z in nitrogen_z]


def main() -> None:
	parser = argparse.ArgumentParser(description=__doc__)
	parser.add_argument(
		"--data-dir",
		type=Path,
		default=Path(__file__).parent,
		help="Directory containing final_data_run_*.data files.",
	)
	args = parser.parse_args()

	data_files = sorted(args.data_dir.glob("final_data_run_*.data"))
	if not data_files:
		print("No final_data_run_*.data files found.")
		return

	all_depths: list[float] = []
	for data_file in data_files:
		depths = compute_depths(data_file)
		if not depths:
			print(f"Skipping {data_file.name}: could not extract nitrogen depth.")
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
