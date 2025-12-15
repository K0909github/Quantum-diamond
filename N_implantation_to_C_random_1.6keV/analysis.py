"""Compute average implantation depth of nitrogen from final_data_run_*.data files."""

from __future__ import annotations

import argparse
from pathlib import Path
from statistics import mean


def read_atoms(path: Path) -> tuple[list[tuple[int, float]], float | None]:
	"""Parse atoms (type, z) and return along with box length in z."""

	z_length: float | None = None
	atoms: list[tuple[int, float]] = []
	reading_atoms = False
	seen_atoms = False
	atom_style: str | None = None

	with path.open("r", encoding="utf-8", errors="ignore") as handle:
		for raw_line in handle:
			line = raw_line.strip()
			if not line:
				if reading_atoms and seen_atoms:
					break
				continue

			tokens = line.split()
			lower_line = line.lower()

			if len(tokens) >= 4 and tokens[2].lower() == "zlo" and tokens[3].lower() == "zhi":
				z_length = float(tokens[1]) - float(tokens[0])
				continue

			if lower_line.startswith("atoms"):
				reading_atoms = True
				seen_atoms = False
				atom_style = None
				if "#" in line:
					atom_style = line.split("#", 1)[1].strip()
				continue

			if reading_atoms:
				if line[0].isalpha():
					break
				atom = _parse_atom_row(tokens, atom_style, z_length)
				if atom is None:
					continue
				atoms.append(atom)
				seen_atoms = True

	return atoms, z_length


def _parse_atom_row(
	tokens: list[str], atom_style: str | None, z_length: float | None
) -> tuple[int, float] | None:
	"""Extract (type, z) while handling optional image flags and styles."""

	style = (atom_style or "").lower()
	count = len(tokens)
	if count < 5:
		return None

	has_images = count >= 8 and all(_looks_like_int(tok) for tok in tokens[-3:])
	ix = iy = iz = 0
	if has_images:
		ix, iy, iz = map(int, tokens[-3:])

	try:
		if style.startswith("atomic") or not style:
			atom_type = int(tokens[1])
			z_coord = float(tokens[4])
		elif style.startswith("charge") or style.startswith("electron"):
			atom_type = int(tokens[1])
			z_coord = float(tokens[5])
		elif style.startswith("molecular"):
			if count < 6:
				return None
			atom_type = int(tokens[2])
			z_coord = float(tokens[5])
		elif style.startswith("full"):
			if count < 7:
				return None
			atom_type = int(tokens[2])
			z_coord = float(tokens[6])
		else:
			atom_type = int(tokens[1])
			z_coord = float(tokens[-4] if has_images else tokens[-1])
	except ValueError:
		return None

	if has_images and z_length is not None:
		z_coord += iz * z_length

	return atom_type, z_coord


def _looks_like_int(token: str) -> bool:
	try:
		int(token)
		return True
	except ValueError:
		return False


def compute_depths(data_path: Path) -> list[float]:
	"""Compute injection depths for all nitrogen atoms in a data file."""

	atoms, _ = read_atoms(data_path)
	if not atoms:
		return []

	substrate_z = [z for atom_type, z in atoms if atom_type in {1, 2}]
	nitrogen_z = [z for atom_type, z in atoms if atom_type == 3]
	if not substrate_z or not nitrogen_z:
		return []

	surface_z = _estimate_surface_level(substrate_z)
	return [surface_z - z for z in nitrogen_z]


def _estimate_surface_level(substrate_z: list[float]) -> float:
	sorted_z = sorted(substrate_z)
	count = len(sorted_z)
	if count <= 10:
		return max(sorted_z)

	top_n = min(count, max(50, count // 1000))
	return mean(sorted_z[-top_n:])


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
