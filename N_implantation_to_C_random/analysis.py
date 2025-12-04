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
	for line in lines:
		stripped = line.strip()
		if not stripped:
			if reading:
				break  # blank line marks end of Atoms block
			continue

		if stripped.lower().startswith("atoms"):
			reading = True
			continue

		if reading:
			if stripped[0].isalpha():  # hit next section
				break
			parts = stripped.split()
			if len(parts) < 5:
				continue
			atom_type = int(parts[1])
			z_coord = float(parts[4])
			atoms.append((atom_type, z_coord))

	return atoms


def compute_depths(data_path: Path) -> list[float]:
	"""Compute injection depths for all nitrogen atoms in a data file."""

	atoms = parse_atoms_section(data_path.read_text().splitlines())
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
