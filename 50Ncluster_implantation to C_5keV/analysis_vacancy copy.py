"""Plot vacancy distribution and compute mean vacancy distance.

This script reads a vacancy coordinate file exported from OVITO and:

1) outputs a bar-graph (histogram) of vacancy distribution along z
2) computes the mean distance between vacancies (mean of all unique pairs)

Supported input formats (best-effort):
- XYZ-like (first line = atom count, second = comment)
- CSV/TSV with columns like x,y,z (case-insensitive)
- whitespace-separated rows containing at least 3 numeric columns (uses last 3)
"""

from __future__ import annotations

import argparse
import csv
import math
from pathlib import Path


def _is_comment_or_empty(line: str) -> bool:
	stripped = line.strip()
	return not stripped or stripped.startswith("#")


def _to_float(token: str) -> float | None:
	try:
		return float(token)
	except ValueError:
		return None


def _extract_xyz_from_tokens(tokens: list[str]) -> tuple[float, float, float] | None:
	if len(tokens) < 3:
		return None

	# OVITO の vacancy_list などでよくある形式:
	#   ParticleType  x  y  z  Occupancy
	# 例: 1 -4.45875 -20.5102 74.0153 0
	# これを "last 3" で取ると (y, z, occupancy) になってしまうので補正する。
	if len(tokens) >= 5:
		t0 = _to_float(tokens[0])
		last = _to_float(tokens[-1])
		if t0 is not None and last is not None:
			is_int_like = abs(t0 - round(t0)) < 1e-9
			is_occupancy_like = last in {0.0, 1.0}
			if is_int_like and is_occupancy_like:
				x = _to_float(tokens[1])
				y = _to_float(tokens[2])
				z = _to_float(tokens[3])
				if x is not None and y is not None and z is not None:
					return float(x), float(y), float(z)

	# 列数4で occupancy が無い場合: type x y z
	if len(tokens) >= 4:
		t0 = _to_float(tokens[0])
		if t0 is not None and abs(t0 - round(t0)) < 1e-9:
			x = _to_float(tokens[1])
			y = _to_float(tokens[2])
			z = _to_float(tokens[3])
			if x is not None and y is not None and z is not None:
				return float(x), float(y), float(z)

	# フォールバック: 行末の3つを座標とみなす（XYZ など）
	last3 = tokens[-3:]
	vals = [_to_float(t) for t in last3]
	if any(v is None for v in vals):
		return None
	x, y, z = (vals[0], vals[1], vals[2])
	return float(x), float(y), float(z)


def _read_xyz_like(path: Path) -> list[tuple[float, float, float]]:
	coords: list[tuple[float, float, float]] = []
	with path.open("r", encoding="utf-8", errors="ignore") as handle:
		lines = [ln.rstrip("\n") for ln in handle]

	if not lines:
		return coords

	# XYZ: first line is an integer atom count
	try:
		int(lines[0].strip())
		start = 2
	except ValueError:
		start = 0

	for raw_line in lines[start:]:
		if _is_comment_or_empty(raw_line):
			continue
		tokens = raw_line.strip().split()
		xyz = _extract_xyz_from_tokens(tokens)
		if xyz is None:
			continue
		coords.append(xyz)

	return coords


def _read_csv_like(path: Path) -> list[tuple[float, float, float]]:
	coords: list[tuple[float, float, float]] = []
	with path.open("r", encoding="utf-8", errors="ignore", newline="") as handle:
		sample = handle.read(4096)
		handle.seek(0)

		dialect = csv.excel
		try:
			dialect = csv.Sniffer().sniff(sample, delimiters=",\t; ")
		except csv.Error:
			# Fall back to comma.
			dialect = csv.excel

		reader = csv.reader(handle, dialect)
		rows = list(reader)

	if not rows:
		return coords

	# Detect header.
	header = [h.strip().lower() for h in rows[0]]
	x_idx = y_idx = z_idx = None
	for i, name in enumerate(header):
		if name in {"x", "posx", "position.x", "position_x"}:
			x_idx = i
		elif name in {"y", "posy", "position.y", "position_y"}:
			y_idx = i
		elif name in {"z", "posz", "position.z", "position_z"}:
			z_idx = i

	data_rows = rows
	if x_idx is not None and y_idx is not None and z_idx is not None:
		data_rows = rows[1:]
		for r in data_rows:
			if len(r) <= max(x_idx, y_idx, z_idx):
				continue
			x = _to_float(r[x_idx].strip())
			y = _to_float(r[y_idx].strip())
			z = _to_float(r[z_idx].strip())
			if x is None or y is None or z is None:
				continue
			coords.append((x, y, z))
		return coords

	# No recognizable header: try treating as whitespace/CSV mixed and parse last 3 numeric.
	for r in rows:
		if not r:
			continue
		flat = " ".join(part for part in r if part is not None)
		if _is_comment_or_empty(flat):
			continue
		tokens = flat.strip().split()
		xyz = _extract_xyz_from_tokens(tokens)
		if xyz is None:
			continue
		coords.append(xyz)

	return coords


def read_vacancy_coordinates(path: Path) -> list[tuple[float, float, float]]:
	# Prefer CSV parser if extension suggests it; otherwise XYZ-like first.
	ext = path.suffix.lower()
	if ext in {".csv", ".tsv"}:
		coords = _read_csv_like(path)
		if coords:
			return coords
		return _read_xyz_like(path)

	coords = _read_xyz_like(path)
	if coords:
		return coords
	return _read_csv_like(path)


def mean_pairwise_distance(coords: list[tuple[float, float, float]]) -> float | None:
	n = len(coords)
	if n < 2:
		return None

	try:
		import numpy as np
	except ModuleNotFoundError:
		total = 0.0
		count = 0
		for i in range(n):
			xi, yi, zi = coords[i]
			for j in range(i + 1, n):
				xj, yj, zj = coords[j]
				dx = xi - xj
				dy = yi - yj
				dz = zi - zj
				total += math.sqrt(dx * dx + dy * dy + dz * dz)
				count += 1
		return total / count if count else None

	arr = np.asarray(coords, dtype=float)
	total = 0.0
	count = 0
	# O(N^2) but without allocating an NxN matrix.
	for i in range(n - 1):
		diffs = arr[i + 1 :] - arr[i]
		dists = np.sqrt(np.sum(diffs * diffs, axis=1))
		total += float(np.sum(dists))
		count += dists.size
	return total / count if count else None


def plot_z_histogram(
	z_values: list[float],
	out_path: Path,
	bin_width: float,
	surface_z: float,
	show: bool,
) -> None:
	try:
		import matplotlib.pyplot as plt
	except ModuleNotFoundError as exc:
		raise ModuleNotFoundError(
			"matplotlib が見つかりません。グラフ出力には matplotlib が必要です。\n"
			"例: pip install matplotlib"
		) from exc

	if not z_values:
		raise ValueError("vacancy 座標が空です。")

	# 位置計算（coords）はそのままに、描画だけ深さ座標に変換する。
	# depth = surface_z - z  (例: surface_z=125 Å)
	x_values = [surface_z - z for z in z_values]

	x_min = min(x_values)
	x_max = max(x_values)
	if bin_width <= 0:
		raise ValueError("bin_width は正の値にしてください。")

	start = math.floor(x_min / bin_width) * bin_width
	end = math.ceil(x_max / bin_width) * bin_width
	if end <= start:
		end = start + bin_width

	edges: list[float] = []
	x = start
	# include end edge
	while x <= end + 1e-9:
		edges.append(float(x))
		x += bin_width

	plt.figure(figsize=(7, 4))
	plt.hist(x_values, bins=edges)
	plt.xlabel(f"depth (Å) = {surface_z:g} - z")
	plt.ylabel("count")
	plt.title("Vacancy distribution (histogram along depth)")
	plt.grid(True, alpha=0.3)
	plt.margins(x=0)
	out_path.parent.mkdir(parents=True, exist_ok=True)
	plt.tight_layout()
	plt.savefig(out_path, dpi=200)
	print(f"Saved plot: {out_path}")
	if show:
		plt.show()


def main() -> None:
	parser = argparse.ArgumentParser(description=__doc__)
	parser.add_argument(
		"vacancy_file",
		type=Path,
		nargs="?",
		default=None,
		help="OVITO から出力した vacancy 座標ファイル (xyz/csv/空白区切りを想定)",
	)
	parser.add_argument(
		"--bin-width",
		type=float,
		default=10.0,
		help="ヒストグラムのビン幅 (Å) (default: 10)",
	)
	parser.add_argument(
		"--out",
		type=Path,
		default=None,
		help="出力画像パス (default: <input>/vacancy_distribution_depth.png)",
	)
	parser.add_argument(
		"--surface-z",
		type=float,
		default=125.0,
		help="深さ座標の基準面 z (Å)。depth = surface_z - z (default: 125)",
	)
	parser.add_argument(
		"--show",
		action="store_true",
		help="保存後にウィンドウ表示 (GUI backend が必要)",
	)
	args = parser.parse_args()

	input_path: Path
	if args.vacancy_file is None:
		# 引数なしで起動した場合は、スクリプトと同じフォルダにある vacancy_list を読む。
		default_path = Path(__file__).resolve().parent / "vacancy_list"
		if not default_path.exists():
			parser.error(
				"vacancy_file が指定されていません。\n"
				"引数なしで使う場合は、このスクリプトと同じフォルダに 'vacancy_list' を置いてください。\n"
				"例: python analysis_vacancy.py vacancy_list"
			)
		input_path = default_path
	else:
		input_path = args.vacancy_file

	coords = read_vacancy_coordinates(input_path)
	if not coords:
		print("vacancy 座標を読み取れませんでした。ファイル形式（xyz/csv）を確認してください。")
		return

	z_values = [z for _, _, z in coords]
	out_path = args.out or (input_path.parent / "vacancy_distribution_depth.png")
	plot_z_histogram(
		z_values,
		out_path=out_path,
		bin_width=args.bin_width,
		surface_z=args.surface_z,
		show=args.show,
	)

	mean_dist = mean_pairwise_distance(coords)
	print("-" * 40)
	print(f"vacancy count: {len(coords)}")
	if mean_dist is None:
		print("mean vacancy distance: N/A (vacancy が 2 個未満)")
	else:
		print(f"mean vacancy distance (all pairs): {mean_dist:.6f} Å")


if __name__ == "__main__":
	main()
