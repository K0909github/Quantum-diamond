"""LAMMPS トラジェクトリから「図の状態」の N 個数を数える。

ユーザー要件に合わせて vacancy は考慮せず、N 原子の空間配置のみで判定する。

ここでの「図の状態」は次の条件を満たす N を指す:

1) N 原子（このデータでは dump の type=3）
2) 表面からの深さ depth = surface_z - z が一定以上（デフォルト 5 nm = 50 Å）
3) その N が、別の N と 5–15 nm の距離にある（磁気双極子相互作用を仮定）
4) 上の(3)を満たす 2つの N が、それ以外の N から 15 nm 以上離れている

入力:
- dump_run_1_min0K.lammpstrj (LAMMPS dump 形式)

例:
  python goodanalysis.py \
	--dump "10Ncluster_implantation to C_0.5keV(final)\\dump_run_1_min0K.lammpstrj" \
	--verbose
"""

from __future__ import annotations

import argparse
import math
from pathlib import Path


ANGSTROM_PER_NM = 10.0


def _to_float(token: str) -> float | None:
	try:
		return float(token)
	except ValueError:
		return None


def _pick_idx(cols_lower: list[str], candidates: list[str]) -> int | None:
	for name in candidates:
		if name in cols_lower:
			return cols_lower.index(name)
	return None


def read_last_dump_n_positions(
	dump_path: Path,
	nitrogen_type: int,
) -> list[tuple[float, float, float]]:
	"""LAMMPS dump から N(type) の (x,y,z) を読む（最後のスナップショット）。

	このリポジトリの dump は通常 `ITEM: ATOMS id type x y z` 形式なので、
	PBC/minimum image や scaled 座標の変換は行わない（ユーザー要望で簡略化）。
	"""

	positions: list[tuple[float, float, float]] = []
	current: list[tuple[float, float, float]] = []
	reading_atoms = False
	x_idx = y_idx = z_idx = None
	type_idx: int | None = None

	with dump_path.open("r", encoding="utf-8", errors="ignore") as handle:
		for raw in handle:
			line = raw.strip()
			if not line:
				continue

			if line.startswith("ITEM:"):
				if line.startswith("ITEM: ATOMS"):
					cols = line.split()[2:]
					cols_l = [c.lower() for c in cols]
					type_idx = cols_l.index("type") if "type" in cols_l else None
					x_idx = _pick_idx(cols_l, ["x", "xu"])
					y_idx = _pick_idx(cols_l, ["y", "yu"])
					z_idx = _pick_idx(cols_l, ["z", "zu"])
					current = []
					reading_atoms = True
					continue
				# 別の ITEM に入ったら、その時点の ATOMS を最後として採用
				if reading_atoms:
					positions = current
				reading_atoms = False
				continue

			if not reading_atoms or type_idx is None or x_idx is None or y_idx is None or z_idx is None:
				continue

			parts = line.split()
			if len(parts) <= max(type_idx, x_idx, y_idx, z_idx):
				continue

			t_val = _to_float(parts[type_idx])
			if t_val is None or int(round(t_val)) != nitrogen_type:
				continue
			x = _to_float(parts[x_idx])
			y = _to_float(parts[y_idx])
			z = _to_float(parts[z_idx])
			if x is None or y is None or z is None:
				continue
			current.append((x, y, z))

	# EOF: ATOMS の最中なら最後を採用
	if reading_atoms:
		positions = current

	return positions


def distance(a: tuple[float, float, float], b: tuple[float, float, float]) -> float:
	dx = a[0] - b[0]
	dy = a[1] - b[1]
	dz = a[2] - b[2]
	return math.sqrt(dx * dx + dy * dy + dz * dz)


def count_n_in_image_state(
	dump_path: Path,
	*,
	nitrogen_type: int,
	surface_z: float,
	min_depth_a: float,
	min_pair_dist_a: float,
	max_pair_dist_a: float,
	isolation_dist_a: float,
) -> tuple[int, dict[str, int]]:
	n_positions = read_last_dump_n_positions(dump_path, nitrogen_type)

	# 深さフィルタ（ペア探索対象）: depth >= min_depth
	deep_indices: list[int] = []
	for idx, (_, _, z) in enumerate(n_positions):
		depth = surface_z - z
		if depth >= min_depth_a:
			deep_indices.append(idx)

	# (3) 5–15 nm の距離にある N ペアを列挙
	candidate_pairs: list[tuple[int, int]] = []
	for a_i in range(len(deep_indices)):
		i = deep_indices[a_i]
		for a_j in range(a_i + 1, len(deep_indices)):
			j = deep_indices[a_j]
			d = distance(n_positions[i], n_positions[j])
			if min_pair_dist_a <= d <= max_pair_dist_a:
				candidate_pairs.append((i, j))

	# (4) その2つが「それ以外のN」から 15nm 以上離れているか
	qualifying: set[int] = set()
	for i, j in candidate_pairs:
		isolated = True
		for k in range(len(n_positions)):
			if k == i or k == j:
				continue
			dik = distance(n_positions[i], n_positions[k])
			if dik < isolation_dist_a:
				isolated = False
				break
			djk = distance(n_positions[j], n_positions[k])
			if djk < isolation_dist_a:
				isolated = False
				break
		if isolated:
			qualifying.add(i)
			qualifying.add(j)

	stats = {
		"N_total_in_dump": len(n_positions),
		"N_deep": len(deep_indices),
		"pair_candidates_5to15nm": len(candidate_pairs),
		"N_in_image_state": len(qualifying),
	}
	return len(qualifying), stats


def main() -> None:
	parser = argparse.ArgumentParser(description="Count nitrogen atoms that match the figure-like conditions.")
	parser.add_argument(
		"--dump",
		type=Path,
		default=Path(
			"10Ncluster_implantation to C_0.5keV(final)/dump_run_1_min0K.lammpstrj"
		),
		help="LAMMPS dump (.lammpstrj) path",
	)
	parser.add_argument("--nitrogen-type", type=int, default=3, help="Atom type ID for nitrogen in dump (default: 3)")
	parser.add_argument(
		"--surface-z",
		type=float,
		default=125.0,
		help="Surface z in Å used as depth = surface_z - z (default: 125.0)",
	)
	parser.add_argument(
		"--min-depth-nm",
		type=float,
		default=5.0,
		help="Minimum depth from surface in nm (default: 5.0)",
	)
	parser.add_argument(
		"--pair-min-nm",
		type=float,
		default=5.0,
		help="N-N 相互作用距離の下限 (nm) (default: 5.0)",
	)
	parser.add_argument(
		"--pair-max-nm",
		type=float,
		default=15.0,
		help="N-N 相互作用距離の上限 (nm) (default: 15.0)",
	)
	parser.add_argument(
		"--isolation-nm",
		type=float,
		default=15.0,
		help="(4) 2つの相互作用Nが他のNから離れている距離の下限 (nm) (default: 15.0)",
	)
	parser.add_argument(
		"--verbose",
		action="store_true",
		help="詳細な内訳も表示",
	)

	args = parser.parse_args()
	dump_path: Path = args.dump
	if not dump_path.exists():
		raise SystemExit(f"dump が見つかりません: {dump_path}")

	min_depth_a = float(args.min_depth_nm) * ANGSTROM_PER_NM
	min_pair_dist_a = float(args.pair_min_nm) * ANGSTROM_PER_NM
	max_pair_dist_a = float(args.pair_max_nm) * ANGSTROM_PER_NM
	isolation_dist_a = float(args.isolation_nm) * ANGSTROM_PER_NM

	n_count, stats = count_n_in_image_state(
		dump_path,
		nitrogen_type=int(args.nitrogen_type),
		surface_z=float(args.surface_z),
		min_depth_a=min_depth_a,
		min_pair_dist_a=min_pair_dist_a,
		max_pair_dist_a=max_pair_dist_a,
		isolation_dist_a=isolation_dist_a,
	)

	print(n_count)
	if args.verbose:
		for k, v in stats.items():
			print(f"{k}: {v}")


if __name__ == "__main__":
	main()

