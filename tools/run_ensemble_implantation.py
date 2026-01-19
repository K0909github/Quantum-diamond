"""Generate an ensemble of LAMMPS input folders with randomized implantation positions.

目的
- 既存の LAMMPS 入力（例: 10atoms_5keV_...txt）をテンプレートとして
  x_pos / y_pos / seed を run ごとに変えたフォルダ(run_01, run_02, ...)を作る。
- LAMMPS を「自動実行」も可能（--lammps を指定した場合）。

想定
- テンプレ入力に以下の行が存在する（少なくとも seed/x_pos/y_pos）：
    variable        seed equal 12345
    variable        x_pos equal 0.0
    variable        y_pos equal 0.0

使い方例（フォルダ生成だけ）
  python tools/run_ensemble_implantation.py \
    --template-dir "10Ncluster_implantation to C_5keV" \
    --input "10atoms_5keV_N_implantation_to_C_ZBL_potential_filedata.txt" \
    --out "10Ncluster_implantation to C_5keV/runs" \
    --runs 10 \
    --seed 12345 \
    --x-range -20 20 --y-range -20 20

使い方例（LAMMPS も回す）
  python tools/run_ensemble_implantation.py ... \
    --lammps "lmp -in in.lmp"

Windows PowerShell ではダブルクォートで囲むのが無難です。
"""

from __future__ import annotations

import argparse
import random
import re
import shlex
import subprocess
from dataclasses import dataclass
from pathlib import Path
from typing import Iterable


@dataclass(frozen=True)
class RunSpec:
    index: int
    seed: int
    x_pos: float
    y_pos: float


_VAR_RE = re.compile(r"^(?P<indent>\s*)variable\s+(?P<name>seed|x_pos|y_pos)\s+equal\s+.*$", re.IGNORECASE)


def _read_text(path: Path) -> str:
    return path.read_text(encoding="utf-8", errors="ignore")


def _write_text(path: Path, text: str) -> None:
    path.write_text(text, encoding="utf-8")


def _patch_lammps_input(template_text: str, run: RunSpec) -> str:
    replaced = {"seed": False, "x_pos": False, "y_pos": False}

    out_lines: list[str] = []
    for raw_line in template_text.splitlines(keepends=False):
        m = _VAR_RE.match(raw_line)
        if not m:
            out_lines.append(raw_line)
            continue

        indent = m.group("indent")
        name = m.group("name").lower()
        if name == "seed":
            out_lines.append(f"{indent}variable        seed equal {run.seed}")
            replaced["seed"] = True
        elif name == "x_pos":
            out_lines.append(f"{indent}variable        x_pos equal {run.x_pos:.6f}")
            replaced["x_pos"] = True
        elif name == "y_pos":
            out_lines.append(f"{indent}variable        y_pos equal {run.y_pos:.6f}")
            replaced["y_pos"] = True
        else:
            out_lines.append(raw_line)

    missing = [k for k, v in replaced.items() if not v]
    if missing:
        raise ValueError(
            "テンプレ入力に variable 行が見つかりません: " + ", ".join(missing)
        )

    return "\n".join(out_lines) + "\n"


def _copy_files(src_dir: Path, dst_dir: Path, patterns: Iterable[str]) -> None:
    dst_dir.mkdir(parents=True, exist_ok=True)
    for pat in patterns:
        for p in src_dir.glob(pat):
            if p.is_file():
                (dst_dir / p.name).write_bytes(p.read_bytes())


def _parse_lammps_cmd(cmd: str) -> list[str]:
    # Windows でも動くよう shlex を posix=False
    return shlex.split(cmd, posix=False)


def main() -> None:
    parser = argparse.ArgumentParser(
        description="LAMMPS入力を複製し、注入位置(x,y)とseedをランダムに変えたrunフォルダを作る"
    )
    parser.add_argument("--template-dir", required=True, help="テンプレがあるフォルダ")
    parser.add_argument("--input", required=True, help="テンプレ入力ファイル名（template-dir内）")
    parser.add_argument("--out", required=True, help="出力先フォルダ（run_01等を作る）")
    parser.add_argument("--runs", type=int, default=10, help="試行回数")
    parser.add_argument("--seed", type=int, default=12345, help="乱数の基準seed")
    parser.add_argument("--x-range", type=float, nargs=2, default=[-20.0, 20.0], metavar=("XMIN", "XMAX"))
    parser.add_argument("--y-range", type=float, nargs=2, default=[-20.0, 20.0], metavar=("YMIN", "YMAX"))
    parser.add_argument(
        "--copy",
        nargs="*",
        default=["*.data", "*.zbl", "*.tersoff*"],
        help="各runフォルダにコピーするファイルのglob（デフォルト: *.data, *.zbl, *.tersoff*）",
    )
    parser.add_argument(
        "--lammps",
        type=str,
        default=None,
        help="指定すると各runでLAMMPSを実行する（例: 'lmp -in in.lmp'）",
    )
    args = parser.parse_args()

    template_dir = Path(args.template_dir)
    if not template_dir.exists():
        raise FileNotFoundError(f"テンプレフォルダが見つかりません: {template_dir}")

    input_arg = Path(args.input)
    template_input = input_arg if input_arg.is_absolute() else (template_dir / input_arg)
    if not template_input.exists():
        # 似た名前のファイルを候補として出す（ユーザーの指定ミスを直しやすくする）
        txt_candidates = sorted(template_dir.glob("*.txt"))
        data_candidates = sorted(template_dir.glob("*.data"))
        other_candidates = sorted(template_dir.glob("*.lmp"))
        candidates = [*txt_candidates, *data_candidates, *other_candidates]
        shown = "\n".join(f"  - {p.name}" for p in candidates[:30])
        if not shown:
            shown = "  (候補ファイルが見つかりませんでした)"

        raise FileNotFoundError(
            "テンプレ入力が見つかりません。指定値を確認してください。\n"
            f"  template-dir: {template_dir}\n"
            f"  input:        {args.input}\n"
            f"  tried:        {template_input}\n"
            "候補（template-dir内）:\n"
            f"{shown}"
        )

    out_root = Path(args.out)
    out_root.mkdir(parents=True, exist_ok=True)

    rng = random.Random(args.seed)
    xmin, xmax = float(args.x_range[0]), float(args.x_range[1])
    ymin, ymax = float(args.y_range[0]), float(args.y_range[1])

    template_text = _read_text(template_input)

    runs: list[RunSpec] = []
    for i in range(1, int(args.runs) + 1):
        x = rng.uniform(xmin, xmax)
        y = rng.uniform(ymin, ymax)
        seed_i = int(args.seed) + i
        runs.append(RunSpec(index=i, seed=seed_i, x_pos=x, y_pos=y))

    lammps_cmd: list[str] | None = None
    if args.lammps:
        lammps_cmd = _parse_lammps_cmd(args.lammps)

    for run in runs:
        run_dir = out_root / f"run_{run.index:02d}"
        run_dir.mkdir(parents=True, exist_ok=True)

        # 必要ファイルをコピー
        _copy_files(template_dir, run_dir, args.copy)

        # 入力ファイルを書き出し
        patched = _patch_lammps_input(template_text, run)
        in_path = run_dir / "in.lmp"
        _write_text(in_path, patched)

        # メモも残す
        _write_text(
            run_dir / "run_params.txt",
            f"index={run.index}\nseed={run.seed}\nx_pos={run.x_pos}\ny_pos={run.y_pos}\n",
        )

        print(f"Prepared: {run_dir} (seed={run.seed}, x={run.x_pos:.3f}, y={run.y_pos:.3f})")

        if lammps_cmd is not None:
            print(f"Running LAMMPS in {run_dir} ...")
            subprocess.run(lammps_cmd, cwd=str(run_dir), check=True)

    print("Done.")


if __name__ == "__main__":
    main()
