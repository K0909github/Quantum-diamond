# 10回ランダム注入（ensemble）で分布を作る手順

このリポジトリは「1回の注入 → vacancy_list / N_list → ヒストグラム」になっていたので、
**注入位置 (x,y) を run ごとに変えた 10回分の結果をまとめて**ヒストグラム化できるようにしました。

## 何を変えたか

- [tools/run_ensemble_implantation.py](tools/run_ensemble_implantation.py)
  - 既存の LAMMPS 入力ファイル（テンプレ）を複製し、run ごとに
    - `variable seed equal ...`
    - `variable x_pos equal ...`
    - `variable y_pos equal ...`
    を書き換えた `run_01 .. run_10` フォルダを作ります。
  - さらに `--lammps "..."` を指定すれば各 run を連続実行もできます。

- [10Ncluster_implantation to C_0.5keV(final)/analysis3_graph_vacancy_list.py](10Ncluster_implantation%20to%20C_0.5keV(final)/analysis3_graph_vacancy_list.py)
  - 複数の `vacancy_list`（フォルダ指定/ワイルドカード指定）をまとめ読みして 1枚の分布ヒストグラムを出せます。

- [10Ncluster_implantation to C_0.5keV(final)/analysis2_graph_N_list3.py](10Ncluster_implantation%20to%20C_0.5keV(final)/analysis2_graph_N_list3.py)
  - 複数入力をまとめ読みできます。
  - LAMMPS dump を直接渡す場合に `--atom-type 3` のように type フィルタできるようにしました（dump に type 列がある前提）。

## 前提

- Python でグラフを保存するために `matplotlib` が必要です。
  - 例: `python -m pip install matplotlib`

## 1) 10回分の run フォルダを作る（注入位置をランダム化）

リポジトリルートで:

```powershell
python tools/run_ensemble_implantation.py `
  --template-dir "10Ncluster_implantation to C_5keV" `
  --input "10atoms_5keV_N_implantation_to_C_ZBL_potential_filedata.txt" `
  --out "10Ncluster_implantation to C_5keV/runs" `
  --runs 10 `
  --seed 12345 `
  --x-range -20 20 `
  --y-range -20 20
```

これで `10Ncluster_implantation to C_5keV/runs/run_01/in.lmp` のように生成されます。

### そのまま LAMMPS も10回回す場合（任意）

LAMMPS のコマンド名は環境により違います（`lmp`, `lmp_mpi`, `lammps` など）。例:

```powershell
python tools/run_ensemble_implantation.py `
  --template-dir "10Ncluster_implantation to C_5keV" `
  --input "10atoms_5keV_N_implantation_to_C_ZBL_potential_filedata.txt" `
  --out "10Ncluster_implantation to C_5keV/runs" `
  --runs 10 `
  --seed 12345 `
  --x-range -20 20 `
  --y-range -20 20 `
  --lammps "lmp -in in.lmp"
```

## 2) vacancy の深さ分布（10回分をまとめて）

`run_*` をフォルダ指定で渡せます:

```powershell
python "10Ncluster_implantation to C_0.5keV(final)/analysis3_graph_vacancy_list.py" `
  "10Ncluster_implantation to C_5keV/runs/run_*" `
  --out "vacancy_depths_ensemble.png" `
  --surface-z 125 `
  --bin-width 1
```

## 3) N の深さ分布（10回分をまとめて）

`dump_run_1_min0K.lammpstrj` を直接読む例（N が type=3 の場合）:

```powershell
python "10Ncluster_implantation to C_0.5keV(final)/analysis2_graph_N_list3.py" `
  "10Ncluster_implantation to C_5keV/runs/run_*/dump_run_1_min0K.lammpstrj" `
  --atom-type 3 `
  --out "nitrogen_depths_ensemble.png" `
  --surface-z 125 `
  --bin-width 5
```

## よくある調整ポイント

- `--x-range/--y-range` は「基板の有効領域（表面の中央付近）」に合わせて調整してください。
- `--surface-z` は基板モデルの表面 z に合わせてください（いまは既存スクリプトと同じ `125 Å` をデフォルトにしています）。
