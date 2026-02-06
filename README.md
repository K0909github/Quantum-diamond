# Quantum-diamond

LAMMPS による分子動力学（MD）シミュレーションを用いて、ダイヤモンド中の欠陥（vacancy）生成や浅い窒素（N）/シリコン（Si）注入を扱う研究用リポジトリです。

このREADMEは「どのフォルダから何を実行し、どう解析するか」を最短で辿れるようにまとめています（詳細手順は別ドキュメントも参照）。

## できること（ざっくり）

- ダイヤモンド/シリコン基板に対する注入（単発・連続注入・クラスター注入など）の条件別シミュレーション
- 複数run（ensemble）を自動生成し、runごとに注入位置・乱数seedを変えて統計分布を作る
- 解析（例）
  - `vacancy_list` の深さ分布ヒストグラム
  - `N_list` または LAMMPS dump から N の深さ分布ヒストグラム
  - dump から「図の状態」に該当する N の個数カウント

## 前提（環境）

- LAMMPS（`lmp`, `lmp_mpi`, `lammps` など。環境でコマンド名が違います）
- （任意）MPI 実行環境（`mpirun`）
- （任意）OVITO（可視化・vacancy抽出など）
- Python 3.10+（解析スクリプトが型ヒント `|` を使います）
  - グラフ保存に `matplotlib` が必要です: `python -m pip install matplotlib`

注意: このリポジトリは OneDrive 配下など「パスにスペースを含む」環境で使うことが多いです。コマンド引数のパスは必ずダブルクォートで囲むか、先に対象フォルダへ `cd` してから実行してください。

## フォルダ構造（主要なもの）

- `diamond_substrate/`, `diamond_substrate_250Å/`: ダイヤモンド基板関連
- `silicon_substrate/`: Si 基板関連
- `Ar_implantation to Si/`: Ar→Si 注入系
- `N_implantation_to_C(final)/`, `N_implantation_to_C_random_*/`: N→C（ダイヤモンド）注入系（ランダム注入/連続注入の派生あり）
- `10Ncluster_implantation to C_*/`, `50Ncluster_implantation to C_5keV/`, `Ncluster_implantation to C_5keV/`: Nクラスター注入系
- `tools/`: ensemble 用の補助スクリプト（runフォルダ生成・自動実行）
- `Tutorials/`: 参考資料
- `*.pdf`: 関連論文・メモ（研究ノート用途）

（フォルダ名に `final` が付くものは、研究の途中経過を整理した版を含みます）

## 最短クイックスタート（ensembleで分布を作る）

詳しい説明は [RUNS_ENSEMBLE_HOWTO.md](RUNS_ENSEMBLE_HOWTO.md) にあります。ここでは最短の流れだけ示します。

### 1) runフォルダを作る（注入位置をrunごとにランダム化）

例（PowerShell）: 既存テンプレを複製し、`run_01`〜`run_10` を作ります。

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

このスクリプトはテンプレ入力内の以下の行を探して run ごとに上書きします（simpleモード）:

- `variable seed equal ...`
- `variable x_pos equal ...`
- `variable y_pos equal ...`

連続注入（noreset）系のテンプレで、入力内に `random(...)` があるタイプは `--template-style loop-random-xy` を使います（例はHOWTO参照）。

### 2) LAMMPS を実行する

runフォルダに入って実行するのが安全です（相対パス参照のため）。例:

```powershell
cd "10Ncluster_implantation to C_5keV/runs/run_01"
mpirun -np 8 lmp -in in.lmp
```

`run_ensemble_implantation.py` に `--lammps "..."` を渡すと、生成後に各runを連続実行もできます。

## 解析（深さ分布など）

### Vacancy の深さ分布（複数runをまとめて）

`vacancy_list`（OVITO出力など）を run ごとに用意できている前提で、まとめてヒストグラム化します。

```powershell
python "10Ncluster_implantation to C_0.5keV(final)/analysis3_graph_vacancy_list.py" `
  "10Ncluster_implantation to C_5keV/runs/run_*" `
  --out "vacancy_depths_ensemble.png" `
  --surface-z 125 `
  --bin-width 1
```

### N の深さ分布（dumpから直接読む）

`analysis2_graph_N_list3.py` は、`N_list` だけでなく LAMMPS dump（`*.lammpstrj`）も読めます。dump側に `type` 列がある場合、`--atom-type 3` のように N だけ抽出できます。

```powershell
python "10Ncluster_implantation to C_0.5keV(final)/analysis2_graph_N_list3.py" `
  "10Ncluster_implantation to C_0.5keV(final)/dump_run_1_min0K.lammpstrj" `
  --atom-type 3 `
  --out "nitrogen_depths.png" `
  --surface-z 125 `
  --bin-width 5
```

### 「図の状態」の N 個数カウント

dump の最後のスナップショットを読み、条件（深さ・N-N距離・孤立性）でカウントします。

```powershell
python goodanalysis.py --dump "10Ncluster_implantation to C_0.5keV(final)/dump_run_1_min0K.lammpstrj" --verbose
```

## パラメータの目安（surface_z など）

- 深さは `depth = surface_z - z`（Å）で定義しています。
- 既定の `surface_z=125 Å` は、`diamond_substrate_250Å` 系の「上面」を想定した値です。
  - 基板厚や座標原点が違うデータでは、解析時に `--surface-z` を調整してください。

## よくあるトラブル

- パスにスペースがあって失敗する: パスを必ずダブルクォートで囲む／先に対象フォルダへ `cd`
- ポテンシャルファイルが見つからない: 入力が相対パス参照のことが多いので、テンプレと同じフォルダで実行する
- `matplotlib` が無い: `python -m pip install matplotlib`

## 参考リンク

- LAMMPS チュートリアル（外部）: https://github.com/mrkllntschpp/lammps-tutorials
- LAMMPS公式ポテンシャルファイル集: https://github.com/lammps/lammps/tree/develop/potentials
