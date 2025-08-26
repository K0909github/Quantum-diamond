# NVセンターを用いた量子コンピュータ研究

これは卒業研究で行った LAMMPS (Large-scale Atomic/Molecular Massively Parallel Simulator) による分子動力学（MD）シミュレーションの記録です。

## 概要

本プロジェクトでは、シリコン（Si）および炭素（C, ダイヤモンド）基板の作成と、それらの基板へのアルゴン（Ar）および窒素（N）原子の照射・欠陥形成プロセスをシミュレートします。NV センターに関連する初期・最終座標も含まれます。

## ディレクトリ構造（要点）

- `si/`: 純Si基板の作成・緩和（`in.silicon_substrate2`, `in.silicon_substrate3`、ポテンシャル`Si.sw`）
- `C/`: ダイヤモンド(C)基板の作成（`in.carbon`、ポテンシャル`SiC.tersoff`、可視化`C.ovito`）
- `C+N/`: C基板へのN照射/堆積（`in.carbon_nitride`, `in.cn_irradiation` ほか、ポテンシャル`BNC.tersoff`）
- `si+Ar/`: Si+Ar照射の結果（`dump.si_ar.atom`, `si+ar.txt`, `si+ar.ovito`）
- `NV/`: NVセンター関連の座標とメモ（`dump.initial_coords`, `dump.final_nv_centers`, `NV.txt`）
- `irradiation/`: 照射条件やログの整理（`irradiation.txt`, `log.lammps` など）
- `Si2/`: 追加の計算ログ（`si2.txt` など）

## 主な実行方法（例）

前提:
- LAMMPS（MPI 版）がインストールされていること（コマンド名は環境により `lmp`, `lmp_mpi`, `lammps` など）。
- MPI 実行環境（Open MPI など）。
- 可能なら OVITO（可視化）。

注意: 本リポジトリはパスにスペースを含みます（OneDrive配下）。`-in` に絶対パスを渡す場合は必ずクォートするか、先に対象ディレクトリへ `cd` してください。

### 1) ディレクトリへ移動して実行（推奨）

Si 基板（例）:
```bash
cd "si"
mpirun -np 4 lmp -in in.silicon_substrate2
# もしくは
mpirun -np 4 lmp -in in.silicon_substrate3
```

ダイヤモンド基板（例）:
```bash
cd "C"
mpirun -np 4 lmp -in in.carbon
```

C+N（窒素照射/堆積, 例）:
```bash
cd "C+N"
mpirun -np 4 lmp -in in.carbon_nitride
# もしくは
mpirun -np 4 lmp -in in.cn_irradiation
```

### 2) フルパス指定で実行（パスをクォート）

```bash
mpirun -np 4 lmp -in \
"/mnt/c/Users/kazuyuki/OneDrive - 兵庫県立大学/Quantum-diamond/si/in.silicon_substrate2"
```

コマンド名が `lmp` でない場合は、環境のコマンドに置き換えてください（例: `lmp_mpi` や `lammps`）。

## 主な出力ファイル

| ファイル/拡張子 | 中身 | 主な用途 |
| --- | --- | --- |
| `log.lammps` | 実行ログ（熱力学量、タイムステップなど） | 実行の記録・確認 |
| `dump.*.atom` など | スナップショット（原子座標/速度 等） | OVITO/解析に使用 |
| `*.ovito` | OVITO セッション | 可視化設定の再現 |
| `*.txt` | 条件メモ/解析結果 | 実験条件の記録 |
| `dump.initial_coords`, `dump.final_nv_centers` | NV 関連の初期/最終座標 | 欠陥解析 |

OVITO での可視化: `dump.*` もしくは `*.ovito` を開いて確認します。

## 可視化手順（OVITO）

### 共通の流れ

1. OVITO を起動し、`File > Import...` から LAMMPS dump（`dump.*.atom` など）または `*.ovito` を開く。
2. タイムステップが複数ある場合は下部のタイムラインで再生・移動。
3. 基本的な表示設定:
	 - `Add modifier > Color coding` で `Particle Type` による色分け（元素ごとに確認しやすい）。
	 - `Add modifier > Slice` で断面表示（厚みを少し持たせると見やすい）。
	 - `Add modifier > Coordination analysis` で配位数分布（sp2/sp3 の傾向を見る際に有用）。
	 - 必要に応じて `Add modifier > Displacement vectors`（初期フレームを参照に原子の変位を可視化）。
4. 高解像度で保存: `Viewport > Save snapshot...`。

注: dump に原子種（type）が含まれていれば `Particle Type` で色分け可能です。原子ごとのエネルギー等は dump に出力されていないと使えません。

### si/（Si 基板）

- 推奨ファイル: `si/dump.atom`（存在する場合）、または `si/Si.ovito`（設定済みセッション）。
- 表示例:
	- `Color coding`（by Particle Type）: Si のみであれば単色、欠陥や混入があればタイプごとに色分け。
	- `Slice` で層構造を観察。
	- 欠陥観察には `Coordination analysis`（理想 Si: 配位数4 付近）。

### C/（ダイヤモンド基板）

- 推奨ファイル: `C/dump.carbon.atom`、または `C/C.ovito`。
- 表示例:
	- `Color coding`（by Particle Type）でCのみ表示。異種原子が無ければ単色でOK。
	- `Coordination analysis` で配位数が3/4に分かれるかを確認（sp2/sp3 の目安）。
	- `Slice` で格子方向を見やすく。

### C+N/（C基板へのN照射/堆積）

- 推奨ファイル: 実行後に生成された `dump.*`（入力により名前が異なる）。
- 例のタイプ割当（必要に応じて調整）:
	- type 1: C、type 2: N（dump のヘッダや入力スクリプトの `mass`/`create_atoms` で確認）。
- 表示例:
	- `Color coding`（by Particle Type）で C と N を別色に。
	- `Displacement vectors` で照射後の変位を可視化。
	- `Coordination analysis` で欠陥周辺の配位変化を確認。

### si+Ar/（Si への Ar 照射）

- 推奨ファイル: `si+Ar/dump.si_ar.atom`、または `si+Ar/si+ar.ovito`。
- 例のタイプ割当（必要に応じて調整）:
	- type 1: Si、type 2: Ar。
- 表示例:
	- `Color coding`（by Particle Type）で Si と Ar を分離表示。
	- `Slice` で侵入深さ・損傷領域を観察。
	- `Displacement vectors` で損傷カスケードの可視化。

### NV/（NV センター関連）

- 推奨ファイル: `NV/dump.initial_coords`, `NV/dump.final_nv_centers`。
- 表示例:
	- `Color coding`（by Particle Type）で C と N を色分け。
	- NV 近傍を `Slice` や `Selection`（球選択）で拡大して局所構造を観察。
	- 初期/最終フレームを並べて変化箇所を比較。

### 既存の OVITO セッションを使う

- `si/Si.ovito`, `C/C.ovito`, `si+Ar/si+ar.ovito` は表示設定済みのテンプレートです。ダブルクリックまたは `File > Open...` で開けば同じ可視化パイプラインを再現できます。
- 別の dump に差し替える場合は、`Pipeline > File source` を新しいファイルに変更してください。

## よくあるトラブルと対処

- パスにスペースが含まれていて失敗する: `-in ".../Quantum-diamond/..."` のように必ずクォートするか、先に対象フォルダへ `cd` してください。
- コマンドが見つからない: 環境での LAMMPS コマンド名（`lmp_mpi` 等）を使用。PATH が通っていない場合はフルパスで指定。
- ポテンシャルファイルが見つからない: 入力スクリプトが相対パスで参照するため、該当ディレクトリで実行するか、入力内のパスを修正。
- 並列数の調整: `-np` をマシンのコア数/キューに合わせて調整。

## メモ

- 本リポジトリの一部フォルダ（`si+Ar/`, `NV/`, `irradiation/`, `Si2/`）は主に結果整理用で、直接実行用の `in.*` がない場合があります。
- 乱数シードやバージョン差により微小差が出ることがあります。再現のために LAMMPS のバージョンを記録しておくと良いです。
