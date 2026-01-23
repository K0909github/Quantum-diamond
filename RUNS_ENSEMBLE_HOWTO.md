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

### WSL/Linux(bash) で実行する場合（重要）

- bash では **改行継続にバッククォート ` は使いません**（PowerShell 用です）。
- bash の改行継続は **バックスラッシュ `\`** です。
- 迷ったら、まずは「1行」で実行してください（スペース入りパスは必ずダブルクォート）。

1行版（bash）:

```bash
python3 tools/run_ensemble_implantation.py --template-dir "10Ncluster_implantation to C_5keV" --input "10atoms_5keV_N_implantation_to_C_ZBL_potential_filedata.txt" --out "10Ncluster_implantation to C_5keV/runs" --runs 10 --seed 12345 --x-range -20 20 --y-range -20 20
```

複数行版（bash）:

```bash
python3 tools/run_ensemble_implantation.py \
  --template-dir "10Ncluster_implantation to C_5keV" \
  --input "10atoms_5keV_N_implantation_to_C_ZBL_potential_filedata.txt" \
  --out "10Ncluster_implantation to C_5keV/runs" \
  --runs 10 \
  --seed 12345 \
  --x-range -20 20 \
  --y-range -20 20
```

### 例: C_7keV(final) に runs を作って、そのまま10回自動実行（WSL/bash）

```bash
python3 tools/run_ensemble_implantation.py --template-dir "10Ncluster_implantation to C_7keV(final)" --input "10atoms_7keV_N_implantation_to_C_ZBL_potential_filedata.txt" --out "10Ncluster_implantation to C_7keV(final)/runs" --runs 10 --seed 12345 --x-range -20 20 --y-range -20 20 --lammps "mpirun -np 8 lmp -in in.lmp"
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

## 1b) 連続注入（noreset, 10回注入）をさらに ensemble で10回回す

例: `N_implantation_to_C_random_5keV_noreset_10_10(final)/1atoms_5keV_N_implantation_to_C_ZBL_potential_filedata.txt`

このタイプのテンプレは **LAMMPS入力の中で `variable m loop 10` のように「10回の連続注入」を自前で回す**作りになっています。
なので「10回連続注入」を **さらに** ensemble で10回（= 10注入×10run）やるには、
テンプレを `run_01 .. run_10` に複製して **各runで乱数の基準seed（例: 12345/67890）だけを変えて** LAMMPS を実行するのが簡単です。

ポイント:

- このテンプレは `x_pos/y_pos` を入力内で `random(...)` で決めます（Python側で `x-range/y-range` を渡す方式ではありません）。
- `12345+${m}` / `67890+${m}` の「12345/67890」を run ごとに変えると、10回注入の乱数系列が run ごとに変わります。

### まずはこれ（おすすめ）: run_ensemble_implantation.py で一発

`--template-style loop-random-xy` を使うと、連続注入テンプレ内の

- `rnd_seed` / `rnd_seed_y` の基準値
- `x_pos` / `y_pos` の `random(xmin,xmax,...)` の範囲

を run ごとに差し替えた `run_01..run_10` を作り、そのまま LAMMPS まで実行できます。

```bash
python3 tools/run_ensemble_implantation.py \
  --template-style loop-random-xy \
  --template-dir "N_implantation_to_C_random_5keV_noreset_10_10(final)" \
  --input "1atoms_5keV_N_implantation_to_C_ZBL_potential_filedata.txt" \
  --out "N_implantation_to_C_random_5keV_noreset_10_10(final)/runs_ensemble" \
  --runs 10 \
  --seed 12345 \
  --x-range -4.6 4.6 \
  --y-range -4.6 4.6 \
  --lammps "mpirun -np 8 lmp -in in.lmp"
```

### PowerShell 例（run_01..run_10 を作って順に実行）

```powershell
$templateDir = "N_implantation_to_C_random_5keV_noreset_10_10(final)"
$templateIn  = "1atoms_5keV_N_implantation_to_C_ZBL_potential_filedata.txt"
$outRoot     = "N_implantation_to_C_random_5keV_noreset_10_10(final)/runs_ensemble"

New-Item -ItemType Directory -Force $outRoot | Out-Null

for ($i=1; $i -le 10; $i++) {
  $runDir = Join-Path $outRoot ("run_{0:D2}" -f $i)
  New-Item -ItemType Directory -Force $runDir | Out-Null

  # テンプレフォルダから必要ファイルをコピー（状況に応じて追加）
  Copy-Item -Force (Join-Path $templateDir "*.data")   $runDir -ErrorAction SilentlyContinue
  Copy-Item -Force (Join-Path $templateDir "*.zbl")    $runDir -ErrorAction SilentlyContinue
  Copy-Item -Force (Join-Path $templateDir "*.tersoff*") $runDir -ErrorAction SilentlyContinue

  # 入力ファイルを in.lmp として複製し、seed基準を run ごとに変更
  $text = Get-Content (Join-Path $templateDir $templateIn) -Raw
  $seedBaseX = 12345 + 1000*$i
  $seedBaseY = 67890 + 1000*$i
  $text = $text.Replace("variable        rnd_seed equal 12345+${m}",  "variable        rnd_seed equal $seedBaseX+${m}")
  $text = $text.Replace("variable        rnd_seed_y equal 67890+${m}","variable        rnd_seed_y equal $seedBaseY+${m}")
  Set-Content -Path (Join-Path $runDir "in.lmp") -Value $text -Encoding utf8

  # 実行（LAMMPSコマンドは環境に合わせて）
  Push-Location $runDir
  # 例: mpirun / lmp は環境で調整
  mpirun -np 8 lmp -in in.lmp
  Pop-Location
}
```

### WSL/Linux(bash) 例（run_01..run_10 を作って順に実行）

```bash
template_dir="N_implantation_to_C_random_5keV_noreset_10_10(final)"
template_in="1atoms_5keV_N_implantation_to_C_ZBL_potential_filedata.txt"
out_root="N_implantation_to_C_random_5keV_noreset_10_10(final)/runs_ensemble"

mkdir -p "$out_root"

for i in $(seq -w 1 10); do
  run_dir="$out_root/run_$i"
  mkdir -p "$run_dir"

  cp -f "$template_dir"/*.data "$run_dir" 2>/dev/null || true
  cp -f "$template_dir"/*.zbl "$run_dir" 2>/dev/null || true
  cp -f "$template_dir"/*.tersoff* "$run_dir" 2>/dev/null || true

  seed_base_x=$((12345 + 1000*10#$i))
  seed_base_y=$((67890 + 1000*10#$i))

  sed \
    -e "s/variable        rnd_seed equal 12345+\${m}/variable        rnd_seed equal ${seed_base_x}+\${m}/" \
    -e "s/variable        rnd_seed_y equal 67890+\${m}/variable        rnd_seed_y equal ${seed_base_y}+\${m}/" \
    "$template_dir/$template_in" > "$run_dir/in.lmp"

  (cd "$run_dir" && mpirun -np 8 lmp -in in.lmp)
done
```

調整したい場合:

- 位置の範囲はテンプレ内の `random(-4.6,4.6,...)` を基板サイズに合わせて調整してください。
- 各runが `dump_all.xyz` や `final_data_run_*.data` を出します（runフォルダ分けしているので上書きはしません）。

## 2) vacancy の深さ分布（10回分をまとめて）

`run_*` をフォルダ指定で渡せます:

```powershell
python "10Ncluster_implantation to C_0.5keV(final)/analysis3_graph_vacancy_list.py" `
  "10Ncluster_implantation to C_5keV/runs/run_*" `
  --out "vacancy_depths_ensemble.png" `
  --surface-z 125 `
  --bin-width 1
```

Ubuntu/WSL(bash) 例（C_7keV(final) の runs から出す）:

```bash
python3 "10Ncluster_implantation to C_0.5keV(final)/analysis3_graph_vacancy_list.py" \
  "10Ncluster_implantation to C_7keV(final)/runs/run_*" \
  --out "10Ncluster_implantation to C_7keV(final)/runs/vacancy_depths_ensemble.png" \
  --surface-z 125 \
  --bin-width 1
```

## 3) N の深さ分布（10回分をまとめて）

注意:

- このスクリプトは **引数なし**で実行すると「スクリプトと同じフォルダ」にある `N_list` を探します。
- runs 配下のデータで分布を作る場合は、必ず `run_*/...` の **入力パターンを引数として渡してください**。

`dump_run_1_min0K.lammpstrj` を直接読む例（N が type=3 の場合）:

```powershell
python "10Ncluster_implantation to C_0.5keV(final)/analysis2_graph_N_list3.py" `
  "10Ncluster_implantation to C_5keV/runs/run_*/dump_run_1_min0K.lammpstrj" `
  --atom-type 3 `
  --out "nitrogen_depths_ensemble.png" `
  --surface-z 125 `
  --bin-width 5 `
  --max-depth 250
```

Ubuntu/WSL(bash) 例（C_7keV(final) の runs から出す）:

```bash
python3 "10Ncluster_implantation to C_0.5keV(final)/analysis2_graph_N_list3.py" \
  "10Ncluster_implantation to C_7keV(final)/runs/run_*/dump_run_1_min0K.lammpstrj" \
  --atom-type 3 \
  --out "10Ncluster_implantation to C_7keV(final)/runs/nitrogen_depths_ensemble.png" \
  --surface-z 125 \
  --bin-width 5 \
  --max-depth 250
```

補足:

- パスにスペースがあるので、`run_*/...` の引数はダブルクォートで囲んでOKです。
  その場合でもスクリプト側でワイルドカードを展開して読み込みます。
- `--max-depth 250` を付けると、基板厚（例: 250Å）より深い値（底から抜けた原子など）を描画から除外できます。

事前確認（ファイルが本当にあるか）:

```bash
ls "10Ncluster_implantation to C_7keV(final)/runs/run_01"
ls "10Ncluster_implantation to C_7keV(final)/runs/run_01/dump_run_1_min0K.lammpstrj"
```

## よくある調整ポイント

- `--x-range/--y-range` は「基板の有効領域（表面の中央付近）」に合わせて調整してください。
- `--surface-z` は基板モデルの表面 z に合わせてください（いまは既存スクリプトと同じ `125 Å` をデフォルトにしています）。

## エラー対処

### FileNotFoundError: テンプレ入力が見つかりません

原因はほぼ次のどれかです:

- `--template-dir` のフォルダ名が実際と違う（例: `...final` が付いている等）
- `--input` のファイル名が実際と違う（keV値や `_filedata` の有無など）
- 実行しているカレントディレクトリがリポジトリルートではない

bash/WSL なら、まず次で「存在する名前」を確認してコピペしてください:

```bash
pwd
ls
ls "10Ncluster_implantation to C_5keV"
```

その上で、`--template-dir` と `--input` を **実在する名前** に合わせます。

### 分布（グラフ出力）でエラーが出る

まず、どのエラーでも共通で次を確認してください:

```bash
# リポジトリルートにいるか確認
pwd
ls

# runフォルダに必要ファイルがあるか確認（例: run_01）
ls "10Ncluster_implantation to C_7keV(final)/runs/run_01"
```

よくある原因:

- `ModuleNotFoundError: No module named 'matplotlib'`
  - 対処: `python3 -m pip install --user matplotlib`
- `...入力が見つかりません`（vacancy_list / dump が無い）
  - N分布: 各runに `dump_run_1_min0K.lammpstrj` があるか確認
  - vacancy分布: 各runに `vacancy_list`（または vacancy_list.xyz 等）が必要です（OVITO等でrunごとに出力）
- パスにスペースがある
  - 対処: 文字列は必ずダブルクォートで囲む（例: `"10Ncluster_implantation to C_7keV(final)/..."`）

こちらで原因を特定するため、エラーが出たときは「実行したコマンド」と「Traceback全文」をそのまま貼ってください。
