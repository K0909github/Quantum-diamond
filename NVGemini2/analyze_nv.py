import numpy as np
import sys

# --- パラメータ設定 (ご自身の環境に合わせて変更してください) ---

# LAMMPSが出力するdumpファイル名
DUMP_FILE = "final_coords.dump"

# シミュレーションボックスのサイズ [x, y, z] (Å)
# in.nv_injection.lmpの box block 0 40 0 40 0 30 と lattice diamond 3.57 から算出
BOX_DIMS = np.array([40 * 3.57, 40 * 3.57, 30 * 3.57])

# 原子タイプの定義
N_TYPE = 2
C_TYPE = 1

# 最近接原子を数えるための距離のカットオフ (Å)
# ダイヤモンドの結合長(約1.54Å)より少し大きい値
BOND_CUTOFF = 1.8

# NVセンターの窒素原子が持つべき配位数（結合数）
NV_COORD_NUM = 3

# ----------------------------------------------------------------

def read_lammps_dump(filename):
    """LAMMPSのdumpファイルを読み込み、原子データを返す"""
    try:
        with open(filename, 'r') as f:
            lines = f.readlines()
    except FileNotFoundError:
        print(f"エラー: dumpファイル '{filename}' が見つかりません。")
        print("LAMMPSの入力スクリプトにdumpコマンドが正しく設定されているか確認してください。")
        sys.exit(1)

    # ヘッダーから原子数を取得
    num_atoms = int(lines[3])
    # 原子データの開始行を探す ( "ITEM: ATOMS" の次の行)
    data_start_line = lines.index("ITEM: ATOMS id type x y z\n") + 1
    
    atoms = []
    for i in range(num_atoms):
        parts = lines[data_start_line + i].split()
        atoms.append([
            int(parts[0]),    # id
            int(parts[1]),    # type
            float(parts[2]),  # x
            float(parts[3]),  # y
            float(parts[4])   # z
        ])
    return np.array(atoms)

def get_distance(coord1, coord2, box_dims):
    """周期境界条件を考慮して2原子間の距離を計算する"""
    delta = np.abs(coord1 - coord2)
    delta = np.where(delta > 0.5 * box_dims, delta - box_dims, delta)
    return np.sqrt((delta**2).sum(axis=-1))

def main():
    """メインの解析処理"""
    # 1. dumpファイルを読み込み、原子を種類別に分ける
    all_atoms = read_lammps_dump(DUMP_FILE)
    n_atoms = all_atoms[all_atoms[:, 1] == N_TYPE]
    c_atoms = all_atoms[all_atoms[:, 1] == C_TYPE]
    
    if len(n_atoms) == 0:
        print("窒素原子が見つかりませんでした。")
        return

    print(f"窒素原子を {len(n_atoms)} 個検出しました。NVセンターを探索します...")

    # 2. 各窒素原子の配位数を計算し、NVセンターを特定する
    nv_centers_coords = []
    for n_atom in n_atoms:
        n_coord = n_atom[2:5]
        neighbor_count = 0
        for c_atom in c_atoms:
            c_coord = c_atom[2:5]
            distance = get_distance(n_coord, c_coord, BOX_DIMS)
            if distance < BOND_CUTOFF:
                neighbor_count += 1
        
        # 配位数が3の窒素原子をNVセンターとして記録
        if neighbor_count == NV_COORD_NUM:
            nv_centers_coords.append(n_coord)

    num_nv_found = len(nv_centers_coords)
    print(f"結果: {num_nv_found} 個のNVセンターが見つかりました。")

    # 3. NVセンター間の距離を計算して表示
    if num_nv_found > 1:
        print("\n--- NVセンター間の距離 ---")
        nv_centers = np.array(nv_centers_coords)
        for i in range(num_nv_found):
            for j in range(i + 1, num_nv_found):
                dist = get_distance(nv_centers[i], nv_centers[j], BOX_DIMS)
                print(f"NVセンター {i+1} と NVセンター {j+1} の距離: {dist:.2f} Å")

if __name__ == "__main__":
    main()