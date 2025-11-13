import sys

def find_nitrogen_z_coordinate(data_file_path):
    """
    LAMMPS data ファイルを読み込み、type 3 の原子の Z 座標を探します。
    """
    in_atoms_section = False
    z_coordinate = None

    try:
        with open(data_file_path, 'r') as f:
            for line in f:
                line = line.strip()

                # "Atoms" という行を見つけたら、そこからが原子データセクション
                if line.startswith('Atoms'):
                    in_atoms_section = True
                    continue  # "Atoms"の次の行から処理開始
                
                if not line:  # 空行はスキップ
                    continue

                if in_atoms_section:
                    # 行をスペースで分割
                    parts = line.split()
                    
                    # 有効なデータ行か確認（最低でもID, type, x, y, z があるはず）
                    # また、コメント行（#）は除外
                    if len(parts) >= 5 and not parts[0].startswith('#'):
                        try:
                            # 2番目の要素（インデックス 1）が atom-type
                            atom_type = int(parts[1])
                            
                            # atom-type が 3 (窒素) なら
                            if atom_type == 3:
                                # 5番目の要素（インデックス 4）が Z 座標
                                z_coordinate = float(parts[4])
                                # 窒素原子は1つだけのはずなので、見つかったらループを抜ける
                                break
                        except (ValueError, IndexError):
                            # ヘッダー後の空行など、数値に変換できない行は無視
                            continue
                            
        if z_coordinate is not None:
            print(f"ファイル '{data_file_path}' 内の Type 3 (窒素) 原子の最終 Z 座標:")
            print(z_coordinate)
        else:
            print(f"ファイル '{data_file_path}' 内に Type 3 の原子が見つかりませんでした。")

    except FileNotFoundError:
        print(f"エラー: ファイル '{data_file_path}' が見つかりません。", file=sys.stderr)
    except Exception as e:
        print(f"ファイルの読み取り中にエラーが発生しました: {e}", file=sys.stderr)

# --- 実行 ---
# 読み込みたいファイル名を指定してください
file_to_read = '1atom_1.6keV_N_to_C.data' 

find_nitrogen_z_coordinate(file_to_read)