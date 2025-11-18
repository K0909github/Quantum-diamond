import sys


def find_nitrogen_positions_with_line_numbers(data_file_path):
    """
    LAMMPS data ファイルを読み込み、type 3 の原子をすべて見つけて
    各原子の (x, y, z) と元ファイルの行番号（1始まり）を出力する。

    汎用性のため、行の中から3つの連続した浮動小数点値を x,y,z として抽出する。
    """
    in_atoms_section = False
    found = []  # (line_number, x, y, z)

    try:
        with open(data_file_path, 'r') as f:
            for idx, raw_line in enumerate(f, start=1):
                line = raw_line.strip()

                # セクション開始判定
                if line.startswith('Atoms'):
                    in_atoms_section = True
                    continue  # 次行から原子データの可能性がある

                if not line:
                    continue

                # セクション見出しが出てきたら Atoms セクションを抜ける
                # たとえば "Velocities" や "Bonds" など
                if in_atoms_section and line[0].isalpha() and not line[0].isdigit():
                    # ただし、コメント行やヘッダ行の可能性もあるため
                    # 最初のトークンが整数に変換できるかを試す
                    tokens = line.split()
                    try:
                        int(tokens[0])
                    except Exception:
                        # 整数でない先頭トークンはセクションヘッダとみなす
                        in_atoms_section = False
                        continue

                if in_atoms_section:
                    parts = line.split()
                    if len(parts) < 3:
                        continue

                    # まず行先頭が atom id であるか試す
                    try:
                        _ = int(parts[0])
                    except ValueError:
                        # 先頭が整数でなければおそらくヘッダ
                        continue

                    # atom_type を推定: 通常は parts[1] が type
                    atom_type = None
                    for type_idx in (1, 2):
                        if type_idx < len(parts):
                            try:
                                atom_type = int(parts[type_idx])
                                type_token_idx = type_idx
                                break
                            except Exception:
                                continue

                    if atom_type is None:
                        continue

                    if atom_type == 3:
                        # type トークンの後から3つ連続する float を探して x,y,z とする
                        coords = None
                        for i in range(type_token_idx + 1, len(parts) - 1):
                            try:
                                x = float(parts[i])
                                y = float(parts[i+1])
                                z = float(parts[i+2])
                                coords = (x, y, z)
                                break
                            except Exception:
                                continue

                        if coords is None:
                            # 最後の手段: 行中で見つかる最初の3つの float を取る
                            floats = []
                            for tok in parts:
                                try:
                                    fval = float(tok)
                                    floats.append(fval)
                                except Exception:
                                    continue
                                if len(floats) >= 3:
                                    break
                            if len(floats) >= 3:
                                coords = (floats[0], floats[1], floats[2])

                        if coords is not None:
                            found.append((idx, coords[0], coords[1], coords[2]))

        if found:
            print(f"ファイル '{data_file_path}' 内の Type 3 (窒素) 原子一覧:")
            for line_no, x, y, z in found:
                print(f"行 {line_no}: x={x}, y={y}, z={z}")
        else:
            print(f"ファイル '{data_file_path}' 内に Type 3 の原子が見つかりませんでした。")

    except FileNotFoundError:
        print(f"エラー: ファイル '{data_file_path}' が見つかりません。", file=sys.stderr)
    except Exception as e:
        print(f"ファイルの読み取り中にエラーが発生しました: {e}", file=sys.stderr)


if __name__ == '__main__':
    # コマンドライン引数でファイル名を受け取る。指定がなければデフォルト名を使用。
    default_file = '1atom_1.6keV_N_to_C.data'
    file_to_read = sys.argv[1] if len(sys.argv) > 1 else default_file
    find_nitrogen_positions_with_line_numbers(file_to_read)