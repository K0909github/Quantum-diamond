import matplotlib.pyplot as plt
import numpy as np

# --- 設定 ---
log_filename = 'log.lammps'
output_filename = 'final_temp_vs_energy_graph_1600K.png'
# --- 設定ここまで ---

print(f"'{log_filename}' を読み込んでいます...")

temps, pot_engs = [], []
is_data_section = False

try:
    with open(log_filename, 'r') as f:
        for line in f:
            if line.strip().startswith('Step'):
                is_data_section = True
                # ヘッダーから 'Temp' と 'PotEng' が何番目の列かを探す
                header = line.strip().split()
                try:
                    temp_col_index = header.index('Temp')
                    pe_col_index = header.index('PotEng')
                except ValueError:
                    is_data_section = False # 探している列がなければリセット
                continue
            
            if line.strip().startswith('Loop time'):
                is_data_section = False
                break

            if is_data_section:
                try:
                    parts = line.split()
                    pe_value = float(parts[pe_col_index])
                    temp_value = float(parts[temp_col_index])
                    pot_engs.append(pe_value)
                    temps.append(temp_value)
                except (ValueError, IndexError):
                    continue

except FileNotFoundError:
    print(f"エラー: '{log_filename}' が見つかりません。")
    print("このスクリプトと同じフォルダに log.lamps ファイルを置いてください。")
    exit()

if not temps:
    print("エラー: ログファイルからデータを読み込めませんでした。")
    exit()

print(f"データ読み込み完了({len(temps)}点)。グラフを作成します...")

# グラフを作成 (散布図を使用)
plt.figure(figsize=(10, 7))
plt.scatter(temps, pot_engs, s=5, label='Data', alpha=0.5)

# グラフの装飾
plt.title('Potential Energy vs. Temperature during Si Melting', fontsize=16)
plt.xlabel('Temperature (K)', fontsize=14)
plt.ylabel('Potential Energy (eV)', fontsize=14)
plt.grid(True, which='both', linestyle='--', linewidth=0.5)

# ★★★ 赤線の位置を 1750 K に修正 ★★★
melting_temp = 1750
plt.axvline(x=melting_temp, color='r', linestyle='--', linewidth=1.5, label=f'Melting Point (~{melting_temp} K)')

plt.legend(fontsize=12)
plt.tight_layout()

plt.savefig(output_filename)
print(f"成功！グラフを '{output_filename}' として保存しました。")