import matplotlib.pyplot as plt
import re

def plot_lammps_log(log_filename, output_filename):
    """
    LAMMPSのログファイルを解析し、温度とポテンシャルエネルギーの関係をプロットして保存する関数。

    Args:
        log_filename (str): LAMMPSのログファイル名。
        output_filename (str): 保存するグラフのファイル名。
    """
    temperatures = []
    potential_energies = []

    data_regex = re.compile(r"^\s*\d+\s+.*")
    start_reading = False
    
    # 【修正点1】文字コードを 'utf-8' に指定し、エラーを無視するオプションを追加
    with open(log_filename, 'r', encoding='utf-8', errors='ignore') as f:
        for line in f:
            if "Step          Time" in line:
                start_reading = True
                continue
            
            if "Loop time" in line:
                start_reading = False
                break
            
            if start_reading and data_regex.match(line):
                try:
                    parts = line.split()
                    temp = float(parts[5])
                    pot_eng = float(parts[3])
                    temperatures.append(temp)
                    potential_energies.append(pot_eng)
                except (ValueError, IndexError):
                    continue

    plt.style.use('default')
    plt.figure(figsize=(12, 8))

    plt.scatter(
        temperatures,
        potential_energies,
        s=10,
        alpha=0.6,
        label='Data'
    )

    melting_point_temp = 1750
    plt.axvline(
        x=melting_point_temp,
        color='red',
        linestyle='--',
        linewidth=2,
        label=f'Melting Point (~{melting_point_temp} K)'
    )

    plt.title('Potential Energy vs. Temperature during Si Melting', fontsize=16)
    plt.xlabel('Temperature (K)', fontsize=14)
    plt.ylabel('Potential Energy (eV)', fontsize=14)
    plt.legend(fontsize=12)
    plt.grid(True, which='both', linestyle='--', linewidth=0.5)
    
    plt.xticks(fontsize=12)
    plt.yticks(fontsize=12)

    plt.tight_layout()

    # 【修正点2】グラフを画面に表示する代わりにファイルに保存
    plt.savefig(output_filename)
    print(f"グラフを {output_filename} として保存しました。")


# --- メインの実行部分 ---
if __name__ == "__main__":
    log_file = 'log.lammps'
    output_file = 'final_temp_vs_energy_graph_1600K.png'
    plot_lammps_log(log_file, output_file)