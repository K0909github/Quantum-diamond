import matplotlib.pyplot as plt
import numpy as np

# データを格納するためのリストを準備
distances = []
g_r = []

# out_rdf.txt ファイルを読み込む
try:
    with open("out_rdf.txt", "r") as f:
        for line in f:
            # コメント行はスキップ
            if line.startswith("#"):
                continue
            
            parts = line.split()
            
            # 3列以上のデータ行のみを処理
            if len(parts) >= 3:
                try:
                    # データをリストに追加
                    # parts[1] が距離 (r)
                    # parts[2] が g(r) の値
                    distances.append(float(parts[1]))
                    g_r.append(float(parts[2]))
                except ValueError:
                    # 数値に変換できない行は無視
                    continue

    # Matplotlib を使ってプロット
    plt.figure(figsize=(8, 6))
    plt.plot(distances, g_r, 'o-', markersize=4, label='Si-Si')
    plt.xlabel("Distance (r) [Å]")
    plt.ylabel("g(r)")
    plt.title("Radial Distribution Function (RDF) of Silicon at 2500 K")
    plt.grid(True)
    plt.axhline(y=1, color='gray', linestyle='--') # g(r)=1 の基準線
    plt.legend()

    # グラフを画像ファイルとして保存
    plt.savefig("rdf_plot.png")

    print("グラフが 'rdf_plot.png' として保存されました。")

except FileNotFoundError:
    print("エラー: 'out_rdf.txt' が見つかりません。")
    print("このスクリプトを 'out_rdf.txt' と同じフォルダに置いて実行してください。")