#!/usr/bin/env python3
N_ids = ids[N_mask]
N_pos = pos[N_mask]


# 近接数（Nの近傍C原子数）をカウント
C_mask = (ty == 1)
C_pos = pos[C_mask]


# 距離計算（単純O(N^2)でOK: 数万以下想定）
NV_idx = []
for i, pN in enumerate(N_pos):
dr = C_pos - pN
dr = delta_pbc(dr)
dist = np.linalg.norm(dr, axis=1)
nC = int(np.sum(dist < RCUT_CN))
if nC == 3:
NV_idx.append(i)


NV_pos = N_pos[NV_idx]
print(f"Detected NV-candidate N atoms: {len(NV_pos)}")


# NV–NV距離の全ペアを計算
pairs = []
M = len(NV_pos)
for i in range(M):
for j in range(i+1, M):
dr = NV_pos[j] - NV_pos[i]
dr = delta_pbc(dr.reshape(1,3)).ravel()
d = float(np.linalg.norm(dr))
pairs.append((i, j, d))


pairs = np.array(pairs, dtype=float)
if pairs.size == 0:
print("No NV–NV pairs to analyze.")
# それでもNV候補数だけは報告
with open('nv_summary.txt','w') as fw:
fw.write(f"NV-candidate count: {M}\n")
raise SystemExit


# 統計
d = pairs[:,2]
stats = {
'n_nv': M,
'n_pairs': len(d),
'mean': float(np.mean(d)),
'median': float(np.median(d)),
'p10': float(np.percentile(d,10)),
'p25': float(np.percentile(d,25)),
'p75': float(np.percentile(d,75)),
'p90': float(np.percentile(d,90)),
}


# 出力
np.savetxt('nv_pairs.csv', pairs, delimiter=',', header='i,j,dist_A', comments='')
with open('nv_summary.txt','w') as fw:
for k,v in stats.items():
fw.write(f"{k}: {v}\n")


print("Saved nv_pairs.csv and nv_summary.txt")