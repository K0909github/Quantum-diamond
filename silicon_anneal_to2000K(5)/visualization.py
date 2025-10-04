import numpy as np
import matplotlib.pyplot as plt

# log.lammpsを読み込む（自分のファイル名に合わせて）
temps, energies = [], []
with open("log.lammps") as f:
    read = False
    for line in f:
        if "Step" in line and "Temp" in line and "TotEng" in line:
            read = True
            continue
        if read:
            if line.strip() == "" or "Loop time" in line:
                read = False
                continue
            parts = line.split()
            temps.append(float(parts[2]))     # Temp (正しい2番目の列)
            energies.append(float(parts[1]))  # TotEng (正しい1番目の列)
            
plt.plot(temps, energies, 'o-', markersize=3)
plt.xlabel("Temperature (K)")
plt.ylabel("Total Energy (eV)")
plt.title("Melting Curve of Silicon")
plt.grid(True)
plt.show()
