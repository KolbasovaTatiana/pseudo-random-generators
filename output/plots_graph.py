import os
import pandas as pd
import matplotlib.pyplot as plt

stats_dir = "output"
speed_file = os.path.join(stats_dir, "speed_results.csv")

if not os.path.exists(speed_file):
    raise FileNotFoundError(f"Не найден файл {speed_file}")

df = pd.read_csv(speed_file)

plt.style.use("seaborn-v0_8-whitegrid")
plt.figure(figsize=(8, 5))

for col in df.columns[1:]:
    label = col.replace("_ms", "")
    plt.plot(df["size"], df[col], marker="o", linewidth=2, label=label)

plt.title("Сравнение скоростей генерации")
plt.xlabel("Размер выборки")
plt.ylabel("Время, мс")
plt.xscale("log")
plt.legend()
plt.tight_layout()

os.makedirs("plots", exist_ok=True)
plt.savefig(os.path.join("plots", "speed_comparison.png"), dpi=300)
plt.close()

print("График скорости сохранён в plots/speed_comparison.png")