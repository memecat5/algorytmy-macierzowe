import pandas as pd
import matplotlib.pyplot as plt

def plot_matrix_benchmark(csv_file="matrix_benchmark.csv"):
    df = pd.read_csv(csv_file)

    n = df["n"]
    time_ms = df["time_ms"]
    flops = df["FLOPs"]
    memory = df["memory_MB"]
    flops_per_sec = df["FLOPs_per_sec"]

    # -----------------------------
    # 1. Time vs n
    # -----------------------------
    plt.figure(figsize=(10, 6))
    plt.plot(n, time_ms)
    plt.xlabel("Rozmiar macierzy n")
    plt.ylabel("Czas (ms)")
    plt.title("Czas obliczeń w zależności od rozmiaru macierzy")
    plt.grid(True)
    plt.savefig("benchmark_time.png", dpi=200)
    plt.close()

    # -----------------------------
    # 2. FLOPs vs n
    # -----------------------------
    plt.figure(figsize=(10, 6))
    plt.plot(n, flops)
    plt.xlabel("Rozmiar macierzy n")
    plt.ylabel("Liczba FLOPs")
    plt.title("Liczba operacji zmiennoprzecinkowych")
    plt.grid(True)
    plt.savefig("benchmark_flops.png", dpi=200)
    plt.close()

    # -----------------------------
    # 3. Memory usage vs n
    # -----------------------------
    plt.figure(figsize=(10, 6))
    plt.plot(n, memory)
    plt.xlabel("Rozmiar macierzy n")
    plt.ylabel("Zużycie pamięci (MB)")
    plt.title("Szacowane zużycie pamięci")
    plt.grid(True)
    plt.savefig("benchmark_memory.png", dpi=200)
    plt.close()

    # -----------------------------
    # 4. FLOPs/sec vs n
    # -----------------------------
    plt.figure(figsize=(10, 6))
    plt.plot(n, flops_per_sec)
    plt.xlabel("Rozmiar macierzy n")
    plt.ylabel("FLOPs / sekundę")
    plt.title("Wydajność obliczeniowa (FLOPs/s)")
    plt.grid(True)
    plt.savefig("benchmark_flops_per_sec.png", dpi=200)
    plt.close()

    print("Wygenerowano wykresy:")
    print("  benchmark_time.png")
    print("  benchmark_flops.png")
    print("  benchmark_memory.png")
    print("  benchmark_flops_per_sec.png")


if __name__ == "__main__":
    plot_matrix_benchmark()
