import pandas as pd
import matplotlib.pyplot as plt

def plot_inversion_benchmark(csv_file="inversion_benchmark.csv"):
    # Load CSV
    df = pd.read_csv(csv_file)

    n = df["n"]

    # Data columns
    rec_time = df["rec_time_ms"]
    gauss_time = df["gauss_time_ms"]
    lu_time = df["lu_time_ms"]

    rec_flops = df["rec_flops"]
    gauss_flops = df["gauss_flops"]
    lu_flops = df["lu_flops"]

    rec_mem = df["rec_mem_MB"]
    gauss_mem = df["gauss_mem_MB"]
    lu_mem = df["lu_mem_MB"]

    rec_fps = df["rec_flops_per_sec"]
    gauss_fps = df["gauss_flops_per_sec"]
    lu_fps = df["lu_flops_per_sec"]

    # -------------------------------------------------------------------------
    # 1. Time plot
    # -------------------------------------------------------------------------
    plt.figure(figsize=(10, 6))
    plt.plot(n, rec_time, label="Rekurencyjna")
    plt.plot(n, gauss_time, label="Gauss")
    plt.plot(n, lu_time, label="LU")
    plt.xlabel("Rozmiar macierzy n")
    plt.ylabel("Czas (ms)")
    plt.title("Czas odwracania macierzy")
    plt.grid(True)
    plt.legend()
    plt.savefig("inversion_time.png", dpi=200)
    plt.close()

    # -------------------------------------------------------------------------
    # 2. FLOPs plot
    # -------------------------------------------------------------------------
    plt.figure(figsize=(10, 6))
    plt.plot(n, rec_flops, label="Rekurencyjna")
    plt.plot(n, gauss_flops, label="Gauss")
    plt.plot(n, lu_flops, label="LU")
    plt.xlabel("Rozmiar macierzy n")
    plt.ylabel("Liczba FLOPs")
    plt.title("Liczba operacji zmiennoprzecinkowych (FLOPs)")
    plt.grid(True)
    plt.legend()
    plt.savefig("inversion_flops.png", dpi=200)
    plt.close()

    # -------------------------------------------------------------------------
    # 3. Memory usage plot
    # -------------------------------------------------------------------------
    plt.figure(figsize=(10, 6))
    plt.plot(n, rec_mem, label="Rekurencyjna")
    plt.plot(n, gauss_mem, label="Gauss")
    plt.plot(n, lu_mem, label="LU")
    plt.xlabel("Rozmiar macierzy n")
    plt.ylabel("Pamięć (MB)")
    plt.title("Zapotrzebowanie pamięci")
    plt.grid(True)
    plt.legend()
    plt.savefig("inversion_memory.png", dpi=200)
    plt.close()

    # -------------------------------------------------------------------------
    # 4. FLOPs per second
    # -------------------------------------------------------------------------
    plt.figure(figsize=(10, 6))
    plt.plot(n, rec_fps, label="Rekurencyjna")
    plt.plot(n, gauss_fps, label="Gauss")
    plt.plot(n, lu_fps, label="LU")
    plt.xlabel("Rozmiar macierzy n")
    plt.ylabel("FLOPs / sekundę")
    plt.title("Wydajność obliczeniowa (FLOPs/s)")
    plt.grid(True)
    plt.legend()
    plt.savefig("inversion_flops_per_sec.png", dpi=200)
    plt.close()

    print("Wygenerowane wykresy:")
    print("  inversion_time.png")
    print("  inversion_flops.png")
    print("  inversion_memory.png")
    print("  inversion_flops_per_sec.png")


if __name__ == "__main__":
    plot_inversion_benchmark()
