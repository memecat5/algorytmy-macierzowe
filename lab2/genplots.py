import pandas as pd
import matplotlib.pyplot as plt

# Use a clean style
plt.style.use("seaborn-v0_8-whitegrid")

# === Multiplication comparison ===
def plot_multiplication_comparison(csv_file="multiplication_comparison.csv"):
    try:
        df = pd.read_csv(csv_file)
    except FileNotFoundError:
        print(f"File {csv_file} not found! Not making plot")
        return
    n = df["n"]

    plt.figure(figsize=(12, 8))

    # --- Operations ---
    plt.subplot(2, 1, 1)
    plt.plot(n, df["classic_ops"], "o-b", label="Classic O(n³)")
    plt.plot(n, df["recursive_ops"], "s-g", label="Recursive O(n³)")
    plt.plot(n, df["strassen_ops"], "^r", label="Strassen O(n²․⁸⁰⁷)")
    plt.xlabel("Matrix size (n)")
    plt.ylabel("Floating-point operations")
    plt.title("Matrix Multiplication: Operation Count Comparison")
    plt.legend()
    plt.grid(True)

    # --- Time ---
    plt.subplot(2, 1, 2)
    plt.plot(n, df["classic_time_ms"], "o-b", label="Classic")
    plt.plot(n, df["recursive_time_ms"], "s-g", label="Recursive")
    plt.plot(n, df["strassen_time_ms"], "^r", label="Strassen")
    plt.xlabel("Matrix size (n)")
    plt.ylabel("Execution time (ms)")
    plt.title("Matrix Multiplication: Execution Time")
    plt.legend()
    plt.grid(True)

    plt.tight_layout()
    plt.savefig("multiplication_comparison_plot.png", dpi=200)
    plt.show()


# === All operations comparison ===
def plot_all_operations(csv_file="all_operations_comparison.csv"):
    try:
        df = pd.read_csv(csv_file)
    except FileNotFoundError:
        print(f"File {csv_file} not found! Not making plot")
        return
    n = df["n"]

    plt.figure(figsize=(10, 6))
    # plt.plot(n, df["multiply_ops"], "o-b", label="Multiplication (Strassen)")
    plt.plot(n, df["invert_ops"], "s-g", label="Inversion")
    plt.plot(n, df["gauss_ops"], "^r", label="Gaussian Elimination")
    plt.plot(n, df["lu_ops"], "d-m", label="LU + Determinant")
    plt.xlabel("Matrix size (n)")
    plt.ylabel("Total floating-point operations")
    plt.title("Comparison of All Matrix Operations")
    plt.legend()
    plt.grid(True)
    plt.tight_layout()
    plt.savefig("all_operations_comparison_plot.png", dpi=200)
    plt.show()


# === Detailed multiplication breakdown ===
def plot_detailed_multiplication(csv_file="detailed_multiplication_analysis.csv"):
    try:
        df = pd.read_csv(csv_file)
    except FileNotFoundError:
        print(f"File {csv_file} not found! Not making plot")
        return
    n = df["n"]

    plt.figure(figsize=(14, 8))

    # Additions
    plt.subplot(3, 1, 1)
    plt.plot(n, df["classic_add"], "o-b", label="Classic")
    plt.plot(n, df["recursive_add"], "s-g", label="Recursive")
    plt.plot(n, df["strassen_add"], "^r", label="Strassen")
    plt.ylabel("Additions/Subtractions")
    plt.title("Addition Operations per Algorithm")
    plt.legend()
    plt.grid(True)

    # Multiplications
    plt.subplot(3, 1, 2)
    plt.plot(n, df["classic_mul"], "o-b", label="Classic")
    plt.plot(n, df["recursive_mul"], "s-g", label="Recursive")
    plt.plot(n, df["strassen_mul"], "^r", label="Strassen")
    plt.ylabel("Multiplications")
    plt.title("Multiplication Operations per Algorithm")
    plt.legend()
    plt.grid(True)

    # Divisions
    plt.subplot(3, 1, 3)
    plt.plot(n, df["classic_div"], "o-b", label="Classic")
    plt.plot(n, df["recursive_div"], "s-g", label="Recursive")
    plt.plot(n, df["strassen_div"], "^r", label="Strassen")
    plt.xlabel("Matrix size (n)")
    plt.ylabel("Divisions")
    plt.title("Division Operations per Algorithm")
    plt.legend()
    plt.grid(True)

    plt.tight_layout()
    plt.savefig("detailed_multiplication_analysis_plot.png", dpi=200)
    plt.show()

def plot_benchmarks(csv_files = ["inversion_benchmark.csv", "gauss_benchmark.csv", "lu_benchmark.csv"]):
    # ======== INVERSION ========
    try:
        df_inv = pd.read_csv(csv_files[0])
        df_gauss = pd.read_csv(csv_files[1])
        df_lu = pd.read_csv(csv_files[2])
    except FileNotFoundError:
        print(f"Some of the files {csv_files} not found! Not making plot")
        return
    n = df_inv["n"]

    plt.figure(figsize=(14, 8)) # time_ms,FLOPs,peak_memory_MB,FLOPs_per_sec

    # plt.subplot(4, 1, 1)
    plt.plot(n, df_inv["time_ms"], label="Inversion")
    plt.plot(n, df_gauss["time_ms"], label="Gaussian Elimination")
    plt.plot(n, df_lu["time_ms"], label="LU factorisation")
    plt.xlabel("Matrix size (n)")
    plt.ylabel("Time [ms]")
    plt.title("Runtime per algorithm")
    plt.legend()
    plt.grid(visible=True)
    plt.savefig("benchark_time.png")

    plt.figure(figsize=(14, 8))
    #plt.subplot(4, 1, 2)
    plt.plot(n, df_inv["FLOPs"], label="Inversion")
    plt.plot(n, df_gauss["FLOPs"], label="Gaussian Elimination")
    plt.plot(n, df_lu["FLOPs"], label="LU factorisation")
    plt.xlabel("Matrix size (n)")
    plt.ylabel("Floating point operations")
    plt.title("Total floating point operations per algorithm")
    plt.legend()
    plt.grid(visible=True)
    plt.savefig("benchark_FLOPs.png")

    plt.figure(figsize=(14, 8))
    # plt.subplot(4, 1, 3)
    plt.plot(n, df_inv["peak_memory_MB"], label="Inversion")
    plt.plot(n, df_gauss["peak_memory_MB"], label="Gaussian Elimination")
    plt.plot(n, df_lu["peak_memory_MB"], label="LU factorisation")
    plt.xlabel("Matrix size (n)")
    plt.ylabel("Memory usage")
    plt.title("Memory usage per algorithm")
    plt.legend()
    plt.grid(visible=True)
    plt.savefig("benchark_memory.png")

    plt.figure(figsize=(14, 8))
    #plt.subplot(4, 1, 4)
    plt.plot(n, df_inv["FLOPs_per_sec"], label="Inversion")
    plt.plot(n, df_gauss["FLOPs_per_sec"], label="Gaussian Elimination")
    plt.plot(n, df_lu["FLOPs_per_sec"], label="LU factorisation")
    plt.xlabel("Matrix size (n)")
    plt.ylabel("FLOPS")
    plt.title("FLOPS per algorithm")
    plt.legend()
    plt.grid(visible=True)
    plt.savefig("benchark_FLOPsPERs.png")

    #plt.tight_layout()
    #plt.savefig("benchmarks_plot.png", dpi=200)
    #plt.show()


# === Run all ===
if __name__ == "__main__":
    print("Generating all plots...")
    # plot_multiplication_comparison()
    plot_all_operations()
    plot_benchmarks()
    # plot_detailed_multiplication()
    print("Plots saved as PNG files.")