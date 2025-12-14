import pandas as pd
import matplotlib.pyplot as plt
import numpy as np

def plot_results():
    # Wczytanie danych
    try:
        df = pd.read_csv('results.csv')
    except FileNotFoundError:
        print("Nie znaleziono pliku results.csv. Uruchom najpierw program w C++.")
        return

    # Ustawienie stylu wykresów
    plt.style.use('seaborn-v0_8-whitegrid')

    # --- Wykres 1: Czas wykonania (Time) ---
    plt.figure(figsize=(8, 6))
    plt.plot(df['N'], df['Time_Compression'], 'o-', label='Kompresja (Compression)', linewidth=2)
    plt.plot(df['N'], df['Time_MVM'], 's-', label='Mnożenie MVM (MVM Time)', linewidth=2)
    
    plt.xscale('log')
    plt.yscale('log')
    plt.xlabel('Rozmiar macierzy N (liczba wierszy)', fontsize=12)
    plt.ylabel('Czas [s]', fontsize=12)
    plt.title('Czas wykonania vs Rozmiar Problemu', fontsize=14)
    plt.legend()
    plt.grid(True, which="both", ls="-", alpha=0.4)
    
    plt.tight_layout()
    plt.savefig('time_complexity.png')
    print("Zapisano wykres czasu do: time_complexity.png")
    plt.close() # Zamknięcie figury, aby nie nadpisywać

    # --- Wykres 2: Złożoność Obliczeniowa (FLOPs) ---
    plt.figure(figsize=(8, 6))
    plt.plot(df['N'], df['MVM_FLOPs'], 'o-', color='green', label='MVM FLOPs (Teoretyczne)', linewidth=2)
    plt.plot(df['N'], df['MMM_FLOPs'], 'd-', color='red', label='MMM A*A FLOPs (Teoretyczne)', linewidth=2)

    # Opcjonalne linie referencyjne O(N)
    ref_n = df['N'] * (df['MVM_FLOPs'].iloc[0] / df['N'].iloc[0])
    plt.plot(df['N'], ref_n, '--', color='gray', alpha=0.5, label='O(N) reference')

    plt.xscale('log')
    plt.yscale('log')
    plt.xlabel('Rozmiar macierzy N', fontsize=12)
    plt.ylabel('Liczba operacji (FLOPs)', fontsize=12)
    plt.title('Złożoność Obliczeniowa (FLOPs)', fontsize=14)
    plt.legend()
    plt.grid(True, which="both", ls="-", alpha=0.4)

    plt.tight_layout()
    plt.savefig('flops_complexity.png')
    print("Zapisano wykres FLOPs do: flops_complexity.png")
    plt.close()

if __name__ == "__main__":
    plot_results()