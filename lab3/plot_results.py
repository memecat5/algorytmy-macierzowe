import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
import os

def plot_results():
    # Sprawdzenie czy plik istnieje
    if not os.path.exists('results.csv'):
        print("Błąd: Nie znaleziono pliku 'results.csv'.")
        print("Uruchom najpierw program w C++, aby wygenerować dane.")
        return

    # Wczytanie danych
    try:
        df = pd.read_csv('results.csv')
    except Exception as e:
        print(f"Błąd podczas wczytywania CSV: {e}")
        return

    # Ustawienie stylu
    plt.style.use('seaborn-v0_8-whitegrid')

    # ---------------------------------------------------------
    # Wykres 1: Czas wykonania (Time Complexity)
    # ---------------------------------------------------------
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
    print("Zapisano: time_complexity.png")
    plt.close()

    # ---------------------------------------------------------
    # Wykres 2: Złożoność Obliczeniowa (FLOPs)
    # ---------------------------------------------------------
    plt.figure(figsize=(8, 6))
    plt.plot(df['N'], df['MVM_FLOPs'], 'o-', color='green', label='MVM FLOPs (Teoretyczne)', linewidth=2)
    plt.plot(df['N'], df['MMM_FLOPs'], 'd-', color='red', label='MMM A*A FLOPs (Teoretyczne)', linewidth=2)

    # Linie referencyjne O(N)
    ref_n = df['N'] * (df['MVM_FLOPs'].iloc[0] / df['N'].iloc[0])
    plt.plot(df['N'], ref_n, '--', color='gray', alpha=0.5, label='Referencja O(N)')

    plt.xscale('log')
    plt.yscale('log')
    plt.xlabel('Rozmiar macierzy N', fontsize=12)
    plt.ylabel('Liczba operacji (FLOPs)', fontsize=12)
    plt.title('Złożoność Obliczeniowa (FLOPs)', fontsize=14)
    plt.legend()
    plt.grid(True, which="both", ls="-", alpha=0.4)

    plt.tight_layout()
    plt.savefig('flops_complexity.png')
    print("Zapisano: flops_complexity.png")
    plt.close()

    # ---------------------------------------------------------
    # Wykres 3: Analiza Błędów (Error Analysis)
    # ---------------------------------------------------------
    # Sprawdzamy czy nowe kolumny istnieją (dla kompatybilności wstecznej)
    if 'SSE_Error' in df.columns and 'Matrix_Diff_Norm' in df.columns:
        plt.figure(figsize=(8, 6))
        
        # SSE Error (Suma Kwadratów Różnic dla wektora wynikowego)
        plt.plot(df['N'], df['SSE_Error'], 'x-', color='purple', label='Błąd MVM (SSE)', linewidth=2)
        
        # Matrix Diff Norm (Norma Frobeniusa różnicy macierzy)
        plt.plot(df['N'], df['Matrix_Diff_Norm'], '^-', color='orange', label='Błąd Rekonstrukcji (Norma Frobeniusa)', linewidth=2)

        plt.xscale('log')
        plt.yscale('log') # Logarytmiczna, bo błędy mogą się różnić rzędami wielkości
        plt.xlabel('Rozmiar macierzy N', fontsize=12)
        plt.ylabel('Wartość błędu', fontsize=12)
        plt.title('Analiza Dokładności Aproksymacji', fontsize=14)
        plt.legend()
        plt.grid(True, which="both", ls="-", alpha=0.4)

        plt.tight_layout()
        plt.savefig('error_analysis.png')
        print("Zapisano: error_analysis.png")
        plt.close()
    else:
        print("Pominięto wykres błędów (brak kolumn SSE_Error/Matrix_Diff_Norm w CSV).")

if __name__ == "__main__":
    plot_results()