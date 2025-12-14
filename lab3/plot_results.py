import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
import os
from scipy.optimize import curve_fit

# Definicja funkcji modelowej: alpha * N^beta
def power_law(x, alpha, beta):
    return alpha * np.power(x, beta)

def fit_and_plot(ax, x_data, y_data, color, label_prefix):
    """
    Funkcja pomocnicza: liczy dopasowanie i rysuje linię trendu
    """
    try:
        # Dopasowanie krzywej
        popt, pcov = curve_fit(power_law, x_data, y_data)
        alpha, beta = popt
        
        # Generowanie punktów do gładkiej linii
        x_fit = np.linspace(min(x_data), max(x_data), 100)
        y_fit = power_law(x_fit, alpha, beta)
        
        # Rysowanie linii dopasowania
        ax.plot(x_fit, y_fit, '--', color=color, alpha=0.8, 
                label=f'{label_prefix} Fit: $\\beta={beta:.3f}$')
        
        return alpha, beta
    except Exception as e:
        print(f"Nie udało się dopasować funkcji dla {label_prefix}: {e}")
        return None, None

def plot_results():
    if not os.path.exists('results.csv'):
        print("Błąd: Nie znaleziono pliku 'results.csv'.")
        return

    try:
        df = pd.read_csv('results.csv')
        df = df.sort_values('N') # Sortowanie dla poprawnego rysowania linii
    except Exception as e:
        print(f"Błąd wczytywania CSV: {e}")
        return

    plt.style.use('seaborn-v0_8-whitegrid')

    # ---------------------------------------------------------
    # Wykres 1: Czas wykonania z dopasowaniem (Time Complexity)
    # ---------------------------------------------------------
    fig, ax = plt.subplots(figsize=(8, 6))
    
    # Dane pomiarowe
    ax.plot(df['N'], df['Time_Compression'], 'o', label='Kompresja (pomiary)', alpha=0.6)
    ax.plot(df['N'], df['Time_MVM'], 's', color='orange', label='MVM (pomiary)', alpha=0.6)
    
    # Dopasowanie funkcji dla MVM Time
    fit_and_plot(ax, df['N'], df['Time_MVM'], 'darkorange', 'MVM Time')
    
    # Dopasowanie funkcji dla Kompresji (opcjonalnie)
    fit_and_plot(ax, df['N'], df['Time_Compression'], 'blue', 'Kompresja')

    ax.set_xscale('log')
    ax.set_yscale('log')
    ax.set_xlabel('Rozmiar macierzy N', fontsize=12)
    ax.set_ylabel('Czas [s]', fontsize=12)
    ax.set_title('Czas wykonania vs Rozmiar Problemu', fontsize=14)
    ax.legend()
    ax.grid(True, which="both", ls="-", alpha=0.4)
    
    plt.tight_layout()
    plt.savefig('time_complexity.png')
    print("Zapisano: time_complexity.png")
    plt.close()

    # ---------------------------------------------------------
    # Wykres 2: Złożoność FLOPs z dopasowaniem
    # ---------------------------------------------------------
    fig, ax = plt.subplots(figsize=(8, 6))
    
    # Dane teoretyczne
    ax.plot(df['N'], df['MVM_FLOPs'], 'o', color='green', label='MVM FLOPs (Dane)', alpha=0.6)
    
    # Dopasowanie funkcji dla MVM FLOPs
    alpha, beta = fit_and_plot(ax, df['N'], df['MVM_FLOPs'], 'darkgreen', 'MVM FLOPs')
    
    #MMM (dla porównania)
    ax.plot(df['N'], df['MMM_FLOPs'], 'd', color='red', label='MMM FLOPs (Dane)', alpha=0.6)
    fit_and_plot(ax, df['N'], df['MMM_FLOPs'], 'darkred', 'MMM FLOPs')

    ax.set_xscale('log')
    ax.set_yscale('log')
    ax.set_xlabel('Rozmiar macierzy N', fontsize=12)
    ax.set_ylabel('Liczba operacji (FLOPs)', fontsize=12)
    ax.set_title(f'Złożoność Obliczeniowa\nDla MVM: $\\alpha \\approx {alpha:.2e}, \\beta \\approx {beta:.3f}$', fontsize=14)
    ax.legend()
    ax.grid(True, which="both", ls="-", alpha=0.4)

    plt.tight_layout()
    plt.savefig('flops_complexity.png')
    print("Zapisano: flops_complexity.png")
    plt.close()

    # ---------------------------------------------------------
    # Wykres 3: Analiza Błędów (bez zmian)
    # ---------------------------------------------------------
    if 'SSE_Error' in df.columns:
        plt.figure(figsize=(8, 6))
        plt.plot(df['N'], df['SSE_Error'], 'x-', color='purple', label='Błąd MVM (SSE)')
        if 'Matrix_Diff_Norm' in df.columns:
            plt.plot(df['N'], df['Matrix_Diff_Norm'], '^-', color='orange', label='Błąd Rekonstrukcji')
        
        plt.xscale('log')
        plt.yscale('log')
        plt.xlabel('Rozmiar macierzy N', fontsize=12)
        plt.ylabel('Błąd', fontsize=12)
        plt.title('Analiza Dokładności', fontsize=14)
        plt.legend()
        plt.grid(True, which="both", ls="-", alpha=0.4)
        plt.tight_layout()
        plt.savefig('error_analysis.png')
        print("Zapisano: error_analysis.png")
        plt.close()

if __name__ == "__main__":
    plot_results()