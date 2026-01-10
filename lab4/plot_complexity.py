import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
from scipy.optimize import curve_fit

# Funkcja modelu: alpha * N^beta
def complexity_model(x, alpha, beta):
    return alpha * np.power(x, beta)

# Wczytanie danych
try:
    df = pd.read_csv("results_task4.csv")
except FileNotFoundError:
    print("Błąd: Nie znaleziono pliku results_task4.csv. Uruchom najpierw program C++.")
    exit()

N_values = df["N"].values

# Przygotowanie wykresów
fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(14, 6))

# Wartości startowe dla optymalizatora: alpha mała, beta ~1.0
initial_guess = [1e-5, 1.0]

# --- Wykres 1: Mnożenie Macierz-Wektor ---
times_vec = df["Time_Vec_ms"].values
try:
    # POPRAWKA: Dodano p0 (start) i maxfev (liczba iteracji)
    popt_vec, _ = curve_fit(complexity_model, N_values, times_vec, p0=initial_guess, maxfev=10000)
    alpha_v, beta_v = popt_vec
    
    # Generowanie punktów do krzywej
    x_fit = np.linspace(min(N_values), max(N_values), 100)
    y_fit_vec = complexity_model(x_fit, alpha_v, beta_v)
    
    ax1.plot(x_fit, y_fit_vec, '--', color='blue', label=f'Fit: {alpha_v:.2e} * N^{beta_v:.2f}')
except RuntimeError:
    print("Nie udało się dopasować krzywej dla Mnożenia Wektora (zbyt mało danych lub szum).")
    alpha_v, beta_v = 0, 0

ax1.scatter(N_values, times_vec, color='red', label='Dane eksperymentalne', zorder=5)
ax1.set_title("Złożoność: Mnożenie Macierz-Wektor")
ax1.set_xlabel("Rozmiar macierzy (N)")
ax1.set_ylabel("Czas (ms)")
ax1.legend()
ax1.grid(True)

# --- Wykres 2: Mnożenie Macierz-Macierz ---
times_mat = df["Time_Mat_ms"].values
try:
    # POPRAWKA: Dodano p0 i maxfev
    popt_mat, _ = curve_fit(complexity_model, N_values, times_mat, p0=initial_guess, maxfev=10000)
    alpha_m, beta_m = popt_mat
    y_fit_mat = complexity_model(x_fit, alpha_m, beta_m)
    
    ax2.plot(x_fit, y_fit_mat, '--', color='blue', label=f'Fit: {alpha_m:.2e} * N^{beta_m:.2f}')
except RuntimeError:
    print("Nie udało się dopasować krzywej dla Mnożenia Macierzy.")
    alpha_m, beta_m = 0, 0

ax2.scatter(N_values, times_mat, color='green', label='Dane eksperymentalne', zorder=5)
ax2.set_title("Złożoność: Mnożenie Macierz-Macierz")
ax2.set_xlabel("Rozmiar macierzy (N)")
ax2.set_ylabel("Czas (ms)")
ax2.legend()
ax2.grid(True)

# Zapis i wyświetlenie
plt.tight_layout()
plt.savefig("complexity_plots.png")
print(f"Raport wygenerowany:")
print(f"Complexity Vector: alpha={alpha_v:.2e}, beta={beta_v:.4f}")
print(f"Complexity Matrix: alpha={alpha_m:.2e}, beta={beta_m:.4f}")
plt.show()