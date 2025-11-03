#!/usr/bin/python3

import matplotlib.pyplot as plt
import pandas as pd

csv_data = pd.read_csv('results.csv', sep=',')

n_data = list(csv_data['n'])
time_data = list(csv_data['time_seconds'])
mem_data = list(csv_data['mem_delta_MB'])

# ==============================
# 1. Wykres czasu działania
# ==============================
plt.figure(figsize=(7,5))
plt.plot(n_data, time_data, 'bo-', label='Czas [s]')
plt.xlabel('Rozmiar macierzy n')
plt.ylabel('Czas działania [s]')
plt.title('Czas działania metody rekurencyjnej')
plt.grid(True)
plt.legend()
plt.savefig('czas_dzialania.png', dpi=150)
# plt.show()

# ==============================
# 2. Wykres zużycia pamięci
# ==============================
plt.figure(figsize=(7,5))
plt.plot(n_data, mem_data, 'go-', label='Pamięć [MB]')
plt.xlabel('Rozmiar macierzy n')
plt.ylabel('Zużycie pamięci [MB]')
plt.title('Zużycie pamięci metody rekurencyjnej')
plt.grid(True)
plt.legend()
plt.savefig('zuzycie_pamieci.png', dpi=150)
# plt.show()

fadd_data = list(csv_data['fadd_count'])
fmult_data = list(csv_data['fmult_count'])
flop_data = [sum(x) for x in zip(fadd_data, fmult_data)]

# ==============================
# 3. Wykres liczby operacji zmiennoprzecinkowych
# ==============================
plt.figure(figsize=(7,5))
plt.plot(n_data, fadd_data, 'ro-', label='Operacje dodawania')
plt.plot(n_data, fmult_data, 'ro-', label='Operacje mnożenia')
plt.plot(n_data, flop_data, 'ro-', label='Wszystkie operacje')
plt.xlabel('Rozmiar macierzy n')
plt.ylabel('Liczba operacji')
plt.title('Liczba operacji metody rekurencyjnej')
plt.grid(True)
plt.legend()
plt.savefig('liczba_operacji.png', dpi=150)
# plt.show()
