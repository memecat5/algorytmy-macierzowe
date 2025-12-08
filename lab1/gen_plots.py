#!/usr/bin/python3

import matplotlib.pyplot as plt
import pandas as pd

csv_data = pd.read_csv('results.csv', sep=',')
csv_data_strassen = pd.read_csv('results_strassen.csv', sep=',')

n_data = list(csv_data['n'])
time_data = list(csv_data['time_seconds'])
mem_data = list(csv_data['mem_delta_MB'])

n_data_strassen = list(csv_data_strassen['n'])
time_data_strassen = list(csv_data_strassen['time_seconds'])
mem_data_strassen = list(csv_data_strassen['mem_delta_MB'])

# ==============================
# 1. Wykres czasu działania
# ==============================
plt.figure(figsize=(7,5))
#plt.plot(n_data, time_data, 'bo-', label='Czas [s] (Binet)')
plt.plot(n_data_strassen, time_data_strassen, 'ro-', label='Czas [s] (Strassen)')
plt.xlabel('Rozmiar macierzy n')
plt.ylabel('Czas działania [s]')
plt.title('Czas działania metody Strassena')
plt.grid(True)
plt.legend()
plt.savefig('czas_dzialania.png', dpi=150)
# plt.show()

# ==============================
# 2. Wykres zużycia pamięci
# ==============================
plt.figure(figsize=(7,5))
#plt.plot(n_data, mem_data, 'bo-', label='Pamięć [MB] (Binet)')
plt.plot(n_data_strassen, mem_data_strassen, 'ro-', label='Pamięć [MB] (Strassen)')
plt.xlabel('Rozmiar macierzy n')
plt.ylabel('Zużycie pamięci [MB]')
plt.title('Zużycie pamięci metody Strassena')
plt.grid(True)
plt.legend()
plt.savefig('zuzycie_pamieci.png', dpi=150)
# plt.show()

fadd_data = list(csv_data_strassen['fadd_count'])
fmult_data = list(csv_data_strassen['fmult_count'])
flop_data = [sum(x) for x in zip(fadd_data, fmult_data)]

# ==============================
# 3. Wykres liczby operacji zmiennoprzecinkowych
# ==============================
plt.figure(figsize=(7,5))
plt.plot(n_data_strassen, fadd_data, 'o-', label='Operacje dodawania')
plt.plot(n_data_strassen, fmult_data, 'o-', label='Operacje mnożenia')
plt.plot(n_data_strassen, flop_data, 'o-', label='Wszystkie operacje')
plt.xlabel('Rozmiar macierzy n')
plt.ylabel('Liczba operacji')
plt.title('Liczba operacji metody Strassena')
plt.grid(True)
plt.legend()
plt.savefig('liczba_operacji.png', dpi=150)
# plt.show()
