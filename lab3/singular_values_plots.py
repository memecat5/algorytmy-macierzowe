import numpy as np
import matplotlib.pyplot as plt
from PIL import Image
import sys

def analyze_image_svd(image_path):
    print(f"--- Analiza SVD dla: {image_path} ---")

    try:
        # 1. Wczytanie obrazu i konwersja na macierz
        img = Image.open(image_path).convert('RGB')
        img_arr = np.array(img)
        
        height, width, _ = img_arr.shape
        print(f"Wymiary obrazu: {width}x{height}")

        # Rozdzielamy kanały (normalizujemy 0-255 -> double)
        # R = img_arr[:, :, 0]
        channels = {
            'Red': img_arr[:, :, 0],
            'Green': img_arr[:, :, 1],
            'Blue': img_arr[:, :, 2]
        }

        plt.figure(figsize=(10, 6))
        colors = {'Red': 'red', 'Green': 'green', 'Blue': 'blue'}

        # 2. Obliczenie SVD dla każdego kanału
        for name, matrix in channels.items():
            print(f"\nLiczenie SVD dla kanału {name}...")
            
            # compute_uv=False oznacza, że liczymy tylko Sigmę (wartości),
            # bez macierzy U i V. To jest BARDZO szybkie.
            singular_values = np.linalg.svd(matrix, compute_uv=False)
            
            # Wypisanie statystyk
            print(f"Top 5 wartości: {singular_values[:5]}")
            print(f"Ostatnie 5 wartości: {singular_values[-5:]}")
            print(f"Suma wszystkich wartości (Energia): {np.sum(singular_values):.2f}")
            
            # Rysowanie wykresu (skala logarytmiczna lepiej pokazuje spadek)
            plt.plot(singular_values, label=f'Kanał {name}', color=colors[name], linewidth=1.5)

        # Konfiguracja wykresu
        plt.title(f'Rozkład wartości osobliwych (SVD)')
        plt.ylabel('Wartość osobliwa (skala log)')
        plt.xlabel('Indeks wartości (k)')
        plt.yscale('log') # Skala logarytmiczna - kluczowa przy SVD!
        plt.grid(True, which="both", ls="-", alpha=0.5)
        plt.legend()
        
        print("\nWyświetlanie wykresu...")
        plt.show()

    except FileNotFoundError:
        print("Błąd: Nie znaleziono pliku. Sprawdź ścieżkę.")
    except Exception as e:
        print(f"Wystąpił błąd: {e}")

if __name__ == "__main__":
    # Domyślnie szuka 'original.ppm', ale możesz podać inny plik jako argument
    filename = sys.argv[1] if len(sys.argv) > 1 else "original.ppm"
    analyze_image_svd(filename)