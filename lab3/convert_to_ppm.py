from PIL import Image
import sys

def convert_image(input_path, output_name="original.ppm", size=(512, 512)):
    try:
        # Otwieramy obraz (JPG, PNG, cokolwiek)
        img = Image.open(input_path)
        
        # Konwertujemy na RGB (zeby usunac przezroczystosc jesli jest)
        img = img.convert("RGB")
        
        # Skalujemy do wymaganego rozmiaru (np. 512x512)
        img = img.resize(size, Image.Resampling.LANCZOS)
        
        # Zapisujemy jako PPM (Format P6 jest domyślny w PIL dla save ppm)
        img.save(output_name)
        print(f"Sukces! Zapisano jako {output_name} ({size[0]}x{size[1]})")
        
    except Exception as e:
        print(f"Błąd: {e}")

if __name__ == "__main__":
    if len(sys.argv) < 2:
        print("Użycie: python convert_to_ppm.py <twoje_zdjecie.jpg>")
    else:
        convert_image(sys.argv[1])