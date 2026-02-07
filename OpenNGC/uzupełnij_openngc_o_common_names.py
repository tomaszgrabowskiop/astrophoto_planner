import pandas as pd
import numpy as np

# Ścieżka do pliku
input_file = 'OpenNGC/NGC.csv'
output_file = 'OpenNGC/NGC_updated.csv'

# Mapa common names
common_names_map = {
    'NGC1746': 'Cluster of Clusters',
    'NGC2168': 'Shoe-Buckle Cluster',
    'NGC7092': 'Pyramid Cluster',
    'NGC7380': 'Wizard Nebula',
    'NGC1333': 'Embryo Nebula',
    'NGC2170': 'Angel Nebula',
    'IC0447': "Dreyer's Nebula",
    'IC1396': "Elephant's Trunk Nebula",
    'NGC7822': "Question Mark Nebula",
}

# Wczytaj CSV – NIC nie kombinujemy z sep/kodowaniem poza tym co już było
df = pd.read_csv(input_file, sep=';', encoding='utf-8', low_memory=False)

# --- KLUCZOWA CZĘŚĆ: napraw kolumnę z Messierem po wczytaniu ---

# Załóżmy, że kolumna z numerem Messiera nazywa się np. 'M' (dostosuj NAZWĘ!)
# Jeśli nazywa się inaczej (np. 'Messier', 'M_ID'), popraw nazwę niżej.
M_COL = 'M'   # <- tu wpisz prawdziwą nazwę kolumny z numerem Messiera

if M_COL in df.columns:
    def normalize_m_id(x):
        if x is None:
            return ""
        s = str(x).strip()
        if not s:
            return ""
        try:
            # Obsłuży '24', '24.0', ' 24.0 ' -> '24'
            return str(int(float(s)))
        except ValueError:
            # Jak są śmieci, to lepiej je wyzerować niż rozwalić 0_opracuj
            return ""

    df[M_COL] = df[M_COL].apply(normalize_m_id)

# --- Uzupełnij kolumnę "Common names" jak wcześniej ---

df['Common names'] = df['Name'].map(common_names_map).fillna(df['Common names'])

# Zapisz – ważne: niech separator i encoding pozostaną takie, jak wymaga Twój pipeline
df.to_csv(output_file, sep=';', index=False, encoding='utf-8')

print("Zaktualizowano:", output_file)
print("Zmienione:")
print(df[df['Name'].isin(common_names_map.keys())][['Name', 'Common names']])
