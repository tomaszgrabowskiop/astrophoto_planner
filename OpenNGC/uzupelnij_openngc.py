import pandas as pd
import numpy as np

# Ścieżka do pliku
input_file = 'NGC.csv'
output_file = 'NGC_updated.csv'

# Mapa common names
common_names_map = {
    'NGC0188': 'North Celestial Pole Cluster',
    'NGC1333': 'Embryo Nebula',
    'NGC1746': 'Cluster of Clusters',
    'NGC1788': 'Cosmic Bat Nebula',
    'NGC1960': 'Pinwheel Cluster',
    'NGC2168': 'Shoe-Buckle Cluster',
    'NGC2170': 'Angel Nebula',
    'NGC2682': 'King Cobra Cluster',
    'NGC5457': 'Pinwheel Galaxy',
    'NGC7078': 'Great Pegasus Cluster',
    'NGC7092': 'Pyramid Cluster',
    'NGC7380': 'Wizard Nebula',
    'NGC7822': "Question Mark Nebula",
    'IC0423': 'Tear Drop Nebula',
    'IC0446': 'Coyote Cloud',
    'IC0447': "Dreyer's Nebula",
    'IC1396': "Elephant's Trunk Nebula",
    'IC1613': 'Cetus Dwarf Galaxy',
    'IC4665': 'Summer Beehive Cluster',
    'IC4756': 'Graff’s Cluster',
}
# Mapa numerów Messiera
messier_map = {
    "NGC1432": "45"
}
# Wczytaj CSV – NIC nie kombinujemy z sep/kodowaniem poza tym co już było
df = pd.read_csv(input_file, sep=';', encoding='utf-8', low_memory=False)

# --- KLUCZOWA CZĘŚĆ: napraw kolumnę z Messierem po wczytaniu ---

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
df[M_COL] = df['Name'].map(messier_map).fillna(df[M_COL])
df['Common names'] = df['Name'].map(common_names_map).fillna(df['Common names'])

# Zapisz
df.to_csv(output_file, sep=';', index=False, encoding='utf-8')

print("Zaktualizowano:", output_file)
print("\nZmienione common names:")
print(df[df['Name'].isin(common_names_map.keys())][['Name', 'Common names']])
print("\nZmienione numery Messiera:")
print(df[df['Name'].isin(messier_map.keys())][['Name', M_COL]])