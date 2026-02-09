import re
import warnings
import numpy as np
import pandas as pd
import networkx as nx
import astropy.units as u
from astroquery.vizier import Vizier
from astropy.coordinates import SkyCoord, search_around_sky
from analiza_katalog import analyze_catalog

warnings.filterwarnings('ignore')

# === SSL FIX (macOS + oficjalny Python nie ma certyfikatów) ===
import ssl
import certifi
ssl._create_default_https_context = lambda: ssl.create_default_context(cafile=certifi.where())

# === KONFIGURACJA DOMYŚLNA ===
RA_DEC_TOLERANCE_ARCMIN = 1.2  # Domyślnie 0.02 stopnia
OUTPUT_FILENAME = "katalog_astro_full.csv"
MIN_SIZE_ARCMIN = 5.0
MAX_MAG = 14.0

PRIORITIES = {'ngc': 1, 'sh2': 2, 'rcw': 3, 'lbn': 4, 'ldn': 5, 'barn': 6, 'ced': 7, 'pgc': 8}
CATALOG_SCORES = {'ngc': 0, 'sh2': 30, 'rcw': 25, 'ced': 20, 'barn': -20, 'ldn': -25, 'lbn': 15, 'pgc': 10}

SH2_COMMON_NAMES = {
    "Sh2-101": "Tulip Nebula", "Sh2-103": "Loop", "Sh2-105": "Crescent Nebula",
    "Sh2-108": "Sadr Region", "Sh2-11": "War and Peace Nebula", "Sh2-117": "N American & Pelican Nebula",
    "Sh2-125": "Cocoon Nebula", "Sh2-129": "Flying Bat Nebula", "Sh2-155": "Cave Nebula",
    "Sh2-162": "Bubble Nebula", "Sh2-184": "Pac Man Nebula", "Sh2-190": "Heart Nebula",
    "Sh2-197": "Maffei 2", "Sh2-199": "Soul Nebula", "Sh2-220": "California Nebula",
    "Sh2-229": "Flaming Star Nebula", "Sh2-234": "Spider Nebula", "Sh2-237": "Fly Nebula",
    "Sh2-238": "Hind's Variable Nebula", "Sh2-244": "Crab Nebula", "Sh2-245": "Fishhook Nebula",
    "Sh2-248": "Jellyfish Nebula", "Sh2-25": "Lagoon Nebula", "Sh2-252": "Monkey Head Nebula",
    "Sh2-261": "Lower's Nebula", "Sh2-264": "Angelfish Nebula", "Sh2-273": "Fox Fur Nebula",
    "Sh2-274": "Medusa Nebula", "Sh2-275": "Rosette Nebula", "Sh2-276": "Barnard's Loop",
    "Sh2-277": "Flame Nebula", "Sh2-279": "Running Man Nebula", "Sh2-281": "Orion Nebula",
    "Sh2-292": "Seagull Nebula head", "Sh2-296": "Seagull Nebula wings", "Sh2-298": "Thor's Helmet",
    "Sh2-30": "Trifid Nebula", "Sh2-311": "Skull & Crossbones Nebula", "Sh2-45": "Omega Nebula",
    "Sh2-49": "Eagle N	ebula", "Sh2-54": "Cauda", "Sh2-6": "Bug Nebula", "Sh2-8": "Cat's Paw Nebula",
}

def print_step(msg):
    print(f"\n[INFO] {msg}")

def fmt(n: int) -> str:
    """Format liczby z separatorem tysięcy (spacja)."""
    return f"{int(n):_}".replace("_", " ")

def convert_ra_dec(df, ra_col, dec_col, unit_type='sexagesimal'):
    if df.empty or ra_col not in df.columns or dec_col not in df.columns:
        return df
    df = df.dropna(subset=[ra_col, dec_col])
    try:
        unit = (u.hourangle, u.deg) if unit_type == 'sexagesimal' else u.deg
        c = SkyCoord(ra=df[ra_col].astype(str).values,
                     dec=df[dec_col].astype(str).values,
                     unit=unit, frame='icrs')
        df['ra'], df['dec'] = np.round(c.ra.degree, 5), np.round(c.dec.degree, 5)
    except Exception as e:
        print(f"  [!] Błąd konwersji RA/DEC: {e}")
    return df

# === 1. POBIERANIE DANYCH ===
def fetch_data():
    print_step("Pobieranie katalogów...")
    
    # NGC/IC z lokalnego pliku
    print("  > NGC/IC z updated_ngc.csv...")
    df_ngc = pd.read_csv("updated_ngc.csv")
    print(f"    {fmt(len(df_ngc))} wierszy")
    
    # Konfiguracja Vizier
    v_std = Vizier(row_limit=-1, columns=['*', '_RAJ2000', '_DEJ2000'])
    
    def get_v(cat_id, name):
        print(f"  > {name} ({cat_id})...")
        try:
            cats = v_std.get_catalogs(cat_id)
            df = cats[0].to_pandas() if cats else pd.DataFrame()
            print(f"    {fmt(len(df))} wierszy")
            return df
        except Exception:
            print("    [!] Błąd połączenia.")
            return pd.DataFrame()

    data = {
        "ngc": df_ngc,
        "sh2": get_v('VII/20/catalog', 'Sharpless'),
        "barn": get_v('VII/220A', 'Barnard'),
        "rcw": get_v('VII/216', 'RCW'),
        "pgc": get_v('VII/119', 'PGC'),
        "ldn": get_v('VII/7A', 'Lynds Dark'),
        "lbn": get_v('VII/9', 'Lynds Bright'),
        "ced": get_v('VII/231', 'Cederblad'),
    }
    
    print_step("PODSUMOWANIE POBIERANIA:")
    total_rows = 0
    for k, v in data.items():
        cnt = len(v)
        print(f"  - {k.upper():8}: {fmt(cnt):>8} wierszy")
        total_rows += cnt
    print("=" * 31)
    print(f"  Razem: {fmt(total_rows):>12} wierszy")
    
    return data

# === 2. PEŁNA NORMALIZACJA ===
def normalize_all(raw):
    print_step("Normalizacja katalogów (unifikacja kolumn)...")
    results = []

    # NGC/IC - Bez zmian w extra_info
    df = raw['ngc'].copy()
    print("  > Normalizacja NGC/IC...")
    df['ra'] = pd.to_numeric(df['ra'], errors='coerce')
    df['dec'] = pd.to_numeric(df['dec'], errors='coerce')
    df['catalog'], df['priority'] = 'ngc', PRIORITIES['ngc']
    results.append(df)

    # SHARPLESS
    df = raw['sh2']
    if not df.empty:
        print("  > Normalizacja Sharpless...")
        df = convert_ra_dec(df, '_RAJ2000', '_DEJ2000', 'deg')
        df['id'] = df['Sh2'].apply(lambda x: f"Sh2-{int(x)}" if pd.notnull(x) else "")
        df['size'] = pd.to_numeric(df['Diam'], errors='coerce')
        df['mag'], df['type'], df['extra_info'] = np.nan, 'HII', ""
        df['common_names'] = df['Sh2'].apply(
            lambda x: SH2_COMMON_NAMES.get(f"Sh2-{int(x)}", "") if pd.notnull(x) else ""
        )
        df['catalog'], df['priority'] = 'sh2', PRIORITIES['sh2']
        results.append(df)

    # BARNARD
    df = raw['barn']
    if not df.empty:
        print("  > Normalizacja Barnard...")
        df = convert_ra_dec(df, '_RAJ2000', '_DEJ2000', 'deg')
        df['id'] = df['Barn'].apply(lambda x: f"B{str(x).strip()}")
        df['size'] = pd.to_numeric(df['Diam'], errors='coerce') if 'Diam' in df.columns else np.nan
        df['mag'], df['type'], df['extra_info'] = np.nan, 'DN', ""
        df['catalog'], df['priority'] = 'barn', PRIORITIES['barn']
        results.append(df)

    # RCW
    df = raw['rcw']
    if not df.empty:
        print("  > Normalizacja RCW...")
        df = convert_ra_dec(df, '_RAJ2000', '_DEJ2000', 'deg')
        df['id'] = df['RCW'].apply(lambda x: f"RCW{str(x).strip()}")
        df['size'], df['mag'], df['type'], df['extra_info'] = np.nan, np.nan, 'HII', ""
        df['catalog'], df['priority'] = 'rcw', PRIORITIES['rcw']
        results.append(df)

    # CEDERBLAD
    df = raw['ced']
    if not df.empty:
        print("  > Normalizacja Cederblad...")
        df = convert_ra_dec(df, '_RAJ2000', '_DEJ2000', 'deg')
        ids = []
        for _, row in df.iterrows():
            num = str(row['Ced']).strip()
            sub = str(row['m_Ced']).strip() if 'm_Ced' in row and pd.notnull(row['m_Ced']) and str(row['m_Ced']) != 'nan' else ""
            ids.append(f"Ced{num}{sub}")
        df['id'] = ids
        df['size'], df['mag'], df['type'], df['extra_info'] = np.nan, np.nan, 'NB', ""
        df['catalog'], df['priority'] = 'ced', PRIORITIES['ced']
        results.append(df)

    # LBN
    df = raw['lbn']
    if not df.empty:
        print("  > Normalizacja LBN...")
        df = convert_ra_dec(df, '_RAJ2000', '_DEJ2000', 'deg')
        df['id'] = df['Seq'].apply(lambda x: f"LBN{x}")
        df['size'], df['mag'], df['type'], df['extra_info'] = np.nan, np.nan, 'NB', ""
        df['catalog'], df['priority'] = 'lbn', PRIORITIES['lbn']
        results.append(df)

    # LDN
    df = raw['ldn']
    if not df.empty:
        print("  > Normalizacja LDN...")
        df = convert_ra_dec(df, '_RAJ2000', '_DEJ2000', 'deg')
        df['id'] = df['LDN'].apply(lambda x: f"LDN{x}")
        df['size'] = np.sqrt(pd.to_numeric(df['Area'], errors='coerce')) * 60
        df['mag'], df['type'], df['extra_info'] = np.nan, 'DN', ""
        df['catalog'], df['priority'] = 'ldn', PRIORITIES['ldn']
        results.append(df)

    # PGC
    df = raw['pgc']
    if not df.empty:
        print("  > Normalizacja PGC...")
        df = convert_ra_dec(df, '_RAJ2000', '_DEJ2000', 'deg')
        df['id'] = df['PGC'].apply(lambda x: f"PGC{int(x)}" if pd.notnull(x) else "")
        df['size'] = pd.to_numeric(df['MajAxis'], errors='coerce') / 60
        df['mag'] = pd.to_numeric(df['Btot'], errors='coerce')
        df['type'], df['extra_info'] = 'Gx', ""
        df['catalog'], df['priority'] = 'pgc', PRIORITIES['pgc']
        results.append(df)

    final_cols = ['id', 'ra', 'dec', 'size', 'mag', 'type', 'extra_info',
                  'common_names', 'priority', 'catalog']
    output = []
    for d in results:
        for c in final_cols:
            if c not in d.columns:
                d[c] = np.nan
        output.append(d[final_cols])
    return pd.concat(output, ignore_index=True)

# === 3. SMART MERGE (LOGIKA GRAFOWA) ===
def smart_merge(df, tolerance_deg):
    df = df.dropna(subset=['ra', 'dec']).reset_index(drop=True)
    print_step(f"Konsolidacja Smart Merge (tolerancja {tolerance_deg*60:.2f} arcmin)...")
    
    coords = SkyCoord(ra=df['ra'].values * u.deg, dec=df['dec'].values * u.deg)
    idx1, idx2, _, _ = search_around_sky(coords, coords, tolerance_deg * u.deg)
    
    g = nx.Graph()
    g.add_nodes_from(range(len(df)))
    for i1, i2 in zip(idx1, idx2):
        if i1 != i2:
            g.add_edge(i1, i2)
    
    clusters = list(nx.connected_components(g))
    print(f"       Znalazłem {fmt(len(clusters))} grup obiektów.")

    merged_rows = []
    for cluster in clusters:
        subset = df.iloc[list(cluster)]
        
        # Jeśli klaster ma tylko jeden obiekt, dodajemy go bez zmian
        if len(subset) == 1:
            merged_rows.append(subset.iloc[0].to_dict())
            continue

        # Sortujemy wg priorytetu (NGC=1, PGC=8) oraz rozmiaru
        subset_sorted = subset.sort_values(by=['priority', 'size'], ascending=[True, False])
        master = subset_sorted.iloc[0].copy()

        # 1. WSPÓŁRZĘDNE: Uśrednianie środka układu
        master['ra'] = np.round(subset['ra'].mean(), 5)
        master['dec'] = np.round(subset['dec'].mean(), 5)

        # 2. JASNOŚĆ (MAG): najlepsza dostępna
        mags_in_hierarchy = subset_sorted['mag'].dropna()
        master['mag'] = mags_in_hierarchy.iloc[0] if not mags_in_hierarchy.empty else np.nan

        # 3. ROZMIAR (SIZE): najlepszy dostępny
        sizes_in_hierarchy = subset_sorted['size'].dropna()
        master['size'] = sizes_in_hierarchy.iloc[0] if not sizes_in_hierarchy.empty else np.nan

        # 4. TEKSTY
        def combine_text(col):
            # Pobierz wartości, usuń rzeczywiste NaN (float)
            items = subset[col].dropna().astype(str)
            
            # Filtruj napisy, które wyglądają jak "nan" lub są puste
            clean_items = [
                x.strip() for x in items 
                if x.strip().lower() not in ['nan', 'none', '']
            ]
            
            # Rozdziel po przecinkach, jeśli w jednej komórce było wiele nazw
            all_parts = []
            for item in clean_items:
                all_parts.extend([p.strip() for p in item.split(',') if p.strip()])
                
            return ", ".join(sorted(set(all_parts)))

        # extra_info + ID klastrowe
        master_id = str(master['id'])
        combined_extra = combine_text('extra_info')
        all_ids_in_cluster = [str(i) for i in subset['id'].unique() if str(i) != master_id]
        final_extra_set = set([p.strip() for p in combined_extra.split(',') if p.strip()])
        final_extra_set.update(all_ids_in_cluster)
        master['extra_info'] = ", ".join(sorted(final_extra_set))
        
        # common_names
        master['common_names'] = combine_text('common_names')

        merged_rows.append(master.to_dict())

    return pd.DataFrame(merged_rows)

def main():
    print("=== GENERATOR KATALOGU ASTRONOMICZNEGO (ENTITY RESOLUTION) ===\n")
    
    raw = fetch_data()
    full = normalize_all(raw)
    
    print_step("DEBUG: priorytety po normalizacji")
    print(
        full.groupby('catalog')['priority']
            .agg(['min', 'max', 'nunique'])
            .reset_index()
            .sort_values('min')
    )
    
    print_step("KONFIGURACJA")
    global RA_DEC_TOLERANCE_ARCMIN, MIN_SIZE_ARCMIN, MAX_MAG
    
    # 1. Tolerancja w minutach z walidacją
    ans = input(
        f"Tolerancja łączenia obiektów w MINUTACH kątowych (arcmin) "
        f"[{RA_DEC_TOLERANCE_ARCMIN}]: "
    ).strip()
    if ans:
        val = float(ans)
        if val > 60:
            print(" [!] OSTRZEŻENIE: Tolerancja > 60' jest zbyt duża. Zredukowano do 60'.")
            val = 60.0
        RA_DEC_TOLERANCE_ARCMIN = val
    
    tol_deg = RA_DEC_TOLERANCE_ARCMIN / 60.0

    # 2. Reszta parametrów
    ans = input(f"Minimalny rozmiar MIN_SIZE_ARCMIN [{MIN_SIZE_ARCMIN}]: ").strip()
    if ans:
        MIN_SIZE_ARCMIN = float(ans)
    
    ans = input(f"Maksymalna jasność MAX_MAG [{MAX_MAG}]: ").strip()
    if ans:
        MAX_MAG = float(ans)

    # Filtrowanie typów PRZED Smart Merge
    print_step("Filtrowanie typów przed Smart Merge")
    pre_merge = full.copy()
    start_count = len(pre_merge)

    forbidden = ['*', '**', '*Ass', 'Dup', 'NonEx', 'SNR', 'PN', 'Nova', 'GPair', 'EmN']

    after_forbidden = pre_merge[~pre_merge['type'].isin(forbidden)]
    removed_forbidden = len(pre_merge) - len(after_forbidden)
    pre_merge = after_forbidden
    print(f"Usunąłem {fmt(removed_forbidden)} obiektów na podstawie typu (forbidden).")
    print(f"Do Smart Merge trafi {fmt(len(pre_merge))} obiektów po odrzuceniu typów forbidden.")

    # Smart Merge na danych po forbidden
    print_step("Smart Merge")
    final = smart_merge(pre_merge, tol_deg)
    print(f"Po Smart Merge: {fmt(len(final))} obiektów (po połączeniu duplikatów w obrębie tolerancji).")

    # Filtrowanie rozmiaru i jasności PO Smart Merge
    print_step("Filtrowanie rozmiaru i jasności po Smart Merge")
    start_count_post_merge = len(final)

    # Rozmiar
    before_size = len(final)
    after_size = final[(final['size'].isna()) | (final['size'] >= MIN_SIZE_ARCMIN)]
    removed_size = before_size - len(after_size)
    final = after_size
    print(
        f"Usunąłem {fmt(removed_size)} obiektów mniejszych niż {MIN_SIZE_ARCMIN} arcmin "
        f"(z wyjątkiem tych bez rozmiaru)."
    )

    # Jasność
    before_mag = len(final)
    after_mag = final[(final['mag'].isna()) | (final['mag'] <= MAX_MAG)]
    removed_mag = before_mag - len(after_mag)
    final = after_mag
    print(
        f"Usunąłem {fmt(removed_mag)} obiektów jaśniejszych niż {MAX_MAG} mag "
        f"(z wyjątkiem tych bez magnitudo)."
    )

    print(
        f"Łącznie usunąłem {fmt(start_count_post_merge - len(final))} obiektów "
        f"po Smart Merge w filtrach rozmiaru i jasności."
    )

    cols = ['id', 'extra_info', 'type', 'ra', 'dec', 'mag', 'size', 'common_names']
    final[cols].to_csv(OUTPUT_FILENAME, index=False)
    
    print_step(f"ZAPISAŁEM {fmt(len(final))} OBIEKTÓW.")
    
    answer = input('\nUruchomić pełną analizę katalogu? [y/n] (domyślnie y): ').strip().lower()
    if answer in ["y", ""]:
        analyze_catalog()
    else:
        print("Pominięto analizę.")

if __name__ == "__main__":
    main()