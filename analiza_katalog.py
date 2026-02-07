#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
ANALIZA KATALOGU ASTRONOMICZNEGO - NAPRAWIONA WERSJA
Pełna analiza kompletności danych w katalog_astro_full.csv
"""

import pandas as pd
import numpy as np

# ==================== KONFIGURACJA ====================
CATALOG_FILE = "katalog_astro_full.csv"

# ==================== FUNKCJE POMOCNICZE ====================

def is_messier_token(token):
    """Sprawdza czy token to Messier (M1-M110)"""
    token = token.strip()
    if not token.startswith('M'):
        return False
    rest = token[1:]
    return rest.isdigit() and 1 <= int(rest) <= 110

def is_caldwell_token(token):
    """Sprawdza czy token to Caldwell (C1-C109)"""
    token = token.strip()
    if not token.startswith('C'):
        return False
    rest = token[1:]
    return rest.isdigit() and 1 <= int(rest) <= 109

def is_herschel_token(token):
    """Sprawdza czy token to Herschel (H)"""
    token = token.strip()
    if not token.startswith('H'):
        return False
    rest = token[1:]
    return rest.isdigit() and 1 <= int(rest) <= 400

def parse_extra_info(extra_info_str):
    """
    Parsuje extra_info i zwraca dict z flagami:
    {
        'has_messier': bool,
        'has_caldwell': bool,
        'has_herschel': bool,
        'messier_num': int lub None,
        'caldwell_num': int lub None
    }
    """
    result = {
        'has_messier': False,
        'has_caldwell': False,
        'has_herschel': False,
        'messier_num': None,
        'caldwell_num': None
    }
    
    if pd.isna(extra_info_str) or extra_info_str == "":
        return result
    
    # Tokenizacja po przecinkach
    tokens = [t.strip() for t in str(extra_info_str).split(',')]
    
    for token in tokens:
        if is_messier_token(token):
            result['has_messier'] = True
            result['messier_num'] = int(token[1:])
        elif is_caldwell_token(token):
            result['has_caldwell'] = True
            result['caldwell_num'] = int(token[1:])
        elif is_herschel_token(token):
            result['has_herschel'] = True
    
    return result

# ==================== GŁÓWNA ANALIZA ====================

def analyze_catalog():
    print("=" * 80)
    print("ANALIZA KATALOGU: katalog_astro_full.csv")
    print("=" * 80)
    print()
    
    # Wczytaj katalog
    try:
        df = pd.read_csv(CATALOG_FILE)
    except FileNotFoundError:
        print(f"BŁĄD: Nie znaleziono pliku {CATALOG_FILE}")
        return
    
    total = len(df)
    
    # [1] OGÓLNE STATYSTYKI
    print("[1] OGÓLNE STATYSTYKI")
    print(f"    Liczba obiektów:     {total:,}")
    print(f"    Liczba kolumn:       {len(df.columns)}")
    print(f"    Kolumny:             {list(df.columns)}")
    print()
    
    # [2] ANALIZA KOLUMNY 'mag' (Jasność całkowita)
    print("[2] ANALIZA KOLUMNY 'mag' (Jasność całkowita)")
    mag_present = df['mag'].notna().sum()
    mag_missing = df['mag'].isna().sum()
    print(f"    Obecne wartości:     {mag_present:,} ({mag_present/total*100:.1f}%)")
    print(f"    Brakujące (NaN):     {mag_missing:,} ({mag_missing/total*100:.1f}%)")
    
    if mag_present > 0:
        print(f"    Min. jasność (mag):  {df['mag'].min():.2f}")
        print(f"    Max. jasność (mag):  {df['mag'].max():.2f}")
        print(f"    Średnia:             {df['mag'].mean():.2f}")
        print(f"    Mediana:             {df['mag'].median():.2f}")
        print()
        print("    Rozkład jasności:")
        bright = (df['mag'] < 5).sum()
        medium_bright = ((df['mag'] >= 5) & (df['mag'] < 8)).sum()
        medium = ((df['mag'] >= 8) & (df['mag'] < 12)).sum()
        faint = (df['mag'] >= 12).sum()
        print(f"      mag < 5:        {bright:5} ({bright/mag_present*100:5.1f}%)  - Bardzo jasne")
        print(f"      5 ≤ mag < 8:     {medium_bright:5} ({medium_bright/mag_present*100:5.1f}%)  - Jasne")
        print(f"      8 ≤ mag < 12:   {medium:5} ({medium/mag_present*100:5.1f}%)  - Średnie")
        print(f"      mag ≥ 12:    {faint:5} ({faint/mag_present*100:5.1f}%)  - Słabe")
    print()
    
    # [3] ANALIZA KOLUMNY 'size' (Rozmiar kątowy w minutach)
    print("[3] ANALIZA KOLUMNY 'size' (Rozmiar kątowy w minutach)")
    size_present = df['size'].notna().sum()
    size_missing = df['size'].isna().sum()
    print(f"    Obecne wartości:     {size_present:,} ({size_present/total*100:.1f}%)")
    print(f"    Brakujące (NaN):     {size_missing:,} ({size_missing/total*100:.1f}%)")
    
    if size_present > 0:
        print(f"    Min. rozmiar (arcmin): {df['size'].min():.2f}'")
        print(f"    Max. rozmiar (arcmin): {df['size'].max():.2f}'")
        print(f"    Średnia:               {df['size'].mean():.2f}'")
        print(f"    Mediana:               {df['size'].median():.2f}'")
        print()
        print("    Rozkład rozmiaru (dla RedCat 61 mm FOV ~3°):")
        very_small = (df['size'] < 5).sum()
        small = ((df['size'] >= 5) & (df['size'] < 20)).sum()
        medium = ((df['size'] >= 20) & (df['size'] < 60)).sum()
        large = ((df['size'] >= 60) & (df['size'] < 180)).sum()
        very_large = (df['size'] >= 180).sum()
        print(f"      < 5' (bardzo małe):     {very_small:5} ({very_small/size_present*100:5.1f}%)  [crop]")
        print(f"      5-20':                   {small:4} ({small/size_present*100:5.1f}%)  [średnie]")
        print(f"      20-60':                   {medium:3} ({medium/size_present*100:5.1f}%)  [dobre]")
        print(f"      60-180' (>1°):            {large:3} ({large/size_present*100:5.1f}%)  [duże]")
        print(f"      ≥ 180' (wielkie):          {very_large:2} ({very_large/size_present*100:5.1f}%)  [bardzo duże]")
    print()
    
    # [4] ANALIZA KOLUMNY 'type' (Typ obiektu)
    print("[4] ANALIZA KOLUMNY 'type' (Typ obiektu)")
    type_present = df['type'].notna().sum()
    type_missing = df['type'].isna().sum()
    print(f"    Obecne wartości:     {type_present:,} ({type_present/total*100:.1f}%)")
    print(f"    Brakujące (NaN):     {type_missing:,} ({type_missing/total*100:.1f}%)")
    print()
    
    if type_present > 0:
        type_counts = df['type'].value_counts().head(20)
        print(f"    Typy obiektów:")
        for i, (obj_type, count) in enumerate(type_counts.items(), 1):
            print(f"      {i:2}. {obj_type:20} {count:5} ({count/total*100:5.1f}%)")
    print()
    
    # [5] ANALIZA KOLUMNY 'catalog' (Katalog źródłowy)
    if 'catalog' in df.columns:
        print("[5] ANALIZA KOLUMNY 'catalog' (Katalog źródłowy)")
        catalog_counts = df['catalog'].value_counts()
        for cat, count in catalog_counts.items():
            print(f"    {cat:10} {count:6} ({count/total*100:5.1f}%)")
        print()
    
    # [6] ANALIZA 'extra_info' (Messier/Caldwell/Herschel) - POPRAWIONA
    print("[6] ANALIZA 'extra_info' (Messier/Caldwell/Herschel)")
    
    # Parsuj wszystkie extra_info
    parsed = df['extra_info'].apply(parse_extra_info)
    
    messier_count = sum(p['has_messier'] for p in parsed)
    caldwell_count = sum(p['has_caldwell'] for p in parsed)
    herschel_count = sum(p['has_herschel'] for p in parsed)
    
    print(f"    Messier (M):            {messier_count} ({messier_count/total*100:5.1f}%)")
    print(f"    Caldwell (C):           {caldwell_count} ({caldwell_count/total*100:5.1f}%)")
    print(f"    Herschel (H):            {herschel_count} ({herschel_count/total*100:5.1f}%)")
    print()
    
    # [7] MACIERZ KOMPLETNOŚCI DANYCH
    print("[7] MACIERZ KOMPLETNOŚCI DANYCH")
    print("    (jakie kolumny zawsze razem się pojawiają)")
    print()
    
    complete = df[df['ra'].notna() & df['dec'].notna() & df['type'].notna() & 
                  df['mag'].notna() & df['size'].notna()]
    print(f"    Pełne dane (ra + dec + type + mag + size):")
    print(f"      {len(complete):,} obiektów ({len(complete)/total*100:.1f}%)")
    print()
    
    pos_only = df[df['ra'].notna() & df['dec'].notna() & df['type'].notna()]
    print(f"    Tylko poz. (ra + dec + type):")
    print(f"      {len(pos_only):,} obiektów ({len(pos_only)/total*100:.1f}%)")
    print()
    
    mag_no_size = df[df['mag'].notna() & df['size'].isna()]
    print(f"    Bez rozmiar (mag, ale brak size):")
    print(f"      {len(mag_no_size):,} obiektów ({len(mag_no_size)/total*100:.1f}%)")
    print()
    
    size_no_mag = df[df['size'].notna() & df['mag'].isna()]
    print(f"    Bez jasności (size, ale brak mag):")
    print(f"      {len(size_no_mag):,} obiektów ({len(size_no_mag)/total*100:.1f}%)")
    print()
    
    # [8] ANALIZA PROBLEMÓW
    print("[8] ANALIZA PROBLEMÓW")
    large_no_mag = df[(df['size'] > 10) & df['mag'].isna()]
    print(f"    Obiekty >10' bez mag: {len(large_no_mag):,} ({len(large_no_mag)/total*100:.1f}%)")
    
    bright_no_size = df[(df['mag'] < 10) & df['size'].isna()]
    print(f"    Obiekty mag<10 bez size: {len(bright_no_size):,} ({len(bright_no_size)/total*100:.1f}%)")
    
    no_type = df[df['type'].isna()]
    print(f"    Bez type (typ obiektu): {len(no_type):,} ({len(no_type)/total*100:.1f}%)")
    print()
    
    # [9] REKOMENDACJE DLA FILTRU SIZE
    print("[9] REKOMENDACJE DLA FILTRU SIZE")
    if size_present > 0:
        print("    Percentyle rozmiarów (dla tych, które go mają):")
        print(f"      Q1 (25%):            {df['size'].quantile(0.25):.1f}'")
        print(f"      Q2 (50%, mediana):   {df['size'].quantile(0.50):.1f}'")
        print(f"      Q3 (75%):            {df['size'].quantile(0.75):.1f}'")
        print(f"      Q4 (90%):            {df['size'].quantile(0.90):.1f}'")
        print()
        print("    Dla RedCat (FOV ~3°):")
        fov_deg = 3.0
        fov_arcmin = fov_deg * 60
        min_strict = fov_arcmin * 0.1
        elimination_pct = (df['size'] < min_strict).sum() / size_present * 100
        print(f"      - Min rozmiar strict:  10% FOV = {min_strict:.0f}' (eliminuje ~{elimination_pct:.1f}%)")
        print(f"      - Min rozmiar soft:     5' (elastycznie, dla crop)")
        print(f"      - Preferowany:         20-120' (duże, bez cropowania)")
    print()
    
    print("=" * 80)

if __name__ == "__main__":
    analyze_catalog()
