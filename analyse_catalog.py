#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
ANALIZA KATALOGU ASTRONOMICZNEGO
Pełna analiza kompletności danych w katalog_astro_full.csv
"""

import pandas as pd
import numpy as np
from shared import CATALOG_PRIORITY, CAT_ORDER, PATHS, fmt, print_step, print_green
# ═══════════════════════════════════════════════════
# KONFIGURACJA
# ═══════════════════════════════════════════════════
CATALOG_FILE = PATHS.catalog_full
SEP_MAIN = "═" * 119
SEP_SUB  = "─" * 114
# ═══════════════════════════════════════════════════
# SEKCJE RAPORTU
# ═══════════════════════════════════════════════════
def _section(num: int, title: str) -> None:
    print_green(f"\n[{num:02d}] {title}")
    print("     " + SEP_SUB)

def _row(label: str, n: int, total: int, width: int = 46) -> None:
    print(f"     {label:<{width}} {n:6,} / {total:,}  ({n / total * 100:5.1f}%)")
# ═══════════════════════════════════════════════════
# BLOK A – podstawowe statystyki kolumn
# ═══════════════════════════════════════════════════
def sec_basic(df: pd.DataFrame, total: int) -> None:
    """[01] Ogólne statystyki"""
    _section(1, "OGÓLNE STATYSTYKI")
    print(f"     Liczba obiektów : {total:,}")
    print(f"     Liczba kolumn   : {len(df.columns)}")
#    print(f"     Kolumny         : {list(df.columns)}")


def sec_magnitude(df: pd.DataFrame, total: int) -> None:
    """[02] Jasność całkowita (mag)"""
    _section(2, "JASNOŚĆ CAŁKOWITA – kolumna 'mag'")
    col = df["mag"]
    present = col.notna().sum()
    missing  = col.isna().sum()

    print(f"     Obecne   : {present:,} ({present / total * 100:.1f}%)")
    print(f"     Brakujące: {missing:,} ({missing / total * 100:.1f}%)")

    if not present:
        return

    b1 = (col < 5).sum()
    b2 = ((col >= 5)  & (col < 8)).sum()
    b3 = ((col >= 8)  & (col < 12)).sum()
    b4 = (col >= 12).sum()

    print("\n     Rozkład jasności:")
    print(f"       mag < 5        : {b1:5,} ({b1 / present * 100:5.1f}%)  – bardzo jasne")
    print(f"       5  ≤ mag < 8   : {b2:5,} ({b2 / present * 100:5.1f}%)  – jasne")
    print(f"       8  ≤ mag < 12  : {b3:5,} ({b3 / present * 100:5.1f}%)  – średnie")
    print(f"       mag ≥ 12       : {b4:5,} ({b4 / present * 100:5.1f}%)  – słabe")


def sec_size(df: pd.DataFrame, total: int) -> None:
    """[03] Rozmiar kątowy (size) + percentyle + rekomendacje FOV"""
    _section(3, "ROZMIAR KĄTOWY – kolumna 'size' [arcmin]")
    col = df["size"]
    present = col.notna().sum()
    missing  = col.isna().sum()

    print(f"     Obecne   : {present:,} ({present / total * 100:.1f}%)")
    print(f"     Brakujące: {missing:,} ({missing / total * 100:.1f}%)")

    if not present:
        return

    s1 = (col < 5).sum()
    s2 = ((col >= 5)   & (col < 20)).sum()
    s3 = ((col >= 20)  & (col < 60)).sum()
    s4 = ((col >= 60)  & (col < 180)).sum()
    s5 = (col >= 180).sum()

    print("\n     Rozkład rozmiaru:")
    print(f"       < 5'             : {s1:5,} ({s1 / present * 100:5.1f}%)  – bardzo małe [crop]")
    print(f"       5' – 20'         : {s2:5,} ({s2 / present * 100:5.1f}%)  – małe/średnie")
    print(f"       20' – 60'        : {s3:5,} ({s3 / present * 100:5.1f}%)  – dobre")
    print(f"       60' – 180'  (>1°): {s4:5,} ({s4 / present * 100:5.1f}%)  – duże")
    print(f"       ≥ 180'           : {s5:5,} ({s5 / present * 100:5.1f}%)  – bardzo duże")

def sec_source_catalog(df: pd.DataFrame, total: int) -> None:
    """[05] Katalog źródłowy (catalog)"""
    if "catalog" not in df.columns:
        return
    _section(5, "KATALOG ŹRÓDŁOWY – kolumna 'catalog'")
    for cat, cnt in df["catalog"].value_counts().items():
        print(f"     {cat:<14} {cnt:6,} ({cnt / total * 100:5.1f}%)")

# ═══════════════════════════════════════════════════
# BLOK B – przynależność do katalogów popularnych i nazwy własne
# ═══════════════════════════════════════════════════
def sec_named_catalogs(df: pd.DataFrame, total: int) -> None:
    """[06] Przynależność do katalogów Messier / Caldwell / Herschel"""
    _section(6, "PRZYNALEŻNOŚĆ DO KATALOGÓW – Messier / Caldwell / Herschel")
    m = (df["messier_nr"]  > 0).sum() if "messier_nr"  in df.columns else 0
    c = (df["caldwell_nr"] > 0).sum() if "caldwell_nr" in df.columns else 0
    h = (df["herschel_nr"] > 0).sum() if "herschel_nr" in df.columns else 0
    print(f"     Messier  (M1–M110)  : {m:4}  ({m / total * 100:5.1f}%)")
    print(f"     Caldwell (C1–C109)  : {c:4}  ({c / total * 100:5.1f}%)")
    print(f"     Herschel (H1–H400)  : {h:4}  ({h / total * 100:5.1f}%)")

def sec_common_names(df: pd.DataFrame, total: int) -> None:
    """[07] Popularne nazwy (common_names)"""
    if "common_names" not in df.columns:
        return
    _section(7, "POPULARNE NAZWY – kolumna 'common_names'")
    named_mask = df["common_names"].notna() & (df["common_names"].str.strip() != "")
    named_cnt  = named_mask.sum()
    print(f"     Z nazwą popularną: {named_cnt:,} ({named_cnt / total * 100:.1f}%)")
    if named_cnt:
        print("\n     Przykłady (pierwsze 10):")
        sample = df.loc[named_mask, ["id", "common_names", "type"]].head(10)
        for _, row in sample.iterrows():
            print(f"       {str(row.get('id', '?')):<15}  "
                  f"{str(row['common_names'])[:40]:<40}  [{row['type']}]")

# ═══════════════════════════════════════════════════
# BLOK C – jakość i kompletność danych
# ═══════════════════════════════════════════════════

def sec_completeness(df: pd.DataFrame, total: int) -> None:
    """[08] Macierz kompletności danych"""
    _section(8, "KOMPLETNOŚĆ DANYCH")
    _row("Pełne  (ra + dec + type + mag + size)",
         (df["ra"].notna() & df["dec"].notna() & df["type"].notna() &
          df["mag"].notna() & df["size"].notna()).sum(), total)
    _row("Pozycja + typ + rozmiar  (ra + dec + type + size)",
         (df["ra"].notna() & df["dec"].notna() &
          df["type"].notna() & df["size"].notna()).sum(), total)
    _row("Tylko pozycja  (ra + dec + type)",
         (df["ra"].notna() & df["dec"].notna() & df["type"].notna()).sum(), total)
    _row("Mag ✓  ale brak size",
         (df["mag"].notna() & df["size"].isna()).sum(), total)
    _row("Size ✓  ale brak mag",
         (df["size"].notna() & df["mag"].isna()).sum(), total)

def sec_problems(df: pd.DataFrame, total: int) -> None:
    """[09] Potencjalne problemy w danych"""
    _section(9, "POTENCJALNE PROBLEMY W DANYCH")
    _row("Obiekty > 10'  bez mag",
         ((df["size"] > 10) & df["mag"].isna()).sum(), total)
    _row("Obiekty mag < 10  bez size",
         ((df["mag"] < 10)  & df["size"].isna()).sum(), total)
    _row("Brak pola 'type'",
         df["type"].isna().sum(), total)

# ═══════════════════════════════════════════════════
# BLOK D – szczegóły per klasa i podsumowanie planistyczne
# ═══════════════════════════════════════════════════

def sec_per_class(df: pd.DataFrame, total: int) -> None:
    """[10] Statystyki per klasa obiektu"""
    _section(10, "STATYSTYKI PER KLASA OBIEKTU")

    groups = {
        "Galaktyki (G / GPair / GGroup)" : df["type"].str.startswith("G", na=False),
        "Mgławice Planetarne (PN)"        : df["type"] == "PN",
        "Mgławice Emisyjne (HII / EmN)"   : df["type"].isin(["HII", "EmN"]),
        "Mgławice Refleksyjne (NB / RfN)" : df["type"].isin(["NB", "RfN"]),
        "Ciemne Mgławice (DN / DrkN)"     : df["type"].isin(["DN", "DrkN"]),
        "Pozostałości Supernowej (SNR)"   : df["type"] == "SNR",
    }

    hdr = f"  {'Klasa':<35}  {'N':>6}  {'%':>5}  {'mag śr':>7}  {'mag med':>7}  {'size śr':>8}  {'size med':>9}"
    print(hdr)
    print("  " + "─" * (len(hdr) - 2))

    for label, mask in groups.items():
        sub = df[mask]
        n = len(sub)
        if n == 0:
            continue
        _f = lambda v: f"{v:.1f}" if pd.notna(v) else "–"
        mag_m  = _f(sub["mag"].mean())   if sub["mag"].notna().any()  else "–"
        mag_md = _f(sub["mag"].median()) if sub["mag"].notna().any()  else "–"
        sz_m   = (f"{sub['size'].mean():.1f}'"   if sub["size"].notna().any() else "–")
        sz_md  = (f"{sub['size'].median():.1f}'" if sub["size"].notna().any() else "–")
        print(f"  {label:<35}  {n:>6}  {n / total * 100:>5.1f}%"
              f"  {mag_m:>7}  {mag_md:>7}  {sz_m:>8}  {sz_md:>9}")

# ═══════════════════════════════════════════════════
# GŁÓWNA FUNKCJA
# ═══════════════════════════════════════════════════

def analyze_catalog() -> None:
    print_green(SEP_MAIN)
    print_green(f"       ANALIZA KATALOGU ASTRONOMICZNEGO  ›  {PATHS.catalog_full}")
    print_green(SEP_MAIN)

    try:
        df = pd.read_csv(CATALOG_FILE)
    except FileNotFoundError:
        print(f"\n  BŁĄD: Nie znaleziono pliku '{CATALOG_FILE}'")
        return

    total = len(df)

    # ── Blok A: statystyki poszczególnych kolumn ───────────────────
    sec_basic(df, total)            # [01] liczba obiektów, kolumny
    sec_magnitude(df, total)        # [02] mag – statystyki, rozkład jasności
    sec_size(df, total)             # [03] size – rozkład, 
    sec_completeness(df, total)     # [08] macierz kompletności

    sec_source_catalog(df, total)   # [05] catalog – podział na katalogi źródłowe

    # ── Blok B: przynależność do katalogów popularnych ───────────────
    sec_named_catalogs(df, total)   # [06] Messier / Caldwell / Herschel
    sec_common_names(df, total)     # [07] common_names – nazwy własne

    # ── Blok C: jakość i kompletność danych ─────────────────────

    # ── Blok D: szczegóły per klasa i podsumowanie planistyczne ────────────
    sec_per_class(df, total)        # [10] statystyki per klasa obiektu

    print(f"\n{SEP_MAIN}\n")


if __name__ == "__main__":
    analyze_catalog()
