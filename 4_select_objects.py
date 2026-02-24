#!/usr/bin/env python3
# =============================================================================
# 4_select_objects.py  –  Astrophotography Planner  |  Etap 4: Wybór obiektów i plan roczny
# =============================================================================
"""
Krok 4: Selekcja obiektów i optymalizacja harmonogramu (Scheduler).

Skrypt odpowiada za stworzenie rocznego planu obserwacyjnego, rozwiązując problem przydziału zasobów (nocy).

Logika działania:
1. Podział na grupy: Obiekty dzielone są wg punktacji (powyżej/poniżej mediany).
2. Algorytm przydziału:
   - Przypisuje obiekty do miesięcy, w których mają najlepsze warunki widoczności.
   - Dba o równomierne rozłożenie celów (Sloty A, B, C dla każdego miesiąca).
   - Priorytetyzuje obiekty rzadkie (widoczne krótko w roku) nad obiektami okołobiegunowymi.
3. Generowanie raportu wstępnego: Tworzy wykresy rozkładu obiektów i statystyki sukcesu planowania.

Wyjście:
- Zaktualizowany plik 'vis_data.json' (flaga 'selected' dla wybranych obiektów).
- Raport ze statystykami planu.
- PDF z miesięcznymi układami obiektów.
"""

import json
import pickle
from dataclasses import dataclass
from datetime import datetime, date, timedelta
from typing import Dict, List, Tuple, Optional

import numpy as np
import pandas as pd

from scipy.optimize import linear_sum_assignment

import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
from matplotlib.table import Cell

import astropy.units as u
from astropy.time import Time
from astropy.coordinates import SkyCoord, AltAz, EarthLocation, get_sun, get_body
import calendar
import pytz

"""
	Wyciszamy ostrzeżenia AstroPy. Dotyczą: 
		
	• ErfaWarning: "dubious year (Note X)": ERFA (silnik metryk czasowych w Astropy) 
	oznacza rok jako „wątpliwy”, gdy brakuje dokładnych danych o skokach sekundowych
	 i modelu czasu dla przyszłych lat.
	• Tried to get polar motions for times after IERS data is valid: Astropy nie ma aktualnych 
	tabel IERS (ruch bieguna, UT1–UTC), więc używa średnich 50‑letnich
	 – dokładność spada do poziomu łuku sekundowego.
	
	Możesz zakomentować poniższe linie kodu, żeby widziec ostrzeżenia. 
"""

import warnings
from astropy.utils.exceptions import AstropyWarning
from erfa import ErfaWarning

warnings.filterwarnings("ignore", category=ErfaWarning)
warnings.filterwarnings("ignore", category=AstropyWarning)

# Import shared.py
from shared import (PATHS, 
                    fmt, print_step, print_green,
                    load_vis_data, save_vis_data, 
                    UserConfig, 
                    CATALOG_PRIORITY, MONTH_NAMES, TYPE_NAMES,
                    H_START, H_END, H_RANGE, N_SAMPLES, CROSSING_SAMPLES,
                    get_engine_raw_hash, get_engine_final_hash, load_observing_data,
                    smart_truncate, build_badge)

# ------------------------------------------------------------
# Stałe / konfiguracja
# ------------------------------------------------------------

# Rozmiar strony i marginesy (A4) – w cm
CM_PER_INCH = 2.54
PAGE_W_CM = 21.0
PAGE_H_CM = 29.7

# marginesy w cm
MARGIN_LEFT_CM = 2.0
MARGIN_RIGHT_CM = 2.0
MARGIN_TOP_CM = 2.0
MARGIN_BOTTOM_CM = 2.0

# przelicz na cale (Matplotlib używa cali)
PAGE_W_IN = PAGE_W_CM / CM_PER_INCH
PAGE_H_IN = PAGE_H_CM / CM_PER_INCH

# obszar roboczy (współrzędne w jednostkach figury 0–1)
WORK_LEFT = MARGIN_LEFT_CM / PAGE_W_CM
WORK_RIGHT = 1.0 - MARGIN_RIGHT_CM / PAGE_W_CM
WORK_BOTTOM = MARGIN_BOTTOM_CM / PAGE_H_CM
WORK_TOP = 1.0 - MARGIN_TOP_CM / PAGE_H_CM

# konfiguracja tabeli
Cell.PAD = 0.01


@dataclass
class MonthlyAssignment:
    """Przypisania obiektów do miesięcy dla jednego wariantu."""
    variant_name: str                 # "A", "B", "C", ...
    month_to_objects: Dict[int, List[str]]  # month -> [obj_id, ...]

# ------------------------------------------------------------
# Określanie priorytetu na podstawie przynależności do katalogu
# ------------------------------------------------------------

def get_catalog_weight(catalog_name: str) -> int:
    """
    Zwraca wagę (priorytet) obiektu na podstawie pola 'catalog' z vis_data.json.
    Korzysta z systemowego słownika CATALOG_PRIORITY z shared.py.
    
    Wartości w CATALOG_PRIORITY to np. ngc=9, ic=8, sh2=7. 
    Mnożymy przez 10, by zachować dotychczasową skalę matematyczną (90, 80, 70...)
    używaną w algorytmie Hungarian. Domyślna waga dla nieznanych = 10.
    """
    base_priority = CATALOG_PRIORITY.get(str(catalog_name).lower(), 0)
    return (base_priority * 10) if base_priority > 0 else 10
# ------------------------------------------------------------
# ZAPIS DO vis_data.json – DODANIE FLAGI SELECTED
# ------------------------------------------------------------

def save_selected_to_vis_data(
    vis_data: Dict,
    variants: List[MonthlyAssignment],
    vis_json_path: str = "vis_data.json",
):
    """
    Dodaje do każdego obiektu w vis_data.json flagę "selected":
    {
        "variant": "A",
        "month": 3,
        "assignment_date": "2026-02-21T10:30:00"
    }
    Obiekty nieprzypisane mają "selected": None
    Zapis do pliku obsługiwany jest przez bibliotekę shared.py.
    """
    # 1) Wyczyść wszystkie 'selected'
    for obj in vis_data.get("objects", []):
        obj["selected"] = None

    # 2) Mapuj aktualne przypisania (z algorytmu)
    selected_map = {}
    for variant in variants:
        for month, obj_ids in variant.month_to_objects.items():
            for obj_id in obj_ids:
                selected_map[obj_id] = {
                    "variant": variant.variant_name,
                    "month": month,
                    "assignment_date": datetime.now().isoformat()
                }

    # 3) Wstaw 'selected' tylko dla obiektów z mapy
    for obj in vis_data.get("objects", []):
        obj_id = obj.get("id")
        if obj_id in selected_map:
            obj["selected"] = selected_map[obj_id]
    
    # Zapis z użyciem funkcji z shared.py (dbającej o indent i UTF-8)
    save_vis_data(vis_data, path=vis_json_path)

    assigned_count = len(selected_map)
    total_objects = len(vis_data.get("objects", []))
    
    print(f"[INFO] Zapisano flagi 'selected' do {assigned_count} obiektów w {vis_json_path}")
    print(f"[INFO] Statystyki: {assigned_count}/{total_objects} obiektów wybranych "
          f"(warianty {' '.join(v.variant_name for v in variants)})")
# ------------------------------------------------------------
# Noc nowiu w danym miesiącu
# ------------------------------------------------------------

def get_nm_day(year: int, month: int) -> int:
    """
    Szuka dnia w miesiącu o minimalnej separacji Słońce-Księżyc (noc nowiu).
    Proste przybliżenie: sprawdzamy dni 1..28 w południe.
    """
    best_day, min_sep = 15, 360.0
    for d in range(1, calendar.monthrange(year, month)[1] + 1):
        t = Time(datetime(year, month, d, 12, 0))
        sep = get_sun(t).separation(get_body("moon", t)).deg
        if sep < min_sep:
            min_sep, best_day = sep, d
    return best_day

def compute_night_length_for_date(
    year: int,
    month: int,
    day: int,
    location: EarthLocation,
    sun_alt_limit_deg: float,
) -> Tuple[int, int]:
    """
    Liczy długość nocy (Słońce < sun_alt_limit_deg) dla doby [00:00–24:00] danego dnia.
    Zwraca (godziny, minuty).
    """
    minutes_step = 5
    n_steps = int(24 * 60 / minutes_step)
    dt_hours = minutes_step / 60.0

    t_start = Time(datetime(year, month, day, 0, 0))
    times = t_start + np.arange(n_steps) * dt_hours * u.hour

    frame = AltAz(obstime=times, location=location)
    s_alt = get_sun(times).transform_to(frame).alt.deg

    night_mask = (s_alt < sun_alt_limit_deg)
    night_hours = np.sum(night_mask) * dt_hours

    total_minutes = int(round(night_hours * 60))
    hours = total_minutes // 60
    minutes = total_minutes % 60
    return hours, minutes

# ------------------------------------------------------------------------------------------------
# Miesięczne najlepsze noce i pełna widoczność roczna  z observing_data.pkl
# ------------------------------------------------------------------------------------------------

def compute_monthly_best_q_hours(observing_data):
    records = []
    for obj_id, day_list in observing_data.items():
        for rec in day_list:
            d = rec["day"]  
            if isinstance(d, datetime):
                d = d.date()
            m = d.month
            q_hours = float(rec.get("q_hours", 0.0))
            records.append({"id": obj_id, "month": m, "q_hours": q_hours})
    df = pd.DataFrame(records)
    if df.empty:
        return pd.DataFrame(columns=["id", "month", "best_q_hours"])
    
    grouped = (
        df.groupby(["id", "month"])["q_hours"]
        .max()                        
        .reset_index()
        .rename(columns={"q_hours": "best_q_hours"})
    )
    return grouped

def compute_yearly_annual_vis(observing_data, minhours: float) -> Dict[str, float]:
    """
    Zwraca słownik: obj_id -> annual_vis,
    gdzie annual_vis = liczba nocy w roku z q_hours > minhours.
    """
    annual_vis_map: Dict[str, float] = {}

    for obj_id, day_list in observing_data.items():
        n_nights = 0
        for rec in day_list:
            qh = float(rec.get("q_hours", 0.0))
            if qh > minhours:
                n_nights += 1
        annual_vis_map[obj_id] = float(n_nights)

    return annual_vis_map

# ---------------------------------------------------------------------
# Budowa wariantów A/B/C... – miesięczna dystrybucja
# ---------------------------------------------------------------------

def build_monthly_variants(
    vis_data: Dict,
    monthly_avg: pd.DataFrame,
    min_avg_q_hours: float,
    per_month_capacity: int = 3,
) -> List[MonthlyAssignment]:
    """
    Tworzy warianty A, B, C używając algorytmu optymalizacji globalnej (Hungarian Algorithm).
    
    LOGIKA HYBRYDOWA:
    1. Obiekty Score > Mediana (ELITA):
       - Priorytet: Score^3 * Quality.
       - Mechanizm: Dostają ogromny offset (1e9), dzięki czemu algorytm ZAWSZE wybiera je pierwsze,
         jeśli tylko pasują do slotu.
         
    2. Obiekty Score <= Mediana (RESZTA):
       - Priorytet: Prestiż Katalogu (NGC > Sh2 > ... > LDN).
       - Mechanizm: Score jest ignorowany. Waga katalogu jest mnożona przez stałą (1000),
         a jakość (czas widoczności) służy tylko do rozstrzygania remisów wewnątrz tego samego katalogu.
    """
    print(f"\n[INFO] Uruchamiam algorytm optymalizacji globalnej (Hungarian Algorithm)...")
    
    # 1. Przygotowanie danych
    scores = {obj["id"]: float(obj.get("score", 0.0)) for obj in vis_data.get("objects", [])}
    all_objs_list = list(scores.keys())
    
    # Szybki lookup do pełnych danych obiektu
    obj_data_map = {obj["id"]: obj for obj in vis_data.get("objects", [])}
    
    # Mapa: obj_id -> {month: hours}
    best_map: Dict[str, Dict[int, float]] = {}
    for _, row in monthly_avg.iterrows():
        oid = row["id"]
        m = int(row["month"])
        v = float(row["best_q_hours"])
        best_map.setdefault(oid, {})[m] = v
    
    # Filtrujemy obiekty, które w ogóle mają sens (spełniają min_avg_q_hours w co najmniej 1 miesiącu)
    valid_objects = []
    for oid in all_objs_list:
        month_map = best_map.get(oid, {})
        if any(v >= min_avg_q_hours for v in month_map.values()):
            valid_objects.append(oid)
            
    n_objects = len(valid_objects)
    n_variants = 3  # Stała liczba wariantów: A, B, C
    n_slots_per_variant = 12 * per_month_capacity
    n_total_slots = n_variants * n_slots_per_variant
    
    print(f"[INFO] Obiektów do rozplanowania: {n_objects}")
    print(f"[INFO] Dostępnych slotów: {n_total_slots} (3 warianty x 12 miesięcy x {per_month_capacity})")

    # --- OBLICZENIE MEDIANY (PRÓG PODZIAŁU LOGIKI) ---
    all_scores_vals = list(scores.values())
    median_score = np.median(all_scores_vals) if all_scores_vals else 0.0
    print(f"[INFO] Mediana Score: {median_score:.2f} (Granica logiki Elita vs Reszta)")

    # 2. Budowa macierzy kosztów
    INVALID_COST = 1e15 # Bardzo duża liczba oznaczająca "niemożliwe"
    cost_matrix = np.full((n_objects, n_total_slots), INVALID_COST, dtype=float)
    
    col_idx_to_slot = {}
    variant_names = ["A", "B", "C"]
    
    col_idx = 0
    for vname in variant_names:
        for month in range(1, 13):
            for slot in range(per_month_capacity):
                col_idx_to_slot[col_idx] = (vname, month)
                col_idx += 1
    
    TIER_1_OFFSET = 1_000_000_000.0 

    for i, oid in enumerate(valid_objects):
        month_map = best_map.get(oid, {})
        obj_score = float(scores.get(oid, 0.0))
        obj_data = obj_data_map.get(oid, {})
        
        # Znajdź globalne maksimum godzin dla tego obiektu (do obliczenia quality ratio)
        if not month_map: continue
        max_possible_hours = max(month_map.values())
        if max_possible_hours <= 0: continue

        # --- LOGIKA HYBRYDOWA ---
        is_elite = obj_score > median_score
        
        catalog_weight = 0
        if not is_elite:
            catalog_name = obj_data.get("catalog", "unknown")
            catalog_weight = get_catalog_weight(catalog_name)

        for c in range(n_total_slots):
            vname, month = col_idx_to_slot[c]
            hours = month_map.get(month, 0.0)
            
            if hours >= min_avg_q_hours:
                # Quality Ratio: 0.0 - 1.0
                quality_ratio = hours / max_possible_hours
                final_cost = 0.0
                
                if is_elite:
                    weighted_val = (obj_score ** 3) * (quality_ratio ** 1.0)
                    final_cost = -(TIER_1_OFFSET + weighted_val)
                else:
                    base_val = catalog_weight * 1000.0
                    quality_bonus = quality_ratio * 100.0 
                    final_cost = -(base_val + quality_bonus)

                cost_matrix[i, c] = final_cost

    # 3. Rozwiązanie (Hungarian Algorithm)
    try:
        row_ind, col_ind = linear_sum_assignment(cost_matrix)
    except Exception as e:
        print(f"[ERROR] Błąd optymalizacji: {e}")
        return []

    # 4. Zbieranie wyników
    variants_month_to_objects = {v: {m: [] for m in range(1, 13)} for v in variant_names}
    assigned_objects = {}  # oid -> (variant, month)
    
    for r, c in zip(row_ind, col_ind):
        cost = cost_matrix[r, c]
        if cost >= INVALID_COST / 100: continue 
            
        oid = valid_objects[r]
        vname, month = col_idx_to_slot[c]
        
        variants_month_to_objects[vname][month].append(oid)
        assigned_objects[oid] = (vname, month)

    # ========================================================================
    # RAPORT KOŃCOWY
    # ========================================================================
    
    print("\n" + "=" * 119)
    print(" RAPORT KOŃCOWY PO PRZYDZIALE DO WARIANTÓW (OPTIMAL + PRESTIGE LOGIC)")
    print("=" * 119)
    
    top_ids = [oid for oid in all_objs_list if scores.get(oid, 0) > median_score]
    top_count = len(top_ids)
    top_assigned = [oid for oid in top_ids if oid in assigned_objects]
    top_unassigned = [oid for oid in top_ids if oid not in assigned_objects]

    if top_count > 0:
        pct_assigned = 100.0 * len(top_assigned) / top_count
        pct_unassigned = 100.0 * len(top_unassigned) / top_count
    else:
        pct_assigned = pct_unassigned = 0.0

    print(f"[INFO] Top (score > mediana) – {top_count} obiektów:")
    print(f"       Przypisanych: {len(top_assigned)}/{top_count} ({pct_assigned:.1f}%)")
    print(f"       Odrzuconych:  {len(top_unassigned)}/{top_count} ({pct_unassigned:.1f}%)")

    if top_unassigned:
        print("[INFO] Odrzucone obiekty z Top (score > mediana):")
        top_unassigned.sort(key=lambda x: scores.get(x, 0), reverse=True)
        limit = 15
        to_show = top_unassigned[:limit]
        remaining = len(top_unassigned) - limit
        
        for oid in to_show:
            score_val = scores.get(oid, 0.0)
            month_map = best_map.get(oid, {})
            best_list = sorted(month_map.items(), key=lambda x: x[1], reverse=True)[:3]
            best_str = ", ".join([f"{m:02d} ({v:.1f}h)" for m, v in best_list]) if best_list else "brak danych"
            print(f"       • {oid:<8} Score: {score_val:5.1f} | Najlepsze: {best_str}")
        
        if remaining > 0:
            print(f"       ...i jeszcze {remaining} obiektów.")

    assigned_from_rest = [oid for oid in assigned_objects if scores.get(oid, 0) <= median_score]
    total_assigned_count = len(assigned_objects)
    
    if total_assigned_count == n_total_slots:
        print("[INFO] Wszystkie sloty zostały wypełnione!")
    else:
        print(f"\n[INFO] Sloty NIE zostały wypełnione w całości ({total_assigned_count}/{n_total_slots}).")
            
    print("\n[INFO] Obiekty użyte z puli poniżej mediany score (Sortowane wg Prestiżu):")
    if assigned_from_rest:
        print(f"       Łącznie: {len(assigned_from_rest)} obiektów z tej puli.")
        
        # ZAKTUALIZOWANA FUNKCJA SORTUJĄCA DO RAPORTU
        def sort_rest_key(oid_inner):
            d = obj_data_map.get(oid_inner, {})
            c_name = d.get("catalog", "unknown")
            return get_catalog_weight(c_name)
            
        assigned_from_rest.sort(key=sort_rest_key, reverse=True)
        
        limit = 15
        to_show = assigned_from_rest[:limit]
        remaining_count = len(assigned_from_rest) - limit
        
        for idx, oid in enumerate(to_show, 1):
            score_val = scores.get(oid, 0.0)
            d = obj_data_map.get(oid, {})
            w = get_catalog_weight(d.get("catalog", "unknown"))
            month_map = best_map.get(oid, {})
            best_list = sorted(month_map.items(), key=lambda x: x[1], reverse=True)[:3]
            best_str = ", ".join([f"{m:02d} ({v:.1f}h)" for m, v in best_list])
            print(f"       {idx:2d}. {oid:<12} (Waga: {w:3d}) Score: {score_val:4.1f} | Najlepsze: {best_str}")
            
        if remaining_count > 0:
            print(f"       ...i jeszcze {remaining_count} obiektów.")
    else:
        print("       Brak")

    pct_total = 100.0 * total_assigned_count / len(all_objs_list) if all_objs_list else 0.0
    print_step(f"Łącznie przypisanych obiektów: {total_assigned_count}/{len(all_objs_list)} ({pct_total:.1f}%)")
    print("=" * 119)

    final_variants = []
    for vname in variant_names:
        final_variants.append(
            MonthlyAssignment(
                variant_name=vname,
                month_to_objects=variants_month_to_objects[vname]
            )
        )
    
    # --- RAPORT KOŃCOWY (Quality & Satisfaction) ---
    print(f" PODSUMOWANIE JAKOŚCI PLANU ({n_objects} obiektów / {n_total_slots} slotów)")
    print("=" * 119)

    total_quality_ratio = 0.0
    quality_counts = 0
    month_scores = {m: [] for m in range(1, 13)}
    sacrificed_gems = []

    for oid, (vname, assigned_month) in assigned_objects.items():
        month_map = best_map.get(oid, {})
        if not month_map: continue
        
        max_possible_hours = max(month_map.values())
        actual_hours = month_map.get(assigned_month, 0.0)
        
        if max_possible_hours > 0:
            ratio = actual_hours / max_possible_hours
            total_quality_ratio += ratio
            quality_counts += 1
            
            obj_score = scores.get(oid, 0.0)
            month_scores[assigned_month].append(obj_score)
            
            if obj_score >= 80.0 and ratio < 0.7:
                best_m = max(month_map, key=month_map.get)
                best_h = month_map[best_m]
                loss_pct = (1.0 - ratio) * 100
                sacrificed_gems.append({
                    'oid': oid, 'score': obj_score,
                    'assigned_m': assigned_month, 'assigned_h': actual_hours,
                    'best_m': best_m, 'best_h': best_h,
                    'loss': loss_pct
                })

    avg_satisfaction = (total_quality_ratio / quality_counts * 100) if quality_counts > 0 else 0.0
    
    if avg_satisfaction >= 90: grade = "WYBITNA"
    elif avg_satisfaction >= 80: grade = "BARDZO DOBRA"
    elif avg_satisfaction >= 70: grade = "DOBRA"
    else: grade = "KOMPROMISOWA (duży tłok)"

    print(f"[INFO] Średnia jakość okna obserwacyjnego: {avg_satisfaction:.1f}% ({grade}).")
    print(f"       Średnia jakość okna obserwacyjnego względem najlepszej możliwej w roku: {avg_satisfaction:.0f}%.")

    print("\n[INFO] Obciążenie kalendarza (Średni Score obiektów w miesiącu):")
    print(f"       {'Miesiąc':<10} {'Śr. Score':<10} {'Liczba':<8} {'Status'}")
    print(f"       {'-'*45}")
    
    M_TH = ["Sty", "Lut", "Mar", "Kwi", "Maj", "Cze", "Lip", "Sie", "Wrz", "Paź", "Lis", "Gru"]
    
    for m in range(1, 13):
        m_scores = month_scores[m]
        avg_score = np.mean(m_scores) if m_scores else 0.0
        count = len(m_scores)
        month_name = M_TH[m-1]
        
        if avg_score >= 80: status = "🔥 ELITA (Top Obiekty)"
        elif avg_score >= 60: status = "✨ DOBRE (Solidne)"
        elif count > 0: status = "☁️  WYPEŁNIACZE (Słabsze)"
        else: status = "⚪ PUSTE"
        
        if count > 0:
            print(f"       {month_name:<10} {avg_score:>6.1f}     {count:>2d} szt.   {status}")
        else:
            print(f"       {month_name:<10} {'-':>6}     {0:>2d} szt.   {status}")

    if sacrificed_gems:
        print("\n[WARN] Kompromisy (Top Obiekty przesunięte do gorszych miesięcy):")
        sacrificed_gems.sort(key=lambda x: x['score'], reverse=True)
        for gem in sacrificed_gems:
            print(f"       • {gem['oid']:<8} (Score {gem['score']:.0f}): "
                  f"Miesiąc {gem['assigned_m']:02d} ({gem['assigned_h']:.1f}h) "
                  f"-> Zamiast {gem['best_m']:02d} ({gem['best_h']:.1f}h). "
                  f"Strata: -{gem['loss']:.0f}%")
    else:
        print("\n[INFO] Brak bolesnych kompromisów (wszystkie Top Obiekty mają dobre warunki).")

    print("=" * 119 + "\n")
    
    return final_variants
# ------------------------------------------------------------
# Rysowanie wykresów – jedna strona na miesiąc, warianty A/B/C
# ------------------------------------------------------------

def plot_month_variant(
    ax,
    year: int,
    month: int,
    variant: MonthlyAssignment,
    vis_data: Dict,
    location: EarthLocation,
    min_alt: float,
    sun_alt_limit_deg: float,
    tz,
):
    objs = variant.month_to_objects.get(month, [])

    nm_day = get_nm_day(year, month)

    t_start_dt = tz.localize(datetime(year, month, nm_day, int(H_START), 0))
    t_start = Time(t_start_dt)
    
    night_h, night_m = compute_night_length_for_date(
        year, month, nm_day, location, sun_alt_limit_deg,
    )

    # Używamy ujednoliconych stałych z shared.py
    h_rel = np.linspace(0, H_RANGE, N_SAMPLES)
    t_utc = t_start + h_rel * u.hour
    altaz = AltAz(obstime=t_utc, location=location)

    # Słońce / Księżyc
    s_alt = get_sun(t_utc).transform_to(altaz).alt.deg
    m_alt = get_body("moon", t_utc).transform_to(altaz).alt.deg

    # Tło wg wysokości Słońca
    for j in range(len(h_rel) - 1):
        s = s_alt[j]
        if s > 0:
            c = "#78909C"
        elif s > -6:
            c = "#90A4AE"
        elif s > -12:
            c = "#B0BEC5"
        elif s > -18:
            c = "#CFD8DC"
        else:
            c = "#ECEFF1"
        ax.axvspan(h_rel[j], h_rel[j + 1], color=c, lw=0)

    # Księżyc
    ax.plot(h_rel, m_alt, color="#003333", lw=1, ls="--", alpha=0.7)

    # Jeśli nie ma obiektów – tylko napis, bez wyłączania osi
    if not objs:
        ax.text(
            0.5, 0.5,
            f"[Brak obiektów w wariancie {variant.variant_name}]",
            ha="center", va="center", transform=ax.transAxes,
        )
    else:
        # MAPOWANIE: Tworzymy jedną mapę dostępnych obiektów
        obj_map = {obj["id"]: obj for obj in vis_data.get("objects", [])}
        cmap = plt.get_cmap("tab20b")

        for idx, oid in enumerate(objs):
            obj = obj_map.get(oid)
            if not obj: continue

            ra, dec = obj.get("ra"), obj.get("dec")
            if ra is None or dec is None: continue

            coord = SkyCoord(ra * u.deg, dec * u.deg)
            o_alt = coord.transform_to(altaz).alt.deg
            
            # 1. Budujemy dynamicznie tekst "sławy" z intów
            famous_labels = []
            if obj.get("messier_nr", 0) > 0:  famous_labels.append(f"M{obj['messier_nr']}")
            if obj.get("caldwell_nr", 0) > 0: famous_labels.append(f"C{obj['caldwell_nr']}")
            if obj.get("herschel_nr", 0) > 0: famous_labels.append(f"H{obj['herschel_nr']}")
            
            famous_str = f" ({', '.join(famous_labels)})" if famous_labels else ""
            
            # 2. Nazwa zwyczajowa
            cname = obj.get("common_names", "")
            short_cname = f" [{smart_truncate(cname, 16)}]" if cname else ""
            
            # 3. Złożenie kompletnej etykiety do legendy (np. "NGC 7000 (M33) [North Americ...]")
            label = f"{oid}{famous_str}{short_cname}"
            
            ax.plot(
                h_rel, o_alt,
                lw=2, color=cmap(idx % 10), label=label,
            )
            
        if objs:
            ax.legend(fontsize=6, loc="upper right")

    # Linia minimalnej wysokości
    ax.axhline(min_alt, color="#880E4F", ls=":", lw=1)

    # Formatowanie osi
    ax.set_ylim(0, 90)
    ax.set_xlim(0, H_RANGE)
    ax.set_xticks(np.arange(0, H_RANGE, 1))
    
    # Etykiety osi X (godziny)
    ax.set_xticklabels([
        (t_start_dt + timedelta(hours=float(h))).strftime("%-H")
        for h in np.arange(0, H_RANGE, 1)
    ])
    ax.axhline(y=0.0, color="black", lw=0.5)
    ax.set_ylabel("Wysokość [deg°]")

    ax.set_title(
        f"Wariant {variant.variant_name}",
        loc="left", fontsize=7,
    )

def generate_monthly_pdf(
    output_path: str,
    year: int,
    vis_data: Dict,
    variants: List[MonthlyAssignment],
    location: EarthLocation,
    min_alt: float,
    sun_limit: float,
    tz,
    starting_page: int = 1,
):
    num_pages = 12
    n_var = len(variants)
    
    with PdfPages(output_path) as pdf:
        # --- SPIS OBIEKTÓW ---
        generate_summary_page(pdf, vis_data, variants)
        # --- KOLEJNE STRONY OBIEKTÓW ---
        for page_idx, month in enumerate(range(1, 13), start=0):
            # --- TU liczymy rzeczy zależne od miesiąca ---
            nm_day = get_nm_day(year, month)
            night_h, night_m = compute_night_length_for_date(
                year,
                month,
                nm_day,
                location,
                sun_limit,
            )
            plt.rc("xtick", labelsize=7)
            plt.rc("ytick", labelsize=7)

            fig, axes = plt.subplots(
                n_var, 1,
                figsize=(PAGE_W_IN, PAGE_H_IN),
                sharex=True,
            )
            if n_var == 1:
                axes = [axes]

            fig.subplots_adjust(
                left=WORK_LEFT,
                right=WORK_RIGHT,
                bottom=WORK_BOTTOM,
                top=0.84,
                hspace=0.2,
            )

            # tytuł / nagłówek strony
            month_name = MONTH_NAMES[month]
            title_line1 = f"Wysokość obiektów – {month_name} {year}"
            title_line2 = f"Noc nowiu {nm_day:02d}/{month:02d}"
            title_line3 = f"Długość nocy: {night_h}h {night_m:02d}m"
            
            fig.text(0.5, 0.97, title_line1,
                     ha="center", va="top", fontsize=16, weight="bold")
            
            fig.text(0.5, 0.94, title_line2,
                     ha="center", va="top", fontsize=12)
            
            fig.text(0.5, 0.92, title_line3,
                     ha="center", va="top", fontsize=12)

            for ax, variant in zip(axes, variants):
                plot_month_variant(
                    ax,
                    year,
                    month,
                    variant,
                    vis_data,
                    location,
                    min_alt=min_alt,
                    sun_alt_limit_deg=sun_limit,
                    tz=tz,
                )
            # Wymuś etykiety X na wszystkich subplotach
            for ax in axes:
                   ax.tick_params(labelbottom=True)

            page_no = starting_page + page_idx
            fig.text(
                   0.5, 0.02, f"{page_no}",
                   ha="center", va="center",
                   fontsize=10, color="gray",
               )

            pdf.savefig(fig)
            plt.close(fig)
            
def generate_summary_page(pdf, vis_data: Dict, variants: List[MonthlyAssignment]):
    """
    Generuje stronę ze spisem w dwóch kolumnach:
    lewa połowa strony, prawa połowa strony.
    Dynamicznie odtwarza flagi sławy (M, C, H) na podstawie pól numerycznych z JSON.
    """
    # Zbierz dane
    table_data = []
    
    for obj in vis_data.get("objects", []):
        sel = obj.get("selected")
        if not sel:
            continue
        
        # 1. Pobranie podstawowych nazw z JSON
        obj_id = str(obj.get("id", "Unknown")).strip()
        c_names = str(obj.get("common_names", "")).strip()
        
        # 2. Dynamiczna rekonstrukcja "odznak sławy" z pól numerycznych
        badge_str = build_badge(obj)
       
        # 3. Budowa pełnej nazwy wyświetlanej w tabeli (zabezpieczenie długości)
        if c_names:
            # Skracamy nazwę zwyczajową do ~22 znaków, żeby nie rozwaliła wąskiej kolumny PDF
            short_cname = smart_truncate(c_names, 32)
            display_name = f"{obj_id} {short_cname}"
        else:
            display_name = f"{obj_id}{badge_str}"
        
        variant = sel.get("variant", "?")
        month_chosen = sel.get("month", 0)
        month_name = f"{month_chosen:02d}"
        table_data.append([display_name, month_name, f"{variant}"])
    
    # Sortuj alfabetycznie po utworzonej nazwie (gwarantuje ładne ułożenie katalogów)
    table_data.sort(key=lambda x: x[0].lower())
    
    # Podziel na dwie części (lewa i prawa kolumna)
    n_total = len(table_data)
    n_half = (n_total + 1) // 2  # zaokrąglenie w górę
    
    left_data = table_data[:n_half]
    right_data = table_data[n_half:]
    
    # Nagłówki kolumn
    col_labels = ["Nazwa", "M", "W"]
    
    # Stwórz figurę
    fig, ax = plt.subplots(figsize=(PAGE_W_IN, PAGE_H_IN))
    ax.axis("off")
    
    # Tytuł strony
    fig.text(0.5, 0.97, "Spis obiektów", 
             ha="center", va="top", fontsize=16, weight="bold")
    
    # ==========================================
    # LEWA TABELA
    # ==========================================
    table_left = ax.table(
        cellText=left_data,
        colLabels=col_labels,
        loc="center",
        bbox=[0, 0, 0.49, 1],  # lewa połowa: x=0, szerokość=0.49
        colWidths=[0.41, 0.04, 0.04],
        cellLoc="left",
    )
    table_left.auto_set_font_size(False)
    table_left.set_fontsize(8)
    table_left.scale(1.1, 1.4)
    
    # Formatowanie komórek lewej tabeli
    for (row, col), cell in table_left.get_celld().items():
        cell.set_linewidth(0)  # Ukryte linie obramowania
        if col in [1, 2]:  # Kolumna "Wariant" i "Miesiąc"
            cell.set_text_props(ha="center")
            
        # Kolor nagłówka
        if row == 0:
            cell.set_facecolor("#90A4AE")
            cell.set_text_props(weight="bold", color="white", size=8, ha="center")
    
    # ==========================================
    # PRAWA TABELA (tylko jeśli są dane)
    # ==========================================
    if right_data:
        table_right = ax.table(
            cellText=right_data,
            colLabels=col_labels,
            loc="center",
            bbox=[0.5, 0, 0.49, 1],  # prawa połowa: x=0.5, szerokość=0.49
            colWidths=[0.41, 0.04, 0.04],
            cellLoc="left",
        )
        table_right.auto_set_font_size(False)
        table_right.set_fontsize(8)
        table_right.scale(1.1, 1.4)
        
        # Formatowanie komórek prawej tabeli
        for (row, col), cell in table_right.get_celld().items():
            cell.set_linewidth(0)  # Ukryte linie obramowania
            if col in [1, 2]:  # Kolumna "Wariant" i „Miesiąc”
                cell.set_text_props(ha="center")
                
            # Kolor nagłówka
            if row == 0:
                cell.set_facecolor("#90A4AE")
                cell.set_text_props(weight="bold", color="white", size=8, ha="center")
    
    pdf.savefig(fig)
    plt.close(fig)

# ------------------------------------------------------------
# Główna funkcja
# ------------------------------------------------------------

def main():
    vis_json_path = PATHS.vis_data
    observing_pkl_path = PATHS.observing_final
    output_pdf_path: Optional[str] = None

    # 1. Wczytanie danych przez współdzielone funkcje z shared.py
    vis = load_vis_data(vis_json_path)
    observing_data = load_observing_data(observing_pkl_path)

    year = vis["year"]
    
    # Dynamiczne ustalenie nazwy PDF, jeśli nie została podana twardo
    if output_pdf_path is None:
        output_pdf_path = f"Astrophotography_Planner_{year}_1.pdf"

    # 2. Wyciągnięcie parametrów konfiguracyjnych (zgodne z nazewnictwem z shared.py)
    params = vis.get("parameters", {})
    min_alt = float(params.get("min_alt", 25.0))
    min_hours = float(params.get("min_hours", 3.0))
    sun_limit = float(params.get("sun_limit", -12.0))

    # Konfiguracja lokalizacji
    lat = vis["location"]["lat"]
    lon = vis["location"]["lon"]
    location = EarthLocation(lat=lat * u.deg, lon=lon * u.deg)
    city = vis.get("location", {}).get("name", "nieznana")
    tz_name = vis.get("location", {}).get("tz", "Europe/Warsaw")
    tz = pytz.timezone(tz_name)

    print_green("\n" + "=" * 119)
    print_green("       OPTYMALIZACJA MIESIĘCZNA I GENEROWANIE PLANERA")
    print_green("=" * 119)
    print(f"[INFO] Rok obserwacji: {year}")
    print(f"[INFO] Lokalizacja: {city} ({lat:.2f}°N, {lon:.2f}°E)")
    print(f"[INFO] Strefa czasowa: {tz_name}")
    print(f"[INFO] Minimalna wysokość nad horyzontem: {min_alt}°")
    print(f"[INFO] Wymagane okno widoczności: {min_hours}h")
    print(f"[INFO] Limit ciemności nieba (Słońce): {sun_limit}°")
    print("=" * 119)
    # 3. Przygotowanie statystyk z observing_data.pkl
    # Używamy najlepszej nocy w danym miesiącu zamiast średniej (best > avg)
    monthly_best = compute_monthly_best_q_hours(observing_data)
    print(f"[INFO] Miesięczna macierz najlepszych nocy: {len(monthly_best)} rekordów.")
    
    annual_vis_map = compute_yearly_annual_vis(observing_data, min_hours)
    print(f"[INFO] Macierz widoczności rocznej: {len(annual_vis_map)} zliczonych obiektów.")

    if monthly_best.empty:
        print("[WARN] Brak danych miesięcznych najlepszych nocy (q_hours). Przerywam.")
        return

    total_objects = len(vis.get("objects", []))
    print(f"[INFO] Całkowita pula obiektów kandydujących z JSON: {total_objects}")
        
    # 4. Zapytanie do użytkownika
    try:
        user_input = input("\n[USER] Podaj liczbę obiektów w pojedynczym wariancie na miesiąc [Enter = 3]: ")
        p_capacity = int(user_input) if user_input.strip() else 3
    except ValueError:
        print("[WARN] Wpisano niepoprawną wartość. Przyjęto domyślnie 3.")
        p_capacity = 3

    # 5. Algorytm optymalizacji (Hungarian) 
    variants = build_monthly_variants(
        vis_data=vis,
        monthly_avg=monthly_best,      # DataFrame z best_q_hours
        min_avg_q_hours=min_hours,     # Próg na best_q_hours w miesiącu
        per_month_capacity=p_capacity  # Limit slotów w wariancie per miesiąc
    )

    obj_in_month = 3 * p_capacity  # Zakładamy zawsze 3 warianty (A, B, C)
    obj_in_year = obj_in_month * 12
    print(
        f"\n[INFO] Zakończono tworzenie wariantów (max. {obj_in_month} DSO w każdym miesiącu w sumie dla wszystkich wariantów {obj_in_year} obiektów)."
    )

    # 6. Zapis wygenerowanych wariantów do pliku JSON
    save_selected_to_vis_data(vis, variants, vis_json_path)
    
    print_step("Rozkład obiektów w utworzonych wariantach:")
    for v in variants:
        total_assigned = sum(len(objs) for objs in v.month_to_objects.values())
        print(f"       Wariant {v.variant_name}: łącznie przypisano {total_assigned} obiekty.")

    # 7. Renderowanie i zapis raportu PDF
    print("Generowanie PDF.")
    generate_monthly_pdf(
        output_path=output_pdf_path,
        year=year,
        vis_data=vis,
        variants=variants,
        location=location,
        min_alt=min_alt,
        sun_limit=sun_limit,
        tz=tz,
        starting_page=4,
    )
    print_green("=" * 119)
    print_step(f"Planer roczny zapisany jako: {output_pdf_path}")



if __name__ == "__main__":
    main()

