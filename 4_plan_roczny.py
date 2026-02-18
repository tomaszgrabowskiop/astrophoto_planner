#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
Monthly overview planner:
- Liczy średnią miesięczną widoczność obiektów (q_hours) z observing_data.pkl.
- Dzieli obiekty na bloki po 36 (warianty A, B, C, ...).
- W każdym bloku rozdziela obiekty na 12 miesięcy (maks. 3 obiekty na miesiąc).
- Generuje PDF: dla każdego miesiąca wykres(y) wysokości obiektów
  podczas nocy nowiu, osobno dla wariantów A, B, C.
- zapisuje wybór do vis_data.json: dodaje flagę "selected" z wariantem/miesiącem.

Wejścia:
- vis_data.json  (katalog obiektów z polami id, ra, dec, score)
- observing_data.pkl  (słownik: obj_id -> lista rekordów per dzień z q_hours)

Wyjście:
- monthly_overview.pdf
- zmodyfikowany vis_data.json (dodana flaga selected)
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

# ------------------------------------------------------------
# Stałe / konfiguracja
# ------------------------------------------------------------

H_START = 14.0      # 14:00 (początek osi do rysowania)
H_RANGE_VIS = 18.0  # rysujemy od 14:00 do 08:00
N_H_SAMPLES = 150

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

YEAR = 2026
MONTH_NAMES_PL = {
    1: "Styczeń",
    2: "Luty",
    3: "Marzec",
    4: "Kwiecień",
    5: "Maj",
    6: "Czerwiec",
    7: "Lipiec",
    8: "Sierpień",
    9: "Wrzesień",
    10: "Październik",
    11: "Listopad",
    12: "Grudzień",
}

@dataclass
class MonthlyAssignment:
    """Przypisania obiektów do miesięcy dla jednego wariantu (bloku 36 obiektów)."""
    variant_name: str                 # "A", "B", "C", ...
    month_to_objects: Dict[int, List[str]]  # month -> [obj_id, ...]


# ------------------------------------------------------------
# Ładowanie danych
# ------------------------------------------------------------

def load_vis_data(path: str = "vis_data.json") -> Dict:
    with open(path, "r", encoding="utf-8") as f:
        return json.load(f)


def load_observing_data(path: str = "observing_data.pkl") -> Dict[str, List[Dict]]:
    with open(path, "rb") as f:
        return pickle.load(f)

# ------------------------------------------------------------
# Post-processing ID: dopisywanie M/C/H w nawiasie dla selected
# ------------------------------------------------------------
def _extract_famous_labels(extra_info: str) -> List[str]:
    """
    Zwraca listę etykiet typu 'M33', 'C20', 'H400' znalezionych w extra_info.
    Zakłada, że identyfikatory są rozdzielone przecinkami.
    """
    if not extra_info:
        return []

    labels: List[str] = []
    tokens = [t.strip() for t in str(extra_info).split(",") if t.strip()]

    for token in tokens:
        t = token.strip().upper()
        if len(t) < 2:
            continue

        prefix = t[0]
        rest = t[1:].strip()

        if prefix not in ("M", "C", "H"):
            continue
        if not rest.isdigit():
            continue

        # normalizujemy: "M 33" -> "M33"
        label = f"{prefix}{int(rest)}"
        labels.append(label)

    # usuwamy duplikaty, zachowując kolejność
    seen = set()
    unique_labels = []
    for lab in labels:
        if lab not in seen:
            seen.add(lab)
            unique_labels.append(lab)
    return unique_labels

def append_famous_labels_to_ids(vis_data: Dict) -> None:
    """
    Modyfikuje vis_data IN-PLACE:
    dla obiektów z selected != None wyciąga z extra_info M/C/H
    i dopisuje je w nawiasie TYLKO do pola 'name', zostawiając 'id' bez zmian, np.:

        id:   "NGC 7000"
        name: "NGC 7000 (M33, C20)"

    Jeśli nawias już jest w name, dodaje po przecinku.
    """
    for obj in vis_data.get("objects", []):
        selected = obj.get("selected")
        if not selected:
            continue  # pomijamy niewybrane

        extra_info = obj.get("extra_info", "")
        labels = _extract_famous_labels(extra_info)
        if not labels:
            continue

        base_id = str(obj.get("id", "")).strip()
        if not base_id:
            continue

        # Pracujemy tylko na name, id zostaje nietknięte
        current_name = str(obj.get("name", base_id)).strip()
        
        # jeśli name już ma nawias, dopisujemy wewnątrz
        if "(" in current_name and current_name.endswith(")"):
            main_part, rest = current_name.split("(", 1)
            rest_inner = rest.rstrip(")")
            existing = [x.strip() for x in rest_inner.split(",") if x.strip()]
            # dodaj nowe, unikając duplikatów
            for lab in labels:
                if lab not in existing:
                    existing.append(lab)
            new_inner = ", ".join(existing)
            new_name = f"{main_part.strip()} ({new_inner})"
        else:
            new_name = f"{current_name} ({', '.join(labels)})"

        obj["name"] = new_name

# ------------------------------------------------------------
# Określanie priorytetu na podstawie common names
# ------------------------------------------------------------

def get_catalog_weight(obj_id: str, common_names_str: str) -> int:
    """
    Zwraca wagę (priorytet) obiektu na podstawie przynależności do katalogów.
    Skanuje zarówno ID, jak i common_names.
    
    Zwraca wagę rosnącą (im wyższa liczba, tym ważniejszy obiekt):
    NGC/IC = 90
    Sh2    = 80
    RCW    = 70
    LBN    = 60
    CED    = 50
    PGC    = 40
    BARN   = 30
    LDN    = 20
    Inne   = 10
    """
    # Hierarchia (odwrócona na wagi dla łatwiejszej matematyki)
    PRIORITY_MAP = {
        'ngc': 90, 'ic': 90,
        'sh2': 80,
        'rcw': 70,
        'lbn': 60,
        'ced': 50,
        'pgc': 40,
        'barn': 30,
        'ldn': 20
    }
    
    # Normalizacja tekstów do szukania
    candidates = set()
    # 1. Dodaj główne ID
    candidates.add(str(obj_id).lower().strip())
    # 2. Rozbij common_names (zakładamy separator przecinek, ale też spacje)
    if common_names_str:
        # Usuwamy nawiasy jeśli są, dzielimy po przecinkach
        clean_cn = str(common_names_str).replace('(', '').replace(')', '')
        tokens = [t.strip().lower() for t in clean_cn.split(',')]
        candidates.update(tokens)
        
    max_weight = 10 # Domyślna waga dla "innych"
    
    for cand in candidates:
        # Sprawdzamy od najsilniejszych
        for key, weight in PRIORITY_MAP.items():
            # Sprawdzenie: czy 'ngc' jest w 'ngc 7000' lub czy 'sh2' jest w 'sh2-155'
            # Używamy startswith lub in, zależnie jak sformatowane są nazwy
            # Bezpieczniej: sprawdźmy czy token zaczyna się od klucza
            # Usuwamy spacje z kandydata żeby złapać "sh2-155" i "sh 2 155"
            cand_nospace = cand.replace(" ", "").replace("-", "")
            if cand_nospace.startswith(key):
                if weight > max_weight:
                    max_weight = weight
    
    return max_weight

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
        "assignment_date": "2026-02-05T16:07:00"
    }
    Obiekty nieprzypisane mają "selected": null
    """
    # 1) wyczyść wszystkie selected
    for obj in vis_data["objects"]:
        obj["selected"] = None

    # 2) mapuj aktualne przypisania
    selected_map = {}
    for variant in variants:
        for month, obj_ids in variant.month_to_objects.items():
            for obj_id in obj_ids:
                selected_map[obj_id] = {
                    "variant": variant.variant_name,
                    "month": month,
                    "assignment_date": datetime.now().isoformat()
                }

    # 3) wstaw selected tylko dla tych z mapy
    for obj in vis_data["objects"]:
        obj_id = obj["id"]
        if obj_id in selected_map:
            obj["selected"] = selected_map[obj_id]
    
    append_famous_labels_to_ids(vis_data)
    
    with open(vis_json_path, "w", encoding="utf-8") as f:
        json.dump(vis_data, f, indent=2, ensure_ascii=False)

    print(f"[INFO] Zapisano flagi 'selected' do {len(selected_map)} obiektów w {vis_json_path}")

    total_objects = len(vis_data["objects"])
    assigned_count = len(selected_map)
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
    for d in range(1, 29):
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
    """
    Zwraca DataFrame: ["id", "month", "best_q_hours"]
    gdzie best_q_hours to NAJLEPSZA noc (max q_hours) w danym miesiącu.
    """
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
    block_size: int = 36,  # Argument zachowany dla kompatybilności, ale IGNOROWANY
    per_month_capacity: int = 3,
    annual_vis_map: Dict[str, float] = None,
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
    if annual_vis_map is None:
        annual_vis_map = {}

    print(f"\n[INFO] Uruchamiam algorytm optymalizacji globalnej (Hungarian Algorithm)...")
    
    # 1. Przygotowanie danych
    scores = {obj["id"]: float(obj.get("score", 0.0)) for obj in vis_data["objects"]}
    all_objs_list = list(scores.keys())
    
    # Szybki lookup do pełnych danych obiektu (żeby wyciągnąć common_names bez pętli)
    obj_data_map = {obj["id"]: obj for obj in vis_data["objects"]}
    
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
    # Wiersze: obiekty
    # Kolumny: sloty (wariant A m1..m12, wariant B m1..m12, ...)
    INVALID_COST = 1e15 # Bardzo duża liczba oznaczająca "niemożliwe"
    
    # Inicjalizacja dużą wartością dodatnią (algorytm szuka minimum, my będziemy wstawiać ujemne)
    cost_matrix = np.full((n_objects, n_total_slots), INVALID_COST, dtype=float)
    
    col_idx_to_slot = {}
    variant_names = ["A", "B", "C"]
    
    col_idx = 0
    for vname in variant_names:
        for month in range(1, 13):
            for slot in range(per_month_capacity):
                col_idx_to_slot[col_idx] = (vname, month)
                col_idx += 1
    
    # Stała separacji - zapewnia, że najgorszy obiekt z Elity jest lepszy niż najlepszy z Reszty
    # Elita będzie miała wartości rzędu -1,000,000,000
    # Reszta będzie miała wartości rzędu -90,000
    TIER_1_OFFSET = 1_000_000_000.0 

    for i, oid in enumerate(valid_objects):
        month_map = best_map.get(oid, {})
        obj_score = float(scores.get(oid, 0.0))
        
        # Pobieramy common_names z mapy
        obj_data = obj_data_map.get(oid, {})
        common_names = obj_data.get("common_names", "")
        
        # Znajdź globalne maksimum godzin dla tego obiektu (do obliczenia quality ratio)
        if not month_map: continue
        max_possible_hours = max(month_map.values())
        if max_possible_hours <= 0: continue

        # --- LOGIKA HYBRYDOWA ---
        is_elite = obj_score > median_score
        
        # Jeśli obiekt jest z "ogona", policz jego wagę katalogową RAZ
        catalog_weight = 0
        if not is_elite:
            catalog_weight = get_catalog_weight(oid, common_names)

        for c in range(n_total_slots):
            vname, month = col_idx_to_slot[c]
            hours = month_map.get(month, 0.0)
            
            if hours >= min_avg_q_hours:
                # Quality Ratio: 0.0 - 1.0
                quality_ratio = hours / max_possible_hours
                
                final_cost = 0.0
                
                if is_elite:
                    # === LOGIKA A: ELITA (Score > Mediana) ===
                    # Stara logika: Score^3 * Quality
                    # Dodajemy TIER_1_OFFSET, żeby "przebić" punktowo wszystko z grupy niższej.
                    # Znak minus, bo solver szuka MINIMUM.
                    
                    weighted_val = (obj_score ** 3) * (quality_ratio ** 1.0)
                    final_cost = -(TIER_1_OFFSET + weighted_val)
                    
                else:
                    # === LOGIKA B: RESZTA (Score <= Mediana) ===
                    # Logika: Prestiż Katalogu.
                    # Wzór: Waga Katalogu (np. 90) * 1000 + Quality * 100
                    # Max wynik stąd to ok. 90,100. To dużo mniej niż 1,000,000,000.
                    
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
        # Sprawdzamy, czy koszt jest "walidny" (mniejszy niż nasza inicjalizacja błędów)
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
    # ... (Tu wklejasz resztę swojej funkcji raportującej od linii: print(f"[INFO] Mediana score..."))
    
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
        # Sortujemy odrzucone po Score malejąco
        top_unassigned.sort(key=lambda x: scores.get(x, 0), reverse=True)
        for oid in top_unassigned:
            score_val = scores.get(oid, 0.0)
            month_map = best_map.get(oid, {})
            best_list = sorted(month_map.items(), key=lambda x: x[1], reverse=True)[:3]
            if best_list:
                best_str = ", ".join([f"{m:02d} ({v:.1f}h)" for m, v in best_list])
            else:
                best_str = "brak danych"
            print(f"       • {oid:<8} Score: {score_val:5.1f} | Najlepsze: {best_str}")

    assigned_from_rest = [oid for oid in assigned_objects if scores.get(oid, 0) <= median_score]
    total_assigned_count = len(assigned_objects)
    
    if total_assigned_count == n_total_slots:
        print("[INFO] Wszystkie sloty zostały wypełnione!")
    else:
        print(f"\n[INFO] Sloty NIE zostały wypełnione w całości ({total_assigned_count}/{n_total_slots}).")
            
    print("\n[INFO] Obiekty użyte z puli poniżej mediany score (Sortowane wg Prestiżu):")
    if assigned_from_rest:
        print(f"       Łącznie: {len(assigned_from_rest)} obiektów z tej puli.")
        
        # Sortowanie do wyświetlania: Najpierw waga katalogu, potem quality
        def sort_rest_key(oid):
            d = obj_data_map.get(oid, {})
            w = get_catalog_weight(oid, d.get("common_names", ""))
            return w
            
        assigned_from_rest.sort(key=sort_rest_key, reverse=True)
        
        limit = 15
        to_show = assigned_from_rest[:limit]
        remaining_count = len(assigned_from_rest) - limit
        
        for idx, oid in enumerate(to_show, 1):
            score_val = scores.get(oid, 0.0)
            
            # Pobieramy wagę żeby pokazać w raporcie
            d = obj_data_map.get(oid, {})
            w = get_catalog_weight(oid, d.get("common_names", ""))
            
            month_map = best_map.get(oid, {})
            best_list = sorted(month_map.items(), key=lambda x: x[1], reverse=True)[:3]
            best_str = ", ".join([f"{m:02d} ({v:.1f}h)" for m, v in best_list])
            
            print(f"       {idx:2d}. {oid:<12} (Waga: {w}) Score: {score_val:4.1f} | Najlepsze: {best_str}")
            
        if remaining_count > 0:
            print(f"       ...i jeszcze {remaining_count} obiektów.")
    else:
        print("       Brak")

    pct_total = 100.0 * total_assigned_count / len(all_objs_list) if all_objs_list else 0.0
    print(f"\n[INFO] Łącznie przypisanych obiektów: {total_assigned_count}/{len(all_objs_list)} ({pct_total:.1f}%)")
    print("=" * 119)

    # Budowanie obiektów wyjściowych
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
    
    MONTH_NAMES = ["Sty", "Lut", "Mar", "Kwi", "Maj", "Cze", "Lip", "Sie", "Wrz", "Paź", "Lis", "Gru"]
    
    for m in range(1, 13):
        m_scores = month_scores[m]
        avg_score = np.mean(m_scores) if m_scores else 0.0
        count = len(m_scores)
        month_name = MONTH_NAMES[m-1]
        
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
        year,
        month,
        nm_day,
        location,
        sun_alt_limit_deg,
    )

    h_rel = np.linspace(0, H_RANGE_VIS, N_H_SAMPLES)
    t_utc = t_start + h_rel * u.hour
    altaz = AltAz(obstime=t_utc, location=location)

    # Słońce / Księżyc
    s_alt = get_sun(t_utc).transform_to(altaz).alt.deg
    m_alt = get_body("moon", t_utc).transform_to(altaz).alt.deg

    # tło wg Słońca
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

    # jeśli nie ma obiektów – tylko napis, ALE bez wyłączania osi
    if not objs:
        ax.text(
            0.5,
            0.5,
            f"[INFO] Brak obiektów w wariancie {variant.variant_name}",
            ha="center",
            va="center",
            transform=ax.transAxes,
        )
    else:
        # RA/DEC obiektów
        ra_dict = {obj["id"]: obj["ra"] for obj in vis_data["objects"]}
        dec_dict = {obj["id"]: obj["dec"] for obj in vis_data["objects"]}
        common_names_dict = {obj["id"]: obj.get("common_names", "") for obj in vis_data["objects"]}

        cmap = plt.get_cmap("tab20b")
        for idx, oid in enumerate(objs):
            ra = ra_dict.get(oid)
            dec = dec_dict.get(oid)
            if ra is None or dec is None:
                continue

            coord = SkyCoord(ra * u.deg, dec * u.deg)
            o_alt = coord.transform_to(altaz).alt.deg
            cname = common_names_dict.get(oid, "")
            short_cname = f" ({cname[:12]}{'...' if len(cname) > 12 else ''})" if cname else ""
            label = oid + short_cname
            ax.plot(
                h_rel,
                o_alt,
                lw=2,
                color=cmap(idx % 10),
                label=label,
            )
        if objs:
            ax.legend(fontsize=6, loc="upper right")

    ax.axhline(min_alt, color="#880E4F", ls=":", lw=1)

    ax.set_ylim(0, 90)
    ax.set_xlim(0, H_RANGE_VIS)
    ax.set_xticks(np.arange(0, H_RANGE_VIS, 1))
    ax.set_xticklabels(
        [
            (t_start_dt + timedelta(hours=float(h))).strftime("%-H")
            for h in np.arange(0, H_RANGE_VIS, 1)
        ]
    )
    ax.axhline(y=0.0, color="black", lw=0.5)
    ax.set_ylabel("Wysokość [deg°]")

    ax.set_title(
        f"Wariant {variant.variant_name}",
        loc="left",
        fontsize=7,
    )

def generate_monthly_pdf(
    output_path: str,
    year: int,
    vis_data: Dict,
    variants: List[MonthlyAssignment],
    location: EarthLocation,
    min_alt: float,
    sunlimit: float,
    tz,
    starting_page: int = 1,
):
    num_pages = 12
    n_var = len(variants)
    
    with PdfPages(output_path) as pdf:
        # --- SPIS OBIEKTÓW ---
        generate_summary_page(pdf, vis_data, variants)
        # --- KLEJNE STRONY OBIEKTÓW ---
        for page_idx, month in enumerate(range(1, 13), start=0):
            # --- TU liczymy rzeczy zależne od miesiąca ---
            nm_day = get_nm_day(year, month)
            night_h, night_m = compute_night_length_for_date(
                year,
                month,
                nm_day,
                location,
                sunlimit,
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
            month_name = MONTH_NAMES_PL[month]
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
                    sun_alt_limit_deg=sunlimit,
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
    """
    # Zbierz dane
    table_data = []
    
    for obj in vis_data["objects"]:
        sel = obj.get("selected")
        if not sel:
            continue
        
        name = obj.get("name", obj["id"])
        
        variant = sel["variant"]
        month_chosen = sel["month"]
        month_name = MONTH_NAMES_PL.get(month_chosen, f"{month_chosen:02d}")
  
        table_data.append([name, month_name, f"{variant}"])
    
    # Sortuj alfabetycznie po nazwie
    table_data.sort(key=lambda x: x[0].lower())
    
    # Podziel na dwie części (lewa i prawa kolumna)
    n_total = len(table_data)
    n_half = (n_total + 1) // 2  # zaokrąglenie w górę
    
    left_data = table_data[:n_half]
    right_data = table_data[n_half:]
    
    # Nagłówki kolumn
    col_labels = ["Nazwa", "Miesiąc", "Wariant"]
    
    # Stwórz figurę
    fig, ax = plt.subplots(figsize=(PAGE_W_IN, PAGE_H_IN))
    ax.axis("off")
    
    # Tytuł strony
    fig.text(0.5, 0.97, "Spis obiektów", 
             ha="center", va="top", fontsize=16, weight="bold")
    
    # LEWA TABELA (bbox: [left, bottom, width, height] w jednostkach figury)
    table_left = ax.table(
        cellText=left_data,
        colLabels=col_labels,
        loc="center",
        bbox=[0, 0, 0.49, 1],  # lewa połowa: x=0.02, szerokość=0.46
        colWidths=[0.35, 0.13, 0.1],
        cellLoc="left",
    )
    table_left.auto_set_font_size(False)
    table_left.set_fontsize(8)
    table_left.scale(1.1, 1.4)
    # TU ustaw grubość linii lewej tabeli
    for (row, col), cell in table_left.get_celld().items():
        cell.set_linewidth(0)  # np. 0.3 – cieńsze niż domyślne
        if col == 2:      # pomiń nagłówek: and row > 0:
            cell.set_text_props(ha="center")
    cells_left = table_left.get_celld()
    for (row, col), cell in cells_left.items():
        if col == 0 and row > 0:  # kolumna 0 = "Nazwa", pomijamy nagłówek (row 0)
            cell.get_text()
    
    # Kolor nagłówka lewej tabeli
    for (row, col), cell in table_left.get_celld().items():
        if row == 0:
            cell.set_facecolor("#90A4AE")
            cell.set_text_props(weight="bold", color="white", size=8)
    
    # PRAWA TABELA (tylko jeśli są dane)
    if right_data:
        table_right = ax.table(
            cellText=right_data,
            colLabels=col_labels,
            loc="center",
            bbox=[0.5, 0, 0.49, 1],  # prawa połowa: x=0.52, szerokość=0.46
            colWidths=[0.35, 0.13, 0.1],
            cellLoc="left",
        )
        table_right.auto_set_font_size(False)
        table_right.set_fontsize(8)
        table_right.scale(1.1, 1.4)
        # TU ustaw grubość linii lewej tabeli
        for (row, col), cell in table_right.get_celld().items():
            cell.set_linewidth(0)  # np. 0.3 – cieńsze niż domyślne
            if col == 2:      # pomiń nagłówek: and row > 0:
                cell.set_text_props(ha="center")
        cells_right = table_right.get_celld()
        for (row, col), cell in cells_right.items():
            if col == 0 and row > 0:
                cell.get_text()
        
        # Kolor nagłówka prawej tabeli
        for (row, col), cell in table_right.get_celld().items():
            if row == 0:
                cell.set_facecolor("#90A4AE")
                cell.set_text_props(weight="bold", color="white", size=8)
    
    pdf.savefig(fig)
    plt.close(fig)

# ------------------------------------------------------------
# Główna funkcja
# ------------------------------------------------------------

def main(
    vis_json_path: str = "vis_data.json",
    observing_pkl_path: str = "observing_data.pkl",
    output_pdf_path = f"Astrophotography_Planner_{YEAR}_1.pdf",
):
    vis = load_vis_data(vis_json_path)
    observing_data = load_observing_data(observing_pkl_path)

    year = vis["year"]
    params = vis.get("parameters", {})
    minalt = float(params.get("minalt", 25.0))
    # minhours interpretujemy teraz jako próg na "najlepszą noc" w miesiącu
    minhours = float(params.get("minhours", 3.0))

    lat = vis["location"]["lat"]
    lon = vis["location"]["lon"]
    location = EarthLocation(lat=lat * u.deg, lon=lon * u.deg)
    city = vis.get("location", {}).get("name", "nieznana")
    sunlimit = float(params.get("sunlimit", -12.0))
    tz_name = vis.get("location", {}).get("tz", "Europe/Warsaw")
    tz = pytz.timezone(tz_name)
    print(f"[INFO] Rok: {year}")
    print(f"       Lokalizacja: {city} ({lat:.2f}°N, {lon:.2f}°E)")
    print(f"       Strefa czasowa: {tz_name}")
    print(f"       Minimalna wysokość nad horyzontem {minalt}°")
    print(f"       Minimalne okno obserwacji: {minhours}h")

    # UŻYWAMY BEST ZAMIAST AVG
    monthly_best = compute_monthly_best_q_hours(observing_data)
    print(f"       Miesięczna macierz najlepszych nocy: {len(monthly_best)} rekordów.")
    annual_vis_map = compute_yearly_annual_vis(observing_data, minhours)
    print(f"       Macierz widoczności rocznej: {len(annual_vis_map)} zliczonych obiektów.")

    if monthly_best.empty:
        print("[WARN] Brak danych miesięcznych najlepszych nocy (q_hours). Przerywam.")
        return

    total_objects = len(vis["objects"])
        
    try:
        user_input = input("[USER] Podaj liczbę obiektów w wariancie [Enter = 3]: ")
        p_capacity = int(user_input) if user_input.strip() else 3
    except ValueError:
        print("[INFO] Błąd: Wpisano niepoprawną wartość. Przyjęto domyślnie 3.")
        p_capacity = 3

    variants = build_monthly_variants(
        vis_data=vis,
        monthly_avg=monthly_best,      # przekazujemy DataFrame z best_q_hours
        per_month_capacity=p_capacity,
        min_avg_q_hours=minhours,      # próg dot. best_q_hours w miesiącu
        annual_vis_map=annual_vis_map,
    )

    obj_in_month = 3 * p_capacity
    obj_in_year = obj_in_month * 12
    print(
        f"[INFO] Utworzono {len(variants)} warianty "
        f"(max. {obj_in_month} DSO w każdym miesiącu, w sumie max. {obj_in_year})."
    )

    # ZAPIS WYBRANYCH MIESIĘCY/OBIEKTÓW DO vis_data.json
    save_selected_to_vis_data(vis, variants, vis_json_path)
    for v in variants:
        total_assigned = sum(len(objs) for objs in v.month_to_objects.values())
        print(f"[INFO] Wariant {v.variant_name}: łącznie {total_assigned} obiekty.")

    generate_monthly_pdf(
        output_pdf_path,
        year,
        vis,
        variants,
        location,
        min_alt=minalt,
        sunlimit=sunlimit,
        tz=tz,
        starting_page=4,
    )
    print(f"[INFO] Zapisano PDF: {output_pdf_path}")


if __name__ == "__main__":
    main()

