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

PRIORITIES = {
    "ngc": 1, "ic": 1,
    "sh2": 2,
    "rcw": 3,
    "lbn": 4,
    "ced": 5,
    "pgc": 6,
    "barn": 7,
    "ldn": 8,
}

def _get_catalog_priority_for_obj(obj: Dict) -> int:
    """
    Niższa wartość = wyższy priorytet.
    Sprawdza id oraz extra_info pod kątem prefiksów z PRIORITIES.
    """
    oid = str(obj.get("id", "")).lower()
    extra = str(obj.get("extra_info", "")).lower()

    def has_prefix(p: str) -> bool:
        return oid.startswith(p) or p in extra

    for pref, pr in PRIORITIES.items():
        if has_prefix(pref):
            return pr
    return 100  # „reszta świata” – najniższy priorytet
        

def build_monthly_variants(
    vis_data: Dict,
    monthly_avg: pd.DataFrame,
    min_avg_q_hours: float,
    block_size: int = 36,
    per_month_capacity: int = 3,
    annual_vis_map: Dict[str, float] = None,
) -> List[MonthlyAssignment]:
    """
    Tworzy warianty A, B, C...

    Logika:
    - dla każdego obiektu mamy:
        * good_months – miesiące, gdzie najlepsza noc >= min_avg_q_hours
        * best_months_sorted – miesiące posortowane po best_q_hours malejąco
        * annual_vis – liczba nocy w roku z q_hours > minhours
    - sloty: 3 warianty (A,B,C) * per_month_capacity * 12 miesięcy.
    - krok 1: bierzemy obiekty z score > mediana ("Top"):
        sortowanie: annual_vis rosnąco, przy remisie score malejąco;
        każdy obiekt próbujemy wstawić:
           najlepszy miesiąc z good_months, przechodząc warianty A->B->C,
           jeśli wszystkie dobre miesiące pełne – obiekt odpada.
    - krok 2: obiekty z score <= mediana ("ogon"):
        wyliczamy priority z katalogów (PRIORITIES),
        sortowanie: priority rosnąco, annual_vis rosnąco, score malejąco,
        takie same zasady przypisywania jak wyżej.
    """
    if annual_vis_map is None:
        annual_vis_map = {}

    # score per obj
    scores = {obj["id"]: float(obj.get("score", 0.0)) for obj in vis_data["objects"]}
    all_objs_sorted = sorted(scores.keys(), key=lambda oid: scores[oid], reverse=True)

    # best_q_hours[obj_id][month] -> float
    best_map: Dict[str, Dict[int, float]] = {}
    for _, row in monthly_avg.iterrows():
        oid = row["id"]
        m = int(row["month"])
        v = float(row["best_q_hours"])
        best_map.setdefault(oid, {})[m] = v

    # mapa priorytetów katalogowych
    catalog_priority_map: Dict[str, int] = {}
    for obj in vis_data["objects"]:
        oid = obj["id"]
        catalog_priority_map[oid] = _get_catalog_priority_for_obj(obj)

    # pomocnicza struktura o obiektach
    def make_obj_info(oid_list: List[str]) -> Dict[str, Dict]:
        info = {}
        for oid in oid_list:
            month_map = best_map.get(oid, {})

            # miesiące z najlepszą nocą >= progu
            good_months = [m for m, v in month_map.items() if v >= min_avg_q_hours]
            n_good = len(good_months)

            # roczna widoczność (z prekomputu)
            annual_vis = float(annual_vis_map.get(oid, 0.0))

            # top miesiące wg best_q_hours
            best_months_sorted = sorted(
                month_map.items(), key=lambda x: x[1], reverse=True
            )

            info[oid] = {
                "id": oid,
                "score": scores.get(oid, 0.0),
                "good_months": good_months,
                "n_good": n_good,
                "annual_vis": annual_vis,
                "best_months_sorted": best_months_sorted,
                "catalog_priority": catalog_priority_map.get(oid, 100),
            }
        return info

    all_obj_ids = set(all_objs_sorted)
    obj_info_map = make_obj_info(list(all_obj_ids))

    # top_months do vis_data (bez zmian względem starej wersji)
    id_to_top_months = {}
    for oid, info in obj_info_map.items():
        top3 = info["best_months_sorted"][:3]
        top3_month_nums = [m for (m, v) in top3]
        if top3_month_nums:
            top3_str = ", ".join(f"{m:02d}" for m in top3_month_nums)
        else:
            top3_str = "—"
        id_to_top_months[oid] = top3_str

    for obj in vis_data["objects"]:
        oid = obj["id"]
        obj["top_months"] = id_to_top_months.get(oid, "—")

    # struktury wariantów
    variant_names = ["A", "B", "C"]
    variants_month_to_objects: Dict[str, Dict[int, List[str]]] = {
        vname: {m: [] for m in range(1, 13)} for vname in variant_names
    }
    assigned_objects: Dict[str, Tuple[str, int]] = {}

    # pomocnicze: aktualnie wolne sloty
    def get_available_slots() -> List[Tuple[str, int]]:
        slots: List[Tuple[str, int]] = []
        for vname in variant_names:
            m_to_objs = variants_month_to_objects[vname]
            for m in range(1, 13):
                if len(m_to_objs[m]) < per_month_capacity:
                    slots.append((vname, m))
        # A,B,C po kolei, w środku miesiące rosnąco
        slots.sort(key=lambda x: (variant_names.index(x[0]), x[1]))
        return slots

    # próba przypisania JEDNEGO obiektu (bez powrotu do niego później)
    def assign_single_object(oid: str) -> bool:
        if oid in assigned_objects:
            return True  # już przydzielony

        info = obj_info_map.get(oid)
        if info is None:
            return False

        good_months = info["good_months"]
        if not good_months:
            return False

        available_slots = get_available_slots()
        if not available_slots:
            return False

        month_to_free_variants: Dict[int, List[str]] = {}
        for vname, m in available_slots:
            month_to_free_variants.setdefault(m, []).append(vname)

        chosen_variant: Optional[str] = None
        chosen_month: Optional[int] = None

        # idziemy po miesiącach od najlepszego (highest best_q_hours)
        for m, _val in info["best_months_sorted"]:
            if m not in good_months:
                continue
            variants_for_month = month_to_free_variants.get(m, [])
            if not variants_for_month:
                continue
            vname = variants_for_month[0]  # A, potem B, C
            chosen_variant = vname
            chosen_month = m
            break

        if chosen_variant is None or chosen_month is None:
            return False

        variants_month_to_objects[chosen_variant][chosen_month].append(oid)
        assigned_objects[oid] = (chosen_variant, chosen_month)
        return True

    # --- KROK 1: Top (score > mediana) ---

    all_scores_list = [scores[oid] for oid in all_objs_sorted]
    median_score = np.median(all_scores_list) if all_scores_list else 0.0

    top_ids = [oid for oid in all_objs_sorted if scores.get(oid, 0.0) > median_score]
    rest_ids = [oid for oid in all_objs_sorted if scores.get(oid, 0.0) <= median_score]

    # sortowanie Top:
    #  - annual_vis rosnąco (rzadkie w roku pierwsze)
    #  - score malejąco
    top_sorted = sorted(
        [oid for oid in top_ids if oid in obj_info_map],
        key=lambda oid: (
            obj_info_map[oid]["annual_vis"],
            -obj_info_map[oid]["score"],
        ),
    )

    for oid in top_sorted:
        if not get_available_slots():
            break
        assign_single_object(oid)

    # --- KROK 2: „ogon” z priorytetami katalogów ---

    rest_sorted = sorted(
        [oid for oid in rest_ids if oid in obj_info_map],
        key=lambda oid: (
            obj_info_map[oid]["catalog_priority"],
            obj_info_map[oid]["annual_vis"],
            -obj_info_map[oid]["score"],
        ),
    )

    for oid in rest_sorted:
        if not get_available_slots():
            break
        assign_single_object(oid)

    # --- RAPORT KOŃCOWY (minimalnie zmodyfikowany, ale sens ten sam) ---

    print("\n" + "=" * 119)
    print(" RAPORT KOŃCOWY PO PRZYDZIALE DO WARIANTÓW")
    print("=" * 119)

    print(f"[INFO] Mediana score w katalogu: {median_score:.2f}")
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
    print(f"       Odrzuconych: {len(top_unassigned)}/{top_count} ({pct_unassigned:.1f}%)")
    # Szczegóły odrzuconych z Top
    if top_unassigned:
        print("\n[INFO] Odrzucone obiekty z Top (score > mediana):")
        for oid in top_unassigned:
            info = obj_info_map.get(oid, {})
            score_val = scores.get(oid, 0.0)
            best_list = info.get("best_months_sorted", [])[:3]
            if best_list:
                best_str_parts = []
                for m, v in best_list:
                    best_str_parts.append(f"{m:02d} ({v:.1f}h)")
                best_str = ", ".join(best_str_parts)
            else:
                best_str = "brak danych"
    
            print(
                f"       • {oid:<8}  Score: {score_val:5.1f} | "
                f"Najlepsze: {best_str}"
            )
            
    free_slots = get_available_slots()
    if free_slots:
        print(f"\n[INFO] Pozostało wolnych slotów: {len(free_slots)}")
        slots_by_var: Dict[str, List[int]] = {}
        for v, m in free_slots:
            slots_by_var.setdefault(v, []).append(m)
        for v in variant_names:
            ms = sorted(slots_by_var.get(v, []))
            if ms:
                ms_str = ", ".join(str(m) for m in ms)
                print(f"       Wariant {v}: miesiące [{ms_str}]")
    else:
        print("\n[INFO] Wszystkie sloty zostały wypełnione!")

    assigned_ids = list(assigned_objects.keys())
    assigned_from_top = [oid for oid in assigned_ids if scores.get(oid, 0.0) > median_score]
    assigned_from_rest = [oid for oid in assigned_ids if scores.get(oid, 0.0) <= median_score]
    print("\n[INFO] Obiekty użyte do uzupełnienia slotów z puli poniżej mediany score:")
    if assigned_from_rest:
        print(f"       Lista wszystkich ({len(assigned_from_rest)}) obiektów z tej puli:")
        # zachowaj kolejność według score malejąco, jak wcześniej
        assigned_from_rest_sorted = sorted(
            assigned_from_rest,
            key=lambda oid: scores.get(oid, 0.0),
            reverse=True,
        )
        for idx, oid in enumerate(assigned_from_rest_sorted, start=1):
            info = obj_info_map.get(oid, {})
            score_val = scores.get(oid, 0.0)
            best_list = info.get("best_months_sorted", [])[:3]
            if best_list:
                best_str_parts = []
                for m, v in best_list:
                    best_str_parts.append(f"{m:02d} ({v:.1f}h)")
                best_str = ", ".join(best_str_parts)
            else:
                best_str = "brak danych"
    
            print(
                f"           {idx:2d}. {oid:<8}  Score: {score_val:5.1f} | "
                f"Najlepsze: {best_str}"
            )
    else:
        print("       (brak obiektów z tej puli użytych do wypełnienia slotów)")
    print("\n[INFO] Obiekty użyte z puli poniżej mediany score:")
    if assigned_from_rest:
        print(f"       Łącznie: {len(assigned_from_rest)} obiektów z tej puli.")

    total_assigned = len(assigned_objects)
    total_objects = len(all_objs_sorted)
    pct_total = 100.0 * total_assigned / total_objects if total_objects > 0 else 0.0
    print(f"\n[INFO] Łącznie przypisanych obiektów: {total_assigned}/{total_objects} ({pct_total:.1f}%)")
    print("=" * 119 + "\n")

    variants: List[MonthlyAssignment] = []
    for vname in variant_names:
        month_to_objects = variants_month_to_objects[vname]
        total_assigned_variant = sum(len(objs) for objs in month_to_objects.values())
        if total_assigned_variant > 0:
            variants.append(
                MonthlyAssignment(
                    variant_name=vname,
                    month_to_objects=month_to_objects,
                )
            )
    return variants

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
        user_input = input(
            f"[USER] Podaj wielkość bloku, na którą będzie dzielona liczba wszystkich "
            f"{total_objects} obiektów [Enter = 36]: "
        )
        chosen_block_size = int(user_input) if user_input.strip() else 36
    except ValueError:
        print("[WARN] Błąd: Wpisano niepoprawną wartość. Przyjęto domyślnie 36.")
        chosen_block_size = 36
    
    try:
        user_input = input("[USER] Podaj liczbę obiektów w wariancie [Enter = 3]: ")
        p_capacity = int(user_input) if user_input.strip() else 3
    except ValueError:
        print("[INFO] Błąd: Wpisano niepoprawną wartość. Przyjęto domyślnie 3.")
        p_capacity = 3

    variants = build_monthly_variants(
        vis_data=vis,
        monthly_avg=monthly_best,      # przekazujemy DataFrame z best_q_hours
        block_size=chosen_block_size,
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

