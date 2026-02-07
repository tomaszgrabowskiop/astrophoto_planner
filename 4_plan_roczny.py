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

import astropy.units as u
from astropy.time import Time
from astropy.coordinates import SkyCoord, AltAz, EarthLocation, get_sun, get_body


# ------------------------------------------------------------
# Stałe / konfiguracja
# ------------------------------------------------------------

H_START = 14.0      # 14:00 (początek osi do rysowania)
H_RANGE_VIS = 18.0  # rysujemy od 14:00 do 08:00
N_H_SAMPLES = 150

# minimalna średnia liczba godzin w miesiącu, żeby uznać, że obiekt "istotnie" widoczny
MIN_AVG_Q_HOURS = 0.5  # możesz dostosować

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


# ------------------------------------------------------------
# Miesięczne średnie q_hours z observing_data.pkl
# ------------------------------------------------------------

def compute_monthly_avg_q_hours(
    observing_data: Dict[str, List[Dict]],
) -> pd.DataFrame:
    """
    Zwraca DataFrame:
    columns = ["id", "month", "avg_q_hours"]
    gdzie avg_q_hours to średnia q_hours w danym miesiącu (ze wszystkich nocy).
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
        return pd.DataFrame(columns=["id", "month", "avg_q_hours"])

    grouped = (
        df.groupby(["id", "month"])["q_hours"]
        .mean()
        .reset_index()
        .rename(columns={"q_hours": "avg_q_hours"})
    )
    return grouped


# ------------------------------------------------------------
# Budowa wariantów A/B/C... – czysto miesięczna dystrybucja
# ------------------------------------------------------------

def build_monthly_variants(
    vis_data: Dict,
    monthly_avg: pd.DataFrame,
    block_size: int = 36,
    per_month_capacity: int = 3,
    min_avg_q_hours: float = MIN_AVG_Q_HOURS,
) -> List[MonthlyAssignment]:
    """
    Tworzy warianty A, B, C... z przepływem obiektów i slotów:
    - obiekty dzielone na bloki po block_size wg score malejąco,
    - dla wariantu A używany jest blok_0,
    - dla wariantu B: najpierw niewykorzystane z A, potem blok_1,
    - dla wariantu C: niewykorzystane z A i B, potem blok_2,
    - wolne sloty z wcześniejszych wariantów (A, B) są dostępne jako miejsca
      do wypełnienia przy kolejnych wariantach (najpierw A, potem B, potem C),
    - po zbudowaniu maksymalnie 3 wariantów próbuje się dopchać wolne sloty
      kolejnymi obiektami spoza pierwszych 3 bloków.
    """
    # score per obj
    scores = {obj["id"]: float(obj.get("score", 0.0)) for obj in vis_data["objects"]}

    # pełna lista obiektów posortowana po score malejąco
    all_objs_sorted = sorted(scores.keys(), key=lambda oid: scores[oid], reverse=True)

    # pomocnicza struktura: avg_q_hours[obj_id][month] -> float
    avg_map: Dict[str, Dict[int, float]] = {}
    for _, row in monthly_avg.iterrows():
        oid = row["id"]
        m = int(row["month"])
        v = float(row["avg_q_hours"])
        avg_map.setdefault(oid, {})[m] = v

    # Funkcja pomocnicza: buduje strukturę info o obiektach
    def make_obj_info(oid_list: List[str]) -> Dict[str, Dict]:
        info = {}
        for oid in oid_list:
            month_map = avg_map.get(oid, {})
            # liczymy tylko miesiące powyżej progu avg_q_hours (min_avg_q_hours = minhours)
            good_months = [m for m, v in month_map.items() if v >= min_avg_q_hours]
            n_good = len(good_months)
            annual_vis = sum(month_map.values()) if month_map else 0.0
            best_months_sorted = sorted(
                month_map.items(),
                key=lambda x: x[1],
                reverse=True
            )
            info[oid] = {
                "id": oid,
                "score": scores.get(oid, 0.0),
                "good_months": good_months,
                "n_good": n_good,
                "annual_vis": annual_vis,
                "best_months_sorted": best_months_sorted,
            }
        return info

    # Podział na bloki po 36 (kandydaci dla A, B, C)
    num_blocks_primary = 3  # A, B, C
    primary_blocks: List[List[str]] = []
    start = 0
    for _ in range(num_blocks_primary):
        end = start + block_size
        block_ids = all_objs_sorted[start:end]
        primary_blocks.append(block_ids)
        start = end

    # Pozostałe obiekty (do dopchnięcia po A/B/C)
    remaining_after_primary = all_objs_sorted[start:]

    # Zbuduj info dla wszystkich obiektów
    all_obj_ids = set(all_objs_sorted)
    obj_info_map = make_obj_info(list(all_obj_ids))
    
    # Zapisz top 3 miesiące do vis_data["objects"] jako pole "top_months"
    id_to_top_months = {}
    for oid, info in obj_info_map.items():
        # best_months_sorted: [(month, avg_q_hours), ...] posortowane malejąco
        top3 = info["best_months_sorted"][:3]
        # weź same numery miesięcy, np. [3, 4, 2]
        top3_month_nums = [m for (m, v) in top3]
        if top3_month_nums:
            # zbuduj string "03, 04, 02"
            top3_str = ", ".join(f"{m:02d}" for m in top3_month_nums)
        else:
            top3_str = "—"
        id_to_top_months[oid] = top3_str

    # wstrzyknij do vis_data["objects"]
    for obj in vis_data["objects"]:
        oid = obj["id"]
        obj["top_months"] = id_to_top_months.get(oid, "—")

    # Przygotuj struktury dla wariantów
    variant_names = ["A", "B", "C"]
    variants_month_to_objects: Dict[str, Dict[int, List[str]]] = {
        vname: {m: [] for m in range(1, 13)} for vname in variant_names
    }

    # Dla śledzenia: które obiekty już zostały przypisane
    assigned_objects: Dict[str, Tuple[str, int]] = {}  # obj_id -> (variant_name, month)

    # Funkcja pomocnicza: lista dostępnych slotów (variant_name, month)
    def get_available_slots() -> List[Tuple[str, int]]:
        slots = []
        for vname in variant_names:
            m_to_objs = variants_month_to_objects[vname]
            for m in range(1, 13):
                if len(m_to_objs[m]) < per_month_capacity:
                    slots.append((vname, m))
        # Priorytet: najpierw wcześniejsze warianty (A przed B, B przed C)
        slots.sort(key=lambda x: (variant_names.index(x[0]), x[1]))
        return slots

    # Funkcja pomocnicza: spróbuj przypisać obiekty do dostępnych slotów
    def assign_objects(candidates: List[str]) -> List[str]:
        """
        candidates: lista obj_id w kolejności rozważania
        Zwraca listę obiektów, których nie udało się przypisać.
        """
        unused = []
        for oid in candidates:
            if oid in assigned_objects:
                continue  # już przypisany wcześniej

            info = obj_info_map.get(oid)
            if info is None:
                unused.append(oid)
                continue

            good_months = info["good_months"]
            if not good_months:
                # obiekt nie ma żadnego miesiąca spełniającego min_avg_q_hours
                unused.append(oid)
                continue

            # odśwież sloty przy każdym obiekcie (bo się zmieniają)
            available_slots = get_available_slots()
            if not available_slots:
                # nie ma żadnych wolnych miejsc, dalsze obiekty też nic nie dostaną
                unused.append(oid)
                continue

            # mapa: month -> lista wariantów z wolnymi slotami (A,B,C w kolejności)
            month_to_free_variants: Dict[int, List[str]] = {}
            for vname, m in available_slots:
                month_to_free_variants.setdefault(m, []).append(vname)

            chosen_variant = None
            chosen_month = None

            # przechodzimy po najlepszych miesiącach obiektu w kolejności malejącej avg_q_hours
            for m, avg_val in info["best_months_sorted"]:
                if m not in good_months:
                    continue
                variants_for_month = month_to_free_variants.get(m, [])
                if not variants_for_month:
                    continue
                # wybierz najwcześniejszy wariant (A, potem B, potem C)
                vname = variants_for_month[0]
                chosen_variant = vname
                chosen_month = m
                break

            if chosen_variant is None:
                # nie udało się znaleźć slotu w żadnym dobrym miesiącu
                unused.append(oid)
                continue

            # dokonujemy przypisania
            variants_month_to_objects[chosen_variant][chosen_month].append(oid)
            assigned_objects[oid] = (chosen_variant, chosen_month)

        return unused

    # --- Główne trzy przebiegi: A, B, C z przepływem obiektów ---

    unused_prev: List[str] = []

    for block_index, vname in enumerate(variant_names):
        block_ids = primary_blocks[block_index] if block_index < len(primary_blocks) else []
        # kandydaci: najpierw niewykorzystane z poprzednich, potem bieżący blok
        # filtr: tylko obiekty, które w ogóle istnieją
        candidates = [oid for oid in unused_prev if oid in obj_info_map] + [
            oid for oid in block_ids if oid in obj_info_map
        ]

        # posortuj kandydatów wg logiki: rzadkość -> całkowita widoczność -> score
        candidates_info = [
            obj_info_map[oid] for oid in candidates
        ]
        candidates_sorted = sorted(
            candidates_info,
            key=lambda o: (
                o["n_good"] if o["n_good"] > 0 else 9999,
                o["annual_vis"],
                -o["score"],
            ),
        )
        candidates_ordered = [o["id"] for o in candidates_sorted]

        # próbujemy przypisać w ramach aktualnej puli slotów (A..aktualny)
        unused_prev = assign_objects(candidates_ordered)

    # --- Dopchanie wolnych slotów obiektami spoza pierwszych 3 bloków ---

    if remaining_after_primary:
        # bierzemy tylko te, które nie są jeszcze przypisane i są w obj_info_map
        remaining_candidates = [
            oid for oid in remaining_after_primary
            if (oid not in assigned_objects and oid in obj_info_map)
        ]
        if remaining_candidates:
            remaining_info = [obj_info_map[oid] for oid in remaining_candidates]
            remaining_sorted = sorted(
                remaining_info,
                key=lambda o: (
                    o["n_good"] if o["n_good"] > 0 else 9999,
                    o["annual_vis"],
                    -o["score"],
                ),
            )
            remaining_ordered = [o["id"] for o in remaining_sorted]
            _unused_final = assign_objects(remaining_ordered)

    # Zbuduj listę wynikowych wariantów (tylko te, które mają cokolwiek)
    variants: List[MonthlyAssignment] = []
    for vname in variant_names:
        month_to_objects = variants_month_to_objects[vname]
        # sprawdź, czy wariant ma choć jeden obiekt
        total_assigned = sum(len(objs) for objs in month_to_objects.values())
        if total_assigned > 0:
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
):
    objs = variant.month_to_objects.get(month, [])

    nm_day = get_nm_day(year, month)

    t_start_dt = datetime(year, month, nm_day, int(H_START), 0)
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
            f"Brak obiektów w wariancie {variant.variant_name}",
            ha="center",
            va="center",
            transform=ax.transAxes,
        )
    else:
        # RA/DEC obiektów
        ra_dict = {obj["id"]: obj["ra"] for obj in vis_data["objects"]}
        dec_dict = {obj["id"]: obj["dec"] for obj in vis_data["objects"]}

        cmap = plt.get_cmap("tab20b")
        for idx, oid in enumerate(objs):
            ra = ra_dict.get(oid)
            dec = dec_dict.get(oid)
            if ra is None or dec is None:
                continue

            coord = SkyCoord(ra * u.deg, dec * u.deg)
            o_alt = coord.transform_to(altaz).alt.deg
            ax.plot(
                h_rel,
                o_alt,
                lw=2,
                color=cmap(idx % 10),
                label=oid,
            )

        ax.legend(fontsize=6, loc="upper right", ncol=2)

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
):
    num_pages = 12
    n_var = len(variants)
    
    with PdfPages(output_path) as pdf:
        # --- SPIS OBIEKTÓW ---
        generate_summary_page(pdf, vis_data, variants)
        # --- KLEJNE STRONY OBIEKTÓW ---
        for month in range(1, 13):
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
                )
            # Wymuś etykiety X na wszystkich subplotach
            for ax in axes:
                   ax.tick_params(labelbottom=True)

            page_no = month
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
        
        # Możliwości: top 3 miesiące (placeholder – policz wcześniej lub tutaj)
        # Przykład: weź z obj lub wstaw "—"
        possibilities = obj.get("top_months", "—")
        
        table_data.append([name, month_name, f"Wariant {variant}", possibilities])
    
    # Sortuj alfabetycznie po nazwie
    table_data.sort(key=lambda x: x[0].lower())
    
    # Podziel na dwie części (lewa i prawa kolumna)
    n_total = len(table_data)
    n_half = (n_total + 1) // 2  # zaokrąglenie w górę
    
    left_data = table_data[:n_half]
    right_data = table_data[n_half:]
    
    # Nagłówki kolumn
    col_labels = ["Nazwa", "Miesiąc", "Wariant", "Możliwości"]
    
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
        bbox=[0, 0, 0.48, 1],  # lewa połowa: x=0.02, szerokość=0.46
        cellLoc="left",
    )
    table_left.auto_set_font_size(False)
    table_left.set_fontsize(8)
    table_left.scale(1.1, 1.4)
    # TU ustaw grubość linii lewej tabeli
    for (row, col), cell in table_left.get_celld().items():
        cell.set_linewidth(0.1)  # np. 0.3 – cieńsze niż domyślne
    cells_left = table_left.get_celld()
    for (row, col), cell in cells_left.items():
        if col == 0 and row > 0:  # kolumna 0 = "Nazwa", pomijamy nagłówek (row 0)
            cell.get_text().set_weight("bold")
    
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
            bbox=[0.5, 0, 0.48, 1],  # prawa połowa: x=0.52, szerokość=0.46
            cellLoc="left",
        )
        table_right.auto_set_font_size(False)
        table_right.set_fontsize(8)
        table_right.scale(1.1, 1.4)
        # TU ustaw grubość linii lewej tabeli
        for (row, col), cell in table_right.get_celld().items():
            cell.set_linewidth(0.1)  # np. 0.3 – cieńsze niż domyślne
        cells_right = table_right.get_celld()
        for (row, col), cell in cells_right.items():
            if col == 0 and row > 0:
                cell.get_text().set_weight("bold")
        
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
    minhours = float(params.get("minhours", 3.0))  # na przyszłość, do filtrowania avg_q_hours

    lat = vis["location"]["lat"]
    lon = vis["location"]["lon"]
    location = EarthLocation(lat=lat * u.deg, lon=lon * u.deg)
    
    sunlimit = float(params.get("sunlimit", -12.0))

    print(f"[INFO] Rok: {year}, minalt={minalt}, minhours={minhours}")
    print(f"[INFO] Lokalizacja: lat={lat}, lon={lon}")

    monthly_avg = compute_monthly_avg_q_hours(observing_data)
    print(f"[INFO] Miesięczna macierz widoczności: {len(monthly_avg)} rekordów.")

    if monthly_avg.empty:
        print("[WARN] Brak danych miesięcznych q_hours. Przerywam.")
        return

    variants = build_monthly_variants(
        vis_data=vis,
        monthly_avg=monthly_avg,
        block_size=36,
        per_month_capacity=3,
        min_avg_q_hours=minhours,  # użyj progu z vis_data.json
    )
    print(f"[INFO] Utworzono {len(variants)} wariantów (bloki po 36 obiektów).")
    # NOWE: ZAPIS DO vis_data.json
    save_selected_to_vis_data(vis, variants, vis_json_path)
    for v in variants:
        total_assigned = sum(len(objs) for objs in v.month_to_objects.values())
        print(f"[INFO] Wariant {v.variant_name}: łącznie {total_assigned} przypisań "
              f"(max {3*12} przy capacity=3).")

    generate_monthly_pdf(
        output_pdf_path,
        year,
        vis,
        variants,
        location,
        min_alt=minalt,
        sunlimit=sunlimit,
    )
    print(f"[INFO] Zapisano PDF: {output_pdf_path}")

if __name__ == "__main__":
    main()
