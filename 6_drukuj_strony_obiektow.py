import pandas as pd
import numpy as np
import pickle
from datetime import datetime, timedelta
import os
import json
import math

from astropy.coordinates import SkyCoord, AltAz, EarthLocation, get_sun, get_body
from astropy.time import Time
import astropy.units as u
import pytz

import matplotlib.pyplot as plt
import matplotlib.dates as mdates
from matplotlib.backends.backend_pdf import PdfPages
from matplotlib.gridspec import GridSpec
import matplotlib.ticker as ticker
from matplotlib.dates import DateFormatter
import textwrap
from tqdm import tqdm

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

# --- KONFIGURACJA ---

PKL_PATH = "observing_data.pkl"
JSON_PATH = "vis_data.json"

# poniższe wartości mogą zostać nadpisane danymi z JSON
LAT, LON = 52.4064, 16.9252
LOCATION = EarthLocation(lat=LAT * u.deg, lon=LON * u.deg)
YEAR = 2026

H_START = 14.0      # początek zakresu godzin (oś czasu na wykresach)
H_RANGE = 17.0      # długość zakresu godzin (np. 14:00–07:00)
TZ_OFFSET_REF = 1.0  # CET jako punkt odniesienia dla wykresu
MIN_ALTITUDE = 30   # minimalna wysokość obiektu na wykresach

MONTH_NAMES = {
    1: "STYCZEŃ",
    2: "LUTY",
    3: "MARZEC",
    4: "KWIECIEŃ",
    5: "MAJ",
    6: "CZERWIEC",
    7: "LIPIEC",
    8: "SIERPIEŃ",
    9: "WRZESIEŃ",
    10: "PAŹDZIERNIK",
    11: "LISTOPAD",
    12: "GRUDZIEŃ",
}

TYPE_NAMES = {
            "Gx": "Galaktyka",
            "G": "Galaktyka",
            "S": "Galaktyka spiralna",
            "SB": "Galaktyka spiralna z poprzeczką",
            "E": "Galaktyka eliptyczna",
            "OC": "Gromada Otwarta",
            "OCl": "Gromada Otwarta",
            "Cl+N": "Gromada + Mgławica",
            "GC": "Gromada Kulista",
            "HII": "Region HII",
            "NB": "Mgławica",
            "Neb": "Mgławica",
            "PN": "Mgławica Planetarna",
            "DN": "Ciemna Mgławica",
            "SNR": "Pozostałość Supernowej",
            "Other": "Inny",
        }

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

# --- FUNKCJE POMOCNICZE ---

def load_vis_data(path=JSON_PATH):
    with open(path, "r", encoding="utf-8") as f:
        data = json.load(f)
    return data

def apply_config_from_vis_data(vis_data):
    global LAT, LON, LOCATION, YEAR, MIN_ALTITUDE

    # Lokalizacja
    LAT = vis_data["location"]["lat"]
    LON = vis_data["location"]["lon"]
    LOCATION = EarthLocation(lat=LAT * u.deg, lon=LON * u.deg)
    
    # Rok
    YEAR = vis_data["year"]

    params = vis_data.get("parameters", {})
    MIN_ALTITUDE = params.get("minalt", MIN_ALTITUDE)

def compute_diagonal_fov_arcmin(camera):
    f = camera["lens_focal_length"]
    w = camera["sensor_width"]
    h = camera["sensor_height"]
    d = math.hypot(w, h)
    fov_rad = 2 * math.atan(d / (2 * f))
    fov_arcmin = fov_rad * 3437.75
    fov_deg = fov_arcmin / 60.0
    return fov_arcmin, fov_deg

def get_image_extent(img):
    h, w = img.shape[:2]
    return [0, w / h, 0, 1]

def format_indeksy(extra_info: str, max_items: int = 7, max_chars: int = 45) -> str:
    """Jedna linia: bierze max_items, skraca do max_chars."""
    if not extra_info:
        return 'brak'
    
    items = [item.strip() for item in extra_info.split(',') if item.strip()]
    short_items = items[:max_items]  # Tylko pierwsze max_items
    
    line = ', '.join(short_items)
    if len(line) > max_chars:
        line = line[:max_chars-3] + '...'
    
    return line


# --- RYSOWANIE STRONY OBIEKTU ---
  
def draw_object_page(pdf, oid, month, nm_day, row, all_data, camera, page_num, tz, best_year_date, best_year_m_hours,):
    days_data = all_data[oid]
    days = [d["day"] for d in days_data]
    
    obj_type_code = row.get("type", "")
    type_pl = TYPE_NAMES.get(obj_type_code, obj_type_code or "Inny")
    
    max_len = 20
    raw_name = row.get('common_names') or 'brak'
    short_name = (raw_name[:max_len - 3] + '...') if len(raw_name) > max_len else raw_name
    line = f"| Nazwa zwyczajowa: {short_name} |\n"

    indeksy_short = format_indeksy(row.get('extra_info', ''), max_items=7, max_chars=45)
        
    ra_rounded = round(float(row['ra']), 2)  
    dec_rounded = round(float(row['dec']), 2) 
    size_rounded = round(float(row['size']),2)
    mag_rounded = round(float(row['mag']),2)
    
    
    fig = plt.figure(figsize=(PAGE_W_IN, PAGE_H_IN))
    fig.subplots_adjust(
    left=WORK_LEFT, right=WORK_RIGHT,
    bottom=WORK_BOTTOM, top=WORK_TOP
    )
    gs = GridSpec(
        4,
        2,
        figure=fig,
        height_ratios=[0.6, 0.9, 1, 1],
        width_ratios=[0.95, 1.05],
        hspace=0.3,
        wspace=0.05,
    )

    # 1. Nagłówek i opis
    ax_txt = fig.add_subplot(gs[0, :])
    ax_txt.axis("off")

    header = (
        f"{oid}\n"
        f"|Nazwa: {short_name} | {type_pl} ({obj_type_code})|\n"
        )

    ax_txt.text(0, 1.0, header, fontsize=16, fontweight="bold", va="top")
    lorem = (
        f"|Typ: {row.get('type', '')} | RA: {ra_rounded:.2f} | Dec: {dec_rounded:.2f}|\n"
        f"|Rozmiar: {size_rounded:.2f}' | Mag.: {mag_rounded:.2f} |\n"
        f"|Indeksy: {indeksy_short}|"
    )
    ax_txt.text(
        0,
        0.45,
        lorem,
        fontfamily="serif",
        fontsize=12,
        linespacing=1.2,
        va="top",
    )
       
    # 2. Noc nowiu
    ax_nm = fig.add_subplot(gs[1, 0])

    t_start = tz.localize(datetime(YEAR, month, nm_day, 14, 0))
    h_rel = np.linspace(0, 18, 150)
    t_utc = Time(t_start) + h_rel * u.hour
    altaz = AltAz(obstime=t_utc, location=LOCATION)

    o_alt = SkyCoord(row["ra"] * u.deg, row["dec"] * u.deg).transform_to(altaz).alt.deg
    s_alt = get_sun(t_utc).transform_to(altaz).alt.deg
    m_alt = get_body("moon", t_utc).transform_to(altaz).alt.deg

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
        ax_nm.axvspan(h_rel[j], h_rel[j + 1], color=c, lw=0)

    ax_nm.plot(h_rel, o_alt, color="#4A148C", lw=2, label="Obiekt")
    ax_nm.plot(h_rel, m_alt, color="#003333", lw=1, ls="--", alpha=0.7)
    ax_nm.axhline(MIN_ALTITUDE, color="#880E4F", ls=":", lw=1)

    ax_nm.set_ylim(0, 90)
    ax_nm.set_xlim(0, 17)
    ax_nm.set_xticks(np.arange(0, 18, 1))
    ax_nm.set_xticklabels(
        [
            (t_start + timedelta(hours=float(h))).strftime("%-H")
            for h in np.arange(0, 18, 1)
        ]
    )
    ax_nm.tick_params(axis='both', which='major', labelsize=8)

    ax_nm.set_title(
            (f"Wysokość obiektu {nm_day:02d}.{month:02d}. Max: {best_year_date.day:02d}.{best_year_date.month:02d}/{best_year_m_hours:.1f}h"
            ),
            loc="left",
            fontsize=10,
        )
    ax_nm.set_ylabel("Wysokość [°]")

    # 3. Kadr FOV
    fov_arcmin, fov_deg = compute_diagonal_fov_arcmin(camera)
    name = row["name"]
    img_path = f"starplots/{oid}.png"
    img = plt.imread(img_path)
    extent = get_image_extent(img)
    ax_img = fig.add_subplot(gs[1, 1])
    try:
        img = plt.imread(img_path)
        ax_img.imshow(img,extent=extent, aspect='equal')
    except FileNotFoundError:
        ax_img.text(
            0.5,
            0.5,
            f"Brak: {img_path}",
            ha="center",
            va="center",
            transform=ax_img.transAxes,
            color="red",
        )
    ax_img.axis("off")
    ax_img.set_title(
        f"Kadr FOV, przekątna {fov_arcmin:.1f}' ({fov_deg:.2f}°).",
        fontsize=10,
        loc="left",
    )

    # 4. Visibility Plot
    if not isinstance(days, pd.DatetimeIndex):
        days_pd = pd.to_datetime(days)
    else:
        days_pd = days

    ax_vis = fig.add_subplot(gs[2, :])
    ax_vis.set_facecolor("#eeeeee")

    s_set, s_rise = [], []
    for d_e in days_data:
        pts = d_e["sun_pts"]
        tz_offset = d_e.get("tz_offset", TZ_OFFSET_REF)
        delta = tz_offset - TZ_OFFSET_REF  # +1h dla CEST względem CET
    
        if len(pts) >= 2:
            t0 = pts[0].hour + pts[0].minute / 60.0
            t1 = pts[1].hour + pts[1].minute / 60.0
    
            # przesunięcie o delta godzin
            t0_shifted = t0 + delta
            t1_shifted = t1 + delta
    
            v0 = (t0_shifted + 24 - H_START) if t0_shifted < 12 else (t0_shifted - H_START)
            v1 = (t1_shifted + 24 - H_START) if t1_shifted < 12 else (t1_shifted - H_START)
    
            s_set.append(max(0, v0))
            s_rise.append(min(H_RANGE, v1))
        else:
            s_set.append(0)
            s_rise.append(H_RANGE)

    ax_vis.fill_between(
        days_pd, s_set, s_rise, color="#6497B1", lw=0, alpha=0.5, zorder=2
    )

    for d_entry in days_data:
        d = d_entry["day"]
        d_pd = pd.to_datetime(d)
    
        tz_offset = d_entry.get("tz_offset", TZ_OFFSET_REF)
        delta = tz_offset - TZ_OFFSET_REF
    
        for s_rel, e_rel in d_entry["qual_segments"]:
            s_rel_shifted = s_rel + delta
            e_rel_shifted = e_rel + delta
    
            # odetnij segmenty wychodzące poza zakres wykresu
            if e_rel_shifted <= 0 or s_rel_shifted >= H_RANGE:
                continue
    
            ax_vis.fill_between(
                [d_pd, d_pd + timedelta(days=1)],
                [s_rel_shifted, s_rel_shifted],
                [e_rel_shifted, e_rel_shifted],
                color="#011F4B",
                lw=0,
                zorder=1,
            )

    ax_vis.set_ylim(0, H_RANGE)
    y_ticks = np.arange(0, H_RANGE + 0.1, 1)
    ax_vis.set_yticks(y_ticks)
    labels_y = [f"{(int(H_START) + int(h)) % 24:02d}:00" for h in y_ticks]
    ax_vis.set_yticklabels(labels_y)

    for i, label in enumerate(ax_vis.get_yticklabels()):
        label.set_fontsize(8 if i % 2 == 0 else 6)

    grid_days = [5, 10, 15, 20, 25]
    ax_vis.xaxis.set_major_locator(mdates.DayLocator(bymonthday=grid_days))
    ax_vis.xaxis.set_major_formatter(mdates.DateFormatter("%-d"))
    ax_vis.xaxis.set_minor_locator(mdates.MonthLocator())
    ax_vis.xaxis.set_minor_formatter(mdates.DateFormatter("%b"))
    ax_vis.xaxis_date()

    ax_vis.grid(
        which="major", axis="x", color="white", linewidth=0.7, alpha=0.9, zorder=3
    )
    ax_vis.grid(
        which="minor", axis="x", color="white", linewidth=1, alpha=1, zorder=3
    )

    ax_vis.tick_params(
        axis="x", which="minor", length=8, pad=5, labelsize=8, color="gray"
    )
    ax_vis.tick_params(
        axis="x", which="major", length=4, pad=2, labelsize=5, color="gray"
    )

    ax_vis.set_xlim(days_pd[0], days_pd[-1])
    ax_vis.yaxis.set_minor_locator(ticker.MultipleLocator(0.25))
    ax_vis.grid(
        which="major", axis="y", color="white", linewidth=0.7, alpha=1, zorder=3
    )
    ax_vis.grid(
        which="minor", axis="y", color="white", linewidth=0.2, alpha=0.9, zorder=3
    )
    ax_vis.set_title(f"Widoczność obiektu w ciągu roku (strefa czasowa: {tz.zone})", fontsize=10, loc="left")

    # 5. Quality Hours
    ax_hrs = fig.add_subplot(gs[3, :])
    ax_hrs.set_facecolor("#eeeeee")

    q_hrs = [d["q_hours"] for d in days_data]
    m_hrs = [d["m_hours"] for d in days_data]

    ax_hrs.set_title(
        "Suma godzin jakościowych (jasny: z Księżycem, ciemny: bez Księżyca)",
        fontsize=10,
        loc="left",
    )
    ax_hrs.fill_between(
        days_pd, q_hrs, color="#6497B1", label="Total Quality Hours", lw=0, zorder=1
    )
    ax_hrs.fill_between(
        days_pd, m_hrs, color="#011F4B", alpha=0.8, label="Moonless Hours", lw=0, zorder=1
    )

    ax_hrs.set_ylim(0, 14)
    ax_hrs.yaxis.set_major_locator(ticker.MultipleLocator(1))
    ax_hrs.yaxis.set_major_formatter(ticker.StrMethodFormatter("{x:.0f}"))
    ax_hrs.yaxis.set_minor_locator(ticker.MultipleLocator(0.5))

    for i, label in enumerate(ax_hrs.get_yticklabels()):
        label.set_fontsize(8 if i % 2 == 0 else 6)

    grid_days = [5, 10, 15, 20, 25]
    ax_hrs.xaxis.set_major_locator(mdates.DayLocator(bymonthday=grid_days))
    ax_hrs.xaxis.set_major_formatter(mdates.DateFormatter("%-d"))
    ax_hrs.xaxis.set_minor_locator(mdates.MonthLocator())
    ax_hrs.xaxis.set_minor_formatter(mdates.DateFormatter("%b"))
    ax_hrs.xaxis_date()

    ax_hrs.grid(
        which="major", axis="x", color="white", linewidth=0.7, alpha=0.9, zorder=3
    )
    ax_hrs.grid(
        which="minor", axis="x", color="white", linewidth=1, alpha=1, zorder=3
    )

    ax_hrs.tick_params(
        axis="x", which="minor", length=8, pad=5, labelsize=8, color="gray"
    )
    ax_hrs.tick_params(
        axis="x", which="major", length=4, pad=2, labelsize=5, color="gray"
    )

    ax_hrs.set_xlim(days_pd[0], days_pd[-1])
    ax_hrs.grid(
        which="major", axis="y", color="white", linewidth=0.7, alpha=1, zorder=3
    )
    ax_hrs.grid(
        which="minor", axis="y", color="white", linewidth=0.2, alpha=0.9, zorder=3
    )
    # Numer strony na dole, pośrodku
    fig.text(
        0.5,
        0.02,
        f"{page_num}",
        ha="center",
        va="center",
        fontsize=10,
        color="gray",
    )

    pdf.savefig(fig, dpi=600)
    plt.close()


# --- Strony z mapami kontekstowymi ---

def draw_context_page(pdf, oid, page_num):
    fig = plt.figure(figsize=(PAGE_W_IN, PAGE_H_IN))
    fig.subplots_adjust(
        left=WORK_LEFT, right=WORK_RIGHT,
        bottom=WORK_BOTTOM, top=WORK_TOP,
    )
    # tytuł strony u góry
    fig.text(
        0.5, 0.933,                      
        f"Mapa kontekstowa {oid}", 
        ha="left", va="top",
        fontsize=16, fontweight="bold",
    )
    
    img_path = f"starplots/context_{oid}.png"
    try:
        img = plt.imread(img_path)
        ih, iw = img.shape[0], img.shape[1]   # height, width
        aspect = iw / ih

        img_h = 0.75          # 3/4 wysokości figury
        img_w = img_h * aspect

        # wycentrowanie ramki o proporcjach obrazka
        x0 = (1.0 - img_w) / 2.0
        y0 = (1.0 - img_h) / 2.0

        ax = fig.add_axes([x0, y0, img_w, img_h])
        ax.set_axis_off()

        # kluczowe: układ danych = piksele, aspect='equal'
        ax.imshow(img, origin="upper")
        ax.set_xlim(0, iw)
        ax.set_ylim(ih, 0)    # odwrócona oś Y, żeby obraz nie był do góry nogami
        ax.set_aspect("equal")  # 1:1 w przestrzeni danych
    except FileNotFoundError:
        ax = fig.add_axes([0, 0, 1, 1])
        ax.set_axis_off()
        ax.text(
            0.5, 0.5, f"Brak: {img_path}",
            ha="center", va="center",
            transform=ax.transAxes, color="red", fontsize=16,
        )

    fig.text(
        0.5, 0.02, f"{page_num}",
        ha="center", va="center",
        fontsize=10, color="gray",
    )

    pdf.savefig(fig, dpi=600)
    plt.close(fig)

# --- MAIN ---

def main():
    print("[INFO] Inicjalizacja danych.")

    vis_data = load_vis_data(JSON_PATH)
    apply_config_from_vis_data(vis_data)
    camera = vis_data["parameters"]["camera"]
    # Strefa czasowa
    tz_name = vis_data.get("location", {}).get("tz", "Europe/Warsaw")
    tz = pytz.timezone(tz_name)
    with open(PKL_PATH, "rb") as f:
        all_data = pickle.load(f)

    objects = vis_data["objects"]
    df = pd.DataFrame(objects)

    # zostaw tylko obiekty, które mają wypełnione selected
    df = df[df["selected"].notna()]

    # sortowanie alfabetyczne po ID
    df = df.sort_values("id")
    #df = df.head(1)

    output_name = f"Astrophotography_Planner_{YEAR}_2.pdf"
    with PdfPages(output_name) as pdf:
        page_num = 13  # numer pierwszej strony obiektu

        # Tworzenie paska postępu
        pbar = tqdm(df.iterrows(), total=len(df), unit="obiekt", ncols=119)
    
        for _, row in pbar:
            oid = row["id"]
            sel = row["selected"]
            month = int(sel["month"])
            days_data = all_data[oid]
            best_month_entry = None
            for d_entry in days_data:
                d = d_entry["day"]
                if d.month != month:
                    continue
                if best_month_entry is None or d_entry["m_hours"] > best_month_entry["m_hours"]:
                    best_month_entry = d_entry
            assignment_date = datetime.fromisoformat(sel["assignment_date"])
            if best_month_entry is not None:
                best_month_date = best_month_entry["day"]
                nm_day = best_month_date.day
            else:
                best_month_date = assignment_date
                nm_day = assignment_date.day
            best_year_entry = None
            for d_entry in days_data:
                if best_year_entry is None or d_entry["m_hours"] > best_year_entry["m_hours"]:
                    best_year_entry = d_entry
            
                if best_year_entry is not None:
                    best_year_date = best_year_entry["day"]
                    best_year_m_hours = float(best_year_entry["m_hours"])
                else:
                    best_year_date = best_month_date
                    best_year_m_hours = 0.0
            
            assignment_date = datetime.fromisoformat(sel["assignment_date"])
            nm_day = assignment_date.day
    
            # KROK 1: Aktualizacja opisu paska i generowanie sekcji
            pbar.set_description(f" Generowanie FOV: {oid}")
            draw_object_page(pdf, oid, month, nm_day,  row, all_data, camera, page_num,  tz,  best_year_date, best_year_m_hours,)
            page_num += 1
    
            # KROK 2: Aktualizacja opisu paska i generowanie mapy
            pbar.set_description(f"       Generowanie MAP: {oid}")
            draw_context_page(pdf, oid, page_num)
            page_num += 1

    print(f"\n[INFO] Proces zakończony. Plik wynikowy: {output_name}.")

if __name__ == "__main__":
    main()
