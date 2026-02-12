#!/usr/bin/env python3

import numpy as np
import pandas as pd
import pytz
import json
import astropy.units as u
from astropy.coordinates import SkyCoord, EarthLocation
from astropy.time import Time
from astroplan import Observer
from datetime import date, datetime, timedelta, timezone
from pathlib import Path
from typing import Tuple, List, Dict, Any
from timezonefinder import TimezoneFinder
import warnings
from astropy.utils.exceptions import AstropyWarning
from erfa import ErfaWarning
import ssl
import certifi
from tqdm import tqdm

# =====================================================
# KONFIGURACJA OSTRZEŻEŃ I SSL
# =====================================================
tqdm.pandas(ncols=119)
warnings.filterwarnings("ignore", category=ErfaWarning)
warnings.filterwarnings("ignore", category=AstropyWarning)
ssl._create_default_https_context = lambda: ssl.create_default_context(cafile=certifi.where())

# =====================================================
# KONFIGURACJA PLIKÓW I PARAMETRÓW
# =====================================================

class Config:
    OUTPUT_FILE = Path("vis_data.json")
    CATALOG_PATH = Path("katalog_astro_full.csv")

    DEFAULT_FOCAL_LENGTH_MM = 300.0
    DEFAULT_SENSOR_WIDTH_MM = 23.5
    DEFAULT_SENSOR_HEIGHT_MM = 15.7
    DEFAULT_SENSOR_PITCH_UM = 3.76
    DEFAULT_SENSOR_ROWS = 4176
    DEFAULT_SENSOR_COLS = 6248
    
CAT_ORDER = ["ngc", "ic", "sh2", "rcw", "lbn", "ced", "pgc", "barn", "ldn"]

# =====================================================
# TABELE PUNKTACJI (SCORECARD)
# =====================================================

# TABELA 1: WPŁYW FILTRÓW
# Klucz: Typ obiektu (Upper Case)
# Wartość: [Punkty "Z Triband" (No NB), Punkty "Z Narrowband" (Has NB)]
SCORING_FILTER = {
    "GX":     [10, 0],
    "G":      [10, 0],
    "DN":     [0,  0],
    "NB":     [20, 30],
    "OCL":    [20, 0],
    "HII":    [20, 30],
    "*":      [10, 0],
    "OTHER":  [10, 0],
    "DUP":    [10, 0],
    "GCL":    [20, 0],
    "**":     [10, 0],
    "GPAIR":  [10, 0],
    "NEB":    [20, 30],
    "*ASS":   [10, 0],
    "CL+N":   [20, 30],
    "RFN":    [0,  0], 
    "PN":     [20, 30],
    "GTRPL":  [10, 0],
    "SNR":    [20, 30],
    "GGROUP": [10, 0],
    "NOVA":   [10, 0]
}

# TABELA 2: WPŁYW NIEBA (BORTLE)
# Klucz: Typ obiektu
# Wartość: [Punkty "Bortle <= 5", Punkty "Bortle > 5"]
SCORING_BORTLE = {
    "GX":     [10, 10],
    "G":      [10, 10],
    "DN":     [30, 0],  
    "NB":     [20, 20],
    "OCL":    [20, 20],
    "HII":    [20, 20],
    "*":      [10, 10],
    "OTHER":  [10, 10],
    "DUP":    [10, 10],
    "GCL":    [20, 20],
    "**":     [10, 10],
    "GPAIR":  [10, 10],
    "NEB":    [20, 20],
    "*ASS":   [10, 10],
    "CL+N":   [20, 20],
    "RFN":    [30, 0],  
    "PN":     [20, 20],
    "GTRPL":  [20, 20],
    "SNR":    [10, 10],
    "GGROUP": [10, 10],
    "NOVA":   [10, 10]
}

# TABELA 3: BONUSY ZA "SŁAWĘ"
SCORE_FAMOUS = {
    "messier":  25,
    "caldwell": 20,
    "herschel": 15
}

# TABELA 4: BONUSY ZA JAKOŚĆ DANYCH
SCORE_DATA_QUALITY = {
    "measured": 15,
    "estimated": 0
}

# =====================================================
# HELPER FUNCTIONS
# =====================================================
def fmt(n: int) -> str:
    """Format liczby z separatorem tysięcy (spacja)."""
    return f"{int(n):_}".replace("_", " ")

def calculate_surface_brightness(mag: float, size_arcmin: float) -> float:
    if size_arcmin <= 0: return 99.0
    # Area in arcmin^2 (uproszczenie: koło)
    area = np.pi * (size_arcmin / 2.0)**2
    if area <= 0: return 99.0
    return mag + 2.5 * np.log10(area)

def is_messier_token(token: str) -> bool:
    if token is None: return False
    t = str(token).strip().upper()
    if not t.startswith("M") or len(t) <= 1: return False
    rest = t[1:].strip()
    return rest.isdigit() and 1 <= int(rest) <= 110

def is_caldwell_token(token: str) -> bool:
    if token is None: return False
    t = str(token).strip().upper()
    if not t.startswith("C") or len(t) <= 1: return False
    rest = t[1:].strip()
    return rest.isdigit() and 1 <= int(rest) <= 109

def is_herschel_token(token: str) -> bool:
    if token is None: return False
    t = str(token).strip().upper()
    if not t.startswith("H") or len(t) <= 1: return False
    rest = t[1:].strip()
    return rest.isdigit() and 1 <= int(rest) <= 400

def analyze_famous_status(extra_info_str: str) -> pd.Series:
    has_m = False
    has_c = False
    has_h = False
    if pd.notna(extra_info_str):
        s = str(extra_info_str).strip()
        if s:
            tokens = [t.strip() for t in s.split(",") if t.strip()]
            for token in tokens:
                if is_messier_token(token): has_m = True
                elif is_caldwell_token(token): has_c = True
                elif is_herschel_token(token): has_h = True
    return pd.Series([has_m, has_c, has_h])

# =====================================================
# ASTRO MATH
# =====================================================

class VectorAstro:
    @staticmethod
    def get_gmst_vector(jd_array: np.ndarray) -> np.ndarray:
        T = (jd_array - 2451545.0) / 36525.0
        gmst = (280.46061837 + 360.98564736629 * (jd_array - 2451545.0) + 0.000387933 * T**2 - T**3 / 38710000.0)
        return gmst % 360.0

    @staticmethod
    def alt_az_vector(ha_rad: np.ndarray, dec_rad: np.ndarray, lat_rad: float) -> np.ndarray:
        sin_alt = (np.sin(dec_rad) * np.sin(lat_rad) + np.cos(dec_rad) * np.cos(lat_rad) * np.cos(ha_rad))
        return np.degrees(np.arcsin(np.clip(sin_alt, -1.0, 1.0)))

def get_full_year_smart_grid(
    year: int,
    observer_loc: Dict[str, Any],
    sun_limit_deg: float,
    step_minutes: int = 15,
) -> Tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
    """
    Generuje siatkę czasu (JD) tylko dla nocy, używając astroplan.
    
    Returns:
        jd_array: Julian Dates dla punktów pomiarowych
        month_array: Indeks miesiąca (0-11) dla każdego punktu
        lst_array: Local Sidereal Time dla każdego punktu
        night_ids: ID nocy dla każdego punktu (0, 1, 2, ..., ~365)
    """
    print(f"[INFO] Generowanie SMART GRID dla roku {year}, limit {sun_limit_deg}°")
    
    obs = Observer(
        latitude=observer_loc["lat"] * u.deg,
        longitude=observer_loc["lon"] * u.deg,
        elevation=0 * u.m,
    )
    
    t_scan = Time(f"{year}-01-01 12:00:00", scale="utc")
    t_end = Time(f"{year+1}-01-01 12:00:00", scale="utc")
    step_day = step_minutes / (60.0 * 24.0)
    
    all_jds: List[np.ndarray] = []
    all_night_ids: List[np.ndarray] = []
    night_counter = 0
    
    n_days = int(t_end.jd - t_scan.jd)    
    with tqdm(total=n_days, desc="       SmartGrid", unit="day", ncols=119) as pbar:
        while t_scan < t_end:
            try:
                t_set = obs.sun_set_time(t_scan, horizon=sun_limit_deg * u.deg, which="next")
                t_rise = obs.sun_rise_time(t_set, horizon=sun_limit_deg * u.deg, which="next")
                
                if t_rise > t_set:
                    night_grid = np.arange(t_set.jd, t_rise.jd, step_day)
                    all_jds.append(night_grid)
                    
                    night_ids_for_this_night = np.full(len(night_grid), night_counter, dtype=int)
                    all_night_ids.append(night_ids_for_this_night)
                    
                    night_counter += 1
            except Exception:
                pass
            
            t_scan = t_scan + 1.0 * u.day
            pbar.update(1)
    
    if not all_jds:
        print("! Ostrzeżenie: Nie znaleziono żadnych nocy (błąd parametrów?).")
        return np.array([]), np.array([]), np.array([]), np.array([])
    
    jd_array = np.concatenate(all_jds)
    night_id_array = np.concatenate(all_night_ids)
    
    print(
        f"       ✓ Wygenerowano {fmt(len(jd_array))} punktów pomiarowych dla {fmt(night_counter)} nocy "
        f"(~{len(jd_array)*step_minutes/60:.0f}h obserwacji)."
    )
    
    # Oblicz miesiące i LST
    t_all = Time(jd_array, format="jd")
    dates = t_all.to_datetime()
    months_array = np.array([d.month - 1 for d in dates], dtype=int)
    
    lon_deg = observer_loc["lon"]
    gmst = VectorAstro.get_gmst_vector(jd_array)
    lst_array = (gmst + lon_deg) % 360.0
    
    return jd_array, months_array, lst_array, night_id_array


def compute_all_visibilities_hybrid(
    df_objects: pd.DataFrame,
    observer_loc: Dict[str, Any],
    jd_array: np.ndarray,
    night_ids: np.ndarray,
    month_array: np.ndarray,
    lst_array: np.ndarray,
    min_alt: float,
    min_hours: float,
    step_minutes: int,
) -> Tuple[np.ndarray, np.ndarray]:
    """
    Oblicza widoczność obiektów HYBRYDOWO:
    1. Per NOC - do precyzyjnego filtrowania (czy była noc z >= min_hours?)
    2. Per MIESIĄC - do kompatybilności z JSON (suma godzin w miesiącu)
    
    Args:
        df_objects: DataFrame z obiektami (musi mieć 'ra', 'dec')
        observer_loc: Dict z 'lat', 'lon'
        jd_array: Tablica Julian Dates (punkty co 15 min)
        night_ids: ID nocy dla każdego punktu JD
        month_array: Indeks miesiąca (0-11) dla każdego punktu JD
        lst_array: Local Sidereal Time dla każdego punktu JD
        min_alt: Minimalna wysokość nad horyzontem (stopnie)
        min_hours: Minimalna liczba godzin w nocy (do zliczania)
        step_minutes: Krok czasowy siatki (typowo 15)
    
    Returns:
        hours_per_month: (n_objects, 12) - suma godzin widoczności per miesiąc
        nights_above_threshold: (n_objects,) - liczba nocy z >= min_hours
    """
    if df_objects.empty or len(jd_array) == 0:
        return np.zeros((len(df_objects), 12)), np.zeros(len(df_objects), dtype=int)
    
    n_objects = len(df_objects)
    print(f"[INFO] Obliczanie macierzowe dla {fmt(n_objects)} obiektów. (To może potrwać kilka minut.)")
    
    # Przygotuj dane
    lat_rad = np.radians(observer_loc["lat"])
    ra_objs = df_objects["ra"].values
    dec_objs = df_objects["dec"].values
    
    ra_rad_row = np.radians(ra_objs)[np.newaxis, :]
    dec_rad_row = np.radians(dec_objs)[np.newaxis, :]
    lst_rad_col = np.radians(lst_array)[:, np.newaxis]
    
    # Oblicz wysokości dla wszystkich punktów czasowych
    ha_rad_matrix = lst_rad_col - ra_rad_row
    alt_matrix = VectorAstro.alt_az_vector(ha_rad_matrix, dec_rad_row, lat_rad)
    
    # Maska widoczności (True jeśli Alt >= min_alt)
    visible_mask = (alt_matrix >= min_alt).astype(float)
    step_hours = step_minutes / 60.0
    weighted_mask = visible_mask * step_hours  # Każdy punkt = 0.25h
    
    # ============================================================
    # AGREGACJA 1: PER NOC (do filtrowania)
    # ============================================================
    unique_nights = np.unique(night_ids)
    n_nights = len(unique_nights)
    nights_above_threshold = np.zeros(n_objects, dtype=int)
        
    for night_id in unique_nights:
        # Punkty należące do tej nocy
        night_mask = (night_ids == night_id)
        
        # Suma godzin w tej nocy dla każdego obiektu
        hours_this_night = np.sum(weighted_mask[night_mask], axis=0)
        
        # Zlicz: które obiekty miały >= min_hours w TEJ nocy
        nights_above_threshold += (hours_this_night >= min_hours).astype(int)
    
    # ============================================================
    # AGREGACJA 2: PER MIESIĄC (do JSON)
    # ============================================================
    hours_per_month = np.zeros((n_objects, 12))
    
    for m in range(12):
        month_mask = (month_array == m)
        if not np.any(month_mask):
            continue
        hours_per_month[:, m] = np.sum(weighted_mask[month_mask], axis=0)
    
    print(f"       ✓ Obliczono widoczność per noc i per miesiąc.")
    
    return hours_per_month, nights_above_threshold

def compute_max_altitude(lat: float, dec: float) -> float:
    lat_rad = np.radians(lat)
    dec_rad = np.radians(dec)
    sin_alt_max = np.sin(lat_rad) * np.sin(dec_rad) + np.cos(lat_rad) * np.cos(dec_rad)
    return np.degrees(np.arcsin(np.clip(sin_alt_max, -1.0, 1.0)))

def compute_rarity_score(max_altitude: float, months_observable: int) -> float:
    alt_rarity = 0.0 if max_altitude >= 70 else (0.3 if max_altitude >= 50 else (0.7 if max_altitude >= 35 else 1.0))
    month_rarity = 0.0 if months_observable >= 6 else (0.3 if months_observable >= 4 else (0.7 if months_observable >= 2 else 1.0))
    return (alt_rarity + month_rarity) / 2.0

# =====================================================
# IMPUTATION (UZUPEŁNIANIE DANYCH)
# =====================================================

def impute_missing_data(row: pd.Series) -> pd.Series:
    """
    Uzupełnia braki Mag i Size. Ustawia flagi 'measured' lub 'estimated'.
    """
    r = row.copy()
    
    # Normalizacja typu
    otype = str(r.get("type_code", "Other")).strip()
    
    # 1. MAGNITUDE
    mag = r.get("mag")
    if pd.notna(mag) and str(mag).strip() != "":
        r["mag_source"] = "measured"
        r["mag"] = float(mag)
    else:
        r["mag_source"] = "estimated"
        # Defaults based on Type (fallback logic)
        if otype in ["Gx", "G", "GPair", "GTrpl", "GGroup"]: val = 13.0
        elif otype in ["PN"]: val = 11.0
        elif otype in ["OCl", "GCl", "Cl+N"]: val = 8.0
        elif otype in ["HII", "SNR", "NB", "Neb"]: val = 10.0
        elif otype in ["DN", "LDN"]: val = 99.0
        else: val = 12.0
        r["mag"] = val

    # 2. SIZE (Arcmin)
    size = r.get("size")
    if pd.notna(size) and str(size).strip() != "" and float(size) > 0:
        r["size_source"] = "measured"
        r["size"] = float(size)
    else:
        r["size_source"] = "estimated"
        # Defaults based on Type
        if otype in ["Gx", "G"]: val = 5.0
        elif otype in ["GPair", "GTrpl", "GGroup"]: val = 8.0
        elif otype in ["PN"]: val = 1.0
        elif otype in ["OCl", "GCl"]: val = 15.0
        elif otype in ["HII", "SNR"]: val = 45.0
        elif otype in ["DN", "NB", "Neb"]: val = 20.0
        else: val = 10.0
        r["size"] = val
        
    return r

# =====================================================
# SCORING (SYSTEM SCORECARD)
# =====================================================

def calculate_scorecard(row: pd.Series, user_params: Dict[str, Any]) -> Dict[str, Any]:
    """
    Oblicza punktację w systemie Scorecard (Karta Wyników).
    Suma punktów z 4 kategorii: Sława + Typ/Filtr + Typ/Bortle + Jakość Danych.
    """
    
    obj_type = str(row.get("type_code", "Other")).upper().strip()
    # Mapowanie 'OCl' -> 'OCL' etc. jeśli potrzebne (ale zakładam zgodność z CSV)
    if obj_type not in SCORING_FILTER:
        obj_type = "OTHER"
        
    current_score = 0.0
    breakdown = {} 

    # A. SŁAWA (FAMOUS BONUS)
    famous_points = 0.0
    if row.get("is_messier"): famous_points += SCORE_FAMOUS["messier"]
    if row.get("is_caldwell"): famous_points += SCORE_FAMOUS["caldwell"]
    if row.get("is_herschel"): famous_points += SCORE_FAMOUS["herschel"]
    
    current_score += famous_points
    breakdown['famous'] = famous_points

    # B. TYP + FILTR
    # has_narrowband = True -> index 1 (Narrowband)
    # has_narrowband = False -> index 0 (Triband/None)
    filter_idx = 1 if user_params.get("has_narrowband", False) else 0
    filter_points = SCORING_FILTER[obj_type][filter_idx]
    current_score += filter_points
    breakdown['filter_type'] = filter_points

    # C. TYP + BORTLE
    # Bortle 1-5 -> index 0, Bortle > 5 -> index 1
    bortle_val = user_params.get("bortle_range", (6, 7))[1]
    bortle_idx = 0 if bortle_val <= 5 else 1
    bortle_points = SCORING_BORTLE[obj_type][bortle_idx]
    current_score += bortle_points
    breakdown['bortle_type'] = bortle_points

    # D. JAKOŚĆ DANYCH
    quality_points = 0.0
    quality_points += SCORE_DATA_QUALITY.get(row.get("mag_source"), 0)
    quality_points += SCORE_DATA_QUALITY.get(row.get("size_source"), 0)
    
    current_score += quality_points
    breakdown['quality'] = quality_points

    return {
        "final_score": current_score,
        "breakdown": breakdown
    }

# =====================================================
# FOV CALCULATOR
# =====================================================

class FOVCalculator:
    @staticmethod
    def calculate_fov_deg(focal, sw, sh):
        return 57.3 * sw / focal, 57.3 * sh / focal
    @staticmethod
    def min_object_size_arcmin(fov_w, fov_h, percent):
        return min(fov_w, fov_h) * percent / 100.0 * 60.0

# =====================================================
# USER INTERFACE
# =====================================================

class YearManager:
    def select_interactive(self) -> int:
        current_year = date.today().year
        print("=" * 119)
        print("ROK OBSERWACJI")
        print("=" * 119)
        year_str = input(
            f"Enter current year lub wpisz wybrany pomiędzy 2000 a 2100 "
            f"(Enter = {current_year}): "
        ).strip()

        if year_str:
            try:
                year = int(year_str)
                if year < 2000 or year > 2100:
                    raise ValueError
            except Exception:
                print(f"Nieprawidłowy, używam {current_year}")
                year = current_year
        else:
            year = current_year

        print(f"Rok: {year}")
        return year

class LocationManager:
    def select_interactive(self) -> Dict[str, Any]:
        print("=" * 119)
        print("KROK 1: WYBÓR LOKALIZACJI")
        print("=" * 119)
        print("1. Poznań, Polska")
        print("2. Kraków, Polska")
        print("3. Wpisz lokalizację")

        choice = input("[1/2/3, domyślnie 1]: ").strip() or "1"

        if choice == "1":
            loc = {
                "lat": 52.4095,
                "lon": 16.9319,
                "tz": "Europe/Warsaw",
                "name": "Poznań, Polska",
            }
            print(f"Wybrano {loc['name']} ({loc['lat']}N, {loc['lon']}E)")
            return loc

        if choice == "2":
            loc = {
                "lat": 50.0647,
                "lon": 19.9450,
                "tz": "Europe/Warsaw",
                "name": "Kraków, Polska",
            }
            print(f"Wybrano {loc['name']} ({loc['lat']}N, {loc['lon']}E)")
            return loc

        city = input("Podaj lokalizację (pisz wielką literą, działa lepiej, np. 'Toruń, Polska'): ").strip()
        if not city:
            print("Brak nazwy, używam domyślnej: Poznań.")
            return {
                "lat": 52.4095,
                "lon": 16.9319,
                "tz": "Europe/Warsaw",
                "name": "Poznań, Polska",
            }

        try:
            print(f"[INFO] Szukanie współrzędnych dla '{city}'...")
            loc_astro = EarthLocation.of_address(city)
            lat = loc_astro.lat.to(u.deg).value
            lon = loc_astro.lon.to(u.deg).value

            tf = TimezoneFinder()
            tz_name = tf.timezone_at(lng=lon, lat=lat)
            if tz_name is None:
                tz_name = "Europe/Warsaw"

            loc = {
                "lat": lat,
                "lon": lon,
                "tz": tz_name,
                "name": city,
            }
            print(f"[INFO] Znaleziono: {loc['name']} {loc['lat']:.4f}N, {loc['lon']:.4f}E, "
                  f"strefa czasowa: {tz_name}")
            return loc
        except Exception as e:
            print(f"Błąd geokodowania ({e}), używam Poznań.")
            return {
                "lat": 52.4095,
                "lon": 16.9319,
                "tz": "Europe/Warsaw",
                "name": "Poznań, Polska",
            }

def get_user_prefs():
    print("="*119 + "\nKROK 2: PARAMETRY WIDOCZNOŚCI I FOV\n" + "="*119)
    
    min_alt = float(input("A. Minimalna wysokość obiektu nad horyzontem (domyślnie 25°): ") or 25.0)
    min_hours = float(input("\nB. Minimalna liczba godzin, którą obiekt jest widoczny w nocy powyżej progu wysokości (domyślnie 3): ") or 3.0)
    print("\nC. Ciemność nieba – kąt słońca pod horyzontem:")
    print("• zmierzch cywilny  (-6°)")
    print("• zmierzch żeglarski (-12°)")
    print("• zmierzch astronomiczny (-18°)")
    sl_choice = input("Możesz wpisać dowolną sensowną wartość, która oznacza liczbę stopni pod horyzontem (domyślnie 12°): ").strip() or "12"
    sun_limit = -float(sl_choice)

    
    print("\nD. Określenie FOV")
    if (input("Użyć domyślnego setupu RedCat61 + ASI2600MC Pro? [y/n, domyślnie y]: ").lower() or "y") == "y":
        focal, sw, sh = Config.DEFAULT_FOCAL_LENGTH_MM, Config.DEFAULT_SENSOR_WIDTH_MM, Config.DEFAULT_SENSOR_HEIGHT_MM
        pitch, rows, cols = Config.DEFAULT_SENSOR_PITCH_UM, Config.DEFAULT_SENSOR_ROWS, Config.DEFAULT_SENSOR_COLS
    else:
        focal = float(input(" Ogniskowa (mm, domyślnie 300): "))
        sw = float(input(" Szerokość sensora (mm, domyślnie 23.5): "))
        sh = float(input(" Wysokość sensora (mm, domyślnie 15.7): "))
        pitch = 3.76; rows = 4000; cols = 6000 # dummy defaults

    fov_w, fov_h = FOVCalculator.calculate_fov_deg(focal, sw, sh)
    print(f"FOV: {fov_w:.2f}° x {fov_h:.2f}°")
    
    perc = float(input("\nE. Minimalny rozmiar obiektu jako % krótszego boku FOV (domyślnie 10): ") or 10.0)
    min_size = FOVCalculator.min_object_size_arcmin(fov_w, fov_h, perc)
    print(f"Minimalny rozmiar obiektu (arcmin): {min_size:.1f}'")
    
    print("\nF. Skala Bortle – określenie zanieczyszczenia światłem:")
    print("• 1  Bortle 1–3  (wieś, ciemne niebo, pomijalne LP)")
    print("• 2  Bortle 4–5  (przedmieścia, umiarkowane LP)")
    print("• 3  Bortle 6–7  (miasto, silne LP)")
    print("• 4  Bortle 8–9  (centrum miasta, ekstremalne LP)")
    b_choice = input("Wybierz Twój stopień zanieczyszczenia nieba światłem (domylśnie 4): ").strip() or "4"
    b_map = {"1": (1,3), "2": (4,5), "3": (6,7), "4": (8,9)}
    bortle_rng = b_map.get(b_choice, (8,9))
    has_nb = (input("\nG. Czy zamierzasz korzystać z filtrów narrowband (H/O/S)? [y/n, domyślnie n]: ").lower() or "n") == "y"
    prefer_famous = (input("\nH. Czy w wyborze obiektów premiować katalogi Messier/Caldwell/Herschel? [y/n, domyślnie y]: ").lower() or "y") == "y"
    
    return {
        "minalt": min_alt, "sunlimit": sun_limit, "minhours": min_hours,
        "minsizearcmin": min_size, "bortle_range": bortle_rng,
        "has_narrowband": has_nb, "prefer_famous": prefer_famous,
        "camera": {"focal": focal, "sw": sw, "sh": sh}
    }

def print_imputation_stats(df_before: pd.DataFrame, df_after: pd.DataFrame):
    """
    Wyświetla statystyki jakości danych po imputacji.
    Porównuje dane przed i po uzupełnieniu.
    """
    total = len(df_after)
    
    # Policz zmierzone mag i size
    mag_measured = (df_after["mag_source"] == "measured").sum()
    size_measured = (df_after["size_source"] == "measured").sum()
    both_measured = ((df_after["mag_source"] == "measured") & 
                     (df_after["size_source"] == "measured")).sum()
    
    mag_pct = (mag_measured / total) * 100
    size_pct = (size_measured / total) * 100
    both_pct = (both_measured / total) * 100
    
    print(f"       ✓ Zmierzone mag: {fmt(mag_measured)}/{fmt(total)} ({mag_pct:.1f}%)")
    print(f"       ✓ Zmierzone size: {fmt(size_measured)}/{fmt(total)} ({size_pct:.1f}%)")
    print(f"       ✓ Oba zmierzone: {fmt(both_measured)}/{fmt(total)} ({both_pct:.1f}%)")

# =====================================================
# ADAPTACYJNY SOFT-CUT
# =====================================================
def apply_adaptive_soft_cut(
    df: pd.DataFrame,
    base_min_score: float = 10.0,
    min_len_for_cut: int = 500,
    keep_top_percent: float = 70.0,
) -> pd.DataFrame:
    """
    Adaptacyjny soft‑cut:
      1) Odrzuca wszystko poniżej base_min_score (twardy, ale niski próg bezpieczeństwa).
      2) Jeśli po tym zostaje dużo obiektów (>= min_len_for_cut),
         ucina dół rozkładu tak, żeby zostawić tylko górne keep_top_percent%.

    Przykład:
      - base_min_score = 10.0
      - keep_top_percent = 70.0  → zachowaj 70% najlepszych (obcinasz 30% najgorszych)
    """
    if df.empty or "final_score" not in df.columns:
        return df

    # 1) Twardy minimalny próg (to Twój dotychczasowy if final_score <= 10: continue)
    before = len(df)
    df = df[df["final_score"] > base_min_score].reset_index(drop=True)
    after = len(df)
    print(f"[INFO] Zastosowano ostre cięcie dla obiektów z wynikiem mniejszym niż {base_min_score:.1f} punktów.\n[INFO] Pozostało {fmt(after)} obiektów z {fmt(before)}. ") 

    if df.empty:
        return df

    # 2) Adaptacyjny cut względem rozkładu tylko jeśli lista jest wystarczająco duża
    if len(df) < min_len_for_cut:
        print(f"[INFO] Soft-cut 2: lista < {min_len_for_cut} obiektów – pomijam cięcie percentylowe.")
        return df

    scores = df["final_score"].astype(float)
    # Dolny percentyl, który utniesz – np. 30% na dole, zostawiasz górne 70%
    lower_tail_percent = 100.0 - keep_top_percent
    p_cut = np.percentile(scores, lower_tail_percent)

    before2 = len(df)
    df = df[df["final_score"] > p_cut].reset_index(drop=True)
    after2 = len(df)

    print(
        f"[INFO] Wartości wyliczone przez cięcie adaptacyjne: final_score > p{lower_tail_percent:.0f} (≈{p_cut:.1f}).\n"
        f"[INFO] Pozostało {fmt(after2)} obiektów z {fmt(before2)} ."
    )
    return df


# =====================================================
# CLUSTERING (Atlas optimization)
# =====================================================

def choose_best_catalog_id(main_id: str, extra_info: str) -> str:
    def cat_rank(obj_id: str) -> int:
        s = obj_id.lower().strip()
        for idx, c in enumerate(CAT_ORDER):
            if s.startswith(c): return idx
        return 99
    
    candidates = [str(main_id).strip()] + [x.strip() for x in str(extra_info).split(",") if x.strip()]
    return sorted(candidates, key=lambda x: (cat_rank(x), x))[0]

def optimize_atlas_pages(df: pd.DataFrame) -> pd.DataFrame:
    print(f"[INFO] Optymalizacja kadrów (Clustering 2.0°)")
    if df.empty:
        return df

    df_sorted = df.sort_values(by=["final_score", "size"], ascending=[False, False]).reset_index(drop=True)
    coords = SkyCoord(ra=df_sorted["ra"].values*u.deg, dec=df_sorted["dec"].values*u.deg)

    final_objs = []
    covered = set()
    MAX_CLUSTER_SIZE_ARCMIN = 300.0  # np. 5° – możesz dostroić

    for i in range(len(df_sorted)):
        if i in covered:
            continue

        leader = df_sorted.iloc[i].copy()

        # znajdź sąsiadów w promieniu 2°
        seps = coords[i].separation(coords)
        neighbors = np.where(seps < 2.0*u.deg)[0]

        # zbuduj indeksy całego klastra (leader + neighbors)
        cluster_idx = [i]
        for n_idx in neighbors:
            if n_idx == i or n_idx in covered:
                continue
            covered.add(n_idx)
            cluster_idx.append(n_idx)

        subset = df_sorted.iloc[cluster_idx]

        # 1) ROZMIAR: bierz max(size) z klastra, ale z limitem na giganty
        sizes = subset["size"].dropna()
        if not sizes.empty:
            reasonable = sizes[sizes <= MAX_CLUSTER_SIZE_ARCMIN]
            if not reasonable.empty:
                leader["size"] = float(reasonable.max())
            else:
                # wszystko > MAX_CLUSTER_SIZE_ARCMIN – przytnij
                #leader["size"] = float(sizes.max())
                leader["size"] = MAX_CLUSTER_SIZE_ARCMIN

        # 2) ZBIERANIE ID I NAZW:
        cluster_ids = set()
        for _, row in subset.iterrows():
            if pd.notna(row.get("id")):
                cluster_ids.add(str(row["id"]))
            # extra_info z pod-obiektów
            if pd.notna(row.get("extra_info")):
                for x in str(row["extra_info"]).split(","):
                    x = x.strip()
                    if x:
                        cluster_ids.add(x)

        # Zbierz common_names z całego klastra
        common_names = []
        for _, row in subset.iterrows():
            if pd.notna(row.get("common_names")):
                cn = str(row["common_names"]).strip()
                if cn and cn not in common_names:
                    common_names.append(cn)

        valid_ids = [x for x in cluster_ids if x and x.lower() != "nan"]
        
        leader_id = str(leader.get("id", "")).strip()
        leader_id_upper = leader_id.upper()
        
        # Jeśli leader JEST już NGC/IC/Sh2 → NIE ZMIENIAJ ID
        if leader_id_upper.startswith(("NGC", "IC", "SH2-")):
            best_id = leader_id
        else:
            # W przeciwnym razie wybierz najlepsze ID wg CAT_ORDER
            extra_ids_str = ", ".join(x for x in valid_ids if x != leader_id)
            best_id = choose_best_catalog_id(leader_id, extra_ids_str)
        
        # Extra_info = wszystkie pozostałe ID oprócz best_id
        rem_ids = [x for x in valid_ids if x != best_id]
        
        leader["id"] = best_id
        leader["extra_info"] = ", ".join(sorted(set(rem_ids)))

        leader["common_names"] = ", ".join(common_names)

        final_objs.append(leader)

    return pd.DataFrame(final_objs)


# =====================================================
# MAIN
# =====================================================

def main():
    # 1. SETUP
    year = YearManager().select_interactive()
    loc = LocationManager().select_interactive()
    params = get_user_prefs()
    
    # 2. LOAD DATA
    print("\n" + "="*119 + "\nKROK 3: ŁADOWANIE I FILTROWANIE\n" + "="*119)
    if not Config.CATALOG_PATH.exists():
        print("Brak pliku katalogu!"); return
        
    df = pd.read_csv(Config.CATALOG_PATH, low_memory=False)
    
    if "type" in df.columns:
        df.rename(columns={"type": "type_code"}, inplace=True)
    
    # Clean coords
    df = df.dropna(subset=["ra", "dec"])
    print(f"[INFO] Załadowano {fmt(len(df))} obiektów z pełnymi danymi.")

    # Add Famous Flags
    print("[INFO] Dodawanie flag Messier/Caldwell/Herschel")
    tqdm.pandas(ncols=119, desc="       Flagi M/C/H")
    if "extra_info" not in df.columns: df["extra_info"] = ""
    flags = df["extra_info"].progress_apply(analyze_famous_status)
    flags.columns = ["is_messier", "is_caldwell", "is_herschel"]
    df = pd.concat([df, flags], axis=1)
    
    # 3. IMPUTATION (Wypełnianie braków + flagowanie źródła)
    print("[INFO] Uzupełnianie danych o rozmiarze i jasności (Measured vs Estimated)")
    tqdm.pandas(ncols=119, desc="       Imputacja")
    df = df.progress_apply(impute_missing_data, axis=1)
    print_imputation_stats(df, df)

    # 4. FAST GEOMETRIC FILTER
    print(f"[INFO] Filtr geometryczny (Alt > {params['minalt']}°, Visible > {params['minhours']}h)")
    total_objects = len(df)
    lat_rad = np.radians(loc["lat"])
    dec_rad = np.radians(df["dec"].values)
    # Fast Max Alt check
    sin_alt_max = np.sin(lat_rad)*np.sin(dec_rad) + np.cos(lat_rad)*np.cos(dec_rad)
    alt_max = np.degrees(np.arcsin(np.clip(sin_alt_max, -1, 1)))
    df_after_alt = df[alt_max >= params["minalt"]].copy()
    print(f"[INFO] Po filtrze min. wysokości ({params['minalt']}°) pozostało {fmt(len(df_after_alt))} z {fmt(total_objects)} obiektów.")
    df = df_after_alt
    
    if df.empty: print("Brak obiektów."); return

    # 5. SMART GRID & VISIBILITY (HYBRYDOWA)
    step = 15
    jd, months, lst, night_ids = get_full_year_smart_grid(
        year, loc, params["sunlimit"], step
    )
    
    if len(jd) == 0:
        print("Brak nocy do analizy!")
        return
    
    # Oblicz widoczność per noc + per miesiąc
    hours_per_month, nights_above_threshold = compute_all_visibilities_hybrid(
        df, loc, jd, night_ids, months, lst, 
        params["minalt"], params["minhours"], step
    )
    
    # 5.5 FAST FILTER: Odrzuć obiekty bez ŻADNEJ nocy z >= minhours
    print(f"[INFO] Filtr czasu widoczności per noc")
    objects_before_filter = len(df)
    
    mask_has_good_night = nights_above_threshold >= 1  
    df = df[mask_has_good_night].reset_index(drop=True)
    hours_per_month = hours_per_month[mask_has_good_night]
    nights_above_threshold = nights_above_threshold[mask_has_good_night]
    
    print(f"       Spośród {fmt(objects_before_filter)} obiektów, {fmt(len(df))} przynajmniej w jedną noc są widoczne przez {params['minhours']}h.")
    
    if df.empty:
        print("Brak obiektów po filtrze czasu widoczności.")
        return

    
    # 6. SCORING LOOP (THE FUNNEL)
    scored_rows = []
    
    print("[INFO] Scoring i Soft Cut")
    for i, (_, row) in enumerate(df.iterrows()):
        # Pobierz dane widoczności
        hours = hours_per_month[i]  # Godziny per miesiąc (dla JSON)
        
        # ODRZUCENIE "DUCHÓW LBN":
        #sb = calculate_surface_brightness(row['mag'], row['size'])
        # Jeśli SB > 16.0 (bardzo ciemne powierzchniowo) ORAZ jasność była szacowana (niepewna)
        # To prawie na pewno jest wielki, ciemny cirrus galaktyczny z katalogu LBN/LDN.
        #if row['mag_source'] == 'estimated' and sb > 16.0:
        #    continue 
            
        # Oblicz miesiące obserwowalne (dla rarity_score)
        months_obs = int(np.sum(hours >= params["minhours"]))
        
        # Calculate Scorecard
        score_res = calculate_scorecard(row, params)
        final_score = score_res["final_score"]
              
        # Max Alt i rarity
        max_alt = compute_max_altitude(loc["lat"], row["dec"])
        rarity = compute_rarity_score(max_alt, months_obs)
        
        # Zapisz do listy
        r_dict = row.to_dict()
        r_dict.update({
            "final_score": final_score,
            "base_score": final_score,
            "score_breakdown": score_res["breakdown"],
            "hours": hours.tolist(),  # ← Godziny per miesiąc
            "maxalt": max_alt,
            "monthsobservable": months_obs,
            "rarity_score": rarity,
            "nights_above_threshold": int(nights_above_threshold[i]) ,
            "is_messier": bool(row.get("is_messier", False)),
            "is_caldwell": bool(row.get("is_caldwell", False)),
            "is_herschel": bool(row.get("is_herschel", False)),
            "is_catalog_famous": bool(row.get("is_catalog_famous", False) if "is_catalog_famous" in row else (row.get("is_messier") or row.get("is_caldwell") or row.get("is_herschel")))
        })
        scored_rows.append(r_dict)

    df_scored = pd.DataFrame(scored_rows)
    if df_scored.empty:
        print("Brak obiektów po scoringu.")
        return
    
    # ZASTĄP dotychczasowy prosty soft-cut:
    df_scored = apply_adaptive_soft_cut(
        df_scored,
        base_min_score=10.0,     # Twój dotychczasowy próg
        min_len_for_cut=500,     # od jakiej liczności używać percentyla
        keep_top_percent=70.0  # zachowaj górne 70% najlepszych
    )

    if df_scored.empty: print("[INFO] Brak obiektów po filtracji."); return

    # 7. CLUSTERING
    df_clustered = optimize_atlas_pages(df_scored)
    print(f"       ✓ Zredukowano liczbe kadrów z {fmt(len(df_scored))} do {fmt(len(df_clustered))}.")
    
    # 8. HARD CONSTRAINTS (Size Filter)
    # Tutaj używamy wartości (imputowanej lub zmierzonej) do ostatecznej decyzji
    print(f"[INFO] Hard Filter: Size >= {params['minsizearcmin']:.1f}'")
    df_final = df_clustered[df_clustered["size"] >= params["minsizearcmin"]].copy()
    
    # Sortowanie końcowe
    df_final = df_final.sort_values(by="final_score", ascending=False).reset_index(drop=True)
    
    print(f"       ✓ Pozostało {fmt(len(df_final))} obiektów.")
    
    # 9. DISPLAY & SAVE
    MAX_PRINT_ROWS = 100  # Ustaw limit (pasuje do Twojego if i >= 300)
    
    print("="*119)
    print(f"{'Nr':>3} {'ID':<12} {'Extra':<25} {'Common name':<25} {'Typ':<8} "
          f"{'RA':>7} {'Dec':>7} {'Mag':>7} {'Size':>8} {'Score':>6}")
    print("-"*119)
    
    json_objects = []
    table_truncated = False
    
    for i, (idx, row) in enumerate(df_final.iterrows()):
        if i >= MAX_PRINT_ROWS:
            print("=" * 119)
            print(f"{f'... i {len(df_final) - MAX_PRINT_ROWS} więcej obiektów (limit wyświetlania).':>119}")
            table_truncated = True
            break
               
        bd = row["score_breakdown"]
        mag_str = f"{'●' if row['mag_source']=='measured' else '○'}{row['mag']:.1f}"
        size_str = f"{'●' if row['size_source']=='measured' else '○'}{row['size']:.1f}'"
        ra_str = f"{row['ra']:.1f}"
        dec_str = f"{row['dec']:.1f}"
        
        print(f"{i+1:3d}. {str(row['id']):<12} {str(row['extra_info'])[:25]:<25} "
              f"{str(row.get('common_names',''))[:25]:<25} {str(row['type_code']):<8} "
              f"{ra_str:>7} {dec_str:>7} {mag_str:>7} {size_str:>8} {row['final_score']:6.1f}")
    
    if table_truncated:
            print("="*119)
            
    # --- GENEROWANIA JSON 
    json_objects = []
    
    for _, row in df_final.iterrows():
        hours_data = row["hours"] if isinstance(row["hours"], list) else np.array(row["hours"], dtype=float).tolist()
        
        is_fam = bool(row.get("is_messier", False) or row.get("is_caldwell", False) or row.get("is_herschel", False))

        json_objects.append({
            # --- Identyfikacja ---
            "id": str(row["id"]),
            "name": str(row["id"]), 
            "common_names": str(row.get("common_names", "")).strip(),
            "extra_info": str(row.get("extra_info", "")).strip(),
            "type": str(row["type_code"]), 

            # --- Dane astronomiczne ---
            "ra": float(row["ra"]),
            "dec": float(row["dec"]),
            "mag": float(row["mag"]) if pd.notna(row["mag"]) else 0.0,
            "size": float(row["size"]) if pd.notna(row["size"]) else 0.0,

            # --- Punktacja i Widoczność (Mapowanie zmiennych) ---
            "score": float(row["final_score"]),          
            "max_altitude_year": float(row["maxalt"]),   
            "monthly_hours_total": hours_data,           
            "months_observable": int(row["monthsobservable"]), 
            "rarity_score": float(row.get("rarity_score", 0.0)), 

            # --- Flagi (Wymagają rzutowania na bool) ---
            "is_messier": bool(row.get("is_messier", False)),
            "is_caldwell": bool(row.get("is_caldwell", False)),
            "is_herschel": bool(row.get("is_herschel", False)),
            "is_catalog_famous": is_fam
        })

    # Save JSON
    out_data = {
        "location": loc,
        "year": year,
        "parameters": params,
        "objects": json_objects
    }
    
    with open(Config.OUTPUT_FILE, "w", encoding="utf-8") as f:
        json.dump(out_data, f, indent=2, ensure_ascii=False)
        
    print(f"[INFO] Zapisano {fmt(len(json_objects))} obiektów do {Config.OUTPUT_FILE}")

if __name__ == "__main__":
    main()