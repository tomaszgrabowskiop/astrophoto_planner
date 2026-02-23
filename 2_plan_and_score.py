#!/usr/bin/env python3
# =============================================================================
# 2_plan_and_score.py  –  Astrophotography Planner  |  Etap 2: Planowanie i scoring
# =============================================================================
"""
Krok 2: Konfiguracja sesji i system punktacji (Scoring).

Ten skrypt pełni dwie kluczowe funkcje: interfejsu konfiguracyjnego oraz silnika oceniającego atrakcyjność obiektów.

Działanie:
1. Konfiguracja (Interaktywna):
   - Użytkownik podaje rok planowania, lokalizację (City/Lat/Lon).
   - Definiuje sprzęt (Kamera/Teleskop) - kluczowe dla obliczeń FOV.
   - Określa warunki brzegowe: minimalna wysokość obiektu, poziom zanieczyszczenia światłem (Bortle), filtry.

2. Scoring (Punktacja):
   - Każdy obiekt z katalogu otrzymuje punkty w oparciu o algorytm uwzględniający:
     * Typ obiektu vs Filtry (np. mgławice emisyjne premiowane przy filtrach narrowband).
     * Jasność powierzchniową i rozmiar.
     * Przynależność do list "The Best of" (Messier, Caldwell, Herschel 400).
     * Warunki lokalne (Bortle).

Wyjście:
- Plik: vis_data.json (zawiera obiekty z przypisanymi punktami i pełną konfigurację użytkownika).
"""
import numpy as np
import pandas as pd
import pytz
import json
import re 
from tqdm import tqdm
from pathlib import Path
from typing import Tuple, List, Dict, Any
import hashlib
import ssl
import certifi

import astropy.units as u
from astropy.coordinates import SkyCoord, EarthLocation
from astropy.time import Time
from astropy.utils.exceptions import AstropyWarning
from astroplan import Observer
from datetime import date, datetime, timedelta, timezone
from timezonefinder import TimezoneFinder

import warnings
from erfa import ErfaWarning

from shared import CATALOG_PRIORITY, CAT_ORDER, PATHS, fmt, print_step, print_green, CameraConfig

# =====================================================
# KONFIGURACJA PASKA, OSTRZEŻEŃ I SSL
# =====================================================

tqdm.pandas(colour="green",ncols=119)
warnings.filterwarnings("ignore", category=ErfaWarning)
warnings.filterwarnings("ignore", category=AstropyWarning)
ssl._create_default_https_context = lambda: ssl.create_default_context(cafile=certifi.where())

# =====================================================
# TABELE PUNKTACJI (SCORECARD)
# =====================================================

# TABELA 1: WPŁYW FILTRÓW
# Wartość: [bez narrowband (triband/brak), z narrowband (H/O/S)]
SCORING_FILTER = {
    "GX":     [15,  0],  "G":      [15,  0],
    "DN":     [ 0,  0],  "NB":     [20, 30],
    "OCL":    [20,  0],  "HII":    [20, 30],
    "*":      [10,  0],  "OTHER":  [10,  0],
    "DUP":    [10,  0],  "GCL":    [20,  0],
    "**":     [10,  0],  "GPAIR":  [15,  0],
    "NEB":    [20, 30],  "*ASS":   [10,  0],
    "CL+N":   [20, 30],  "RFN":    [ 0,  0],
    "PN":     [20, 30],  "GTRPL":  [15,  0],
    "SNR":    [20, 30],  "GGROUP": [15,  0],
    "NOVA":   [10,  0],
}

# TABELA 2: WPŁYW NIEBA (BORTLE)
# Wartość: [Bortle ≤5 (ciemne niebo), Bortle >5 (miasto/LP)]
SCORING_BORTLE = {
    "GX":     [15, 15],  "G":      [15, 15],
    "DN":     [30,  0],  "NB":     [20, 20],
    "OCL":    [20, 20],  "HII":    [20, 20],
    "*":      [10, 10],  "OTHER":  [10, 10],
    "DUP":    [10, 10],  "GCL":    [20, 20],
    "**":     [10, 10],  "GPAIR":  [15, 15],
    "NEB":    [20, 20],  "*ASS":   [10, 10],
    "CL+N":   [20, 20],  "RFN":    [30,  0],
    "PN":     [20, 20],  "GTRPL":  [15, 15],
    "SNR":    [10, 10],  "GGROUP": [15, 15],
    "NOVA":   [10, 10],
}

# TABELA 3: BONUSY ZA "SŁAWĘ" (Famous)
# Naliczamy tylko jedną nagrodę: najwyższa kategoria wygrywa (hierarchia M > C > H)
SCORE_FAMOUS = {
    "messier":  25,
    "caldwell": 20,
    "herschel": 15,
}

# TABELA 4: NAGRODY ZA JAKOŚĆ DANYCH
# Nagroda jest przyznawana za pole, które NIE było zgadywane (imputed == False)
SCORE_DATA_QUALITY = {
    "measured": 15,
}

# =====================================================
# FUNKCJE POMOCNICZE
# =====================================================

def calculate_surface_brightness(mag: float, size_arcmin: float) -> float:
    """Wylicza jasność powierzchniową [mag/arcmin²]. Zwraca 99.0 dla braku danych."""
    if size_arcmin <= 0:
        return 99.0
    area = np.pi * (size_arcmin / 2.0) ** 2
    if area <= 0:
        return 99.0
    return mag + 2.5 * np.log10(area)

def add_parameters_hash_to_output(out_data: Dict[str, Any], params: Dict[str, Any]) -> None:
    """
    Dodaje hash parametrów do słownika output (używany przez krok 3_ do cache).
    Hash obejmuje tylko min_alt i sun_limit — jedyne parametry wpływające na maskę widoczności.
    """
    obj_min_alt   = params.get("min_alt", 20.0)
    sun_alt_limit = params.get("sun_limit", -6.0)
    s = f"{obj_min_alt:.6f}|{sun_alt_limit:.6f}"
    out_data["parameters_hash"] = hashlib.md5(s.encode()).hexdigest()

def _as_int0(x: Any) -> int:
    """Bezpieczna konwersja do int (np. z NaN) w celu ochrony przed błędem rzutowania."""
    return 0 if pd.isna(x) else int(x)

def _as_bool(x: Any) -> bool:
    """Bezpieczna konwersja do bool (zapobiega potraktowaniu stringa 'False' jako True)."""
    if pd.isna(x):
        return False
    if isinstance(x, bool):
        return x
    if isinstance(x, (int, float)):
        return bool(int(x))
    return str(x).strip().lower() in ("1", "true", "t", "yes", "y")

def format_indeksy(extra_info: str) -> str:
    """
    Sortuje zawartość extra_info według priorytetów katalogów z CAT_ORDER z shared.py.
    Wyodrębnia prefiksy (np. NGC, Sh2) i w obrębie tego samego katalogu sortuje numerycznie.
    """
    if pd.isna(extra_info) or str(extra_info).strip().lower() in ('nan', 'none', ''):
        return ""

    items = [item.strip() for item in str(extra_info).split(',') if item.strip()]
    if not items:
        return ""

    def cat_rank(item: str) -> tuple:
        key = item.lower()
        
        # 1. Wyodrębnienie prefiksu katalogu
        if key.startswith("sh2"):
            prefix = "sh2"
        elif re.match(r'^b\d+', key):
            prefix = "barn"  # Mapowanie katalogu Barnarda (np. B33 -> barn)
        else:
            # Dowolny inny katalog: wyciągamy same litery z początku (np. 'ic' z 'ic1795')
            m = re.match(r"^[a-z]+", key)
            prefix = m.group(0) if m else key[:3]
            
        # 2. Ustalenie priorytetu na podstawie CAT_ORDER
        try:
            cat_idx = CAT_ORDER.index(prefix)
        except ValueError:
            cat_idx = len(CAT_ORDER)  # Nieznany katalog spada na sam koniec
            
        # 3. Wyciągnięcie numeru do sortowania wewnętrznego (np. NGC20 -> 20)
        match = re.search(r'\d+', key)
        num = int(match.group()) if match else 999999
        
        # Sortowanie: 1. Priorytet z CAT_ORDER, 2. Numer obiektu, 3. Oryginalny string (fallback)
        return (cat_idx, num, key)

    items.sort(key=cat_rank)
    return ', '.join(items)

def calculate_rich_field_bonus(extra_info: str) -> float:
    """
    Oblicza Bonus Bogactwa Kadru na podstawie obiektów wchłoniętych (w extra_info).
    Punkty wg wagi katalogu:
      +3.0 pkt: NGC, IC
      +2.0 pkt: Sh2, RCW, LBN, Ced
      +1.0 pkt: Barnard (B)
      +0.5 pkt: PGC, LDN, reszta
    Maksymalny bonus wynosi +15.0 punktów.
    """
    if not extra_info or str(extra_info).strip().lower() in ('nan', 'none', ''):
        return 0.0

    bonus = 0.0
    items = [item.strip().lower() for item in str(extra_info).split(',') if item.strip()]
    
    for item in items:
        if item.startswith(('ngc', 'ic')):
            bonus += 3.0
        elif item.startswith(('sh2', 'rcw', 'lbn', 'ced')):
            bonus += 2.0
        elif re.match(r'^b\d+', item):  # Barnard
            bonus += 1.0
        else:  # PGC, LDN i wszystko inne
            bonus += 0.5
            
    # Zwracamy bonus, ale nie więcej niż "Capping" = 15.0
    return min(bonus, 15.0)

# =====================================================
# SCORING
# =====================================================

def calculate_scorecard(row: pd.Series, user_params: Dict[str, Any]) -> Dict[str, Any]:
    """
    Oblicza punktację w systemie Scorecard (Karta Wyników).
    Suma punktów z 4 kategorii: Sława + Typ/Filtr + Typ/Bortle + Jakość Danych.
    """
    obj_type = str(row.get("type", "OTHER")).upper().strip()
    if obj_type not in SCORING_FILTER:
        obj_type = "OTHER"

    current_score = 0.0
    breakdown = {}

    # A. SŁAWA (FAMOUS BONUS)
    m = _as_int0(row.get("messier_nr"))
    c = _as_int0(row.get("caldwell_nr"))
    h = _as_int0(row.get("herschel_nr"))
    
    famous_points = 0.0
    if user_params.get("prefer_famous", True):
        if m > 0:
            famous_points = SCORE_FAMOUS["messier"]
        elif c > 0:
            famous_points = SCORE_FAMOUS["caldwell"]
        elif h > 0:
            famous_points = SCORE_FAMOUS["herschel"]
        
    current_score += famous_points
    breakdown["famous"] = famous_points

    # B. TYP + FILTR 
    filter_idx = 1 if user_params.get("has_narrowband", False) else 0
    filter_points = SCORING_FILTER[obj_type][filter_idx]
    current_score += filter_points
    breakdown["filter_type"] = filter_points

    # C. TYP + BORTLE 
    bortle_val = user_params.get("bortle_range", (6, 7))[1]
    bortle_idx = 0 if bortle_val <= 5 else 1
    bortle_points = SCORING_BORTLE[obj_type][bortle_idx]
    current_score += bortle_points
    breakdown["bortle_type"] = bortle_points

    # D. JAKOŚĆ DANYCH
    mag_imp  = _as_bool(row.get("mag_imputed", False))
    size_imp = _as_bool(row.get("size_imputed", False))

    quality_points = 0.0
    quality_points += 0.0 if mag_imp else SCORE_DATA_QUALITY["measured"]
    quality_points += 0.0 if size_imp else SCORE_DATA_QUALITY["measured"]
    current_score += quality_points
    breakdown["quality"] = quality_points

    return {
        "final_score": current_score,
        "breakdown":   breakdown,
    }
# =====================================================
# ASTRO MATH & SMART GRID (Widoczność Hybrydowa)
# =====================================================

class VectorAstro:
    @staticmethod
    def get_gmst_vector(jd_array: np.ndarray) -> np.ndarray:
        T = (jd_array - 2451545.0) / 36525.0
        gmst = (280.46061837 + 360.98564736629 * (jd_array - 2451545.0) +
                0.000387933 * T**2 - T**3 / 38710000.0)
        return gmst % 360.0

    @staticmethod
    def alt_az_vector(ha_rad: np.ndarray, dec_rad: np.ndarray, lat_rad: float) -> np.ndarray:
        sin_alt = (np.sin(dec_rad) * np.sin(lat_rad) +
                   np.cos(dec_rad) * np.cos(lat_rad) * np.cos(ha_rad))
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
    print_step(f"Generowanie SMART GRID dla roku {year}, limit słońca: {sun_limit_deg}°")

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
    with tqdm(total=n_days, desc=" SmartGrid", unit="day", ncols=119, colour="green") as pbar:
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
        print(" [!] Ostrzeżenie: Nie znaleziono żadnych nocy (błąd parametrów?).")
        return np.array([]), np.array([]), np.array([]), np.array([])

    jd_array = np.concatenate(all_jds)
    night_id_array = np.concatenate(all_night_ids)

    print(f"       ✓ Wygenerowano {fmt(len(jd_array))} punktów pomiarowych dla {fmt(night_counter)} nocy "
          f"(~{len(jd_array)*step_minutes/60.0:.0f}h obserwacji).")

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
    """
    if df_objects.empty or len(jd_array) == 0:
        return np.zeros((len(df_objects), 12)), np.zeros(len(df_objects), dtype=int)

    n_objects = len(df_objects)
    print_step(f"Obliczanie macierzowe dla {fmt(n_objects)} obiektów. (To może potrwać kilkanaście sekund.)")

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
        night_mask = (night_ids == night_id)
        hours_this_night = np.sum(weighted_mask[night_mask], axis=0)
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

    print("       ✓ Obliczono widoczność per noc i per miesiąc.")
    return hours_per_month, nights_above_threshold

def compute_max_altitude(lat: float, dec: float) -> float:
    lat_rad = np.radians(lat)
    dec_rad = np.radians(dec)
    sin_alt_max = np.sin(lat_rad) * np.sin(dec_rad) + np.cos(lat_rad) * np.cos(dec_rad)
    return np.degrees(np.arcsin(np.clip(sin_alt_max, -1.0, 1.0)))

# =====================================================
# USER INTERFACE / MANAGERS
# =====================================================
class SessionConfigManager:
    def select_interactive(self) -> Tuple[int, Dict[str, Any]]:
        print_green("\n" + "=" * 119)
        print_green("KROK 1: WYBÓR ROKU I LOKALIZACJI")
        print_green("=" * 119)
        current_year = date.today().year
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

        print_step("Wybierz opcję lub wpisz własne miejsce:")
        print("1. Poznań, Polska")
        print("2. Kraków, Polska")
        print("3. Wpisz lokalizację")
        choice = input("[1/2/3, domyślnie 1]: ").strip() or "1"
        
        if choice == "1":
            loc = {"lat": 52.4095, "lon": 16.9319, "tz": "Europe/Warsaw", "name": "Poznań, Polska"}
            print(f"Wybrano {loc['name']} ({loc['lat']}N, {loc['lon']}E)")
            return year, loc
            
        if choice == "2":
            loc = {"lat": 50.0647, "lon": 19.9450, "tz": "Europe/Warsaw", "name": "Kraków, Polska"}
            print(f"Wybrano {loc['name']} ({loc['lat']}N, {loc['lon']}E)")
            return year, loc
            
        city = input("Podaj lokalizację (pisz wielką literą, działa lepiej, np. 'Toruń, Polska'): ").strip()
        if not city:
            print("Brak nazwy, używam domyślnej: Poznań.")
            loc = {"lat": 52.4095, "lon": 16.9319, "tz": "Europe/Warsaw", "name": "Poznań, Polska"}
            return year, loc
        try:
            print(f"[INFO] Szukanie współrzędnych dla '{city}'...")
            loc_astro = EarthLocation.of_address(city)
            lat = loc_astro.lat.to(u.deg).value
            lon = loc_astro.lon.to(u.deg).value
            
            tf = TimezoneFinder()
            tz_name = tf.timezone_at(lng=lon, lat=lat)
            if tz_name is None:
                tz_name = "Europe/Warsaw"
                
            loc = {"lat": lat, "lon": lon, "tz": tz_name, "name": city}
            print(f"[INFO] Znaleziono: {loc['name']} {loc['lat']:.4f}N, {loc['lon']:.4f}E, strefa czasowa: {tz_name}")
            return year, loc
            
        except Exception as e:
            print(f"Błąd geokodowania ({e}), używam Poznań.")
            loc = {"lat": 52.4095, "lon": 16.9319, "tz": "Europe/Warsaw", "name": "Poznań, Polska"}
            return year, loc

def get_user_prefs() -> Dict[str, Any]:
    cam = CameraConfig()
    print_green("\n" + "=" * 119)
    print_green("KROK 2: PARAMETRY WIDOCZNOŚCI I FOV")
    print_green("=" * 119)
    
    try:
        min_alt = float(input("A. Minimalna wysokość obiektu nad horyzontem [domyślnie 25°]: ") or 25.0)
    except Exception: min_alt = 25.0
    
    try:
        min_hours = float(input("\nB. Minimalna liczba godzin, którą obiekt jest widoczny w nocy powyżej progu wysokości [domyślnie 3]: ") or 3.0)
    except Exception: min_hours = 3.0

    print("\nC. Ciemność nieba – kąt słońca pod horyzontem:")
    print("   • zmierzch cywilny (-6°)")
    print("   • zmierzch żeglarski (-12°)")
    print("   • zmierzch astronomiczny (-18°)")
    sl_choice = input("   Możesz wpisać dowolną sensowną wartość (liczba stopni pod horyzontem) [domyślnie 12]: ").strip() or "12"
    try: sun_limit = -abs(float(sl_choice))
    except Exception: sun_limit = -12.0

    print("\nD. Określenie FOV")
    use_default = input("   Użyć domyślnego setupu RedCat61 + ASI2600MC Pro? [t/n, domyślnie t]: ").strip().lower() or "t"
    
    if use_default != "t":
        try:
            cam.lens_focal_length = float(input("   • Ogniskowa [mm, domyślnie 300]: ").strip())
        except Exception:
            pass  
        try:
            cam.sensor_width = float(input("   • Szerokość sensora [mm, domyślnie 23.5]: ").strip())
            cam.sensor_height = float(input("   • Wysokość sensora [mm, domyślnie 15.7]: ").strip())
        except Exception:
            pass
        try:
            cam.sensor_pitch = float(
                input("   • Wielkość piksela [µm, domyślnie 3.76]: ").strip())
        except Exception:
            pass
        try:
            cam.sensor_rows = int(
                input("   • Liczba wierszy sensora [domyślnie 4176]: ").strip())
        except Exception:
            pass
        try:
            cam.sensor_cols = int(
                input("   • Liczba kolumn sensora [domyślnie 6248]: ").strip())
        except Exception:
            pass
        print(f"[INFO] Ogniskowa: {cam.lens_focal_length}mm," 
            f"sensor: {cam.sensor_width}x{cam.sensor_height}mm, "
            f"pitch: {cam.sensor_pitch}µm, "
            f"{cam.sensor_cols}x{cam.sensor_rows}px"
            )

    fov_w, fov_h, _ = cam.calculate_fov()
    print_step(f"Szacowany FOV: {fov_w:.2f}° x {fov_h:.2f}°")

    try:
        percent_fov = float(input("\nE. Minimalny rozmiar obiektu jako % krótszego boku FOV [domyślnie 10]: ") or 10.0)
    except Exception:
        percent_fov = 10.0
        
    min_size_arcmin = cam.get_min_match_size_arcmin(percent_fov)
    print(f"[INFO] Minimalny rozmiar obiektu (arcmin): {min_size_arcmin:.1f}'.")

    print("\nF. Skala Bortle – określenie zanieczyszczenia światłem:")
    print("   1) Bortle 1–3 (wieś, ciemne niebo, pomijalne LP)")
    print("   2) Bortle 4–5 (przedmieścia, umiarkowane LP)")
    print("   3) Bortle 6–7 (miasto, silne LP)")
    print("   4) Bortle 8–9 (centrum miasta, ekstremalne LP)")
    b_choice = input("   Wybierz Twój stopień zanieczyszczenia nieba światłem [domyślnie 4]: ").strip() or "4"
    b_map = {"1": (1, 3), "2": (4, 5), "3": (6, 7), "4": (8, 9)}
    bortle_rng = b_map.get(b_choice, (8, 9))

    has_nb = (input("\nG. Czy zamierzasz korzystać z filtrów narrowband (H/O/S)? [t/n, domyślnie n]: ").strip().lower() or "n") == "t"
    prefer_famous = (input("\nH. Czy w wyborze obiektów premiować katalogi Messier/Caldwell/Herschel? [t/n, domyślnie t]: ").strip().lower() or "t") == "t"

    return {
        "min_alt": min_alt,
        "sun_limit": sun_limit,
        "min_hours": min_hours,
        "min_size_arcmin": min_size_arcmin,
        "bortle_range": bortle_rng,
        "has_narrowband": has_nb,
        "prefer_famous": prefer_famous,
        "camera": cam.to_dict(),
    }
# =====================================================
# OPTYMALIZACJA ATLASU (CLUSTERING)
# =====================================================

def choose_best_catalog_id(main_id: str, extra_info: str) -> str:
    """
    Wybiera najważniejsze ID na podstawie katalogu (np. NGC > IC).
    """
    def cat_rank(obj_id: str) -> int:
        s = obj_id.lower().strip()
        for idx, c in enumerate(CAT_ORDER):
            if s.startswith(c):
                return idx
        return 99

    candidates = [str(main_id).strip()] + [x.strip() for x in str(extra_info).split(",") if x.strip()]
    return sorted(candidates, key=lambda x: (cat_rank(x), x))[0]

def optimize_atlas_pages(df: pd.DataFrame, user_params: Dict[str, Any]) -> pd.DataFrame:
    """
    Zastępuje statyczny graf. Używa dynamicznego FOV użytkownika (80% krótszego boku).
    Grupuje bliskie cele w jeden kadr dla atlasu PDF, chroniąc wielkie obiekty przed wchłanianiem.
    """
    print_step("Optymalizacja kadrów (Clustering z dynamicznym FOV)")
    
    if df.empty:
        return df

    # Wyciągnij parametry aparatu do FOV
    c = user_params.get("camera", {})
    cam = CameraConfig.from_dict(c)
    fov_w, fov_h, _ = cam.calculate_fov()
    
    # 80% krótszego boku to bezpieczny margines — kadr pomieści obiekty w rogach
    clustering_radius_deg = min(fov_w, fov_h) * 0.80
    
    print(f"       • FOV sprzętu: {fov_w:.2f}° x {fov_h:.2f}°")
    print(f"       • Promień klastrowania ustalony na: {clustering_radius_deg:.2f}°")

    df_sorted = df.sort_values(by=["final_score", "size"], ascending=[False, False]).reset_index(drop=True)
    coords = SkyCoord(ra=df_sorted["ra"].values * u.deg, dec=df_sorted["dec"].values * u.deg)

    final_objs = []
    covered = set()
    MAX_CLUSTER_SIZE_ARCMIN = 300.0 
    
    # Próg bezpieczeństwa: obiekt większy niż 50% promienia nie da się wchłonąć
    protection_size_arcmin = (clustering_radius_deg * 60.0) * 1.0

    for i in range(len(df_sorted)):
        if i in covered:
            continue

        leader = df_sorted.iloc[i].copy()
        leader_coord = coords[i]
        
        # Znajdź sąsiadów w promieniu
        seps = coords.separation(leader_coord).deg
        neighbors_idx = np.where((seps <= clustering_radius_deg) & (seps > 0))[0]
        
        # Odrzuć zużyte ORAZ chronione, ogromne obiekty (aby nie zostały pożarte przez lidera punktowego)
        valid_neighbors_idx = []
        leader_size = float(leader.get("size", 0.0))
        
        for n in neighbors_idx:
            if n not in covered:
                neighbor_size = float(df_sorted.iloc[n].get("size", 0.0))
                
                # Zezwól na połączenie jeśli:
                # 1. Sąsiad jest mniejszy niż próg bezpieczeństwa
                # 2. LUB Lider jest mniejszy niż próg bezpieczeństwa (pozwalamy wchłonąć giganta przez zwarty punkt)
                if (neighbor_size < protection_size_arcmin) or (leader_size < protection_size_arcmin):
                    valid_neighbors_idx.append(n)

        if not valid_neighbors_idx:
            final_objs.append(leader)
            covered.add(i)
            continue

        subset = df_sorted.iloc[valid_neighbors_idx]
        
        # 1. ROZMIAR (limitowany)
        sizes = subset["size"].dropna()
        if not sizes.empty:
            reasonable = sizes[sizes <= MAX_CLUSTER_SIZE_ARCMIN]
            if not reasonable.empty:
                leader["size"] = float(max(leader["size"], reasonable.max()))

        # 2. IDENTYFIKATORY i EXTRA INFO
        cluster_ids = set([str(leader["id"])])
        if pd.notna(leader.get("extra_info")):
            for x in str(leader["extra_info"]).split(","):
                if x.strip(): cluster_ids.add(x.strip())

        for _, row_n in subset.iterrows():
            cluster_ids.add(str(row_n["id"]))
            if pd.notna(row_n.get("extra_info")):
                for x in str(row_n["extra_info"]).split(","):
                    if x.strip(): cluster_ids.add(x.strip())

        cluster_ids = {x for x in cluster_ids if x and x != "nan"}

        # 3. COMMON NAMES
        c_names = set()

        def _add_common_name(val):
            """Dodaje common_name tylko jeśli nie jest NaN/None/pusty/'nan'."""
            if val is None:
                return
            s = str(val).strip()
            if not s:
                return
            if s.lower() in ("nan", "none"):
                return
            c_names.add(s)

        # lider klastra
        _add_common_name(leader.get("common_names"))

        # pozostali z klastra
        for _, row_n in subset.iterrows():
            _add_common_name(row_n.get("common_names"))

        leader["common_names"] = ", ".join(sorted(c_names)) if c_names else ""

        # 4. SŁAWA (MAX, aby zachować flagę na kadrze)
        m_max = _as_int0(leader.get("messier_nr"))
        c_max = _as_int0(leader.get("caldwell_nr"))
        h_max = _as_int0(leader.get("herschel_nr"))
        
        for _, row_n in subset.iterrows():
            m_max = max(m_max, _as_int0(row_n.get("messier_nr")))
            c_max = max(c_max, _as_int0(row_n.get("caldwell_nr")))
            h_max = max(h_max, _as_int0(row_n.get("herschel_nr")))

        leader["messier_nr"] = m_max
        leader["caldwell_nr"] = c_max
        leader["herschel_nr"] = h_max

        # 5. BEST ID (Mądry wybór wg wielkości/jasności w klastrze)
        candidates_for_best = []
        
        # 5a. Budowa puli kandydatów (leader + subset)
        for _, c_row in pd.concat([pd.DataFrame([leader]), subset]).iterrows():
            c_id = str(c_row.get("id", "")).strip()
            if not c_id or c_id.lower() == "nan":
                continue
                
            # Interesują nas tylko "główne" katalogi jako kandydaci na etykietę kadru
            if not c_id.upper().startswith(("NGC", "IC", "SH2-")):
                continue
                
            c_size = float(c_row.get("size", 0.0)) if pd.notna(c_row.get("size")) else 0.0
            c_mag = float(c_row.get("mag", 99.0)) if pd.notna(c_row.get("mag")) else 99.0
            
            candidates_for_best.append({
                "id": c_id,
                "size": c_size,
                "mag": c_mag
                })

        # 5b. Wybór najlepszego ID wg [size -> mag -> choose_best_catalog_id]
        if candidates_for_best:
            # Sortuj malejąco po rozmiarze, a przy remisie rosnąco po jasności (mag)
            candidates_for_best.sort(key=lambda x: (x["size"], -x["mag"]), reverse=True)
            # Bierzemy cały słownik zwycięzcy
            best_id_candidate = candidates_for_best[0]["id"]
            
            # Weryfikacja remisu dla bezpieczeństwa (rzadkie, ale jeśli dwa giganty są identyczne)
            tied_candidates = [c["id"] for c in candidates_for_best 
                               if c["size"] == candidates_for_best[0]["size"] 
                               and c["mag"] == candidates_for_best[0]["mag"]]
            
            if len(tied_candidates) > 1:
                # Jeśli jest idealny remis size i mag, użyj starej funkcji z CAT_ORDER
                combined_tied = ",".join(sorted(tied_candidates[1:]))
                best_id = choose_best_catalog_id(tied_candidates[0], combined_tied)
            else:
                best_id = best_id_candidate
        else:
            # Fallback dla egzotycznych klastrów (bez NGC/IC/SH2) - użyj starej logiki
            combined_extra = ",".join(sorted(cluster_ids))
            best_id = choose_best_catalog_id(str(leader["id"]), combined_extra)

        # 5c. Przypisanie wyniku
        if best_id in cluster_ids:
            cluster_ids.remove(best_id)
        if best_id != str(leader["id"]).strip():
            # Szukamy oryginalnego wiersza zwycięzcy w puli obiektów z tego klastra
            pool = pd.concat([pd.DataFrame([leader]), subset])
            winner_match = pool[pool["id"].astype(str).str.strip() == best_id]
            
            if not winner_match.empty:
                w_row = winner_match.iloc[0]
                leader["type"] = w_row["type"]
                leader["ra"]   = float(w_row["ra"])
                leader["dec"]  = float(w_row["dec"])
                # Jasność nadpisujemy tylko, jeśli obiekt ją posiadał
                if pd.notna(w_row.get("mag")):
                    leader["mag"] = float(w_row["mag"])
                leader["mag_imputed"]  = _as_bool(w_row.get("mag_imputed", False))
                leader["size_imputed"] = _as_bool(w_row.get("size_imputed", False))

        leader["id"] = best_id
        raw_extra = ",".join(cluster_ids) if cluster_ids else ""
        leader["extra_info"] = format_indeksy(raw_extra)
        
        #  6. BONUS BOGACTWA KADRU (Rich Field Bonus)
        rich_bonus = calculate_rich_field_bonus(leader["extra_info"])
        
        if rich_bonus > 0:
            current_score = float(leader.get("final_score", 0.0))
            leader["final_score"] = current_score + rich_bonus
            # Bezpieczna aktualizacja szczegółów punktacji (breakdown) dla JSON/diagnostyki
            breakdown = leader.get("score_breakdown")
            if isinstance(breakdown, dict):
                breakdown["rich_field_bonus"] = rich_bonus
            leader["score_breakdown"] = breakdown
            
        # Oznacz jako załatwione
        covered.add(i)
        for n_idx in valid_neighbors_idx:
            covered.add(n_idx)

        final_objs.append(leader)

    df_clustered = pd.DataFrame(final_objs)
    print(f"       ✓ Zredukowano kadrów z {fmt(len(df))} do {fmt(len(df_clustered))}.")
    
    return df_clustered
# =====================================================
# ADAPTIVE SOFT CUT
# =====================================================

def apply_adaptive_soft_cut(
    df: pd.DataFrame, 
    base_min_score: float = 10.0,
    min_len_for_cut: int = 500,
    keep_top_percent: float = 70.0,
) -> pd.DataFrame:
    """
    Adaptacyjny softcut:
    1) Odrzuca wszystko poniżej base_min_score (twardy, ale niski próg bezpieczeństwa).
    2) Jeśli po tym zostaje dużo obiektów (> min_len_for_cut), ucina dół rozkładu tak, 
       żeby zostawić tylko górne keep_top_percent%.
    """
    if df.empty or "final_score" not in df.columns:
        return df

    # 1. Twardy minimalny próg 
    before = len(df)
    df = df[df["final_score"] >= base_min_score].reset_index(drop=True)
    after = len(df)
    print_step(f"Zastosowano ostre cięcie dla obiektw z wynikiem < {base_min_score:.1f} pkt.")
    print(f"[INFO] Pozostało {fmt(after)} obiektów z {fmt(before)}.")
    
    if df.empty:
        return df

    # 2. Adaptacyjny cut względem rozkładu (tylko jeśli lista jest wystarczająco długa)
    if len(df) < min_len_for_cut:
        print_step(f"Soft-cut: lista < {min_len_for_cut} obiektów — pomijam cięcie percentylowe.")
        return df

    scores = df["final_score"].astype(float)
    lower_tail_percent = 100.0 - keep_top_percent
    
    # Dolny percentyl, który utniesz (np. 30% na dole, zostawiasz górne 70%)
    pcut = np.percentile(scores, lower_tail_percent)
    
    before2 = len(df)
    df = df[df["final_score"] >= pcut].reset_index(drop=True)
    after2 = len(df)
    
    print(f"[INFO] Wartości wyliczone przez cięcie adaptacyjne: "
          f"final_score >= p{lower_tail_percent:.0f} ({pcut:.1f}).")
    print(f"[INFO] Pozostało {fmt(after2)} obiektów z {fmt(before2)}.")

    return df

def _clean_str(val):
    s = str(val).strip()
    return "" if s.lower() in ["nan", "none", ""] else s

# =====================================================
# GŁÓWNA PĘTLA APLIKACJI
# =====================================================

def main():
    # ---------------------------------------------------------
    # 1. SETUP UX / WYBÓR PARAMETRÓW
    # ---------------------------------------------------------
    year, loc = SessionConfigManager().select_interactive()
    params = get_user_prefs()

    # ---------------------------------------------------------
    # 2. LOAD DATA
    # ---------------------------------------------------------
    print_green("\n" + "=" * 119)
    print_green("KROK 3. ŁADOWANIE I FILTROWANIE")
    print_green("=" * 119)
    
    if not PATHS.catalog_full.exists():
        print(f"Brak pliku katalogu: {PATHS.catalog_full}!")
        return

    # Wymuszamy typy bool dla flag wygenerowanych w kroku 1, aby uniknąć pomyłki "False" -> True
    df = pd.read_csv(
        PATHS.catalog_full, 
        low_memory=False,
        converters={
            "size_imputed": _as_bool,
            "mag_imputed": _as_bool
        }
    )

    # ---------------------------------------------------------
    # 3. FAST GEOMETRIC FILTER
    # ---------------------------------------------------------
    print_step(f"Filtr geometryczny (Alt >= {params['min_alt']}°, Visible >= {params['min_hours']}h)")
    total_objects = len(df)

    lat_rad = np.radians(loc["lat"])
    dec_rad = np.radians(df["dec"].values)
    
    # Szybki Max Alt check
    sin_alt_max = np.sin(lat_rad)*np.sin(dec_rad) + np.cos(lat_rad)*np.cos(dec_rad)
    alt_max = np.degrees(np.arcsin(np.clip(sin_alt_max, -1, 1)))

    df_after_alt = df[alt_max >= params["min_alt"]].copy()
    print_step(f"Po filtrze min. wysokości ({params['min_alt']}°) "
          f"pozostało {fmt(len(df_after_alt))} z {fmt(total_objects)} obiektów.")
    
    df = df_after_alt
    if df.empty:
        print("Brak obiektów.")
        return

    # ---------------------------------------------------------
    # 4. SMART GRID & VISIBILITY (HYBRYDOWA)
    # ---------------------------------------------------------
    step = 15
    jd, months, lst, night_ids = get_full_year_smart_grid(
        year, loc, params["sun_limit"], step
    )
    if len(jd) == 0:
        print("Brak nocy do analizy!")
        return

    # Oblicz widoczność per noc / per miesiąc
    hours_per_month, nights_above_threshold = compute_all_visibilities_hybrid(
        df, loc, jd, night_ids, months, lst, 
        params["min_alt"], params["min_hours"], step
    )

    # ---------------------------------------------------------
    # 5. FAST FILTER (Odrzuć obiekty bez ŻADNEJ nocy z min_hours)
    # ---------------------------------------------------------
    objects_before_filter = len(df)
    mask_has_good_night = (nights_above_threshold >= 1)
    
    df = df[mask_has_good_night].reset_index(drop=True)
    hours_per_month = hours_per_month[mask_has_good_night]
    nights_above_threshold = nights_above_threshold[mask_has_good_night]

    print(f"       ✓ Spośród {fmt(objects_before_filter)} obiektów, {fmt(len(df))} przynajmniej w jedną noc "
          f"są widoczne przez >= {params['min_hours']}h.")

    if df.empty:
        print("Brak obiektów po filtrze czasu widoczności.")
        return

    # ---------------------------------------------------------
    # 6. SCORING
    # ---------------------------------------------------------
    print_green("\n" + "=" * 119)
    print_green("KROK 4. PUNKTACJA (SCORING) OBIEKTÓW")
    print_green("=" * 119)
    scored_rows = []
    
    for i, row in df.iterrows():
        # Max Alt (do JSON)
        max_alt = compute_max_altitude(loc["lat"], row["dec"])
        
        # Scoring
        score_res = calculate_scorecard(row, params)
        final_score = score_res["final_score"]
        
        # Pobierz godziny per miesiąc dla JSON
        hours = hours_per_month[i]
        
        r_dict = row.to_dict()
        r_dict.update({
            "final_score": final_score,
            "base_score": final_score,
            "score_breakdown": score_res["breakdown"],
            "max_alt": max_alt,
            "hours_per_month": list(np.round(hours, 2))  # Zapisz widoczność do pliku wyjściowego
        })
        scored_rows.append(r_dict)

    df_scored = pd.DataFrame(scored_rows)
    if df_scored.empty:
        print("Brak obiektów po scoringu.")
        return

    # ---------------------------------------------------------
    # 7. PRE-FILTERING (Tylko najgorsze śmieci, żeby nie blokować tła!)
    # ---------------------------------------------------------
    df_scored = apply_adaptive_soft_cut(
        df_scored, 
        base_min_score=10.0,     # Niska podłoga! (puszczamy plankton do klastrowania)
        min_len_for_cut=5000,   # Celowo za wysoko, żeby wyłączyć percentyl na tym etapie
        keep_top_percent=100.0
    )
    if df_scored.empty:
        print("\n[INFO] Brak obiektów po wstępnej filtracji adaptacyjnej.")
        return

    # ---------------------------------------------------------
    # 8. CLUSTERING (Atlas Optimization)
    # ---------------------------------------------------------
    print_green("\n" + "=" * 119)
    print_green("KROK 5. ŁĄCZENIE OBIEKTÓW W KADRY (CLUSTERING)")
    print_green("=" * 119)
    df_clustered = optimize_atlas_pages(df_scored, params)

    # ---------------------------------------------------------
    # 8b. FINAL ADAPTIVE SOFT CUT (Ostre cięcie gotowych kadrów do Atlasu)
    # ---------------------------------------------------------
    print_step("Ostateczna selekcja punktowa zoptymalizowanych kadrów...")
    df_clustered = apply_adaptive_soft_cut(
        df_clustered, 
        base_min_score=50.0,       # Twarde odcięcie nudnych kadrów
        min_len_for_cut=150,          # Jeśli kandydatów do  jest >150...
        keep_top_percent=70.0     # ...zostaw tylko najlepsze x%.
    )

    # ---------------------------------------------------------
    # 9. HARD FILTER (Size)
    # ---------------------------------------------------------
    print_step(f"Hard Filter Size >= {params['min_size_arcmin']:.1f}'")
    df_final = df_clustered[df_clustered["size"] >= params["min_size_arcmin"]].copy()
    
    # Sortowanie końcowe (najlepsze na górę)
    df_final = df_final.sort_values(by=["final_score"], ascending=[False]).reset_index(drop=True)
    print(f"    ✓ Pozostało {fmt(len(df_final))} finalnych kadrów.")

    # ---------------------------------------------------------
    # 10. ZAPIS DO JSON I WYŚWIETLENIE
    # ---------------------------------------------------------
    print_green("\n" + "=" * 119)
    print_green("KROK 6. ZAPIS DANYCH DO VIS_DATA.JSON")
    print_green("=" * 119)
    MAX_PRINT_ROWS = 200
    
    print("\n" + "=" * 119)
    print(f"{'Nr':>3} {'ID':<12} {'Extra':<25} {'Common name':<25} {'Typ':<8} "
          f"{'RA':>7} {'Dec':>7} {'Mag':>7} {'Size':>8} {'Score':>6}")
    print("-" * 119)

    json_objects = []
    table_truncated = False
    
    for i, (_, row) in enumerate(df_final.iterrows()):
        # WYCZYSZCZENIE NAZW ("nan" -> "")
        c_id = _clean_str(row.get("id", ""))
        c_names = _clean_str(row.get("common_names", ""))
        c_extra = _clean_str(row.get("extra_info", ""))
        c_type = _clean_str(row.get("type", "OTHER"))

        # Budowa obiektu do JSON (zapisujemy WSZYSTKIE kadry)
        obj_dict = {
            "id": _clean_str(row.get("id", "")),
            "common_names": _clean_str(row.get("common_names", "")),
            "extra_info": _clean_str(row.get("extra_info", "")),
            "type": _clean_str(row.get("type", "OTHER")),
            "catalog": str(row.get("catalog", "unknown")),
            "ra": float(row["ra"]),
            "dec": float(row["dec"]),
            "size": float(row["size"]),
            "mag": float(row["mag"]) if pd.notna(row["mag"]) else None,
            "max_alt": float(row["max_alt"]),
            "score": float(row["final_score"]),
            "base_score": float(row["base_score"]),
            "score_breakdown": row["score_breakdown"],
            "hours_per_month": row["hours_per_month"]
        }
        
        # Zapisz konkretne numery z katalogów (0 oznacza brak przynależności)
        obj_dict["messier_nr"] = _as_int0(row.get("messier_nr"))
        obj_dict["caldwell_nr"] = _as_int0(row.get("caldwell_nr"))
        obj_dict["herschel_nr"] = _as_int0(row.get("herschel_nr"))
        
        json_objects.append(obj_dict)

        # Drukowanie (zabezpieczenie nałożone na sam print)
        if i >= MAX_PRINT_ROWS:
            if not table_truncated:
                print("=" * 119)
                msg = f"... i {len(df_final) - MAX_PRINT_ROWS} więcej kadrów (limit wyświetlania)."
                print(f"{msg:>119}")
                print("=" * 119)
                table_truncated = True
            continue

        # Weryfikacja jakości (● = pełne/zmierzone (imputed==False), ○ = niepełne/zgadywane (imputed==True))
        mag_imputed = _as_bool(row.get("mag_imputed", False))
        size_imputed = _as_bool(row.get("size_imputed", False))
        
        mag_mark = '○' if mag_imputed else '●'
        size_mark = '○' if size_imputed else '●'

        # Ochrona przed ewentualnymi NaN na wydruku (teoretycznie nie powinno ich być w Twoim CSV)
        mag_str = f"{mag_mark}{row['mag']:.1f}" if pd.notna(row.get("mag")) else "   None"
        size_str = f"{size_mark}{row['size']:.1f}'" if pd.notna(row.get("size")) else "   None"

        ra_str = f"{row['ra']:.1f}"
        dec_str = f"{row['dec']:.1f}"

        # Wydruk z zachowaniem ucięcia nazw (np. Extra na 25 znaków), aby utrzymać perfekcyjną siatkę!
        print(f"{i+1:3d}. {c_id[:12]:<12} {c_extra[:25]:<25} {c_names[:25]:<25} {c_type[:8]:<8} "
              f"{ra_str:>7} {dec_str:>7} {mag_str:>7} {size_str:>8} {row['final_score']:6.1f}")


    if not table_truncated:
        print("=" * 119)

    # Zapis
    out_data = {
        "location": loc,
        "year": year,
        "parameters": params,
        "objects": json_objects
    }
    
    # Hash do cache'owania w kroku 3
    add_parameters_hash_to_output(out_data, params)

    save_path = PATHS.vis_data
    with open(save_path, "w", encoding="utf-8") as f:
        json.dump(out_data, f, indent=2, ensure_ascii=False)
        
    print_step(f"Zapisano {fmt(len(json_objects))} zoptymalizowanych kadrów do {save_path}\n")

if __name__ == "__main__":
    main()
