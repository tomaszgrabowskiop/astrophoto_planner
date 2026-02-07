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

# =====================================================
# CONFIGURATION
# =====================================================

class Config:
    OUTPUT_FILE = Path("vis_data.json")
    CATALOG_PATH = Path("katalog_astro_full.csv")

    # Domyślny setup: RedCat 61 + ASI2600MC Pro
    DEFAULT_FOCAL_LENGTH_MM = 300.0
    DEFAULT_SENSOR_WIDTH_MM = 23.5
    DEFAULT_SENSOR_HEIGHT_MM = 15.7

    # Dodatkowe domyślne parametry kamery (ASI2600MC Pro)
    DEFAULT_SENSOR_PITCH_UM = 3.76
    DEFAULT_SENSOR_ROWS = 4176
    DEFAULT_SENSOR_COLS = 6248
    
CAT_ORDER = ["ngc", "ic", "sh2", "rcw", "lbn", "ced", "pgc", "barn", "ldn"]

# =====================================================
# ASTRO MATH
# =====================================================

class AstroMath:
    """Fast astronomical calculations"""

    @staticmethod
    def get_julian_date(dt: datetime) -> float:
        if dt.tzinfo is not None:
            dt = dt.astimezone(timezone.utc)

        Y, M, D = dt.year, dt.month, dt.day
        H, MN, S = dt.hour, dt.minute, dt.second

        if M <= 2:
            Y -= 1
            M += 12

        A = Y // 100
        B = 2 - A + (A // 4)

        return (
            int(365.25 * (Y + 4716))
            + int(30.6001 * (M + 1))
            + D
            + B
            - 1524.5
            + (H + MN / 60 + S / 3600) / 24
        )

    @staticmethod
    def get_gmst(jd: float) -> float:
        T = (jd - 2451545.0) / 36525.0
        gmst = (
            280.46061837
            + 360.98564736629 * (jd - 2451545.0)
            + 0.000387933 * T**2
            - T**3 / 38710000.0
        )
        return gmst % 360.0

    @staticmethod
    def get_sun_pos(jd: float) -> Tuple[float, float]:
        T = (jd - 2451545.0) / 36525.0

        L0 = 280.46646 + 36000.76983 * T + 0.0003032 * T**2
        M = 357.52911 + 35999.05029 * T - 0.0001537 * T**2

        C = (
            (1.914602 - 0.004817 * T - 0.000014 * T**2) * np.sin(np.radians(M))
            + (0.019993 - 0.000101 * T) * np.sin(np.radians(2 * M))
            + 0.000289 * np.sin(np.radians(3 * M))
        )

        true_long = L0 + C
        omega = 125.04 - 1934.136 * T
        lambda_sun = true_long - 0.00569 - 0.00478 * np.sin(np.radians(omega))

        epsilon_mean = 23 + 26 / 60 + 21.448 / 3600 - (46.8150 * T) / 3600
        epsilon = epsilon_mean + 0.00256 * np.cos(np.radians(omega))

        rad_lambda = np.radians(lambda_sun)
        rad_eps = np.radians(epsilon)

        alpha = np.arctan2(
            np.cos(rad_eps) * np.sin(rad_lambda),
            np.cos(rad_lambda),
        )
        delta = np.arcsin(np.sin(rad_eps) * np.sin(rad_lambda))

        return np.degrees(alpha) % 360, np.degrees(delta)

    @staticmethod
    def alt_az(ha_deg: float, dec_deg: float, lat_deg: float) -> float:
        ha = np.radians(ha_deg)
        dec = np.radians(dec_deg)
        lat = np.radians(lat_deg)

        sin_alt = (
            np.sin(dec) * np.sin(lat)
            + np.cos(dec) * np.cos(lat) * np.cos(ha)
        )

        return np.degrees(np.arcsin(np.clip(sin_alt, -1.0, 1.0)))

# =====================================================
# VECTORIZED ASTRO
# =====================================================

class VectorAstro:
    """Zwektoryzowane obliczenia astronomiczne dla Numpy"""

    @staticmethod
    def get_gmst_vector(jd_array: np.ndarray) -> np.ndarray:
        """Wektorowe obliczanie GMST dla tablicy Julian Dates"""
        T = (jd_array - 2451545.0) / 36525.0
        gmst = (
            280.46061837
            + 360.98564736629 * (jd_array - 2451545.0)
            + 0.000387933 * T**2
            - T**3 / 38710000.0
        )
        return gmst % 360.0

    @staticmethod
    def alt_az_vector(ha_rad: np.ndarray, dec_rad: np.ndarray, lat_rad: float) -> np.ndarray:
        """Wektorowe obliczanie wysokości."""
        sin_alt = (
            np.sin(dec_rad) * np.sin(lat_rad)
            + np.cos(dec_rad) * np.cos(lat_rad) * np.cos(ha_rad)
        )
        return np.degrees(np.arcsin(np.clip(sin_alt, -1.0, 1.0)))

# =====================================================
# SMART NIGHT GRID
# =====================================================

def get_full_year_smart_grid(
    year: int,
    observer_loc: Dict[str, Any],
    sun_limit_deg: float,
    step_minutes: int = 15,
) -> Tuple[np.ndarray, np.ndarray, np.ndarray]:
    """
    Generuje siatkę czasu (JD) tylko dla nocy, używając astroplan
    do wyznaczenia dokładnych momentów zmierzchu i świtu.
    Zwraca: (jd_array, month_array, lst_array).
    """
    print(f"Generowanie SMART GRID (tylko noce) dla roku {year}, limit {sun_limit_deg}°...")

    obs = Observer(
        latitude=observer_loc["lat"] * u.deg,
        longitude=observer_loc["lon"] * u.deg,
        elevation=0 * u.m,
    )

    t_scan = Time(f"{year}-01-01 12:00:00", scale="utc")
    t_end = Time(f"{year+1}-01-01 12:00:00", scale="utc")

    step_day = step_minutes / (60.0 * 24.0)

    all_jds: List[np.ndarray] = []
    day_count = 0

    while t_scan < t_end:
        try:
            t_set = obs.sun_set_time(t_scan, horizon=sun_limit_deg * u.deg, which="next")
            t_rise = obs.sun_rise_time(t_set, horizon=sun_limit_deg * u.deg, which="next")

            if t_rise > t_set:
                night_grid = np.arange(t_set.jd, t_rise.jd, step_day)
                all_jds.append(night_grid)
        except Exception:
            pass

        t_scan = t_scan + 1.0 * u.day
        day_count += 1

    if not all_jds:
        print("! Ostrzeżenie: Nie znaleziono żadnych nocy (błąd parametrów?).")
        return np.array([]), np.array([]), np.array([])

    jd_array = np.concatenate(all_jds)

    print(
        f"✓ Wygenerowano {len(jd_array)} punktów pomiarowych "
        f"(efektywnie {len(jd_array)*step_minutes/60:.0f}h obserwacji)."
    )

    t_all = Time(jd_array, format="jd")
    dates = t_all.to_datetime()
    months_array = np.array([d.month - 1 for d in dates])

    lon_deg = observer_loc["lon"]
    gmst = VectorAstro.get_gmst_vector(jd_array)
    lst_array = (gmst + lon_deg) % 360.0

    return jd_array, months_array, lst_array

def compute_all_visibilities_matrix(
    df_objects: pd.DataFrame,
    observer_loc: Dict[str, Any],
    jd_array: np.ndarray,
    month_array: np.ndarray,
    lst_array: np.ndarray,
    min_alt: float,
    step_minutes: int,
) -> np.ndarray:
    """
    Oblicza widoczność (Alt > min_alt) dla przygotowanej siatki nocnej,
    zwraca godziny/miesiąc dla każdego obiektu.
    """
    if df_objects.empty or len(jd_array) == 0:
        return np.zeros((len(df_objects), 12))

    print(
        f"Obliczanie macierzowe (Core) dla {len(df_objects)} obiektów "
        f"i {len(jd_array)} punktów czasu..."
    )

    lat_rad = np.radians(observer_loc["lat"])

    ra_objs = df_objects["ra"].values
    dec_objs = df_objects["dec"].values

    ra_rad_row = np.radians(ra_objs)[np.newaxis, :]
    dec_rad_row = np.radians(dec_objs)[np.newaxis, :]

    lst_rad_col = np.radians(lst_array)[:, np.newaxis]
    ha_rad_matrix = lst_rad_col - ra_rad_row

    alt_matrix = VectorAstro.alt_az_vector(ha_rad_matrix, dec_rad_row, lat_rad)

    visible_mask = (alt_matrix >= min_alt).astype(float)

    step_hours = step_minutes / 60.0
    weighted_mask = visible_mask * step_hours

    hours_per_month = np.zeros((len(df_objects), 12))
    unique_months = np.unique(month_array)

    for m in unique_months:
        time_indices = (month_array == m)
        if not np.any(time_indices):
            continue
        hours_per_month[:, m] = np.sum(weighted_mask[time_indices], axis=0)

    return hours_per_month

def compute_max_altitude(lat: float, dec: float) -> float:
    lat_rad = np.radians(lat)
    dec_rad = np.radians(dec)
    sin_alt_max = np.sin(lat_rad) * np.sin(dec_rad) + np.cos(lat_rad) * np.cos(dec_rad)
    return np.degrees(np.arcsin(np.clip(sin_alt_max, -1.0, 1.0)))

def compute_rarity_score(max_altitude: float, months_observable: int) -> float:
    """
    Oblicz "rarity_score" w zakresie 0–1 na podstawie:
    - maksymalnej wysokości w roku,
    - liczby miesięcy, w których obiekt jest sensownie obserwowalny.
    """

    # Składnik od wysokości maksymalnej
    if max_altitude >= 70:
        alt_rarity = 0.0
    elif max_altitude >= 50:
        alt_rarity = 0.3
    elif max_altitude >= 35:
        alt_rarity = 0.7
    else:
        alt_rarity = 1.0

    # Składnik od długości sezonu
    if months_observable >= 6:
        month_rarity = 0.0
    elif months_observable >= 4:
        month_rarity = 0.3
    elif months_observable >= 2:
        month_rarity = 0.7
    else:
        month_rarity = 1.0

    return (alt_rarity + month_rarity) / 2.0

# =====================================================
# IMPUTATION
# =====================================================

def impute_magnitude(row):
    """Estymuj mag jeśli brakuje, na podstawie type_code i catalog"""
    mag = row.get("mag", np.nan)
    if pd.notna(mag):
        return float(mag)

    obj_type_code = str(row.get("type_code", "OTHER")).upper().strip()
    catalog = str(row.get("catalog", "ngc")).lower().strip()

    type_estimates = {
        "HII": 8.5,
        "SNR": 8.0,
        "DN": 28.0,
        "OC": 7.5,
        "GC": 8.5,
        "GX": 11.0,
        "G": 11.0,
        "S": 11.5,
        "SB": 11.5,
        "E": 11.0,
        "PN": 9.5,
        "NB": 12.0,
        "CL+N": 8.0,
        "OCL": 7.5,
        "GCL": 8.5,
    }

    catalog_estimates = {
        "sh2": 8.5,
        "rcw": 8.5,
        "barn": 28.0,
        "ldn": 28.0,
        "lbn": 12.0,
        "pgc": 11.0,
        "ced": 12.0,
        "ngc": 10.0,
    }

    if obj_type_code in type_estimates:
        return type_estimates[obj_type_code]

    if catalog in catalog_estimates:
        return catalog_estimates[catalog]

    return 11.0  # fallback

def impute_size(row):
    """Estymuj size jeśli brakuje, zwraca (size, source_flag)"""
    size = row.get("size", np.nan)
    if pd.notna(size) and float(size) > 0:
        return float(size), "measured"

    obj_type_code = str(row.get("type_code", "OTHER")).upper().strip()
    catalog = str(row.get("catalog", "ngc")).lower().strip()

    type_estimates = {
        "HII": 40.0,
        "SNR": 35.0,
        "DN": 25.0,
        "OC": 20.0,
        "GC": 12.0,
        "GX": 25.0,
        "G": 25.0,
        "S": 25.0,
        "SB": 25.0,
        "E": 20.0,
        "PN": 5.0,
        "NB": 18.0,
        "CL+N": 30.0,
        "OCL": 20.0,
        "GCL": 12.0,
    }

    catalog_estimates = {
        "sh2": 40.0,
        "rcw": 30.0,
        "barn": 20.0,
        "ldn": 15.0,
        "lbn": 15.0,
        "pgc": 20.0,
        "ced": 10.0,
        "ngc": 20.0,
    }

    if obj_type_code in type_estimates:
        return type_estimates[obj_type_code], "estimated"

    if catalog in catalog_estimates:
        return catalog_estimates[catalog], "estimated"

    return 15.0, "estimated"

def load_and_impute_catalog(df_catalog: pd.DataFrame) -> pd.DataFrame:
    """Uzupełnij brakujące dane w katalogu"""
    df = df_catalog.copy()

    df["has_mag_measured"] = df["mag"].notna()
    df["has_size_measured"] = (
        df["size"].notna()
        & (pd.to_numeric(df["size"], errors="coerce") > 0)
    )

    df["mag"] = df.apply(impute_magnitude, axis=1)
    size_imputed = df.apply(impute_size, axis=1)
    df["size"] = size_imputed.apply(lambda x: x[0])
    df["size_source"] = size_imputed.apply(lambda x: x[1])

    df["size"] = pd.to_numeric(df["size"], errors="coerce").fillna(15.0)
    df["mag"] = pd.to_numeric(df["mag"], errors="coerce").fillna(11.0)

    return df

# =====================================================
# SURFACE BRIGHTNESS + SCORING (pełny z v2)
# =====================================================

def calculate_surface_brightness(mag, size_arcmin):
    """Oblicz SB (mag/arcsec²) - kluczowe dla miasta"""
    if pd.isna(mag) or pd.isna(size_arcmin) or size_arcmin <= 0:
        return 99.0

    radius_sec = (size_arcmin * 60.0) / 2.0
    area_arcsec2 = np.pi * (radius_sec**2)
    sb = mag + 2.5 * np.log10(area_arcsec2)
    return sb

def get_data_quality_label(row):
    """Etykieta jakości danych"""
    has_mag = row.get("has_mag_measured", False)
    has_size = row.get("has_size_measured", False)

    if has_mag and has_size:
        return "✓ Complete"
    elif has_mag or has_size:
        return "◐ Partial"
    else:
        return "○ Estimated"

def city_optimized_score(row, user_params):
    """
    Zaawansowany scoring dla miasta + filtry narrowband.
    """
    has_narrowband = user_params.get("has_narrowband", False)

    # NOWE: jeśli nie ma filtrów, stosuj prosty, „niemiejskI” scoring
    if not has_narrowband:
        mag = float(row["mag"]) if pd.notna(row["mag"]) else 99.0
        size = float(row["size"]) if pd.notna(row["size"]) else 0.0

        # bardzo prosty score: jasność + rozmiar, bez Bortle/SB itd.
        score = 0.0
        if mag < 7:
            score += 20
        elif mag < 10:
            score += 10

        if size > 60:
            score += 25
        elif size > 30:
            score += 15
        elif size > 10:
            score += 5

        final_score = max(0.0, score)
        return {
            "base_score": score,
            "quality_factor": 1.0,
            "final_score": final_score,
            "data_quality": get_data_quality_label(row),
            "size_source": row.get("size_source", "unknown"),
            "mag_source": "measured" if row.get("has_mag_measured") else "estimated",
        }

    obj_type_code = str(row.get("type_code", "")).upper().strip()
    mag = float(row["mag"]) if pd.notna(row["mag"]) else 99.0
    size = float(row["size"]) if pd.notna(row["size"]) else 0.0
    catalog = str(row.get("catalog", "")).lower()
    extra = str(row.get("extra_info", ""))

    bortle_min, bortle_max = user_params.get("bortle_range", (6, 7))
    has_narrowband = user_params.get("has_narrowband", False)
    prefer_famous = user_params.get("prefer_famous", True)

    score = 0.0

    is_narrowband_target = (
        obj_type_code in ["HII", "SNR", "PN"]
        or catalog in ["sh2", "rcw"]
    )

    is_cluster = obj_type_code in ["OC", "GC", "CL+N", "OCL", "GCL"]
    is_galaxy = obj_type_code in ["GX", "G", "S", "SB", "E"]
    is_dark = (obj_type_code in ["DN"]) or (catalog in ["barn", "ldn"])

    if is_narrowband_target:
        if has_narrowband:
            if bortle_max >= 7:
                score += 80
            else:
                score += 50
        else:
            if bortle_max >= 7:
                score += 20
            else:
                score += 40

    if is_cluster:
        # Dla gromad SB jest mylące (gwiazdy punktowe, ciemne tło między nimi).
        # Dlatego NIE filtrujemy po SB, tylko patrzymy na zintegrowaną jasność.
        if mag < 7:
            score += 20
        elif mag < 10:
            score += 10
    else:
        # Dla mgławic/galaktyk SB jest krytyczne, zwłaszcza w mieście.
        if size > 0:
            sb = calculate_surface_brightness(mag, size)

            if bortle_max <= 5:
                # Ciemniejsze niebo (Bortle 1–5): można tolerować trochę słabsze SB
                if sb < 22:
                    score += 30
                elif sb < 24:
                    score += 15
            else:
                # Miasto / silne LP (Bortle 6+): trzeba mocno premiować jasne SB
                if sb < 20:
                    score += 40
                elif sb < 22:
                    score += 20
                # Jeśli obiekt ma bardzo słabe SB i NIE jest celem NB – odrzucamy
                elif sb > 23 and not is_narrowband_target:
                    return {
                        "base_score": 0.0,
                        "quality_factor": 1.0,
                        "final_score": 0.0,
                        "data_quality": get_data_quality_label(row),
                        "size_source": row.get("size_source", "unknown"),
                        "mag_source": (
                            "measured" if row.get("has_mag_measured") else "estimated"
                        ),
                    }

    if is_galaxy:
        if bortle_max <= 5:
            score += 20
        elif bortle_max <= 7:
            score += 5
        else:
            score -= 20

    if is_dark:
        score -= 100

    if size > 0:
        sb = calculate_surface_brightness(mag, size)

        if bortle_max <= 5:
            if sb < 22:
                score += 30
            elif sb < 24:
                score += 15
        else:
            if sb < 20:
                score += 40
            elif sb < 22:
                score += 20
            elif sb > 23 and not is_narrowband_target:
                return {
                    "base_score": 0.0,
                    "quality_factor": 1.0,
                    "final_score": 0.0,
                    "data_quality": get_data_quality_label(row),
                    "size_source": row.get("size_source", "unknown"),
                    "mag_source": (
                        "measured" if row.get("has_mag_measured") else "estimated"
                    ),
                }

    if prefer_famous:
        if "M" in extra:
            score += 100
        elif "C" in extra:
            score += 60
        elif "H" in extra:
            score += 30
    else:
        if "M" in extra:
            score += 40
        elif "C" in extra:
            score += 20

    if mag < 7:
        score += 20
    elif mag < 10:
        score += 10

    if size > 60:
        score += 25
    elif size > 30:
        score += 15
    elif size > 10:
        score += 5

    quality = 1.0
    if not row.get("has_mag_measured", True):
        quality *= 0.90
    if not row.get("has_size_measured", True):
        quality *= 0.85

    final_score = max(0.0, score * quality)

    return {
        "base_score": score,
        "quality_factor": quality,
        "final_score": final_score,
        "data_quality": get_data_quality_label(row),
        "size_source": row.get("size_source", "unknown"),
        "mag_source": (
            "measured" if row.get("has_mag_measured") else "estimated"
        ),
    }

# =====================================================
# FOV CALCULATOR
# =====================================================

class FOVCalculator:
    @staticmethod
    def calculate_fov_deg(
        focal_length_mm: float,
        sensor_width_mm: float,
        sensor_height_mm: float,
    ) -> Tuple[float, float]:
        fov_w = 57.3 * sensor_width_mm / focal_length_mm
        fov_h = 57.3 * sensor_height_mm / focal_length_mm
        return fov_w, fov_h

    @staticmethod
    def min_object_size_arcmin(
        fov_w_deg: float, fov_h_deg: float, percent_fov: float
    ) -> float:
        min_fov_deg = min(fov_w_deg, fov_h_deg)
        min_size_deg = min_fov_deg * percent_fov / 100.0
        return min_size_deg * 60.0

# =====================================================
# YEAR  
# =====================================================
class YearManager:
    def select_interactive(self) -> int:
        current_year = date.today().year
        print("=" * 141)
        print("ROK OBSERWACJI")
        print("=" * 141)
        
        year_str = input(f"[Enter={current_year} lub wpisz wybrany pomiędzy 2000 a 2100]: ").strip()
        if year_str:
            try:
                year = int(year_str)
                if year < 2000 or year > 2100:
                    raise ValueError()
            except:
                print(f"Nieprawidłowy, używam {current_year}")
                year = current_year
        else:
            year = current_year
        
        print(f"✓ Rok: {year}")
        return year

# =====================================================
# LOCATION MANAGER
# =====================================================

class LocationManager:
    def select_interactive(self) -> Dict[str, Any]:
        print("=" * 141)
        print("KROK 1: WYBÓR LOKALIZACJI")
        print("=" * 141)
        print("1. Poznań, Polska")
        print("2. Kraków, Polska")
        print("3. Wpisz lokalizację\n")

        choice = input(r"[1/2/3, domyślnie 1]: ").strip() or "1"

        if choice == "1":
            loc = {
                "lat": 52.4095,
                "lon": 16.9319,
                "tz": "Europe/Warsaw",
                "name": "Poznań, Polska",
            }
            print(
                f"\n✓ Wybrano {loc['name']} ({loc['lat']}°N, {loc['lon']}°E)"
            )
            return loc

        if choice == "2":
            loc = {
                "lat": 50.0647,
                "lon": 19.9450,
                "tz": "Europe/Warsaw",
                "name": "Kraków, Polska",
            }
            print(
                f"\n✓ Wybrano {loc['name']} ({loc['lat']}°N, {loc['lon']}°E)"
            )
            return loc

        city = input("Podaj lokalizację (pisz wielką literą, działa lepiej, np. Toruń, Polska): ").strip()
        if not city:
            print("Brak nazwy, używam domyślnej: Poznań.")
            return {
                "lat": 52.4095,
                "lon": 16.9319,
                "tz": "Europe/Warsaw",
                "name": "Poznań, Poland",
            }

        try:
            print(f"Szukanie współrzędnych dla '{city}'...")
            loc_astro = EarthLocation.of_address(city)
            lat = loc_astro.lat.to(u.deg).value
            lon = loc_astro.lon.to(u.deg).value
            loc = {
                "lat": lat,
                "lon": lon,
                "tz": "Europe/Warsaw",
                "name": city,
            }
            print(
                f"\n✓ Znaleziono {loc['name']} ({loc['lat']:.4f}°N, {loc['lon']:.4f}°E)"
            )
            return loc
        except Exception as e:
            print(f"Błąd geokodowania: {e}, używam Poznań.")
            return {
                "lat": 52.4095,
                "lon": 16.9319,
                "tz": "Europe/Warsaw",
                "name": "Poznań, Polska",
            }

# =====================================================
# DSO CATALOG
# =====================================================
def choose_best_catalog_id(main_id: str, extra_info: str) -> str:
    """
    Wybierz najlepszy identyfikator katalogowy do pola 'id'
    na podstawie priorytetu katalogów (NGC/IC > Sh2 > ... > LDN).
    Nie zmienia żadnych innych pól, tylko zwraca nowy string ID.
    """
    def cat_rank(obj_id: str) -> tuple[int, str]:
        s = obj_id.lower().strip()
        if s.startswith("ngc"):
            cat = "ngc"
        elif s.startswith("ic"):
            cat = "ic"
        elif s.startswith("sh2-"):
            cat = "sh2"
        elif s.startswith("rcw"):
            cat = "rcw"
        elif s.startswith("lbn"):
            cat = "lbn"
        elif s.startswith("ced"):
            cat = "ced"
        elif s.startswith("pgc"):
            cat = "pgc"
        elif s.startswith("b"):
            cat = "barn"
        elif s.startswith("ldn"):
            cat = "ldn"
        else:
            cat = "zzz"

        try:
            rank = CAT_ORDER.index(cat)
        except ValueError:
            rank = len(CAT_ORDER)
        return (rank, s)

    # Rozbij extra_info na listę ID
    extra_ids = [
        x.strip()
        for x in str(extra_info).split(",")
        if x and x.strip()
    ]

    # Kandydaci: aktualne id + wszystkie z extra_info
    candidates = [str(main_id).strip()] + extra_ids
    # wybierz najlepszy wg katalogu, potem alfabetycznie
    best = sorted(candidates, key=cat_rank)[0]
    return best
    
class DSOCatalog:
    def __init__(self):
        self.df: pd.DataFrame | None = None

    def load(self):
        if not Config.CATALOG_PATH.exists():
            raise FileNotFoundError(f"Brak pliku katalogu: {Config.CATALOG_PATH}")

        print(f"Ładowanie bazy obiektów: {Config.CATALOG_PATH}...")
        self.df = pd.read_csv(Config.CATALOG_PATH, low_memory=False)
        if "common_names" in self.df.columns:
            self.df["common_names"] = (
                self.df["common_names"]
                .fillna("")           # NaN → ""
                .astype(str)
                .replace("nan", "")   # literalne "nan" → ""
                .str.strip()
            )

        if "type" in self.df.columns:
            self.df.rename(columns={"type": "type_code"}, inplace=True)

        self.df["mag"] = pd.to_numeric(self.df["mag"], errors="coerce")
        self.df["size"] = pd.to_numeric(self.df["size"], errors="coerce")

        if "base_score" not in self.df.columns:
            print("! UWAGA: Brak 'base_score' w CSV. Skrypt będzie liczył dynamicznie.")
            self.df["base_score"] = 1.0
        else:
            self.df["base_score"] = pd.to_numeric(
                self.df["base_score"], errors="coerce"
            )

        self.df = self.df.dropna(subset=["ra", "dec"])
        print(f"✓ Załadowano {len(self.df)} obiektów z pełnymi danymi.")

    def get_display_info(self, row):
        type_map = {
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

        code = str(row.get("type_code", "")).strip()
        obj_type = type_map.get(code, code)

        main_name = str(row["id"])
        extra = (
            str(row.get("extra_info", ""))
            if pd.notna(row.get("extra_info"))
            else ""
        )

        if extra:
            display_name = f"{main_name} ({extra})"
        else:
            display_name = main_name

        return display_name, obj_type

# =====================================================
# GEOMETRY FILTER
# =====================================================

def filter_candidates_geometric(
    df: pd.DataFrame,
    lat: float,
    min_alt: float,
    min_hours: float,
) -> pd.DataFrame:
    print("[INFO] Optymalizacja geometryczna (filtrowanie matematyczne)...")

    lat_rad = np.radians(lat)
    dec_rad_all = np.radians(df["dec"].values)

    sin_alt_max = np.sin(lat_rad) * np.sin(dec_rad_all) + np.cos(lat_rad) * np.cos(
        dec_rad_all
    )
    alt_max_deg = np.degrees(np.arcsin(np.clip(sin_alt_max, -1.0, 1.0)))

    mask_alt = alt_max_deg >= min_alt
    df_alt = df[mask_alt].copy()
    print(f" - Po filtrze min. wysokości ({min_alt}°): {len(df_alt)}/{len(df)}")

    if df_alt.empty:
        return df_alt

    dec_rad = np.radians(df_alt["dec"].values)
    alt_rad = np.radians(min_alt)

    numerator = np.sin(alt_rad) - np.sin(lat_rad) * np.sin(dec_rad)
    denominator = np.cos(lat_rad) * np.cos(dec_rad)

    cos_h = np.divide(
        numerator,
        denominator,
        out=np.zeros_like(numerator),
        where=(denominator != 0.0),
    )

    cos_h = np.clip(cos_h, -1.0, 1.0)

    durations = np.zeros_like(cos_h)

    mask_circumpolar = cos_h <= -1
    durations[mask_circumpolar] = 24.0

    mask_normal = (cos_h > -1) & (cos_h < 1)
    durations[mask_normal] = (
        2.0 * np.degrees(np.arccos(cos_h[mask_normal])) / 15.041
    )

    mask_dur = durations >= min_hours
    df_final = df_alt[mask_dur].copy()

    print(
        f" - Po filtrze czasu nad horyzontem ({min_hours}h): "
        f"{len(df_final)}/{len(df_alt)}"
    )

    return df_final

# =====================================================
# USER PREFERENCES
# =====================================================

def explain_and_get_prefs():
    print("=" * 141)
    print("KROK 2: PARAMETRY WIDOCZNOŚCI I FOV")
    print("=" * 141)

    try:
        min_alt = float(
            input("A. Min. wysokość obiektu nad horyzontem (°, domyślnie 25): ")
            or 25.0
        )
    except Exception:
        min_alt = 25.0

    print("\nB. CIEMNOŚĆ NIEBA (kąt Słońca pod horyzontem):")
    print(" -6° = zmierzch cywilny (6)")
    print(" -12° = zmierzch żeglarski (12)")
    print(" -18° = zmierzch astronomiczny (18)\n")

    choice = input("[6/12/18, domyślnie 12]: ").strip() or "12"
    if choice == "6":
        sun_limit = -6.0
    elif choice == "12":
        sun_limit = -12.0
    else:
        sun_limit = -18.0

    try:
        min_hours = float(
            input("\nC. Min. godzin w nocy (domyślnie 3): ") or 3.0
        )
    except Exception:
        min_hours = 3.0

    print("\nD. FOV (sprzęt):")
    use_default = input("Użyć RedCat61 + ASI2600MC Pro? [y/n, domyślnie y]: ").strip() or "y"
    use_default = use_default.lower()

    if use_default == "y":
        focal = Config.DEFAULT_FOCAL_LENGTH_MM
        sw = Config.DEFAULT_SENSOR_WIDTH_MM
        sh = Config.DEFAULT_SENSOR_HEIGHT_MM
        pitch = Config.DEFAULT_SENSOR_PITCH_UM
        rows = Config.DEFAULT_SENSOR_ROWS
        cols = Config.DEFAULT_SENSOR_COLS
        print(f"✓ Ogniskowa: {focal} mm, sensor: {sw} x {sh} mm")
        print(f" (pitch: {pitch} µm, {cols} x {rows} px)")
    else:
        try:
            focal = float(input(" Ogniskowa (mm): ").strip())
        except Exception:
            focal = Config.DEFAULT_FOCAL_LENGTH_MM

        try:
            sw = float(input(" Szerokość sensora (mm): ").strip())
            sh = float(input(" Wysokość sensora (mm): ").strip())
        except Exception:
            sw = Config.DEFAULT_SENSOR_WIDTH_MM
            sh = Config.DEFAULT_SENSOR_HEIGHT_MM

        # nowy blok – parametry dyskretne matrycy
        try:
            pitch = float(input(" Wielkość piksela (µm, domyślnie 3.76): ").strip() or 3.76)
        except Exception:
            pitch = 3.76

        try:
            rows = int(input(" Liczba wierszy sensora (domyślnie 4176): ").strip() or 4176)
            cols = int(input(" Liczba kolumn sensora (domyślnie 6248): ").strip() or 6248)
        except Exception:
            rows = 4176
            cols = 6248

    fov_w, fov_h = FOVCalculator.calculate_fov_deg(focal, sw, sh)
    print(f"✓ Szacowany FOV: {fov_w:.2f}° x {fov_h:.2f}°")

    try:
        percent_fov = float(
            input(
                "\nE. Minimalny rozmiar obiektu "
                "(% krótszego boku FOV, domyślnie 10): "
            )
            or 10.0
        )
    except Exception:
        percent_fov = 10.0

    min_size_arcmin = FOVCalculator.min_object_size_arcmin(
        fov_w, fov_h, percent_fov
    )

    print(
        f"✓ Min. rozmiar obiektu: {min_size_arcmin:.1f}' = "
        f"{percent_fov:.1f}% krótszego boku FOV"
    )

    print("\nF. BORTLE SCALE (Jasność nieba)")
    print(" 1 = Bortle 1-3 (ciemne niebo, wieś)")
    print(" 2 = Bortle 4-5 (przedmieścia, umiarkowane LP)")
    print(" 3 = Bortle 6-7 (miasto, silne LP)")
    print(" 4 = Bortle 8-9 (centrum miasta, ekstremalne LP)\n")

    bortle_choice = input("[1/2/3/4, domyślnie 4]: ").strip() or "4"
    bortle_map = {
        "1": (1, 3),
        "2": (4, 5),
        "3": (6, 7),
        "4": (8, 9),
    }
    bortle_range = bortle_map.get(bortle_choice, (6, 7))

    print("\nG. FILTRY")
    use_narrowband = (
        input("Masz filtry Ha/OIII/SHO? [y/n, domyślnie n]: ")
        .strip()
        .lower()
        or "n"
    )
    has_narrowband = use_narrowband == "y"

    print("\nH. PREFERENCJE KATALOGÓW")
    prefer_famous_in = (
        input(
            "Premiować obiekty Messier/Caldwell/Herschel? "
            "[y/n, domyślnie y]: "
        )
        .strip()
        .lower()
        or "y"
    )
    prefer_famous = prefer_famous_in == "y"

    return {
        "minalt": min_alt,
        "sunlimit": sun_limit,
        "minhours": min_hours,
        "minsizearcmin": min_size_arcmin,
        "bortle_range": bortle_range,
        "has_narrowband": has_narrowband,
        "prefer_famous": prefer_famous,
        "camera": {
            "lens_focal_length": focal,
            "sensor_pitch": pitch,
            "sensor_rows": rows,
            "sensor_cols": cols,
            "sensor_width": sw,
            "sensor_height": sh,
        },
    }

# =====================================================
# SAFETY CUT 
# =====================================================

def apply_safety_cut(
    canddf: pd.DataFrame,                           # Wejściowy DataFrame z obiektami po clusteringu/FOV i scoringu
    global_min_score: float = 3.0,              # Opcjonalny globalny minimalny final_score; wszystko poniżej jest od razu wyrzucane
    min_len_for_cut: int = 500,                   # Minimalna liczba obiektów, od której włączamy globalne cięcie po medianie
    main_percentile: float = 30.0,              # Percentyl (np. 50 = mediana) do globalnego obcięcia dołu rozkładu final_score
    pri_tail_percentile: float = 10.0,           # Dolny percentyl ogona w GRUPIE PRIORYTETOWEJ do ucięcia (np. 10% najgorszych)
    nonpri_tail_percentile: float = 50.0,             # Dolny percentyl ogona w GRUPIE NIEPRIORYTETOWEJ (mocniejsze cięcie, np. 30%)
    target_max_objects: int | None = None       # Opcjonalny twardy limit liczby obiektów; jeśli None, nie ograniczamy „na sztywno”
) -> pd.DataFrame:                                               # Zwraca przefiltrowany DataFrame po zastosowaniu wszystkich cięć
    df = canddf.copy()

    if "final_score" not in df.columns:
        print("[PLAN] Brak kolumny 'final_score' – safety-cut pomijam.")
        return df

    # Miękki globalny próg
    if global_min_score is not None:
        before = len(df)
        df = df[df["final_score"] >= global_min_score].reset_index(drop=True)
        print(f"[PLAN] Hard-cut global_min_score >= {global_min_score:.1f}: {before} -> {len(df)}")

    if len(df) == 0:
        return df

    # 1) Globalne odcięcie dołu rozkładu po medianie
    if len(df) >= min_len_for_cut:
        scores = df["final_score"].astype(float)
        p_cut = np.percentile(scores, main_percentile)
        before = len(df)
        df = df[df["final_score"] > p_cut].reset_index(drop=True)
        after = len(df)
        print(
            f"[PLAN] Global hard-cut 1: final_score > p{main_percentile:.0f} (≈{p_cut:.1f}). "
            f"Obiekty: {before} -> {after}."
        )
    else:
        print(f"[PLAN] Lista < {min_len_for_cut} obiektów – pomijam globalne cięcie percentylowe.")

    if len(df) == 0:
        return df

    # 2) Odtworzenie flag priorytetów, jeśli nie ma (opiera się na tych samych kolumnach co wcześniej)
    if "data_quality" in df.columns:
        dq = df["data_quality"].fillna("○ Estimated")
        df["dq_rank"] = dq.map({
            "✓ Complete": 0,
            "◐ Partial": 1,
            "○ Estimated": 2,
        }).fillna(2).astype(int)
    else:
        df["dq_rank"] = 1

    if "is_famous" not in df.columns:
        def _is_famous_row(row: pd.Series) -> bool:
            extra = str(row.get("extrainfo", "") or "")
            has_catalog = any(x in extra for x in ("M", "C", "H"))
            cname = row.get("common_names", None)
            has_common = (
                cname is not None
                and not pd.isna(cname)
                and str(cname).strip().lower() != "nan"
                and str(cname).strip() != ""
            )
            return has_catalog or has_common
        df["is_famous"] = df.apply(_is_famous_row, axis=1)

    if "is_fully_measured" not in df.columns:
        df["is_fully_measured"] = (
            (df.get("mag_source", "") == "measured") &
            (df.get("size_source", "") == "measured")
        )

    if "priority_group" not in df.columns:
        df["priority_group"] = df["is_famous"] | df["is_fully_measured"]

    # Dla czytelności: sort jak w głównej logice
    df = df.sort_values(
        by=["priority_group", "final_score"],
        ascending=[False, False],
    ).reset_index(drop=True)

    # 3) Rozbicie na dwie grupy
    df_pri = df[df["priority_group"]].copy()
    df_nonpri = df[~df["priority_group"]].copy()

    print(f"[PLAN] Podział na grupy: priorytetowa={len(df_pri)}, reszta={len(df_nonpri)}")

    def _cut_tail(group: pd.DataFrame, tail_percentile: float, label: str) -> pd.DataFrame:
        if len(group) == 0:
            return group
        scores = group["final_score"].astype(float)
        p_tail = np.percentile(scores, tail_percentile)
        before = len(group)
        group_cut = group[group["final_score"] > p_tail].reset_index(drop=True)
        after = len(group_cut)
        print(
            f"[PLAN] Tail-cut {label}: final_score > p{tail_percentile:.0f} (≈{p_tail:.1f}). "
            f"{before} -> {after}."
        )
        return group_cut

    # priorytetowa – delikatne cięcie ogona
    if len(df_pri) > 100:
            df_pri = _cut_tail(
                df_pri,
                tail_percentile=pri_tail_percentile,
                label="priority_group",
            )
    else:
            print(
                f"[PLAN] Grupa priorytetowa ma tylko {len(df_pri)} obiektów "
                f"– pomijam tail-cut dla priority_group."
            )

    # niepriorytetowa – mocniejsze cięcie
    df_nonpri = _cut_tail(df_nonpri, tail_percentile=nonpri_tail_percentile, label="non-priority")

    df_final = pd.concat([df_pri, df_nonpri], ignore_index=True)

    # 4) Opcjonalny twardy limit liczby obiektów
    if target_max_objects is not None and len(df_final) > target_max_objects:
        before = len(df_final)
        df_final = df_final.sort_values(
            by=["priority_group", "final_score"],
            ascending=[False, False],
        ).head(target_max_objects).reset_index(drop=True)
        print(
            f"[PLAN] Ostateczny limit do {target_max_objects} obiektów "
            f"(top po priority_group + final_score): {before} -> {len(df_final)}."
        )

    df_final = df_final.drop(columns=["dq_rank"], errors="ignore")

    print(f"[PLAN] Obiekty w bazie po safety-cut: {len(df_final)}.")
    return df_final

from typing import Dict, Any, List  # masz już na górze

def ensure_unique_ids_in_output(objects: List[Dict[str, Any]]) -> None:
    """
    Modyfikuje listę obiektów IN PLACE tak, aby pola 'id' (i spójnie 'name')
    były unikalne. Jeśli to samo id pojawia się kilka razy, kolejnym
    wystąpieniom dopisuje '(2)', '(3)', itd.

    Przykład:
      NGC1027      -> NGC1027
      NGC1027      -> NGC1027 (2)
      NGC1027      -> NGC1027 (3)
    """
    seen: Dict[str, int] = {}

    for obj in objects:
        orig_id = str(obj.get("id", "")).strip()
        if not orig_id:
            continue

        count = seen.get(orig_id, 0) + 1
        seen[orig_id] = count

        if count == 1:
            new_id = orig_id
        else:
            new_id = f"{orig_id} ({count})"

        # spójnie podmieniamy zarówno id, jak i name
        obj["id"] = new_id
        obj["name"] = new_id

# =====================================================
# MAIN
# =====================================================

def main():
    year_mgr = YearManager()
    year = year_mgr.select_interactive()

    locmgr = LocationManager()
    location = locmgr.select_interactive()

    user_params = explain_and_get_prefs()

    print("=" * 141)
    print("KROK 3: FILTROWANIE KATALOGU")
    print("=" * 141)

    catalog = DSOCatalog()
    catalog.load()

    print("\n[INFO] Uzupełnianie brakujących danych (mag, size)...")
    df_all = load_and_impute_catalog(catalog.df)
    
    # UJEDNOLICENIE NAZWY: pracujemy na 'extra_info'
    if "extrainfo" in df_all.columns and "extra_info" not in df_all.columns:
        df_all["extra_info"] = df_all["extrainfo"]

    mag_measured = df_all["has_mag_measured"].sum()
    size_measured = df_all["has_size_measured"].sum()
    both = (df_all["has_mag_measured"] & df_all["has_size_measured"]).sum()

    print(
        f" ✓ Zmierzone mag: {mag_measured}/{len(df_all)} "
        f"({mag_measured/len(df_all)*100:.1f}%)"
    )
    print(
        f" ✓ Zmierzone size: {size_measured}/{len(df_all)} "
        f"({size_measured/len(df_all)*100:.1f}%)"
    )
    print(
        f" ✓ Oba zmierzone: {both}/{len(df_all)} "
        f"({both/len(df_all)*100:.1f}%)"
    )

    df_filtered = filter_candidates_geometric(
        df_all,
        location["lat"],
        user_params["minalt"],
        user_params["minhours"],
    )

    if df_filtered.empty:
        print("Brak obiektów spełniających kryteria geometryczne.")
        return

    print(f"\n✓ {len(df_filtered)} kandydatów do analizy.")
    print("Trwa analiza widoczności z uwzględnieniem Słońca (może zająć chwilę)...")

    step_minutes = 15
    jd_arr, month_arr, lst_arr = get_full_year_smart_grid(
        year=year,
        observer_loc=location,
        sun_limit_deg=user_params["sunlimit"],
        step_minutes=step_minutes,
    )

    all_hours_matrix = compute_all_visibilities_matrix(
        df_filtered,
        location,
        jd_arr,
        month_arr,
        lst_arr,
        min_alt=user_params["minalt"],
        step_minutes=step_minutes,
    )

    visible_candidates: List[tuple] = []
    total = len(df_filtered)

    for i, (_, row) in enumerate(df_filtered.iterrows()):
        if (i + 1) % 50 == 0:
            print(f"{i+1}/{total}...", end="", flush=True)

        hours = all_hours_matrix[i]
        if not np.any(hours >= user_params["minhours"]):
            continue

        maxalt = compute_max_altitude(location["lat"], row["dec"])
        score_result = city_optimized_score(row, user_params)
        if score_result["final_score"] <= 0:
            continue

        visible_candidates.append((row, hours, maxalt, score_result))

    print(
        f"\n✓ Analiza zakończona. Znaleziono "
        f"{len(visible_candidates)} w pełni widocznych obiektów."
    )

    if not visible_candidates:
        print("Brak obiektów.")
        return

    rows: List[Dict[str, Any]] = []
    min_hours = user_params["minhours"]  # już masz userparams z explain_and_get_prefs

    for row, hours, maxalt, scoreresult in visible_candidates:
        rawextra = row.get("extra_info", "")
        if pd.isna(rawextra):
            extraclean = ""
        else:
            extraclean = str(rawextra)
        
        # 1) WYBÓR NAJLEPSZEGO ID KATALOGOWEGO
        original_id = str(row["id"]).strip()
        best_id = choose_best_catalog_id(original_id, extraclean)
        
        # 2) JEŚLI best_id ≠ original_id → dopisz stare ID do extra_info
        if best_id != original_id:
            parts = [x.strip() for x in extraclean.split(",") if x.strip()]
            if original_id not in parts:
                parts.append(original_id)
            extraclean = ", ".join(sorted(set(parts)))

        hours_arr = np.array(hours, dtype=float)
        months_obs = int(np.sum(hours_arr >= min_hours))
        rarity = compute_rarity_score(maxalt, months_obs)
        cname = row.get("common_names", None)

        rows.append({
            "id": row["id"],
            "extra_info": extraclean,
            "type_code": row.get("type_code", ""),
            "ra": float(row["ra"]),
            "dec": float(row["dec"]),
            "size": float(row["size"]),
            "mag": float(row["mag"]) if pd.notna(row["mag"]) else None,
            "base_score": scoreresult["base_score"],
            "final_score": scoreresult["final_score"],
            "quality_factor": scoreresult["quality_factor"],
            "data_quality": scoreresult["data_quality"],
            "size_source": scoreresult["size_source"],
            "mag_source": scoreresult["mag_source"],
            "hours": hours,
            "maxalt": maxalt,
            "common_names": row.get("common_names", None),
            "monthsobservable": months_obs,
            "rarity_score": rarity,
        })

    canddf = pd.DataFrame(rows)

    print("=" * 141)
    print("KROK 4: OPTYMALIZACJA STRON ATLASU (CLUSTERING)")
    print("=" * 141)

    # --- LOGIKA SORTOWANIA (Democja słabych danych) ---

    # 1. Sprawdź, czy obiekt jest "Sławny" (M, C, H lub posiada Common Name)
    def check_is_famous(row):
        # Sprawdzamy flagi w extrainfo (Messier, Caldwell, Herschel)
        extra = str(row.get("extrainfo", ""))
        has_catalog = any(x in extra for x in ["M", "C", "H"])

        # Sprawdzamy czy ma nazwę potoczną (Common Name)
        cname = row.get("common_names")
        has_common_name = (cname is not None) and (not pd.isna(cname)) and (str(cname).strip() != "")

        return has_catalog or has_common_name

    canddf["is_famous"] = canddf.apply(check_is_famous, axis=1)

    # 2. Sprawdź jakość danych (Czy Mag i Size są zmierzone, a nie szacowane?)
    canddf["is_fully_measured"] = (
        (canddf["mag_source"] == "measured") &
        (canddf["size_source"] == "measured")
    )

    # 3. Ustal priorytet:
    # Zostaje na górze JEŚLI: (Jest Sławny) LUB (Ma pełne dane pomiarowe)
    # Spada na dół JEŚLI: (Nie jest sławny) I (Ma dane szacunkowe/częściowe)
    canddf["priority_group"] = canddf["is_famous"] | canddf["is_fully_measured"]

    # 4. Sortowanie
    # Najpierw grupa priorytetowa (True > False), potem wynik punktowy
    canddf = canddf.sort_values(
        by=["priority_group", "final_score"],
        ascending=[False, False],
    ).reset_index(drop=True)

    # Wyświetl statystykę dla pewności
    demoted_count = len(canddf) - canddf["priority_group"].sum()
    if demoted_count > 0:
        print(f"[INFO] Przesunięto na koniec listy {demoted_count} obiektów "
              f"z danymi szacunkowymi (bez nazw własnych/M/C/H).")

    # =========================================================================
    # KROK 4: OPTYMALIZACJA STRON ATLASU (CLUSTERING Z MERGE'OWANIEM NAZW)
    # =========================================================================

    ATLAS_FOV_RADIUS = 2.0 * u.deg
    final_atlas_objects: List[pd.Series] = []
    covered_indices: set[int] = set()

    coords = SkyCoord(
        ra=canddf["ra"].values * u.deg,
        dec=canddf["dec"].values * u.deg,
    )
    # Definicja funkcji rankingu lokalnie, aby mieć do niej dostęp w pętli
    def get_catalog_rank(obj_id: str) -> int:
        s = str(obj_id).lower().strip()
        if s.startswith("ngc"): return 0
        if s.startswith("ic"): return 1
        if s.startswith("sh2"): return 2
        if s.startswith("rcw"): return 3
        if s.startswith("lbn"): return 4
        if s.startswith("ced"): return 5
        if s.startswith("pgc"): return 6
        if s.startswith("barn"): return 7 # b lub barn
        if s.startswith("b"): return 7    # obsługa skrótu Barnard
        if s.startswith("ldn"): return 8
        return 99 # inne
        
    print(f"\n[INFO] Rozpoczynam grupowanie kadrów (promień {ATLAS_FOV_RADIUS})...")
    
    for i in range(len(canddf)):
        # Jeśli obiekt został już "wchłonięty", pomijamy go
        if i in covered_indices:
            continue

        # Kopia lidera
        leader = canddf.iloc[i].copy()
        
        # Zbiór wszystkich ID w tej grupie (Lider + Sąsiedzi + ich Extra Info)
        name_pool = set()

        # 1. Dodajemy ID i Extra Info Lidera do puli
        name_pool.add(str(leader['id']).strip())
        if pd.notna(leader.get('extra_info')):
            # Dodajemy wszystko co jest w extra_info lidera
            for x in str(leader['extra_info']).split(','):
                if x.strip(): name_pool.add(x.strip())

        # Szukamy sąsiadów
        separations = coords[i].separation(coords)
        neighbor_indices = np.where(separations < ATLAS_FOV_RADIUS)[0]

        # Listy do zbierania Common Names
        gathered_common_names = []

        # Pobieramy Common Name lidera
        leader_cn = leader.get('common_names')
        if pd.notna(leader_cn) and str(leader_cn).strip() and str(leader_cn).lower() != 'nan':
            gathered_common_names.append(str(leader_cn).strip())

        # PĘTLA PO SĄSIADACH
        for n_idx in neighbor_indices:
            if n_idx == i:
                continue

            # Oznaczamy sąsiada jako obsłużonego
            covered_indices.add(n_idx)
            neighbor_row = canddf.iloc[n_idx]
            
            # 2. Dodajemy ID i Extra Info Sąsiada do puli
            name_pool.add(str(neighbor_row['id']).strip())
            
            n_extra = neighbor_row.get('extra_info')
            if pd.notna(n_extra):
                 for x in str(n_extra).split(','):
                     if x.strip(): name_pool.add(x.strip())

            # Pobieramy Common Name sąsiada
            n_cn = neighbor_row.get('common_names')
            if pd.notna(n_cn) and str(n_cn).strip() and str(n_cn).lower() != 'nan':
                n_cn_str = str(n_cn).strip()
                if n_cn_str not in gathered_common_names:
                    gathered_common_names.append(n_cn_str)

        # --- KONIEC PĘTLI PO SĄSIADACH ---
        # (Kod poniżej musi być cofnięty o jedno wcięcie w lewo względem Twojej wersji)

        # 3. Sortujemy pulę nazw i wybieramy najlepszą DLA CAŁEJ GRUPY
        valid_names = [n for n in name_pool if n and n.lower() != 'nan']
        
        if valid_names:
            sorted_names = sorted(valid_names, key=lambda x: (get_catalog_rank(x), x))
            
            # Zwycięzca zostaje nowym ID Lidera
            best_id = sorted_names[0]
            leader['id'] = best_id
            
            # Reszta trafia do extra_info
            remaining_ids = sorted_names[1:]
            leader['extra_info'] = ", ".join(remaining_ids)

        # 4. Aktualizacja Common Names
        if gathered_common_names:
            leader['common_names'] = ", ".join(gathered_common_names)

        final_atlas_objects.append(leader)

    print(
        f"✓ Zredukowano listę obiektów z {len(canddf)} "
        f"do {len(final_atlas_objects)} unikalnych kadrów."
    )

    # Odtwarzamy DataFrame
    canddf = pd.DataFrame(final_atlas_objects)

    # =========================================================================
    # 1. FILTROWANIE (LOGIKA)
    # =========================================================================
    min_size_arcmin = user_params["minsizearcmin"]
    count_before = len(canddf)

    # FILTR: odsiać obiekty ze znanym size < min_size_arcmin
    canddf = canddf[
        canddf["size"].isna() | (canddf["size"] >= min_size_arcmin)
    ]

    count_after = len(canddf)

    # Wyświetlamy informację o filtracji 
    if count_after < count_before:
        diff = count_before - count_after
        print(f"\n[INFO] Ukryto {diff} obiektów mniejszych niż {min_size_arcmin:.1f}' "
              f"(filtr rozmiaru z kroku 2.E).")
        print(f" Pozostało {count_after} obiektów do wyświetlenia.")

    # Dodanie kategorii rozmiaru (tylko dla pozostałych obiektów)
    def _size_cat(s):
        if pd.isna(s):
            return "? UNKNOWN"
        if s >= min_size_arcmin:
            return "●●● LARGE"
        if s >= min_size_arcmin * 0.4:
            return "●● MEDIUM"
        return "● SMALL (crop)"

    canddf["size_category"] = canddf["size"].apply(_size_cat)

    # *** NOWY KROK: SAFETY CUT PO WSZYSTKIM ***

    canddf = apply_safety_cut(
        canddf,
        global_min_score=0.0,
        min_len_for_cut=500,
        main_percentile=50.0,
        pri_tail_percentile=10.0,
        nonpri_tail_percentile=30.0,
        target_max_objects=None,  # na start bez twardego limitu
    )

    # =========================================================================
    # 2. RYSOWANIE TABELI (PREZENTACJA)
    # =========================================================================
    print("\n" + "=" * 141)
    print(
        f"{'Nr':>3} {'Nazwa':10s} {'Extra':15s} {'Common name':30s} {'Typ':23s} "
        f"{'RA':>7s} {'Dec':>7s} {'Mag':>7s} {'Size':>9s} "
        f"{'Quality':>11s} {'Score':>7s}"
    )
    print("-" * 141)

    choices_map: Dict[int, Dict[str, Any]] = {}
    MAX_PRINT_ROWS = 200
    table_truncated = False  # Inicjalizacja flagi przed pętlą

    for i, (_, crow) in enumerate(canddf.iterrows(), start=1):
        clean_name = str(crow["id"]).strip()
        _, objtypefull = catalog.get_display_info(crow)
        raw_extra = crow.get("extra_info", "")
        if pd.isna(raw_extra) or str(raw_extra).lower() == 'nan':
            raw_extra = ""
        else:
            raw_extra = str(raw_extra)

        extra_str = (raw_extra[:12] + "..." if len(raw_extra) > 12 else raw_extra).ljust(15)

        commonname = crow.get("common_names", None)
        commonname_display = (
            (str(commonname)[:27] + "..." if len(str(commonname)) > 27 else str(commonname)).ljust(30)
            if commonname is not None and not pd.isna(commonname)
            else ""
        )

        ra_val = crow.get("ra", None)
        dec_val = crow.get("dec", None)
        ra_str = f"{float(ra_val):.1f}" if pd.notna(ra_val) else ""
        dec_str = f"{float(dec_val):.1f}" if pd.notna(dec_val) else ""

        mag_val = crow.get("mag", None)
        mag_source = crow.get("mag_source", "unknown")
        mag_icon = "●" if mag_source == "measured" else "○"
        mag_str = (
            f"{mag_icon}{float(mag_val):.1f}" if pd.notna(mag_val) else ""
        )

        size_val = crow.get("size", None)
        size_source = crow.get("size_source", "unknown")
        size_icon = "●" if size_source == "measured" else "○"
        size_cat = crow.get("size_category", "")
        size_str = (
            f"{size_icon}{float(size_val):.1f}'" if pd.notna(size_val) else ""
        )

        data_quality = crow.get("data_quality", "")
        score_val = crow.get("final_score", 0)
        score_str = f"{float(score_val):.0f}"

        # Sprawdzenie limitu wierszy
        if i > MAX_PRINT_ROWS:
            print("=" * 141)
            print(f"{f'... i {len(canddf) - MAX_PRINT_ROWS} więcej obiektów (limit wyświetlania).':>141}")
            table_truncated = True
            break

        print(
            f"{i:3d}. {clean_name:10s} {extra_str:15s} {commonname_display:30s} {objtypefull:23s}"
            f"{ra_str:>7s} {dec_str:>7s} {mag_str:>7s} {size_str:>9s} "
            f"{data_quality:11s} {score_str:>7s} "
        )

        choices_map[i] = {
            "name": clean_name,
            "id": crow["id"],
            "extrainfo": crow["extra_info"],
            "ra": round(float(crow["ra"]), 3),
            "dec": round(float(crow["dec"]), 3),
            "size": (
                round(float(crow["size"]), 3)
                if pd.notna(crow["size"])
                else None
            ),
            "mag": None
            if pd.isna(crow["mag"])
            else round(float(crow["mag"]), 2),
            "type": objtypefull,
            "typecode": crow["type_code"],
            "score": round(float(crow["final_score"]), 2),
            "base_score": round(float(crow["base_score"]), 2),
            "maxalt": round(float(crow["maxalt"]), 2),
            "hours": np.round(crow["hours"], 1).tolist(),
            "commonname": None
            if commonname is None or pd.isna(commonname)
            else str(commonname),
            "data_quality": data_quality,
            "size_category": size_cat,
        }

    # Stopka (tylko jeśli tabela nie została ucięta przez limit)
    if not table_truncated:
        print("=" * 141)
        print(f"{'Koniec listy':>141}")

    # =========================================================================
    # 3. ZAPIS
    # =========================================================================

    outputdata = {
        "location": location,
        "year": year,
        "parameters": {
            "minalt": user_params["minalt"],
            "sunlimit": user_params["sunlimit"],
            "minhours": user_params["minhours"],
            "minsizearcmin": user_params["minsizearcmin"],
            "bortlerange": user_params["bortle_range"],
            "hasnarrowband": user_params["has_narrowband"],
            "preferfamous": user_params["prefer_famous"],
            "camera": user_params["camera"],
        },
        "objects": [],
    }
    
    for _, obj in canddf.iterrows():  
        hours_total = np.array(obj["hours"], dtype=float)
        months_obs = int(np.sum(hours_total >= user_params["minhours"]))
        rarity = compute_rarity_score(obj["maxalt"], months_obs)
        cname = obj.get("common_names", "")

        outputdata["objects"].append({
            "id": obj["id"],
            "name": obj["id"],
            "ra": float(obj["ra"]),
            "dec": float(obj["dec"]),
            "score": float(obj["final_score"]),
            "type": obj.get("type_code", ""),
            "size": float(obj["size"]) if obj["size"] is not None else 0.0,
            "mag": float(obj["mag"]) if obj["mag"] is not None else 0.0,
            "commonname": "" if pd.isna(cname) else str(cname).strip(),
            "extra_info": obj.get("extra_info", ""),
            "monthly_hours_total": hours_total.tolist(),
            "max_altitude_year": round(float(obj["maxalt"]), 1),
            "months_observable": months_obs,
            "rarity_score": round(float(rarity), 2),
        })
    ensure_unique_ids_in_output(outputdata["objects"])
    with open(Config.OUTPUT_FILE, "w", encoding="utf-8") as f:
        json.dump(outputdata, f, indent=2, ensure_ascii=False)

    print(f"\n✓ Zapisano dane do pliku: {Config.OUTPUT_FILE}")
    print(f"✓ Liczba wybranych obiektów: {len(canddf)}")

if __name__ == "__main__":
    main()