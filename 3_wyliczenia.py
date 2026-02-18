import json
import pickle
import hashlib
import os
from datetime import datetime, timedelta
from concurrent.futures import ProcessPoolExecutor, as_completed
from functools import partial
from dataclasses import dataclass
from typing import Dict, List, Any, Tuple

import numpy as np
import pandas as pd
from astropy.coordinates import SkyCoord, EarthLocation, get_sun, get_body, AltAz
from astropy.time import Time
import astropy.units as u
import pytz
from tqdm import tqdm

# =====================================================
# KONFIGURACJA
# =====================================================

H_START = 14.0
H_END = 31.0
H_RANGE = H_END - H_START
N_SAMPLES = int(H_RANGE * 60 / 5)

CROSSING_SAMPLES = 1440

# Pliki cache / stanu
RAW_DATA_PKL = "observing_data_raw.pkl"
FINAL_DATA_PKL = "observing_data.pkl"
FINAL_HASH_FILE = "observing_data_final.hash"
ENGINE_STATE_PKL = "observing_engine_state.pkl"

VERBOSE = True  # globalny przełącznik logów

# =====================================================
# DATACLASS DO PRZECHOWYWANIA SUROWYCH DANYCH
# =====================================================

@dataclass
class RawObjectData:
    """Dane niezależne od parametrów filtrowania (wysokości, crossingi)."""
    obj_id: str
    o_alt_all: np.ndarray  # [ndays, N_SAMPLES]
    sun_pts_list: List[List[datetime]]  # [ndays]
    obj_pts_list: List[List[datetime]]  # [ndays]


# =====================================================
# FUNKCJE POMOCNICZE
# =====================================================

def load_vis_data(path: str = "vis_data.json") -> Dict[str, Any]:
    with open(path, "r", encoding="utf-8") as f:
        return json.load(f)

def save_vis_data(data: Dict[str, Any], path: str = "vis_data.json"):
    with open(path, "w", encoding="utf-8") as f:
        json.dump(data, f, indent=2, ensure_ascii=False)

def get_params_hash(
    obj_min_alt_deg: float,
    sun_alt_limit_deg: float,
    lat: float,
    lon: float,
    year: int,
) -> str:
    """
    Hash parametrów, które wpływają na wynik obliczeń:
    - minimalna wysokość obiektu
    - limit wysokości Słońca
    - lokalizacja (lat, lon)
    - rok
    """
    s = (
        f"{obj_min_alt_deg:.6f}|{sun_alt_limit_deg:.6f}|"
        f"{lat:.6f}|{lon:.6f}|{int(year)}"
    )
    return hashlib.md5(s.encode()).hexdigest()

def get_crossings(target, height, t_noon, location) -> List[Time]:
    times = t_noon + np.linspace(0, 24, CROSSING_SAMPLES) * u.hour
    altaz = target.transform_to(AltAz(obstime=times, location=location))
    alts = altaz.alt.deg
    diff = alts - height
    idx = np.where(np.diff(np.sign(diff)))[0]

    res = []
    for i in idx:
        t_exact = times[i] + (times[i + 1] - times[i]) * (-diff[i] / (diff[i + 1] - diff[i]))
        res.append(t_exact)

    return sorted(res)

def load_engine_state() -> dict:
    if os.path.exists(ENGINE_STATE_PKL):
        try:
            with open(ENGINE_STATE_PKL, "rb") as f:
                return pickle.load(f)
        except Exception:
            return {}
    return {}

def save_engine_state(state: dict):
    with open(ENGINE_STATE_PKL, "wb") as f:
        pickle.dump(state, f)

def mask_to_segments(mask, h_start=H_START, h_end=H_END) -> List[Tuple[float, float]]:
    mask = np.asarray(mask, dtype=bool)
    if mask.size == 0 or not mask.any():
        return []

    diff = np.diff(mask.astype(int))
    starts_idx = list(np.where(diff == 1)[0] + 1)
    ends_idx = list(np.where(diff == -1)[0] + 1)

    if mask[0]:
        starts_idx = [0] + starts_idx
    if mask[-1]:
        ends_idx = ends_idx + [mask.size]

    segments = []
    t_grid = np.linspace(h_start, h_end, mask.size, endpoint=False)

    for s_i, e_i in zip(starts_idx, ends_idx):
        t_start = t_grid[s_i]
        t_end = t_grid[e_i] if e_i < mask.size else h_end

        start_rel = t_start - h_start
        end_rel = t_end - h_start
        if end_rel > start_rel:
            segments.append((start_rel, end_rel))

    return segments

# =====================================================
# ETAP 1: OBLICZANIE DANYCH SUROWYCH (RAW)
# =====================================================

def process_single_object_raw(
    obj_data: Tuple[str, float, float],
    days: pd.DatetimeIndex,
    location: EarthLocation,
    days_noon: Time,
    t_night_offsets_hours: np.ndarray,
    precomputed_sun_pts: List[List[datetime]],  # <--- NOWY ARGUMENT
) -> RawObjectData:
    obj_id, ra, dec = obj_data
    coord = SkyCoord(ra * u.deg, dec * u.deg)
    ndays = len(days)

    # 1. Obliczanie wysokości w siatce godzinowej (wektoryzowane - szybkie)
    t_grid = days_noon[:, None] + t_night_offsets_hours[None, :] * u.hour
    frame = AltAz(obstime=t_grid.reshape(-1), location=location)
    o_alt_all = coord.transform_to(frame).alt.deg.reshape(ndays, N_SAMPLES)

    obj_pts_list = []

    # 2. Pętla po dniach - teraz liczy TYLKO crossingi obiektu
    # Słońce (sun_pts) jest już policzone i przekazane w argumencie
    for day_idx, t_noon in enumerate(days_noon):
        # USUNIĘTO: sun_pts = get_crossings(get_sun(t_noon), ...)
        
        # Liczymy tylko obiekt
        obj_pts = get_crossings(coord, 0, t_noon, location)
        obj_pts_list.append([p.datetime for p in obj_pts])

    return RawObjectData(
        obj_id=obj_id,
        o_alt_all=o_alt_all,
        sun_pts_list=precomputed_sun_pts,  # Zwracamy gotowe dane słońca
        obj_pts_list=obj_pts_list,
    )

def compute_raw_data(
    json_path: str,
    object_limit: int,
    max_workers: int = None,
) -> Dict[str, RawObjectData]:
    """
    Oblicza LUB UZUPEŁNIA dane surowe dla obiektów z vis_data.json
    ograniczonych do object_limit.
    """
    vis = load_vis_data(json_path)

    lat = vis["location"]["lat"]
    lon = vis["location"]["lon"]
    location = EarthLocation(lat=lat * u.deg, lon=lon * u.deg)
    year = vis["year"]

    params = vis.get("parameters", {})
    obj_min_alt_deg = params.get("minalt", 20.0)
    sun_alt_limit_deg = params.get("sunlimit", -6.0)

    objects = vis["objects"][:object_limit]
    objects_data = [(obj["id"], obj["ra"], obj["dec"]) for obj in objects]

    days = pd.date_range(f"{year}-01-01", f"{year}-12-31", freq="D")
    ndays = len(days)

    # === SMART CACHE DLA RAW ===
    current_raw_hash = get_params_hash(
        obj_min_alt_deg, sun_alt_limit_deg, lat, lon, year
    )
    raw_hash_file = "observing_data_raw.hash"

    # === INICJALIZACJA CACHE ===
    all_raw_data: Dict[str, RawObjectData] = {}
    cache_valid = False
    
    # Sprawdź hash + wczytaj cache
    if os.path.exists(RAW_DATA_PKL) and os.path.exists(raw_hash_file):
        try:
            with open(raw_hash_file, "r") as f:
                saved_hash = f.read().strip()
            
            if current_raw_hash == saved_hash:
                with open(RAW_DATA_PKL, "rb") as f:
                    all_raw_data = pickle.load(f)
                cache_valid = True
                if VERBOSE:
                    print(f"[CACHE] ✓ Załadowano {len(all_raw_data)} obiektów z RAW cache (RAW hash: {current_raw_hash[:8]}...).")
            else:
                if VERBOSE:
                    print(f"[RAW CACHE] Invalidated! saved_hash={saved_hash[:8]} ≠ current={current_raw_hash[:8]}")
                os.remove(RAW_DATA_PKL)
                os.remove(raw_hash_file)
        except Exception as e:
            if VERBOSE:
                print(f"[ERROR] Błąd walidacji RAW hash/cache: {e}")
    else:
        if VERBOSE:
            print("[RAW CACHE] Brak hash lub RAW cache - obliczam od nowa.")

    if not cache_valid:
        all_raw_data = {}  # ← KLUCZOWE! Pusty cache
        if VERBOSE:
            print(f"[RAW]        Brak ważnego RAW cache - przygotowuję przeliczenie wszystkich obiektów.")

    # Oblicz brakujące (lub wszystkie)
    existing_ids = set(all_raw_data.keys())
    missing_data = [
        (obj_id, ra, dec)
        for (obj_id, ra, dec) in objects_data
        if obj_id not in existing_ids
    ]

    if VERBOSE:
        if missing_data:
            print(f"[RAW]     Brakujących obiektów w RAW cache: {len(missing_data)}.")
        else:
            print("[RAW]      RAW cache zawiera wszystkie obiekty z aktualnego limitu.")

    if not missing_data:
        # ZAPIS HASH (nawet jeśli nie przeliczaliśmy)
        with open(raw_hash_file, "w") as f:
            f.write(current_raw_hash)
        return all_raw_data

    days_noon = Time([
        datetime.combine(d.date(), datetime.min.time()) + timedelta(hours=12)
        for d in days
    ])

    t_night_offsets_hours = np.linspace(2, 19, N_SAMPLES)

    # === OPTYMALIZACJA: PRE-COMPUTE SUN ===
    # Liczymy dane słoneczne RAZ, zamiast robić to w każdym wątku dla każdego obiektu.
    if VERBOSE:
        print(f"[PRE-COMPUTE] Obliczam dane Słońca dla {ndays} dni.")
    
    precomputed_sun_pts = []
    # To jest szybkie (365 iteracji), więc robimy to w głównym wątku
    for t_noon in tqdm(days_noon, desc="          Sun Calc", ncols=119):
        # 0 oznacza horyzont (wschód/zachód słońca)
        s_pts = get_crossings(get_sun(t_noon), 0, t_noon, location)
        precomputed_sun_pts.append([p.datetime for p in s_pts])

    if VERBOSE:
        print(f"[COMPUTE] Obliczam RAW dla {len(missing_data)} brakujących obiektów...")

    process_raw_func = partial(
        process_single_object_raw,
        days=days,
        location=location,
        days_noon=days_noon,
        t_night_offsets_hours=t_night_offsets_hours,
        precomputed_sun_pts=precomputed_sun_pts,  # <--- Przekazujemy tutaj
    )

    with ProcessPoolExecutor(max_workers=max_workers) as executor:
        future_to_objid = {
            executor.submit(process_raw_func, obj_data): obj_data[0]
            for obj_data in missing_data
        }

        with tqdm(
            total=len(missing_data),
            desc="          Raw compute",
            unit="obj",
            ncols=119,
        ) as pbar:
            for future in as_completed(future_to_objid):
                obj_id = future_to_objid[future]
                try:
                    raw_result = future.result()
                    all_raw_data[raw_result.obj_id] = raw_result
                    pbar.set_postfix_str(f"ostatni={obj_id}")
                except Exception as e:
                    print(f"[ERROR] Błąd dla {obj_id}: {e}")
                pbar.update(1)

    # === ZAPIS RAW + HASH ===
    with open(RAW_DATA_PKL, "wb") as f:
        pickle.dump(all_raw_data, f)
    
    with open(raw_hash_file, "w") as f:
        f.write(current_raw_hash)
    
    if VERBOSE:
        print(f"[SAVE] Raw data zapisane do {RAW_DATA_PKL} ({len(all_raw_data)} obiektów).")
        print(f"[SAVE] RAW hash zapisany: {current_raw_hash[:8]}...")

    return all_raw_data

# =====================================================
# ETAP 2: KONWERSJA RAW -> FINAL
# =====================================================

def process_raw_to_final(
    raw_data: RawObjectData,
    days: pd.DatetimeIndex,
    sun_alt_all: np.ndarray,
    moon_alt_all: np.ndarray,
    sun_alt_limit_deg: float,
    obj_min_alt_deg: float,
    vis: dict,
) -> Tuple[str, List[Dict[str, Any]]]:
    obj_id = raw_data.obj_id
    o_alt_all = raw_data.o_alt_all
    ndays = len(days)
    results = []

    t_step = H_RANGE / (N_SAMPLES - 1)
    tz_name = vis["location"]["tz"]
    tz = pytz.timezone(tz_name)

    for day_idx, d in enumerate(days):
        d_date = d.date() if hasattr(d, "date") else d
        d_midnight = datetime.combine(d_date, datetime.min.time())
        d_midnight_local = tz.localize(d_midnight)
        tz_offset = d_midnight_local.utcoffset().total_seconds() / 3600.0

        o_alt = o_alt_all[day_idx]
        s_alt = sun_alt_all[day_idx]
        m_alt = moon_alt_all[day_idx]

        quality_mask = (o_alt > obj_min_alt_deg) & (s_alt < sun_alt_limit_deg)

        transit_rel = np.argmax(o_alt) * (H_RANGE / N_SAMPLES)
        q_hours = np.sum(quality_mask) * (H_RANGE / N_SAMPLES)
        m_hours = np.sum(quality_mask & (m_alt < 0.0)) * (H_RANGE / N_SAMPLES)
        qual_segments = mask_to_segments(quality_mask)

        new_sun_pts = []
        diff = s_alt - sun_alt_limit_deg
        crossings = np.where(np.diff(np.signbit(diff)))[0]

        if len(crossings) > 0:
            for idx in crossings:
                y0, y1 = diff[idx], diff[idx + 1]
                fraction = -y0 / (y1 - y0)
                h_val = H_START + (idx + fraction) * t_step
                pt_time = d_midnight + timedelta(hours=h_val)
                new_sun_pts.append(pt_time)
        else:
            mean_alt = np.mean(s_alt)
            if mean_alt > sun_alt_limit_deg:
                dummy_t = d_midnight + timedelta(hours=H_START)
                new_sun_pts = [dummy_t, dummy_t]
            else:
                new_sun_pts = []

        results.append({
            "day": d_date,
            "sun_pts": new_sun_pts,
            "obj_pts": raw_data.obj_pts_list[day_idx],
            "transit_rel": transit_rel,
            "q_hours": q_hours,
            "m_hours": m_hours,
            "qual_segments": qual_segments,
            "tz_offset": tz_offset,
        })

    return obj_id, results

def reprocess_to_final(
    json_path: str,
    all_raw_data: Dict[str, RawObjectData],
    max_workers: int = None,  # Dodano parametr
) -> Tuple[Dict[str, List[Dict]], str]:
    vis = load_vis_data(json_path)

    params = vis.get("parameters", {})
    obj_min_alt_deg = params.get("minalt", 20.0)
    sun_alt_limit_deg = params.get("sunlimit", -6.0)

    year = vis["year"]
    lat = vis["location"]["lat"]
    lon = vis["location"]["lon"]
    location = EarthLocation(lat=lat * u.deg, lon=lon * u.deg)

    days = pd.date_range(f"{year}-01-01", f"{year}-12-31", freq="D")
    ndays = len(days)

    current_params_hash = get_params_hash(obj_min_alt_deg, sun_alt_limit_deg, lat, lon, year)

    if VERBOSE:
        print(f"[INFO] Przeliczanie FINAL (Parallel): konwersja {len(all_raw_data)} obiektów.")

    days_noon = Time([
        datetime.combine(d.date(), datetime.min.time()) + timedelta(hours=12)
        for d in days
    ])

    t_night_offsets_hours = np.linspace(2, 19, N_SAMPLES)
    t_grid_all = days_noon[:, None] + t_night_offsets_hours[None, :] * u.hour
    obstime_all = t_grid_all.reshape(-1)

    frame_all = AltAz(obstime=obstime_all, location=location)

    # Obliczenia wspólne (wykonywane raz w głównym procesie)
    sun_icrs_all = get_sun(obstime_all)
    sun_alt_all = sun_icrs_all.transform_to(frame_all).alt.deg.reshape(ndays, N_SAMPLES)

    moon_icrs_all = get_body("moon", obstime_all)
    moon_alt_all = moon_icrs_all.transform_to(frame_all).alt.deg.reshape(ndays, N_SAMPLES)

    all_final_data: Dict[str, List[Dict]] = {}

    # Przygotowanie funkcji partial z argumentami stałymi
    process_final_func = partial(
        process_raw_to_final,
        days=days,
        sun_alt_all=sun_alt_all,
        moon_alt_all=moon_alt_all,
        sun_alt_limit_deg=sun_alt_limit_deg,
        obj_min_alt_deg=obj_min_alt_deg,
        vis=vis,
    )

    # Uruchomienie ProcessPoolExecutor
    # Uwaga: Przekazujemy tylko raw_data jako zmienną, reszta jest w partial
    with ProcessPoolExecutor(max_workers=max_workers) as executor:
        # Tworzymy listę futures
        futures = {
            executor.submit(process_final_func, raw_data): raw_data.obj_id 
            for raw_data in all_raw_data.values()
        }

        for future in tqdm(
            as_completed(futures),
            desc="          Reprocess",
            unit="obj",
            ncols=119,
            total=len(futures)
        ):
            try:
                obj_id_result, results = future.result()
                all_final_data[obj_id_result] = results
            except Exception as e:
                print(f"[ERROR] Błąd w FINAL dla obiektu: {e}")

    if VERBOSE:
        print(f"[REPROCESS] ✓ Przetworzono {len(all_final_data)} obiektów.")

    return all_final_data, current_params_hash

def should_reprocess(json_path: str) -> bool:
    """
    True jeśli trzeba reprocessować FINAL (parametry lub lokalizacja/rok się zmieniły).
    """
    if not os.path.exists(FINAL_HASH_FILE):
        return True

    vis = load_vis_data(json_path)
    params = vis.get("parameters", {})
    obj_min_alt_deg = params.get("minalt", 20.0)
    sun_alt_limit_deg = params.get("sunlimit", -6.0)
    year = vis["year"]
    lat = vis["location"]["lat"]
    lon = vis["location"]["lon"]

    current_hash = get_params_hash(obj_min_alt_deg, sun_alt_limit_deg, lat, lon, year)

    try:
        with open(FINAL_HASH_FILE, "r") as f:
            old_hash = f.read().strip()
        return current_hash != old_hash
    except Exception:
        return True

# =====================================================
# PRZYROSTOWY FINAL CACHE
# =====================================================

def get_target_ids_from_vis(vis: Dict[str, Any], object_limit: int) -> List[str]:
    objects = vis.get("objects", [])
    ids = [o["id"] for o in objects]
    return ids[:object_limit]

def reprocess_missing_final(
    json_path: str,
    all_raw_data: Dict[str, RawObjectData],
    existing_final: Dict[str, List[Dict]],
    target_ids: List[str],
    max_workers: int = None, # Dodano parametr
) -> Tuple[Dict[str, List[Dict]], str]:
    # ... (kod inicjalizacji zmiennych bez zmian: vis, params, year, location, days) ...
    vis = load_vis_data(json_path)
    params = vis.get("parameters", {})
    obj_min_alt_deg = params.get("minalt", 20.0)
    sun_alt_limit_deg = params.get("sunlimit", -6.0)
    year = vis["year"]
    lat = vis["location"]["lat"]
    lon = vis["location"]["lon"]
    location = EarthLocation(lat=lat * u.deg, lon=lon * u.deg)
    
    days = pd.date_range(f"{year}-01-01", f"{year}-12-31", freq="D")
    ndays = len(days)
    current_params_hash = get_params_hash(obj_min_alt_deg, sun_alt_limit_deg, lat, lon, year)

    # Obliczenia wspólne Słońca/Księżyca
    days_noon = Time([
        datetime.combine(d.date(), datetime.min.time()) + timedelta(hours=12)
        for d in days
    ])
    t_night_offsets_hours = np.linspace(2, 19, N_SAMPLES)
    t_grid_all = days_noon[:, None] + t_night_offsets_hours[None, :] * u.hour
    obstime_all = t_grid_all.reshape(-1)
    frame_all = AltAz(obstime=obstime_all, location=location)

    sun_icrs_all = get_sun(obstime_all)
    sun_alt_all = sun_icrs_all.transform_to(frame_all).alt.deg.reshape(ndays, N_SAMPLES)
    moon_icrs_all = get_body("moon", obstime_all)
    moon_alt_all = moon_icrs_all.transform_to(frame_all).alt.deg.reshape(ndays, N_SAMPLES)

    existing_ids = set(existing_final.keys())
    missing_ids = [oid for oid in target_ids if oid not in existing_ids]

    if not missing_ids:
        final_subset = {oid: existing_final[oid] for oid in target_ids if oid in existing_final}
        return final_subset, current_params_hash

    to_process_raw = [
        all_raw_data[oid]
        for oid in missing_ids
        if oid in all_raw_data
    ]

    if VERBOSE:
        print(f"[INFO] Przeliczam FINAL (Parallel) dla {len(to_process_raw)} brakujących obiektów.")

    new_final: Dict[str, List[Dict]] = {}
    
    process_final_func = partial(
        process_raw_to_final,
        days=days,
        sun_alt_all=sun_alt_all,
        moon_alt_all=moon_alt_all,
        sun_alt_limit_deg=sun_alt_limit_deg,
        obj_min_alt_deg=obj_min_alt_deg,
        vis=vis,
    )

    with ProcessPoolExecutor(max_workers=max_workers) as executor:
        futures = {
            executor.submit(process_final_func, raw_data): raw_data.obj_id 
            for raw_data in to_process_raw
        }

        for future in tqdm(
            as_completed(futures),
            desc="       Reprocess",
            unit="obj",
            ncols=119,
            total=len(futures),
        ):
             try:
                obj_id_result, results = future.result()
                new_final[obj_id_result] = results
             except Exception as e:
                 print(f"[ERROR] Błąd w FINAL: {e}")

    merged = dict(existing_final)
    merged.update(new_final)
    final_subset = {oid: merged[oid] for oid in target_ids if oid in merged}

    return final_subset, current_params_hash

# =====================================================
# GŁÓWNA FUNKCJA
# =====================================================

def run_engine_from_vis_json(
    json_path: str = "vis_data.json",
    max_workers: int = None,
    force_all: bool = False,
) -> Dict[str, List[Dict]]:
    vis = load_vis_data(json_path)

    city_name = vis["location"].get("name", "Unknown")
    year = vis["year"]
    lat = vis["location"]["lat"]
    lon = vis["location"]["lon"]
    params = vis.get("parameters", {})
    obj_min_alt_deg = params.get("minalt", 20.0)
    sun_alt_limit_deg = params.get("sunlimit", -6.0)

    prev_state = load_engine_state()
    if VERBOSE:
        if prev_state:
            print(
                f"[INFO] Poprzednia lokalizacja: {prev_state.get('city_name', 'Unknown')}, "
                f"({prev_state.get('lat', float('nan')):.2f}, {prev_state.get('lon', float('nan')):.2f}), rok {prev_state.get('year')}, "
                f"nad horyzontem {prev_state.get('minalt')}, zmierzch {prev_state.get('sunlimit')}."
            )
            print(
                    f"[INFO] Bieżąca lokalizacja:    {city_name}, "
                    f"({lat:.2f}, {lon:.2f}), rok {year}, "
                    f"nad horyzontem {obj_min_alt_deg}, zmierzch {sun_alt_limit_deg}."
                )
        else:
            print("[INFO] Brak poprzedniego stanu silnika (ENGINE_STATE_PKL).")

    vis_objects = vis.get("objects", [])
    total_objects = len(vis_objects)

    if VERBOSE:
        print(f"[INFO] W vis_data.json jest {total_objects} obiektów.")

    default_limit = 108
    max_limit = total_objects

    try:
        user_input = input(
            f"       Ile obiektów przeliczyć? [domyślnie {default_limit}, Enter = wszystkie]: "
        ).strip()

        if user_input == "":
            object_limit = max_limit
        else:
            val = int(user_input)
            if 1 <= val <= max_limit:
                object_limit = val
            else:
                print(f"[WARN] Podano liczbę spoza zakresu 1..{max_limit}, używam wszystkich.")
                object_limit = max_limit
    except Exception:
        print("[WARN] Błąd wejścia, używam wszystkich obiektów.")
        object_limit = max_limit

    if VERBOSE:
        print(f"[INFO] Do przeliczenia: {object_limit} obiektów.")

    if VERBOSE:
        lat_str = f"{lat:.2f}"
        lon_str = f"{lon:.2f}"
        print("=" * 119)
        print("ENGINE: Smart Cache System")
        print("=" * 119)
        print(f"[INFO] Lokalizacja: {city_name} ({lat_str}°, {lon_str}°)")
        print(f"[INFO] Rok: {year}")
        print(f"[INFO] Parametry: minimalna wysokość obiektu: {obj_min_alt_deg}°, "
              f"wysokość słońca: {sun_alt_limit_deg}°")
        print(f"[INFO] Cache RAW:   {RAW_DATA_PKL}")
        print(f"[INFO] Cache FINAL: {FINAL_DATA_PKL}")
        print("=" * 119)

    # --- FORCE-ALL: usunięcie cache ---
    if force_all:
        if os.path.exists(RAW_DATA_PKL):
            os.remove(RAW_DATA_PKL)
            if VERBOSE:
                print("[FORCE] Usunięto RAW cache.")
        if os.path.exists(FINAL_DATA_PKL):
            os.remove(FINAL_DATA_PKL)
            if VERBOSE:
                print("[FORCE] Usunięto FINAL cache.")
        if os.path.exists(FINAL_HASH_FILE):
            os.remove(FINAL_HASH_FILE)
            if VERBOSE:
                print("[FORCE] Usunięto plik hash FINAL.")

    # ========== ETAP 1: RAW ==========
    all_raw_data = compute_raw_data(json_path, object_limit, max_workers=max_workers)

    # ========== ETAP 2: FINAL ==========
    needs_reprocess = should_reprocess(json_path)
    target_ids = get_target_ids_from_vis(vis, object_limit)

    if needs_reprocess or not os.path.exists(FINAL_DATA_PKL):
        if VERBOSE:
            print("[REPROCESS] Parametry/lokalizacja/rok się zmieniły lub brak FINAL cache – pełne przeliczenie FINAL.")
        all_final_data_full, current_hash = reprocess_to_final(json_path, all_raw_data, max_workers=max_workers)
        all_final_data = {oid: all_final_data_full[oid] for oid in target_ids if oid in all_final_data_full}
        with open(FINAL_DATA_PKL, "wb") as f:
            pickle.dump(all_final_data_full, f)
        with open(FINAL_HASH_FILE, "w") as f:
            f.write(current_hash)
        if VERBOSE:
            print(f"[SAVE] Final data zapisane do: {FINAL_DATA_PKL}")
            print(f"[SAVE] Hash zapisany do: {FINAL_HASH_FILE}")
    else:
        if VERBOSE:
            print("[CACHE] Parametry się nie zmieniły, sprawdzam brakujące obiekty w FINAL cache.")
        with open(FINAL_DATA_PKL, "rb") as f:
            final_cache = pickle.load(f)
        all_final_data, current_hash = reprocess_missing_final(
            json_path, all_raw_data, final_cache, target_ids, max_workers=max_workers
        )
        merged_cache = dict(final_cache)
        merged_cache.update({oid: all_final_data[oid] for oid in all_final_data})
        with open(FINAL_DATA_PKL, "wb") as f:
            pickle.dump(merged_cache, f)
        with open(FINAL_HASH_FILE, "w") as f:
            f.write(current_hash)
        if VERBOSE:
            print(f"[CACHE] ✓ Używam FINAL cache (po przyrostowym uzupełnieniu).")

    engine_state = {
        "last_processed": datetime.now().isoformat(),
        "raw_data_objects": len(all_raw_data),
        "final_data_objects": len(all_final_data),
        "lat": lat,
        "lon": lon,
        "year": year,
        "minalt": obj_min_alt_deg,
        "sunlimit": sun_alt_limit_deg,
        "city_name": city_name,
    }
    save_engine_state(engine_state)

    if VERBOSE:
        print("=" * 119)
        print("[INFO] Silnik zakończył pracę.")
        print(f"[INFO] Raw objects: {len(all_raw_data)}. Final objects: {len(all_final_data)}.")
        print("=" * 119)

    return all_final_data

if __name__ == "__main__":
    import sys

    force_all = "--force-all" in sys.argv
    run_engine_from_vis_json(
        "vis_data.json",
        max_workers=None,
        force_all=force_all,
    )
