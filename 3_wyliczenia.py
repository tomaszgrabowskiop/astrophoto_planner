import json
import pickle
import hashlib
import os
from datetime import datetime, timedelta
from concurrent.futures import ProcessPoolExecutor, as_completed
from functools import partial
from dataclasses import dataclass, asdict
from typing import Dict, List, Any, Tuple

import numpy as np
import pandas as pd
from astropy.coordinates import SkyCoord, EarthLocation, get_sun, get_body, AltAz
from astropy.time import Time
import astropy.units as u
from tqdm import tqdm

# =====================================================
# KONFIGURACJA
# =====================================================

H_START = 14.0
H_END = 31.0
H_RANGE = H_END - H_START
N_SAMPLES = int(H_RANGE * 60 / 5)

CROSSING_SAMPLES = 1440

# Pliki cache
RAW_DATA_PKL = "observing_data_raw.pkl"
FINAL_DATA_PKL = "observing_data.pkl"
FINAL_HASH_FILE = "observing_data_final.hash"
CONFIG_FILE = "config_3.json"

# =====================================================
# DATACLASS DO PRZECHOWYWANIA SUROWYCH DANYCH
# =====================================================

@dataclass
class RawObjectData:
    """Dane niezależne od parametrów filtrowania (wysokości, crossingi)."""
    obj_id: str
    o_alt_all: np.ndarray  # shape: [ndays, N_SAMPLES] - wysokości obiektu
    sun_pts_list: List[List[datetime]]  # [ndays] - momenty przejścia słońca (non ISO format)
    obj_pts_list: List[List[datetime]] # [ndays] - momenty przejścia obiektu (nonISO format)
    

@dataclass
class Config:
    """Konfiguracja dla 3_wyliczenia."""
    raw_pkl_path: str = RAW_DATA_PKL
    final_pkl_path: str = FINAL_DATA_PKL
    final_hash_file: str = FINAL_HASH_FILE
    auto_reprocess: bool = True
    verbose: bool = True


# =====================================================
# FUNKCJE POMOCNICZE
# =====================================================

def load_vis_data(path: str = "vis_data.json") -> Dict[str, Any]:
    """Ładuje dane z vis_data.json."""
    with open(path, "r", encoding="utf-8") as f:
        return json.load(f)


def save_vis_data(data: Dict[str, Any], path: str = "vis_data.json"):
    """Zapisuje dane do vis_data.json."""
    with open(path, "w", encoding="utf-8") as f:
        json.dump(data, f, indent=2, ensure_ascii=False)


def get_params_hash(obj_min_alt_deg: float, sun_alt_limit_deg: float) -> str:
    """
    Hash tylko parametrów które wpływają na maskowanie w 3_.
    (min_alt i sun_alt_limit)
    """
    s = f"{obj_min_alt_deg:.6f}|{sun_alt_limit_deg:.6f}"
    return hashlib.md5(s.encode()).hexdigest()


def load_config(path: str = CONFIG_FILE) -> Config:
    """Wczytuje konfigurację z pliku, jeśli istnieje."""
    if os.path.exists(path):
        try:
            with open(path, "r") as f:
                cfg_dict = json.load(f)
            return Config(**cfg_dict)
        except Exception as e:
            print(f"[WARNING] Błąd wczytywania {path}: {e}, używam domyślnej konfiguracji.")
    return Config()


def save_config(cfg: Config, path: str = CONFIG_FILE):
    """Zapisuje konfigurację do pliku."""
    cfg_dict = asdict(cfg)
    with open(path, "w") as f:
        json.dump(cfg_dict, f, indent=2)


def get_crossings(target, height, t_noon, location) -> List[Time]:
    """Szuka precyzyjnych momentów przejścia przez wysokość 'height'."""
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


def mask_to_segments(mask, h_start=H_START, h_end=H_END) -> List[Tuple[float, float]]:
    """Zamienia maskę bool na listę ciągłych odcinków [start_rel, end_rel]."""
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
) -> RawObjectData:
    """
    Oblicza TYLKO dane niezależne od parametrów filtrowania.
    (Wysokości, punkty przejścia)
    
    Wynik: RawObjectData
    """
    obj_id, ra, dec = obj_data
    coord = SkyCoord(ra * u.deg, dec * u.deg)
    ndays = len(days)

    # --- 1. Transformacja AltAz dla obiektu (wszystkie dni/czasy) ---
    t_grid = days_noon[:, None] + t_night_offsets_hours[None, :] * u.hour
    frame = AltAz(obstime=t_grid.reshape(-1), location=location)
    o_alt_all = coord.transform_to(frame).alt.deg.reshape(ndays, N_SAMPLES)

    # --- 2. Punkty przejścia (per dzień) ---
    sun_pts_list = []
    obj_pts_list = []
    
    for day_idx, d in enumerate(days):
        t_noon = days_noon[day_idx]
        
        # Crossingi na 0° (uniwersalny punkt odniesienia)
        sun_pts = get_crossings(get_sun(t_noon), 0, t_noon, location)
        obj_pts = get_crossings(coord, 0, t_noon, location)
        
        # Konwersja do ISO format (do serializacji)
        sun_pts_list.append([p.datetime for p in sun_pts])
        obj_pts_list.append([p.datetime for p in obj_pts])

    return RawObjectData(
        obj_id=obj_id,
        o_alt_all=o_alt_all,
        sun_pts_list=sun_pts_list,
        obj_pts_list=obj_pts_list,
    )


def compute_raw_data(
    json_path: str,
    cfg: Config,
    max_workers: int = None,
) -> Dict[str, RawObjectData]:
    """
    Oblicza lub wczytuje dane surowe dla wszystkich obiektów.
    
    Returns: Dict[obj_id, RawObjectData]
    """
    vis = load_vis_data(json_path)
    
    lat = vis["location"]["lat"]
    lon = vis["location"]["lon"]
    location = EarthLocation(lat=lat * u.deg, lon=lon * u.deg)
    year = vis["year"]
    
    objects = vis["objects"]
    objects_data = [(obj["id"], obj["ra"], obj["dec"]) for obj in objects]
    
    days = pd.date_range(f"{year}-01-01", f"{year}-12-31", freq="D")
    ndays = len(days)
    
    # --- Sprawdzenie cache ---
    if os.path.exists(cfg.raw_pkl_path):
        if cfg.verbose:
            print(f"[CACHE] Wczytuję dane surowe z {cfg.raw_pkl_path}.")
        try:
            with open(cfg.raw_pkl_path, "rb") as f:
                all_raw_data = pickle.load(f)
            
            if cfg.verbose:
                print(f"[CACHE] ✓ Załadowano {len(all_raw_data)} obiektów z cache.")
            return all_raw_data
        except Exception as e:
            if cfg.verbose:
                print(f"[ERROR] Błąd wczytywania raw cache: {e}, obliczam od nowa...")
    
    # --- Obliczenia raw data ---
    if cfg.verbose:
        print(f"[COMPUTE] Obliczam dane surowe dla {len(objects_data)} obiektów podczas {ndays} dni.")
    
    days_noon = Time([
        datetime.combine(d.date(), datetime.min.time()) + timedelta(hours=12)
        for d in days
    ])
    
    t_night_offsets_hours = np.linspace(2, 19, N_SAMPLES)
    
    all_raw_data = {}
    
    process_raw_func = partial(
        process_single_object_raw,
        days=days,
        location=location,
        days_noon=days_noon,
        t_night_offsets_hours=t_night_offsets_hours,
    )
    
    with ProcessPoolExecutor(max_workers=max_workers) as executor:
        future_to_objid = {
            executor.submit(process_raw_func, obj_data): obj_data[0]
            for obj_data in objects_data
        }
        
        with tqdm(total=len(objects_data), desc="Raw compute", unit="obj", ncols=119) as pbar:
            for future in as_completed(future_to_objid):
                obj_id = future_to_objid[future]
                try:
                    raw_result = future.result()
                    all_raw_data[raw_result.obj_id] = raw_result
                    pbar.set_postfix_str(f"ostatni={obj_id}")
                except Exception as e:
                    print(f"[ERROR] Błąd dla {obj_id}: {e}")
                pbar.update(1)
    
    # --- Zapisz raw data ---
    with open(cfg.raw_pkl_path, "wb") as f:
        pickle.dump(all_raw_data, f)
    
    if cfg.verbose:
        print(f"[SAVE] Raw data zapisane do {cfg.raw_pkl_path} ({len(all_raw_data)} obiektów).")
    
    return all_raw_data


# =====================================================
# ETAP 2: KONWERSJA RAW -> FINAL (SZYBKA REPROCESSACJA)
# =====================================================

def process_raw_to_final(
    raw_data: RawObjectData,
    days: pd.DatetimeIndex,
    sun_alt_all: np.ndarray,
    moon_alt_all: np.ndarray,
    sun_alt_limit_deg: float,
    obj_min_alt_deg: float,
) -> Tuple[str, List[Dict[str, Any]]]:
    """
    Konwertuje surowe dane do formatu finalnego z nowymi parametrami.
    Teraz oblicza również poprawne czasy zmierzchu/świtu dla zadanego sun_alt_limit_deg.
    """
    obj_id = raw_data.obj_id
    o_alt_all = raw_data.o_alt_all
    ndays = len(days)
    results = []

    # Parametry siatki czasowej (muszą być zgodne z N_SAMPLES i H_RANGE)
    # H_START = 14.0, H_END = 31.0
    t_step = H_RANGE / (N_SAMPLES - 1)

    for day_idx, d in enumerate(days):
        d_date = d.date() if hasattr(d, "date") else d
        d_midnight = datetime.combine(d_date, datetime.min.time())
        
        o_alt = o_alt_all[day_idx]
        s_alt = sun_alt_all[day_idx]
        m_alt = moon_alt_all[day_idx]

        # --- 1. MASKOWANIE (Jakość) ---
        quality_mask = (o_alt > obj_min_alt_deg) & (s_alt < sun_alt_limit_deg)
        
        transit_rel = np.argmax(o_alt) * (H_RANGE / N_SAMPLES)
        q_hours = np.sum(quality_mask) * (H_RANGE / N_SAMPLES)
        m_hours = np.sum(quality_mask & (m_alt < 0.0)) * (H_RANGE / N_SAMPLES)
        qual_segments = mask_to_segments(quality_mask)

        # --- 2. OBLICZANIE PUNKTÓW ZMIERZCHU (FIX dla 6_drukuj...) ---
        # Zamiast przekazywać raw_data.sun_pts_list (które są dla 0 stopni),
        # interpolujemy czas przejścia przez sun_alt_limit_deg na podstawie s_alt.
        
        new_sun_pts = []
        
        # Wykrywamy momenty przecięcia granicy sun_alt_limit_deg
        # diff zmienia znak w punkcie przecięcia
        diff = s_alt - sun_alt_limit_deg
        crossings = np.where(np.diff(np.signbit(diff)))[0]
        
        if len(crossings) > 0:
            for idx in crossings:
                # Interpolacja liniowa dla dokładniejszego czasu
                y0, y1 = diff[idx], diff[idx+1]
                # x to frakcja między indeksem idx a idx+1, gdzie y=0
                # y = y0 + (y1 - y0) * x  =>  0 = y0 + (y1 - y0) * x  => x = -y0 / (y1 - y0)
                fraction = -y0 / (y1 - y0)
                
                # Czas w godzinach od północy (startujemy od H_START=14.0)
                h_val = H_START + (idx + fraction) * t_step
                
                # Dodajemy do północy
                pt_time = d_midnight + timedelta(hours=h_val)
                new_sun_pts.append(pt_time)
        else:
            # BRAK PRZEJŚĆ (np. białe noce lub noc polarna)
            mean_alt = np.mean(s_alt)
            if mean_alt > sun_alt_limit_deg:
                # Słońce zawsze POWYŻEJ limitu (białe noce, lato) -> Brak nocy
                # Skrypt 6_ traktuje pustą listę jako "cała noc", więc musimy podać
                # punkt zerowy (start = koniec), żeby narysował brak nocy.
                dummy_t = d_midnight + timedelta(hours=H_START)
                new_sun_pts = [dummy_t, dummy_t]
            else:
                # Słońce zawsze PONIŻEJ limitu (zima stulecia, noc polarna) -> Cała noc
                # Pusta lista w skrypcie 6_ powoduje wypełnienie całego wykresu kolorem nocy.
                new_sun_pts = []

        results.append({
            "day": d_date,
            "sun_pts": new_sun_pts,  # <--- Podmieniamy na poprawne punkty
            "obj_pts": raw_data.obj_pts_list[day_idx],
            "transit_rel": transit_rel,
            "q_hours": q_hours,
            "m_hours": m_hours,
            "qual_segments": qual_segments,
        })

    return obj_id, results


def reprocess_to_final(
    json_path: str,
    all_raw_data: Dict[str, RawObjectData],
    cfg: Config,
) -> Tuple[Dict[str, List[Dict]], str]:
    """
    Szybka konwersja raw -> final z nowymi parametrami.
    
    Returns: (all_final_data, current_hash)
    """
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
    
    current_params_hash = get_params_hash(obj_min_alt_deg, sun_alt_limit_deg)
    
    if cfg.verbose:
        print(f"[REPROCESS] Konwersja {len(all_raw_data)} obiektów do ostatecznego formatu.")
        print(f"[PARAMS] Minimalna wysokość obiektu: {obj_min_alt_deg}°, wysokość słońca pod horyzontem: {sun_alt_limit_deg}°.")
        print(f"[HASH] {current_params_hash}.")
    
    # --- Precompute Sun/Moon altitudes dla całego roku ---
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
    
    # --- Konwersja per obiekt ---
    all_final_data = {}
    
    for obj_id, raw_data in tqdm(
        all_raw_data.items(),
        desc="Reprocess",
        unit="obj",
        ncols=119,
        total=len(all_raw_data)
    ):
        obj_id_result, results = process_raw_to_final(
            raw_data, days, sun_alt_all, moon_alt_all,
            sun_alt_limit_deg, obj_min_alt_deg
        )
        all_final_data[obj_id_result] = results
    
    if cfg.verbose:
        print(f"[REPROCESS] ✓ Przetworzono {len(all_final_data)} obiektów.")
    
    return all_final_data, current_params_hash


def should_reprocess(json_path: str, cfg: Config) -> bool:
    """
    Sprawdza czy warto reprocessować (czy parametry się zmieniły).
    
    Returns: True jeśli trzeba reprocessować
    """
    if not os.path.exists(cfg.final_hash_file):
        return True  # Brak hash pliku - zawsze reprocess
    
    vis = load_vis_data(json_path)
    params = vis.get("parameters", {})
    obj_min_alt_deg = params.get("minalt", 20.0)
    sun_alt_limit_deg = params.get("sunlimit", -6.0)
    
    current_hash = get_params_hash(obj_min_alt_deg, sun_alt_limit_deg)
    
    try:
        with open(cfg.final_hash_file, "r") as f:
            old_hash = f.read().strip()
        
        if current_hash == old_hash:
            return False  # Hash się zgadza - skip reprocess
        else:
            return True  # Hash się zmienił - reprocess
    except Exception:
        return True


# =====================================================
# GŁÓWNA FUNKCJA
# =====================================================

def run_engine_from_vis_json(
    json_path: str = "vis_data.json",
    max_workers: int = None,
    force_raw: bool = False,
    force_final: bool = False,
    cfg: Config = None,
) -> Dict[str, List[Dict]]:
    """
    Inteligentny system cachowania:
    
    1. Oblicza/wczytuje raw data (wysokości, crossingi) - POWOLNE
    2. Szybko reprocessuje final data gdy parametry się zmienią - SZYBKIE
    
    Args:
        json_path: ścieżka do vis_data.json
        max_workers: liczba CPU cores dla parallelizacji
        force_raw: wymuś recalc raw data (ignoruj cache)
        force_final: wymuś reprocess final data (ignoruj hash)
        cfg: obiekt Config (jeśli None, wczytaj z pliku)
    
    Returns:
        Dict[obj_id, list_of_results_per_day]
    """
    if cfg is None:
        cfg = load_config()
    
    vis = load_vis_data(json_path)
    
    city_name = vis["location"].get("name", "Unknown")
    year = vis["year"]
    lat = vis["location"]["lat"]
    lon = vis["location"]["lon"]
    
    params = vis.get("parameters", {})
    obj_min_alt_deg = params.get("minalt", 20.0)
    sun_alt_limit_deg = params.get("sunlimit", -6.0)
    
    if cfg.verbose:
        print("=" * 119)
        print(f"ENGINE: Smart Cache System")
        print("=" * 119)
        print(f"[INFO] Lokalizacja: {city_name} ({lat}°, {lon}°)")
        print(f"[INFO] Rok: {year}")
        print(f"[INFO] Parametry: minimalna wysokość obiektu: {obj_min_alt_deg}°, wysokość słońca pod horyzontem: {sun_alt_limit_deg}°")
        print(f"[INFO] Cache: plik z danymi raw - {cfg.raw_pkl_path}, plik z maskami wyboru - {cfg.final_pkl_path}")
        print("=" * 119)
    
    # ========== ETAP 1: Raw Data ==========
    if force_raw and os.path.exists(cfg.raw_pkl_path):
        if cfg.verbose:
            print("[FORCE] Usuwam stary cache raw data.")
        os.remove(cfg.raw_pkl_path)
    
    all_raw_data = compute_raw_data(json_path, cfg, max_workers=max_workers)
    
    # ========== ETAP 2: Final Data (z cachingiem parametrów) ==========
    needs_reprocess = force_final or should_reprocess(json_path, cfg)
    
    if needs_reprocess:
        if cfg.verbose:
            print("[REPROCESS] Parametry się zmieniły, ponowne przeliczenie.")
        
        all_final_data, current_hash = reprocess_to_final(json_path, all_raw_data, cfg)
        
        # Zapisz final data + hash
        with open(cfg.final_pkl_path, "wb") as f:
            pickle.dump(all_final_data, f)
        
        with open(cfg.final_hash_file, "w") as f:
            f.write(current_hash)
        
        if cfg.verbose:
            print(f"[SAVE] Final data zapisane do: {cfg.final_pkl_path}")
            print(f"[SAVE] Hash zapisany do: {cfg.final_hash_file}")
    else:
        if cfg.verbose:
            print(f"[CACHE] Parametry się nie zmieniły, wczytuję final cache.")
        
        try:
            with open(cfg.final_pkl_path, "rb") as f:
                all_final_data = pickle.load(f)
            
            if cfg.verbose:
                print(f"[CACHE] ✓ Załadowano {len(all_final_data)} obiektów z cache.")
        except Exception as e:
            if cfg.verbose:
                print(f"[ERROR] Błąd wczytywania final cache: {e}, przeliczam ponownie.")
            
            all_final_data, current_hash = reprocess_to_final(json_path, all_raw_data, cfg)
            
            with open(cfg.final_pkl_path, "wb") as f:
                pickle.dump(all_final_data, f)
            
            with open(cfg.final_hash_file, "w") as f:
                f.write(current_hash)
    
    # ========== ETAP 3: Zapis parametrów do vis_data.json ==========
    # (dla referencji)
    vis["_engine_state"] = {
        "last_processed": datetime.now().isoformat(),
        "raw_data_objects": len(all_raw_data),
        "final_data_objects": len(all_final_data),
    }
    save_vis_data(vis, json_path)
    
    if cfg.verbose:
        print("=" * 119)
        print(f"[SUCCESS] Silnik zakończył pracę.")
        print(f"[STATS] Raw objects: {len(all_raw_data)}. Final objects: {len(all_final_data)}.")
        print("=" * 119)
    
    return all_final_data


if __name__ == "__main__":
    import sys
    
    # Opcje:
    # python 3_wyliczenia.py                  <- auto (smart cache)
    # python 3_wyliczenia.py --force-raw      <- wymuś recalc raw
    # python 3_wyliczenia.py --force-final    <- wymuś recalc final
    # python 3_wyliczenia.py --force-all      <- wszystko od nowa
    # python 3_wyliczenia.py --verboseoff        <- verbose output is off
    
    force_raw = "--force-raw" in sys.argv or "--force-all" in sys.argv
    force_final = "--force-final" in sys.argv or "--force-all" in sys.argv
    verbose = "--verboseoff" not in sys.argv  # domyślnie verbose
    
    cfg = load_config()
    cfg.verbose = verbose
    
    run_engine_from_vis_json(
        "vis_data.json",
        max_workers=None,
        force_raw=force_raw,
        force_final=force_final,
        cfg=cfg,
    )