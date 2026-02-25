#!/usr/bin/env python3
# =============================================================================
# 3_compute.py  –  Astrophotography Planner  |  Etap 3: Obliczenia astronomiczne
# =============================================================================
"""
Krok 3: Silnik obliczeń astronomicznych (Compute Engine).

Najbardziej czasochłonny etap procesu generowania planera. Wykorzystuje bibliotekę AstroPy
do precyzyjnych obliczeń mechaniki nieba dla Twojej konkretnej lokalizacji.

Główne obliczenia:
1. Wyznaczenie widoczności: Dla każdego obiektu liczona jest jego wysokość nad horyzontem
   dla każdej nocy w roku (z dokładnością minutową lub zgodną z siatką).
2. Analiza okien obserwacyjnych:
   - Uwzględnia zmierzch (Astronomiczny/Nautyczny/Cywilny).
   - Uwzględnia fazę i pozycję Księżyca (wyklucza czas, gdy Księżyc przeszkadza).
   - Uwzględnia limity sprzętowe (minimalna wysokość nad horyzontem).
3. Obliczenie sumarycznych "Jakościowych Godzin" (Imaging Hours) dla każdego obiektu w skali roku.

Optymalizacja:
- Skrypt wykorzystuje Multiprocessing (wszystkie rdzenie CPU).
- Wyniki są cache'owane w pliku 'observing_data.pkl'. Ponowne uruchomienie bez zmiany lokalizacji/roku  jest znacznie szybsze.

Wyjście:
- Aktualizacja cache i danych w pamięci dla kolejnych kroków.
"""
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
from pathlib import Path

from shared import (PATHS, 
                    fmt, print_step, print_green,
                    load_vis_data, save_vis_data, 
                    UserConfig,  
                    H_START, H_END, H_RANGE, N_SAMPLES, CROSSING_SAMPLES,
                    H_NOON_TO_START, H_NOON_TO_END,
                    get_engine_raw_hash, get_engine_final_hash)

# =====================================================
# KONFIGURACJA
# =====================================================

VERBOSE = True  # Jeśli ustawione na True, program "mówi" nam co aktualnie robi (printuje logi)

# =====================================================
# DATACLASS DO PRZECHOWYWANIA SUROWYCH DANYCH
# =====================================================

@dataclass
class RawObjectData:
    """Dane niezależne od parametrów filtrowania (wysokości, crossingi).
       Jest to nasze "surowe pudełko", do którego pakujemy informacje
       przed obrabianiem z ograniczeniami np. zanieczyszczenia światłem.
    """
    obj_id: str
    o_alt_all: np.ndarray  # Tablica wysokości obiektu [ndays, N_SAMPLES]
    sun_pts_list: List[List[datetime]]  # Kiedy słońce przecina horyzont [ndays]
    obj_pts_list: List[List[datetime]]  # Kiedy obiekt przecina horyzont [ndays]

# =====================================================
# FUNKCJE POMOCNICZE
# =====================================================

def get_crossings(target, height, t_noon, location) -> List[Time]:
    """
    Funkcja, która znajduje dokładny moment przecięcia przez obiekt (target)
    konkretnej wysokości (height) na niebie.
    """
    # Sprawdzamy dobę w CROSSING_SAMPLES równych odstępach
    times = t_noon + np.linspace(0, 24, CROSSING_SAMPLES) * u.hour
    altaz = target.transform_to(AltAz(obstime=times, location=location))
    alts = altaz.alt.deg
    
    # diff to różnica między pozycją na niebie a zadaną linią (np. linią horyzontu)
    diff = alts - height
    
    # Tu sprawdzamy, gdzie diff zmieniło znak - czyli obiekt przeszedł przez zadaną wysokość
    idx = np.where(np.diff(np.sign(diff)))[0]

    res = []
    # Dla każdego znalezionego przecięcia robimy dokładniejsze przybliżenie
    for i in idx:
        t_exact = times[i] + (times[i + 1] - times[i]) * (-diff[i] / (diff[i + 1] - diff[i]))
        res.append(t_exact)

    return sorted(res)

def load_engine_state() -> dict:
    """
    Wczytuje stan poprzedniego wywołania algorytmu, żeby wiedzieć,
    czy parametry się zmieniły i czy cache wymaga zresetowania (invalidacji).
    """
    if os.path.exists(PATHS.engine_state):
        try:
            with open(PATHS.engine_state, "rb") as f:
                return pickle.load(f)
        except Exception:
            return {}
    return {}

def save_engine_state(state: dict):
    """
    Zapisuje obecne parametry (jak lokalizacja i rok) jako stan silnika na przyszłość.
    """
    with open(PATHS.engine_state, "wb") as f:
        pickle.dump(state, f)

def mask_to_segments(mask, h_start=H_START, h_end=H_END) -> List[Tuple[float, float]]:
    """
    Zamienia tablicę "mask" (składającą się z Prawda/Fałsz) na odcinki (początek, koniec).
    Używamy tego do określenia dokładnych ram czasowych obserwacji w nocy.
    """
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
    if mask.size < 2:
        return segments
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
    precomputed_sun_pts: List[List[datetime]],  # Z góry obliczone punkty słońca, dla szybszego działania
) -> RawObjectData:
    """
    Liczy "surowe" (bazowe) wznoszenie się wybranego obiektu nad horyzont w każdym dniu roku.
    """
    obj_id, ra, dec = obj_data
    coord = SkyCoord(ra * u.deg, dec * u.deg)
    ndays = len(days)

    # 1. Obliczanie wysokości w siatce godzinowej w sposób "wektoryzowany" 
    t_grid = days_noon[:, None] + t_night_offsets_hours[None, :] * u.hour
    frame = AltAz(obstime=t_grid.reshape(-1), location=location)
    o_alt_all = coord.transform_to(frame).alt.deg.reshape(ndays, N_SAMPLES)

    obj_pts_list = []

    # 2. Pętla po dniach - tutaj szukamy kiedy dokładnie obiekt wschodzi i zachodzi (crossing horyzontu 0)
    for day_idx, t_noon in enumerate(days_noon):
        obj_pts = get_crossings(coord, 0, t_noon, location)
        obj_pts_list.append([p.datetime for p in obj_pts])

    # Zwracamy gotowe "pudełko" z zapisanymi wyliczeniami
    return RawObjectData(
        obj_id=obj_id,
        o_alt_all=o_alt_all,
        sun_pts_list=precomputed_sun_pts,  # Zwracamy gotowe dane słońca
        obj_pts_list=obj_pts_list,
    )

def compute_raw_data(
    json_path: Path | str,
    object_limit: int,
    max_workers: int = None,
) -> Dict[str, RawObjectData]:
    """
    Oblicza LUB UZUPEŁNIA dane surowe dla obiektów z vis_data.json
    ograniczonych do object_limit.
    """
    # 1. Wczytujemy plik z danymi i ustawieniami użytkownika (z shared.py)
    vis = load_vis_data(json_path)
    cfg = UserConfig.from_vis_data(vis)
    
    # Tworzymy obiekt reprezentujący nasze miejsce obserwacji na Ziemi.
    # Musimy użyć konkretnych nazw argumentów 'lat' (szerokość) i 'lon' (długość).
    location = EarthLocation(
        lat=cfg.location.lat * u.deg, 
        lon=cfg.location.lon * u.deg
    )
    
    # Wyciągamy z pliku tylko tyle obiektów, na ile pozwala 'object_limit'
    objects = vis["objects"][:object_limit]
    objects_data = [(obj["id"], obj["ra"], obj["dec"]) for obj in objects]

    # Generujemy listę wszystkich dni w wybranym roku
    days = pd.date_range(f"{cfg.year}-01-01", f"{cfg.year}-12-31", freq="D")
    ndays = len(days)

    # === SMART CACHE DLA RAW ===
    # Obliczamy hash dla obecnych ustawień RAW.
    # Hash zależy tylko od miejsca na Ziemi i roku. Jeśli te rzeczy się nie zmienią, 
    raw_hash = get_engine_raw_hash(cfg.location.lat, cfg.location.lon, cfg.year)
    raw_hash_file = PATHS.observing_raw_hash

    # === INICJALIZACJA CACHE (PAMIĘCI PODRĘCZNEJ) ===
    all_raw_data: Dict[str, RawObjectData] = {}
    cache_valid = False
    
    # Sprawdzamy, czy mamy już zapisane jakieś wcześniejsze obliczenia (plik raw) oraz ich hash
    if os.path.exists(PATHS.observing_raw) and os.path.exists(raw_hash_file):
        try:
            with open(raw_hash_file, "r") as f:
                saved_hash = f.read().strip()
            
            if raw_hash == saved_hash:
                with open(PATHS.observing_raw, "rb") as f:
                    all_raw_data = pickle.load(f)
                cache_valid = True
                if VERBOSE:
                    print(f"[CACHE] Załadowano {len(all_raw_data)} obiektów z RAW cache (RAW hash: {raw_hash[:8]}...).")
            else:
                if VERBOSE:
                    print(f"[CACHE] Invalidated! saved_hash={saved_hash[:8]} ≠ current={raw_hash[:8]}")
                os.remove(PATHS.observing_raw)
                os.remove(raw_hash_file)
        except Exception as e:
            if VERBOSE:
                print(f"[ERROR] Błąd walidacji RAW hash/cache: {e}")
    else:
        if VERBOSE:
            print("[CACHE] Brak hash lub RAW cache - obliczam od nowa.")

    # Jeśli dane w cache były niepoprawne lub ich nie było, upewniamy się, że słownik jest pusty
    if not cache_valid:
        all_raw_data = {}  # ← KLUCZOWE! Pusty cache, zaczynamy z czystą kartą
        if VERBOSE:
            print(f"[RAW]   Brak ważnego RAW cache - przygotowuję przeliczenie wszystkich obiektów.")

    # Obliczamy, których obiektów nam jeszcze brakuje
    existing_ids = set(all_raw_data.keys())
    missing_data = [
        (obj_id, ra, dec)
        for (obj_id, ra, dec) in objects_data
        if obj_id not in existing_ids
    ]

    if VERBOSE:
        if missing_data:
            print(f"[RAW]   Brakujących obiektów w RAW cache: {len(missing_data)}.")
        else:
            print("[RAW]   RAW cache zawiera wszystkie obiekty z aktualnego limitu.")

    # Jeśli nie ma niczego do policzenia (wszystko było w cache), to po prostu zapisujemy hash i kończymy
    if not missing_data:
        # ZAPIS HASH (nawet jeśli nie przeliczaliśmy)
        with open(raw_hash_file, "w") as f:
            f.write(raw_hash)
        return all_raw_data

    # Przygotowujemy "środek dnia" (południe) dla każdego dnia roku jako punkt odniesienia
    days_noon = Time([
        datetime.combine(d.date(), datetime.min.time()) + timedelta(hours=12)
        for d in days
    ])

    # Kiedy prowadzimy obserwacje w nocy (odstępy od południa)
    t_night_offsets_hours = np.linspace(H_NOON_TO_START, H_NOON_TO_END, N_SAMPLES)

    # === OPTYMALIZACJA: WYLICZAMY SŁOŃCE TYLKO RAZ ===
    # Liczymy dane słoneczne RAZ na początku, żeby nie obciążać komputera powtarzaniem tego 
    # dla każdego obiektu
    if VERBOSE:
        print(f"[COMPUTE] Obliczam dane Słońca dla {ndays} dni.")
    
    precomputed_sun_pts = []
    for t_noon in tqdm(days_noon, desc="          Sun Calc", ncols=119, colour="green"):
        # Parametr "0" oznacza wysokość horyzontu. Szukamy kiedy Słońce go przecina.
        s_pts = get_crossings(get_sun(t_noon), 0, t_noon, location)
        precomputed_sun_pts.append([p.datetime for p in s_pts])

    if VERBOSE:
        print(f"[COMPUTE] Obliczam RAW dla {len(missing_data)} brakujących obiektów...")

    # "Pakujemy" funkcję process_single_object_raw ze wszystkimi stałymi danymi, 
    # żeby łatwo było ją wrzucić do wielu rdzeni procesora naraz.
    process_raw_func = partial(
        process_single_object_raw,
        days=days,
        location=location,
        days_noon=days_noon,
        t_night_offsets_hours=t_night_offsets_hours,
        precomputed_sun_pts=precomputed_sun_pts,
    )

    # Rozdzielamy pracę na wiele procesów
    with ProcessPoolExecutor(max_workers=max_workers) as executor:
        # Zlecamy obliczenia dla wszystkich brakujących obiektów
        future_to_objid = {
            executor.submit(process_raw_func, obj_data): obj_data[0]
            for obj_data in missing_data
        }
        # Pasek postępu
        with tqdm(
            total=len(missing_data),
            desc="          Raw compute",
            unit="obj",
            ncols=119,
            colour = "green",
        ) as pbar:
            for future in as_completed(future_to_objid):
                obj_id = future_to_objid[future]
                try:
                    # Wyciągamy wynik pojedynczego obiektu i dodajemy go do "głównego"
                    raw_result = future.result()
                    all_raw_data[raw_result.obj_id] = raw_result
                    pbar.set_postfix_str(f"ostatni={obj_id}")
                except Exception as e:
                    print(f"[ERROR] Błąd dla {obj_id}: {e}")
                pbar.update(1)

    # === ZAPISUJEMY WYNIKI NA DYSK  ===
    with open(PATHS.observing_raw, "wb") as f:
        pickle.dump(all_raw_data, f)
    
    with open(raw_hash_file, "w") as f:
        f.write(raw_hash)
    
    if VERBOSE:
        print(f"[SAVE]  Raw data zapisane do {PATHS.observing_raw} ({len(all_raw_data)} obiektów).")
        print(f"[SAVE]  RAW hash zapisany: {raw_hash[:8]}...")

    # Na koniec  funkcja zwraca wypełniony surowe dane
    return all_raw_data

# =====================================================
# ETAP 2: KONWERSJA RAW -> FINAL
# =====================================================

def process_raw_to_final(
    raw_data: RawObjectData,
    days: pd.DatetimeIndex,
    sun_alt_all: np.ndarray,
    moon_alt_all: np.ndarray,
    sun_limit: float,        # cfg.sun_limit z UserConfig
    min_altitude: float,     # cfg.min_altitude z UserConfig
    vis: dict,
) -> Tuple[str, List[Dict[str, Any]]]:
    """
    Przetwarza surowe dane (RAW) jednego obiektu na ostateczny format (FINAL),
    uwzględniając limity wysokości obiektu i wysokości Słońca.

    Wejście:
    - raw_data: RawObjectData (wysokości obiektu i punkty przecięcia horyzontu)
    - days: index dni roku (pd.DatetimeIndex)
    - sun_alt_all: wysokość Słońca [ndays, N_SAMPLES]
    - moon_alt_all: wysokość Księżyca [ndays, N_SAMPLES]
    - sun_limit: próg wysokości Słońca (np. -12.0 dla zmierzchu żeglarskiego)
    - min_altitude: minimalna wysokość obiektu nad horyzontem (deg)
    - vis: pełna struktura vis_data.json (dla strefy czasowej)

    Wyjście:
    - (obj_id, results), gdzie results to lista słowników per dzień:
      {
          "day": d_date,
          "sun_pts": [...],
          "obj_pts": [...],
          "transit_rel": ...,
          "q_hours": ...,
          "m_hours": ...,
          "qual_segments": [...],
          "tz_offset": ...,
      }
    """
    obj_id = raw_data.obj_id
    o_alt_all = raw_data.o_alt_all
    ndays = len(days)
    results: List[Dict[str, Any]] = []

    # Krok czasowy w godzinach między kolejnymi próbkami w siatce
    t_step = H_RANGE / (N_SAMPLES - 1) 
    # "-1": żeby wizualnie mieściło się w obrębie nocy, nie jest to idealne piec minut, a 5:02, 
    # ale dzięki temu wchodzi w "klepsydrę" i ładnie się prezentuje na wykresie.

    # Strefa czasowa z vis_data.json
    tz_name = vis["location"]["tz"]
    tz = pytz.timezone(tz_name)

    for day_idx, d in enumerate(days):
        # Data kalendarzowa (bez czasu)
        d_date = d.date() if hasattr(d, "date") else d

        # Północ lokalna danego dnia, do której odnosimy godziny H_START..H_END
        d_midnight = datetime.combine(d_date, datetime.min.time())
        d_midnight_local = tz.localize(d_midnight)
        tz_offset = d_midnight_local.utcoffset().total_seconds() / 3600.0

        # Wysokości obiektu, Słońca i Księżyca dla tej nocy
        o_alt = o_alt_all[day_idx]
        s_alt = sun_alt_all[day_idx]
        m_alt = moon_alt_all[day_idx]

        # Maska jakości:
        # - obiekt powyżej min_altitude
        # - Słońce poniżej sun_limit (ciemno)
        quality_mask = (o_alt > min_altitude) & (s_alt < sun_limit)

        # Relatywny czas tranzytu (maksimum wysokości obiektu) w godzinach od H_START
        transit_rel = float(np.argmax(o_alt) * t_step)

        # Przedziały czasowe, w których quality_mask == True
        qual_segments = mask_to_segments(quality_mask)
        q_hours = float(sum(e - s for s, e in qual_segments))
        
        # Przedziały "moonless" (dobrze + Księżyc poniżej horyzontu)
        moonless_mask = quality_mask & (m_alt < 0.0)
        moonless_segments = mask_to_segments(moonless_mask)
        m_hours = float(sum(e - s for s, e in moonless_segments))

        # Wyliczenie punktów przecięcia Słońca z progiem sun_limit
        new_sun_pts: List[datetime] = []
        diff = s_alt - sun_limit
        crossings = np.where(np.diff(np.signbit(diff)))[0]

        if len(crossings) > 0:
            # Dla każdego przecięcia interpolujemy dokładny czas
            for idx in crossings:
                y0, y1 = diff[idx], diff[idx + 1]
                fraction = -y0 / (y1 - y0) if (y1 - y0) != 0 else 0.0
                h_val = H_START + (idx + fraction) * t_step
                pt_time = d_midnight + timedelta(hours=h_val)
                new_sun_pts.append(pt_time)
        else:
            # Brak przecięć – fallback zgodnie z oryginalną logiką
            mean_alt = np.mean(s_alt)
            if mean_alt > sun_limit:
                # Słońce cały czas powyżej progu – dwa identyczne punkty dummy
                dummy_t = d_midnight + timedelta(hours=H_START)
                new_sun_pts = [dummy_t, dummy_t]
            else:
                # Słońce cały czas poniżej progu – brak punktów
                new_sun_pts = []

        # Zapis wyników dla tej doby – dokładnie jak w starej wersji
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
    json_path: Path | str,
    all_raw_data: Dict[str, RawObjectData],
    max_workers: int = None,  # Dodano parametr
) -> Tuple[Dict[str, List[Dict]], str]:
    """
    Kiedy parametry zanieczyszczenia (zmierzch, wys. nad horyzontem) się zmienią,
    ta funkcja przelicza wszystkie surowe dane (RAW) na ostateczne (FINAL).
    """
    vis = load_vis_data(json_path)
    cfg = UserConfig.from_vis_data(vis)
    
    location = EarthLocation(lat=cfg.location.lat * u.deg, lon=cfg.location.lon * u.deg)

    days = pd.date_range(f"{cfg.year}-01-01", f"{cfg.year}-12-31", freq="D")
    ndays = len(days)

    final_hash = get_engine_final_hash(
        cfg.min_altitude, 
        cfg.sun_limit, 
        cfg.location.lat, 
        cfg.location.lon, 
        cfg.year
    )

    if VERBOSE:
        print(f"          Przeliczanie FINAL (Parallel): konwersja {len(all_raw_data)} obiektów.")

    days_noon = Time([
        datetime.combine(d.date(), datetime.min.time()) + timedelta(hours=12)
        for d in days
    ])

    t_night_offsets_hours = np.linspace(H_NOON_TO_START, H_NOON_TO_END, N_SAMPLES)
    t_grid_all = days_noon[:, None] + t_night_offsets_hours[None, :] * u.hour
    obstime_all = t_grid_all.reshape(-1)

    frame_all = AltAz(obstime=obstime_all, location=location)

    # Obliczenia wspólne (wykonywane raz w głównym procesie) dla Słońca i Księżyca
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
        sun_limit=cfg.sun_limit,       
        min_altitude=cfg.min_altitude, 
        vis=vis,
    )

    # Uruchomienie ProcessPoolExecutor
    with ProcessPoolExecutor(max_workers=max_workers) as executor:
        futures = {
            executor.submit(process_final_func, raw_data): raw_data.obj_id 
            for raw_data in all_raw_data.values()
        }

        for future in tqdm(
            as_completed(futures),
            desc="          Reprocess",
            unit="obj",
            ncols=119,
            colour="green",
            total=len(futures)
        ):
            try:
                obj_id_result, results = future.result()
                all_final_data[obj_id_result] = results
            except Exception as e:
                print(f"[ERROR] Błąd w FINAL dla obiektu: {e}")

    if VERBOSE:
        print(f"[COMPUTE] Przetworzono {len(all_final_data)} obiektów.")

    return all_final_data, final_hash

def should_reprocess(json_path: Path | str) -> bool:
    """
    Zwraca True, jeśli trzeba ponownie przeliczyć bazę FINAL 
    (bo np. zmieniono minimalną wysokość lub limit słońca w ustawieniach).
    """
    if not os.path.exists(PATHS.observing_hash):
        return True

    vis = load_vis_data(json_path)
    cfg = UserConfig.from_vis_data(vis)

    current_hash = get_engine_final_hash(
        cfg.min_altitude, 
        cfg.sun_limit, 
        cfg.location.lat, 
        cfg.location.lon, 
        cfg.year
    )

    try:
        with open(PATHS.observing_hash, "r") as f:
            old_hash = f.read().strip()
        # Jeśli nowy hash jest inny niż stary zapisany w pliku - musimy liczyć
        return current_hash != old_hash
    except Exception:
        return True

# =====================================================
# PRZYROSTOWY FINAL CACHE
# =====================================================

def get_target_ids_from_vis(vis: Dict[str, Any], object_limit: int) -> List[str]:
    """Pobiera listę ID z jsona aż do ustalonego limitu."""
    objects = vis.get("objects", [])
    ids = [o["id"] for o in objects]
    return ids[:object_limit]

def reprocess_missing_final(
    json_path: Path | str,
    all_raw_data: Dict[str, RawObjectData],
    existing_final: Dict[str, List[Dict]],
    target_ids: List[str],
    max_workers: int = None,
) -> Tuple[Dict[str, List[Dict]], str]:
    """
    Jeśli parametry (hash) się nie zmieniły, ale dodaliśmy np. nowe obiekty w vis_data.json,
    ta funkcja przeliczy tylko te brakujące.
    """
    vis = load_vis_data(json_path)
    cfg = UserConfig.from_vis_data(vis)
    
    location = EarthLocation(lat=cfg.location.lat * u.deg, lon=cfg.location.lon * u.deg)
    
    days = pd.date_range(f"{cfg.year}-01-01", f"{cfg.year}-12-31", freq="D")
    ndays = len(days)
    
    final_hash = get_engine_final_hash(
        cfg.min_altitude, 
        cfg.sun_limit, 
        cfg.location.lat, 
        cfg.location.lon, 
        cfg.year
    )

    # Obliczenia wspólne Słońca/Księżyca
    days_noon = Time([
        datetime.combine(d.date(), datetime.min.time()) + timedelta(hours=12)
        for d in days
    ])
    t_night_offsets_hours = np.linspace(H_NOON_TO_START, H_NOON_TO_END, N_SAMPLES)
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
        # Jeśli nic nie brakuje, zwracamy od razu stare dane + hash
        final_subset = {oid: existing_final[oid] for oid in target_ids if oid in existing_final}
        return final_subset, final_hash

    to_process_raw = [
        all_raw_data[oid]
        for oid in missing_ids
        if oid in all_raw_data
    ]

    if VERBOSE:
        print(f"[INFO]  Przeliczam FINAL (Parallel) dla {len(to_process_raw)} brakujących obiektów.")

    new_final: Dict[str, List[Dict]] = {}
    
    process_final_func = partial(
        process_raw_to_final,
        days=days,
        sun_alt_all=sun_alt_all,
        moon_alt_all=moon_alt_all,
        sun_limit=cfg.sun_limit,
        min_altitude=cfg.min_altitude,
        vis=vis,
    )

    with ProcessPoolExecutor(max_workers=max_workers) as executor:
        futures = {
            executor.submit(process_final_func, raw_data): raw_data.obj_id 
            for raw_data in to_process_raw
        }

        for future in tqdm(
            as_completed(futures),
            desc="          Reprocess",
            unit="obj",
            ncols=119,
            colour="green",
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

    return final_subset, final_hash

# =====================================================
# GŁÓWNA FUNKCJA
# =====================================================

def run_engine_from_vis_json(
    json_path: Path | str = PATHS.vis_data,
    max_workers: int = None,
    force_all: bool = False,
) -> Dict[str, List[Dict]]:
    
    vis = load_vis_data(json_path)
    cfg = UserConfig.from_vis_data(vis)

    prev_state = load_engine_state()
    if VERBOSE:
        if prev_state:
            # Używamy get() by uniknąć błędu, jeśli jakiegoś klucza by zabrakło
            print(
                f"[INFO] Poprzednia lokalizacja: {prev_state.get('location_name', 'Unknown')}, "
                f"({prev_state.get('lat', float('nan')):.2f}, {prev_state.get('lon', float('nan')):.2f}), rok {prev_state.get('year')}, "
                f"nad horyzontem {prev_state.get('minalt')}, zmierzch {prev_state.get('sunlimit')}."
            )
            print(
                f"[INFO] Bieżąca lokalizacja:    {cfg.location.name}, "
                f"({cfg.location.lat:.2f}, {cfg.location.lon:.2f}), rok {cfg.year}, "
                f"nad horyzontem {cfg.min_altitude}, zmierzch {cfg.sun_limit}."
            )
        else:
            print("[INFO] Brak poprzedniego stanu silnika (PATHS.engine_state).")

    vis_objects = vis.get("objects", [])
    total_objects = len(vis_objects)

    if VERBOSE:
        print(f"[INFO] W vis_data.json jest {total_objects} obiektów.")

    default_limit = 108
    max_limit = total_objects
    object_limit = max_limit

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
        print_green("=" * 119)
        print_green("ENGINE: Smart Cache System")
        print_green("=" * 119)
        print(f"[INFO] Lokalizacja: {cfg.location.name} ({cfg.location.lat:.2f}°, {cfg.location.lon:.2f}°)")
        print(f"[INFO] Rok: {cfg.year}")
        print(f"[INFO] Parametry: minimalna wysokość obiektu: {cfg.min_altitude}°, "
              f"wysokość słońca: {cfg.sun_limit}°")
        print(f"[INFO] Cache RAW:   {PATHS.observing_raw}")
        print(f"[INFO] Cache FINAL: {PATHS.observing_final}")
        print("=" * 119)

    # --- FORCE-ALL: usunięcie cache (Pamięci podręcznej) ---
    if force_all:
        if os.path.exists(PATHS.observing_raw):
            os.remove(PATHS.observing_raw)
            if VERBOSE:
                print("[FORCE] Usunięto RAW cache.")
        if os.path.exists(PATHS.observing_final):
            os.remove(PATHS.observing_final)
            if VERBOSE:
                print("[FORCE] Usunięto FINAL cache.")
        if os.path.exists(PATHS.observing_hash):
            os.remove(PATHS.observing_hash)
            if VERBOSE:
                print("[FORCE] Usunięto plik hash FINAL.")

    # ========== ETAP 1: RAW ==========
    all_raw_data = compute_raw_data(json_path, object_limit, max_workers=max_workers)

    # ========== ETAP 2: FINAL ==========
    needs_reprocess = should_reprocess(json_path)
    target_ids = get_target_ids_from_vis(vis, object_limit)

    if needs_reprocess or not os.path.exists(PATHS.observing_final):
        if VERBOSE:
            print("[COMPUTE] Parametry/lokalizacja/rok się zmieniły lub brak FINAL cache – pełne przeliczenie FINAL.")
        
        # Przeliczamy od zera
        all_final_data_full, current_hash = reprocess_to_final(json_path, all_raw_data, max_workers=max_workers)
        
        # Wybieramy tylko tyle obiektów, o ile prosił użytkownik (object_limit)
        all_final_data = {oid: all_final_data_full[oid] for oid in target_ids if oid in all_final_data_full}
        
        with open(PATHS.observing_final, "wb") as f:
            pickle.dump(all_final_data_full, f)
        with open(PATHS.observing_hash, "w") as f:
            f.write(current_hash)
        
        if VERBOSE:
            print(f"[SAVE] Final data zapisane do: {PATHS.observing_final}")
            print(f"[SAVE] Hash zapisany do: {PATHS.observing_hash}")
    else:
        if VERBOSE:
            print("[CACHE] Parametry się nie zmieniły, sprawdzam brakujące obiekty w FINAL cache.")

        with open(PATHS.observing_final, "rb") as f:
            final_cache = pickle.load(f)
        # fast-path: jeśli wszystkie target_ids już są w cache, zero liczenia
        missing = [oid for oid in target_ids if oid not in final_cache]
        if not missing:
            # tylko wyciągamy final_data, bez żadnych obliczeń
            all_final_data = {oid: final_cache[oid] for oid in target_ids}
            with open(PATHS.observing_hash, "r") as f:
                current_hash = f.read().strip()
        else:
            # dopiero tutaj wolne reprocess_missing_final
            all_final_data, current_hash = reprocess_missing_final(
                json_path, all_raw_data, final_cache, target_ids, max_workers=max_workers
            )
            final_cache.update(all_final_data)
            with open(PATHS.observing_final, "wb") as f:
                pickle.dump(final_cache, f)
            with open(PATHS.observing_hash, "w") as f:
                f.write(current_hash)

        merged_cache = dict(final_cache)
        merged_cache.update({oid: all_final_data[oid] for oid in all_final_data})

        with open(PATHS.observing_final, "wb") as f:
            pickle.dump(merged_cache, f)
        with open(PATHS.observing_hash, "w") as f:
            f.write(current_hash)
            
        if VERBOSE:
            print(f"[CACHE] ✓ Używam FINAL cache (po ewentualnym uzupełnieniu).")

    # Zapisujemy stan algorytmu na przyszłość
    engine_state = {
        "last_processed": datetime.now().isoformat(),
        "raw_data_objects": len(all_raw_data),
        "final_data_objects": len(all_final_data),
        "lat": cfg.location.lat,
        "lon": cfg.location.lon,
        "year": cfg.year,
        "minalt": cfg.min_altitude,
        "sunlimit": cfg.sun_limit,
        "location_name": cfg.location.name,
    }
    save_engine_state(engine_state)

    if VERBOSE:
        print(f"[INFO] Raw objects: {len(all_raw_data)}. Final objects: {len(all_final_data)}.")
        print_green("=" * 119)
        print_green("[INFO] Silnik zakończył pracę.")
        print_green("=" * 119)

    return all_final_data

if __name__ == "__main__":
    import sys

    force_all = "--force-all" in sys.argv
    run_engine_from_vis_json(
        PATHS.vis_data,
        max_workers=None,
        force_all=force_all,
    )
