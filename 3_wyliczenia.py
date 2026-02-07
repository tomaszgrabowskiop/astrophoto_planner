import json
import pickle
from datetime import datetime, timedelta
from concurrent.futures import ProcessPoolExecutor
from functools import partial

import numpy as np
import pandas as pd
from astropy.coordinates import SkyCoord, EarthLocation, get_sun, get_body, AltAz
from astropy.time import Time
import astropy.units as u

from tqdm import tqdm
from concurrent.futures import ProcessPoolExecutor, as_completed

# Stałe siatki (jak w oryginale)
H_START = 14.0  # 14:00
H_END = 31.0    # 07:00
H_RANGE = H_END - H_START  # 14h
N_SAMPLES = int(H_RANGE * 60 / 5)  # 14h * 60 / 15 = 56: siatka obecności na niebie co 15 minut.

# Zredukowana rozdzielczość dla get_crossings (co 5 min)
CROSSING_SAMPLES = 1440 # co minutę


def load_vis_data(path: str = "vis_data.json"):
    """Ładuje dane z pliku vis_data.json i zwraca słownik."""
    with open(path, "r", encoding="utf-8") as f:
        data = json.load(f)
    return data


def get_crossings(target, height, t_noon, location):
    """Szuka precyzyjnych momentów (UTC) przejścia przez wysokość 'height'."""
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


def mask_to_segments(mask, h_start=H_START, h_end=H_END):
    """
    Zamienia maskę bool (długości N_SAMPLES) na listę ciągłych odcinków [start_rel, end_rel]
    w skali 0–14 (relatywnie do 17:00).
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
    t_grid = np.linspace(h_start, h_end, mask.size, endpoint=False)

    for s_i, e_i in zip(starts_idx, ends_idx):
        t_start = t_grid[s_i]
        if e_i < mask.size:
            t_end = t_grid[e_i]
        else:
            t_end = h_end

        start_rel = t_start - h_start
        end_rel = t_end - h_start
        if end_rel > start_rel:
            segments.append((start_rel, end_rel))

    return segments


def process_single_object(
    obj_data,
    days,
    location,
    sun_alt_limit_deg,
    obj_min_alt_deg,
    days_noon,
    t_night_offsets_hours,
    sun_alt_all,
    moon_alt_all,
):
    """
    Obliczenia dla jednego obiektu - wywoływane równolegle.

    obj_data: (obj_id, ra_deg, dec_deg)
    days: lista dat (pandas.DatetimeIndex albo lista datetime.date)
    location: EarthLocation
    sun_alt_limit_deg: float, np. -12.0
    obj_min_alt_deg: float, np. 25.0
    days_noon: Time array [ndays] z południami dni
    t_night_offsets_hours: 1D np.array z godzinami [2..19] (N_SAMPLES)
    sun_alt_all: np.ndarray [ndays, N_SAMPLES] – wysokość Słońca
    moon_alt_all: np.ndarray [ndays, N_SAMPLES] – wysokość Księżyca
    """
    obj_id, ra, dec = obj_data
    coord = SkyCoord(ra * u.deg, dec * u.deg)

    ndays = len(days)
    results = []

    # --- 1. Pełna siatka czasu nocy dla wszystkich dni tego obiektu ---
    # t_noon + offsety godzinowe (2..19h => 17:00–07:00 lokalnie)
    # shape: (ndays, N_SAMPLES)
    t_grid = days_noon[:, None] + t_night_offsets_hours[None, :] * u.hour

    # Jeden AltAz dla wszystkich punktów czasu (flatten, potem reshape)
    frame = AltAz(obstime=t_grid.reshape(-1), location=location)

    # Wektoryzowana transformacja RA/DEC -> AltAz dla wszystkich dni/czasów
    o_alt_all = coord.transform_to(frame).alt.deg.reshape(ndays, N_SAMPLES)

    # --- 2. Dzień po dniu: używamy precomputed Sun/Moon + o_alt_all ---
    for day_idx, d in enumerate(days):
        d_date = d.date() if hasattr(d, "date") else d
        t_noon = days_noon[day_idx]

        # Punkty charakterystyczne (dokładne crossingi) liczymy jak dotąd, per dzień
        sun_pts = get_crossings(get_sun(t_noon), sun_alt_limit_deg, t_noon, location)
        obj_pts = get_crossings(coord, obj_min_alt_deg, t_noon, location)

        # Wysokości dla tego dnia: 1D [N_SAMPLES]
        o_alt = o_alt_all[day_idx]
        s_alt = sun_alt_all[day_idx]
        m_alt = moon_alt_all[day_idx]

        # maska jakościowa: alt > obj_min_alt_deg, Sun < sun_alt_limit_deg
        quality_mask = (o_alt > obj_min_alt_deg) & (s_alt < sun_alt_limit_deg)

        # tranzit względny (maksimum alt w skali 0–14)
        transit_rel = np.argmax(o_alt) * (H_RANGE / N_SAMPLES)

        # sumaryczne godziny jakościowe
        q_hours = np.sum(quality_mask) * (H_RANGE / N_SAMPLES)
        m_hours = np.sum(quality_mask & (m_alt < 0.0)) * (H_RANGE / N_SAMPLES)

        # ciągłe segmenty jakościowe w skali 0–14
        qual_segments = mask_to_segments(quality_mask)

        results.append({
            "day": d_date,
            "sun_pts": [p.datetime for p in sun_pts],
            "obj_pts": [p.datetime for p in obj_pts],
            "transit_rel": transit_rel,
            "q_hours": q_hours,
            "m_hours": m_hours,
            "qual_segments": qual_segments,
        })

    return obj_id, results


def run_engine_from_vis_json(json_path="vis_data.json", max_workers=None):
    """
    Główna funkcja - wczytuje vis_data.json, wyciąga:
    - lokalizację
    - rok
    - listę obiektów
    i uruchamia obliczenia równolegle po obiektach.
    Dodatkowo:
    - prekomputuje wysokości Słońca i Księżyca w siatce [dzień, czas] dla całego roku,
    - wektoryzuje AltAz dla obiektów po czasie, używając broadcastingu z Sun/Moon.
    """
    vis = load_vis_data(json_path)

    # Lokalizacja
    lat = vis["location"]["lat"]
    lon = vis["location"]["lon"]
    location = EarthLocation(lat=lat * u.deg, lon=lon * u.deg)
    tz_name = vis["location"].get("tz", "UTC")
    city_name = vis["location"].get("name", "Unknown")
    # Rok
    year = vis["year"]

    # Parametry ograniczeń
    params = vis.get("parameters", {})
    obj_min_alt_deg = params.get("minalt", 20.0)
    sun_alt_limit_deg = params.get("sunlimit", -6.0)
    min_hours = params.get("minhours", None)
    
    # Lista obiektów
    objects = vis["objects"]
    objects_data = [
        (obj["id"], obj["ra"], obj["dec"])
        for obj in objects
    ]

    # Zakres dni w danym roku
    days = pd.date_range(f"{year}-01-01", f"{year}-12-31", freq="D")
    ndays = len(days)

    print(f"Engine: Start obliczeń dla {len(objects_data)} obiektów, {len(days)} dni.")
    print(f"Engine: Używam równoległości z max_workers={max_workers or 'auto'}.")
    print(f"Engine: Lokalizacja lat={lat}, lon={lon}, miasto={city_name}, strefa_czasowa={tz_name}, rok={year}.")
    print(f"Engine: Parametry: minhours={min_hours}, minalt={obj_min_alt_deg}, sunlimit={sun_alt_limit_deg}.")

    # --- PREKOMPUTACJA SŁOŃCA I KSIĘŻYCA DLA WSZYSTKICH DNI/NOCY ---

    # Południa każdego dnia jako Time (UTC)
    days_noon = Time([
        datetime.combine(d.date(), datetime.min.time()) + timedelta(hours=12)
        for d in days
    ])

    # Odstępy czasowe dla nocy (2..19h od południa = 17:00–07:00)
    t_night_offsets_hours = np.linspace(2, 19, N_SAMPLES)

    # Pełna siatka czasu: shape (ndays, N_SAMPLES)
    t_grid_all = days_noon[:, None] + t_night_offsets_hours[None, :] * u.hour
    obstime_all = t_grid_all.reshape(-1)  # (ndays * N_SAMPLES,)

    # Wspólna rama AltAz dla całego roku i nocy
    frame_all = AltAz(obstime=obstime_all, location=location)

    # Słońce
    sun_icrs_all = get_sun(obstime_all)
    sun_alt_all = sun_icrs_all.transform_to(frame_all).alt.deg.reshape(ndays, N_SAMPLES)

    # Księżyc
    moon_icrs_all = get_body("moon", obstime_all)
    moon_alt_all = moon_icrs_all.transform_to(frame_all).alt.deg.reshape(ndays, N_SAMPLES)

    all_data = {}
    
    # Równoległe przetwarzanie obiektów z paskiem postępu
    with ProcessPoolExecutor(max_workers=max_workers) as executor:
        process_func = partial(
            process_single_object,
            days=days,
            location=location,
            sun_alt_limit_deg=sun_alt_limit_deg,
            obj_min_alt_deg=obj_min_alt_deg,
            days_noon=days_noon,
            t_night_offsets_hours=t_night_offsets_hours,
            sun_alt_all=sun_alt_all,
            moon_alt_all=moon_alt_all,
        )
    
        # uruchamiamy zadania
        future_to_objid = {
            executor.submit(process_func, obj_data): obj_data[0]
            for obj_data in objects_data
        }
    
        # pasek postępu – 1 krok = 1 obiekt
        with tqdm(total=len(objects_data), desc="Engine", unit="obj") as pbar:
            for future in as_completed(future_to_objid):
                obj_id = future_to_objid[future]
                obj_id_result, obj_results = future.result()
                all_data[obj_id_result] = obj_results
    
                pbar.set_postfix_str(f"ostatni={obj_id}")
                pbar.update(1)


    with open("observing_data.pkl", "wb") as f:
        pickle.dump(all_data, f)

    print("Engine: Dane zapisane.")
    return all_data


if __name__ == "__main__":
    # Użyj max_workers=None dla automatycznego doboru (liczba rdzeni)
    run_engine_from_vis_json("vis_data.json", max_workers=None)
