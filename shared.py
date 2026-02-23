# shared.py
"""
Moduł konfiguracji współdzielonej (Shared Configuration).

Ten plik pełni rolę centralnego magazynu ustawień dla całego projektu Astrophotography Planner.
Zawiera:
1. PATHS: Ścieżki do plików wejściowych, wynikowych i tymczasowych.
2. Klasy konfiguracyjne:
   - UserConfig: Przechowuje preferencje użytkownika (lokalizacja, rok).
   - CameraConfig: Przechowuje parametry sprzętu (ogniskowa, rozmiar sensora).
3. Stałe astronomiczne i systemowe: Priorytety katalogów, definicje kolorów, stałe fizyczne.

UWAGA DLA DEVELOPERÓW:
Zmienne zdefiniowane w tym pliku są importowane przez wszystkie pozostałe skrypty (kroki 1-7).
Zmiana nazwy klucza w słowniku PATHS lub atrybutu w klasie Config wymaga aktualizacji w całym projekcie.
"""
from __future__ import annotations
from dataclasses import dataclass, field
from pathlib import Path
from typing import Any, Dict, List, Tuple
import json
import hashlib
import pickle
import math
import re
import pandas as pd

# ═══════════════════════════════════════════
# 1. ŚCIEŻKI PLIKÓW (system config)
# ═══════════════════════════════════════════

@dataclass
class PathConfig:
    base_dir: Path = Path(".")
    
    @property
    def ngc_csv(self) -> Path:  # Lub inna nazwa
        return self.base_dir / "data" / "NGC.csv"
    @property
    def vis_data(self) -> Path:
        return self.base_dir / "data" / "vis_data.json"
    @property
    def catalog_full(self) -> Path:
        return self.base_dir / "data" / "katalog_astro_full.csv"
    @property
    def observing_raw(self) -> Path:
        return self.base_dir / "data"  / "observing_data_raw.pkl"
    @property
    def observing_final(self) -> Path:
        return self.base_dir / "data"  / "observing_data.pkl"
    @property
    def observing_hash(self) -> Path:
        return self.base_dir / "data"  / "observing_data_final.hash"  
    @property
    def observing_raw_hash(self) -> Path:
        return self.base_dir / "data"  / "observing_data_raw.hash"
    @property
    def engine_state(self) -> Path:
        return self.base_dir / "data"  / "observing_engine_state.pkl"
    @property
    def starplots(self) -> Path:
        return self.base_dir / "data"  / "starplots"

# Domyślna instancja — skrypty importują ją i używają od razu
PATHS = PathConfig()


# ═══════════════════════════════════════════
# 2. USER CONFIG (modeluje vis_data.json)
# ═══════════════════════════════════════════

@dataclass
class LocationConfig:
    lat: float
    lon: float
    tz: str = "Europe/Warsaw"
    name: str = "Unknown"

    def to_dict(self) -> dict:
        return {"lat": self.lat, "lon": self.lon,
                "tz": self.tz, "name": self.name}


@dataclass
class CameraConfig:
    lens_focal_length: float = 300.0
    sensor_width: float = 23.5
    sensor_height: float = 15.7
    sensor_pitch: float = 3.76
    sensor_rows: int = 4176
    sensor_cols: int = 6248

    def to_dict(self) -> dict:
        return {
            "lens_focal_length": self.lens_focal_length,
            "sensor_width": self.sensor_width,
            "sensor_height": self.sensor_height,
            "sensor_pitch": self.sensor_pitch,
            "sensor_rows": self.sensor_rows,
            "sensor_cols": self.sensor_cols,
        }
    
    @classmethod
    def from_dict(cls, d: dict) -> "CameraConfig":
        return cls(
            lens_focal_length=d.get("lens_focal_length", 300.0),
            sensor_width=d.get("sensor_width", 23.5),
            sensor_height=d.get("sensor_height", 15.7),
            sensor_pitch=d.get("sensor_pitch", 3.76),
            sensor_rows=int(d.get("sensor_rows", 4176)),
            sensor_cols=int(d.get("sensor_cols", 6248)),
        )
    
    def calculate_fov(self) -> Tuple[float, float, float]:
        """
        Oblicza FOV (width, height, diagonal) w stopniach.
        Używa dokładnego wzoru: 2 * atan(sensor / 2*focal).
        """
        def calc_deg(sensor_mm):
            return 2 * math.degrees(math.atan(sensor_mm / (2 * self.lens_focal_length)))
        w_deg = calc_deg(self.sensor_width)
        h_deg = calc_deg(self.sensor_height)
        d_deg = calc_deg(math.hypot(self.sensor_width, self.sensor_height))
        return w_deg, h_deg, d_deg
    
    def get_min_match_size_arcmin(self, percent: float) -> float:
        """
        Zwraca minimalny rozmiar obiektu (w minutach kątowych), 
        który zajmuje 'percent' krótszego boku kadru.
        Zastępuje: FOVCalculator.min_object_size_arcmin
        """
        w, h, _ = self.calculate_fov()
        min_fov_deg = min(w, h)
        return (min_fov_deg * 60.0) * (percent / 100.0)

@dataclass
class UserConfig:
    year: int
    location: LocationConfig
    min_altitude: float = 25.0       # params["minalt"]
    sun_limit: float = -12.0         # params["sunlimit"]
    min_hours: float = 3.0           # params["minhours"]
    min_size_arcmin: float = 10.0    # params["minsizearcmin"]
    bortle_range: Tuple[int, int] = (8, 9)
    has_narrowband: bool = False
    prefer_famous: bool = True
    camera: CameraConfig = field(default_factory=CameraConfig)

    # ─── Wczytanie z vis_data.json ───────────────────────────
    @classmethod
    def from_vis_data(cls, vis: dict) -> "UserConfig":
        """
            lat = vis["location"]["lat"]
            lon = vis["location"]["lon"]
            year = vis["year"]
            params = vis.get("parameters", {})
            obj_min_alt_deg = params.get("minalt", 20.0)
            ...
        """
        loc = vis["location"]
        p = vis.get("parameters", {})
        cam = p.get("camera", {})
        return cls(
            year=vis["year"],
            location=LocationConfig(
                lat=loc["lat"],
                lon=loc["lon"],
                tz=loc.get("tz", "Europe/Warsaw"),
                name=loc.get("name", "Unknown"),
            ),
            min_altitude=p.get("min_alt", 25.0),
            sun_limit=p.get("sun_limit", -12.0),
            min_hours=p.get("min_hours", 3.0),
            min_size_arcmin=p.get("min_size_arcmin", 10.0),
            bortle_range=tuple(p.get("bortle_range", [8, 9])),
            has_narrowband=p.get("has_narrowband", False),
            prefer_famous=p.get("prefer_famous", True),
            camera=CameraConfig.from_dict(cam) if cam else CameraConfig(),
        )

    # ─── Zapis z powrotem do struktury vis_data.json ─────────
    def to_params_dict(self) -> dict:
        """
        Odwrotność from_vis_data() — zwraca params gotowe do zapisu w JSON.
        Używane w 2_ zamiast ręcznego budowania słownika.
        """
        return {
            "min_alt": self.min_altitude,
            "sun_limit": self.sun_limit,
            "min_hours": self.min_hours,
            "min_size_arcmin": self.min_size_arcmin,
            "bortle_range": list(self.bortle_range),
            "has_narrowband": self.has_narrowband,
            "prefer_famous": self.prefer_famous,
            "camera": self.camera.to_dict(),
        }


# ═══════════════════════════════════════════
# 3. STAŁE 
# ═══════════════════════════════════════════
# Słownik polskich nazw miesięcy
MONTH_NAMES = {
    1: "STYCZEŃ", 2: "LUTY", 3: "MARZEC", 4: "KWIECIEŃ",
    5: "MAJ", 6: "CZERWIEC", 7: "LIPIEC", 8: "SIERPIEŃ",
    9: "WRZESIEŃ", 10: "PAŹDZIERNIK", 11: "LISTOPAD", 12: "GRUDZIEŃ",
    }
MONTH_NAMES_PRINT = {
    1: "stycznia", 2: "lutego", 3: "marca", 4: "kwietnia",
    5: "maja", 6: "czerwca", 7: "lipca", 8: "sierpnia",
    9: "września", 10: "października", 11: "listopada", 12: "grudnia",
}
# Słownik typów i polskie odpowiedniki nazw
TYPE_NAMES = {
            "*": "Gwiazda",
            "**": "Gwiazda Podwójna",
            "*Ass": "Association of stars",
            "Nova": "Nova",
            "Gx": "Galaktyka",
            "G": "Galaktyka",
            "S": "Galaktyka spiralna",
            "SB": "Galaktyka spiralna z poprzeczką",
            "E": "Galaktyka eliptyczna",
            "GPair": "Para Galaktyk",
            "GTrpl": "Triplet Galaktyk",
            "GGroup":"Grupa Galaktyk",
            "OC": "Gromada Otwarta",
            "OCl": "Gromada Otwarta",
            "GCl": "Gromada kulista",
            "Cl+N": "Gromada z Mgławicą",
            "GC": "Gromada Kulista",
            "HII": "Obszar HII",
            "NB": "Mgławica",
            "Neb": "Mgławica",
            "EmN": "Mgławica Emisyjna",
            "PN": "Mgławica Planetarna",
            "RfN": "Mgławica Refleksyjna",
            "DN": "Ciemna Mgławica",
            "DrkN": "Ciemna Mgławica",
            "SNR": "Pozostałość po SuperN",
            "Other": "Inny",
            "NonEx": "Nieistniejący",
    }
# Hierarchia ważności katalogów 
# Wartość: im wyżej tym lepszy katalog do wybrania jako "główny ID"
CATALOG_PRIORITY: Dict[str, int] = {
    "ngc": 9, "ic": 8, "sh2": 7, "rcw": 6,
    "lbn": 5, "ced": 4, "pgc": 3, "barn": 2, "ldn": 1,
}

CAT_ORDER: List[str] = list(CATALOG_PRIORITY.keys())  # zachowana kolejność

# =====================================================
# 4. STAŁE ALGORYTMU OBLICZENIOWEGO
# =====================================================

H_START = 15.0                                              # Początek okna (15:00)
H_END = 33.0                                                  # Koniec okna (9:00 rano następnego dnia)
H_RANGE = H_END - H_START
H_Y_RANGE = int(H_RANGE - 1)             #Zakres osi Y na wykresie widoczności
H_Y_START  = 1                                              # Zacząć od 0 + H_Y_START
N_SAMPLES = int(H_RANGE * 60 / 5)   # Liczba próbek (co 5 minut)
CROSSING_SAMPLES = 1440                   # Gęstość próbkowania dla algorytmu skrzyżowań (crossings)


# ═══════════════════════════════════════════
# 5. UTILITY FUNCTIONS
# ═══════════════════════════════════════════
def fmt(n: int) -> str:
    """Format liczby z separatorem tysięcy (spacja)."""
    return f"{int(n):_}".replace("_", " ")

def print_step(msg: str) -> None:
    print(f"\n[INFO] {msg}")

def print_green(msg: str) -> None:
    GREEN = "\033[92m"
    RESET = "\033[0m"
    print(f"{GREEN}{msg}{RESET}")

def load_vis_data(path: Path | str = "vis_data.json") -> Dict[str, Any]:
    with open(path, "r", encoding="utf-8") as f:
        return json.load(f)

def save_vis_data(data: Dict[str, Any], path: Path | str = "vis_data.json") -> None:
    with open(path, "w", encoding="utf-8") as f:
        json.dump(data, f, indent=2, ensure_ascii=False)

def load_observing_data(path: Path | str) -> Dict[str, List[Dict]]:
    with open(path, "rb") as f:
        return pickle.load(f)

def smart_truncate(text: str, max_len: int = 32, ellipsis: str = "...") -> str:
    """Skraca tekst do max_len, tnąc po spacji + ellipsis."""
    if len(text) <= max_len:
        return text 
    # Znajdź ostatnią spację przed max_len - len(ellipsis)
    cut_pos = text.rfind(" ", 0, max_len - len(ellipsis))
    if cut_pos == -1:
        # Brak spacji – tnij brutalnie
        return text[:max_len - len(ellipsis)] + ellipsis    
    return text[:cut_pos].rstrip() + ellipsis

def extract_famous_labels(extra_info: str) -> List[str]:
    """Zwraca listę etykiet 'M33', 'C20', 'H400' z pola extra_info."""
    if not extra_info:
        return []
    seen, result = set(), []
    for token in str(extra_info).split(","):
        t = token.strip().upper()
        if len(t) < 2:
            continue
        prefix, rest = t[0], t[1:].strip()
        if prefix not in ("M", "C", "H") or not rest.isdigit():
            continue
        label = f"{prefix}{int(rest)}"
        if label not in seen:
            seen.add(label)
            result.append(label)
    return result

def build_badge(row) -> str:
    badges = []

    m = row.get("messier_nr")
    c = row.get("caldwell_nr")
    h = row.get("herschel_nr")

    # zabezpieczenie na None / puste
    m = int(m) if m not in (None, "") else 0
    c = int(c) if c not in (None, "") else 0
    h = int(h) if h not in (None, "") else 0

    if m > 0:
        badges.append(f"M{m}")
    if c > 0:
        badges.append(f"C{c}")
    if h > 0:
        badges.append(f"H{h}")

    return f" [{', '.join(badges)}]" if badges else ""

# Funkcje hashujące

def get_engine_raw_hash(lat: float, lon: float, year: int) -> str:
    """
    Hash parametrów dla danych SUROWYCH (zależny tylko od geometrii i czasu).
    Zmiana min_alt/sun_limit NIE wpływa na ten hash. Używane w 3_compute.
    """
    s = f"{lat:.6f}|{lon:.6f}|{int(year)}"
    return hashlib.md5(s.encode()).hexdigest()

def get_engine_final_hash(min_alt: float, sun_limit: float, lat: float, lon: float, year: int) -> str:
    """
    Hash parametrów, które wpływają na ostateczny wynik obliczeń. Używane w 3_compute.
    Obejmuje limity zanieczyszczenia światłem i wysokości, plus geometrię.
    """
    s = f"{min_alt:.6f}|{sun_limit:.6f}|{lat:.6f}|{lon:.6f}|{int(year)}"
    return hashlib.md5(s.encode()).hexdigest()

