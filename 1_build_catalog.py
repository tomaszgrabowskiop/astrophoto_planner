#!/usr/bin/env python3
"""
Krok 1: Budowa i agregacja katalogu astronomicznego.

Skrypt ten odpowiada za przygotowanie głównej, surowej bazy danych obiektów ('katalog_astro_full.csv').

Główne zadania skryptu:
1. Wczytanie bazowej listy obiektów (domyślnie OpenNGC z pliku data/NGC.csv).
2. Pobranie uzupełniających danych z serwisu VizieR dla katalogów specjalistycznych
   (Sharpless, Barnard, LDN, LBN, RCW, Cederblad itp.).
3. Normalizacja danych: ujednolicenie jednostek, formatów współrzędnych i nazw.
4. Entity Resolution: Inteligentne łączenie duplikatów (np. rozpoznanie, że NGC 1499 i Sh2-220 to ten sam obiekt).
5. Wstępne filtrowanie: Interaktywne odrzucenie obiektów zbyt małych lub zbyt słabych (wg preferencji użytkownika).

Wymagania:
- Aktywne połączenie z internetem (do zapytań VizieR).
- Plik wejściowy: data/NGC.csv

Wyjście:
- Plik: data/katalog_astro_full.csv (baza gotowa do scoringu).
"""
import re                              
import warnings
import numpy as np
import pandas as pd
import networkx as nx
import astropy.units as u
from astroquery.vizier import Vizier
from astropy.coordinates import SkyCoord, search_around_sky
from tqdm import tqdm

warnings.filterwarnings('ignore')

import ssl, certifi
ssl._create_default_https_context = lambda: ssl.create_default_context(cafile=certifi.where())

from shared import CATALOG_PRIORITY, CAT_ORDER, PATHS, fmt, print_step, print_green
from analyse_catalog import analyze_catalog

# === KONFIGURACJA DOMYŚLNA ===

TRASH_MAG   = 24.0   # obiekty ciemniejsze (wyższe) odrzucamy przed merge
TRASH_SIZE  = 0.5    # obiekty mniejsze niż 0.5' odrzucamy przed merge

# Wartości domyślne (użytkownik może je nadpisać interaktywnie)
_DEFAULT_MIN_SIZE  = 5.0
_DEFAULT_MAX_MAG   = 17.0
_DEFAULT_MAX_CLUST = 300.0   # limit rozmiaru klastra po merge

def _ask_float(prompt: str, default: float) -> float:
    """Pyta użytkownika o wartość float; Enter = wartość domyślna."""
    try:
        raw = input(f"  {prompt} [domyślnie: {default}]: ").strip()
        return float(raw) if raw else default
    except ValueError:
        print(f"  ⚠  Nieprawidłowa wartość – używam domyślnej: {default}")
        return default

def _configure_params() -> tuple[float, float, float]:
    print()
    print_green("=" * 119)
    print_green(f"{'=' * 44} PARAMETRY FILTROWANIA KATALOGU {'=' * 43} ")
    print_green("=" * 119)
    print("  (naciśnij Enter, aby zaakceptować wartość domyślną)")
    min_size  = _ask_float("Minimalny rozmiar obiektu  [arcmin]", _DEFAULT_MIN_SIZE)
    max_mag   = _ask_float("MInimalna jasność (granica słaba) [mag]", _DEFAULT_MAX_MAG)
    max_obj = _ask_float("Maksymalny dopuszczalny rozmiar obiektu (cap) [arcmin]", _DEFAULT_MAX_CLUST) # obiekty większe, będą przycięte do tej wielkości
 
    return min_size, max_mag, max_obj

MIN_SIZE_ARCMIN, MAX_MAG, MAX_OBJECT_SIZE_ARCMIN = _configure_params()

SH2_COMMON_NAMES = {
    "Sh2-25": "The Lagoon Nebula",
    "Sh2-101": "Tulip Nebula", "Sh2-103": "Loop", "Sh2-105": "Crescent Nebula",
    "Sh2-108": "Sadr Region", "Sh2-11": "War and Peace Nebula", "Sh2-117": "N American & Pelican Nebula",
    "Sh2-125": "Cocoon Nebula", "Sh2‑126": "Great Lacerta Nebula",  "Sh2-129": "Flying Bat Nebula", 
    "Sh2-136": "Ghost Nebula", "Sh2-155": "Cave Nebula",
    "Sh2-162": "Bubble Nebula", "Sh2-184": "Pac Man Nebula", "Sh2-190": "Heart Nebula",
    "Sh2-197": "Maffei 2", "Sh2-199": "Soul Nebula", "Sh2-220": "California Nebula",
    "Sh2-229": "Flaming Star Nebula", "Sh2-234": "Spider Nebula", "Sh2-237": "Fly Nebula",
    "Sh2-238": "Hind's Variable Nebula", "Sh2-244": "Crab Nebula", "Sh2-245": "Fishhook Nebula",
    "Sh2-248": "Jellyfish Nebula", "Sh2-252": "Monkey Head Nebula",
    "Sh2-261": "Lower's Nebula", "Sh2-264": "Angelfish Nebula", "Sh2-273": "Fox Fur Nebula",
    "Sh2-274": "Medusa Nebula", "Sh2-275": "Rosette Nebula", "Sh2-276": "Barnard's Loop",
    "Sh2-277": "Flame Nebula", "Sh2-279": "Running Man Nebula", "Sh2-281": "Orion Nebula",
    "Sh2-292": "Seagull Nebula head", "Sh2-296": "Seagull Nebula wings", "Sh2-298": "Thor's Helmet",
    "Sh2-30": "Trifid Nebula", "Sh2-311": "Skull & Crossbones Nebula", "Sh2-45": "Omega Nebula",
    "Sh2-49": "Eagle Nebula", "Sh2-54": "Cauda", "Sh2-6": "Bug Nebula", "Sh2-8": "Cat's Paw Nebula",
}

NGC_EXTRA_COMMON_NAMES = {
    'NGC0188': 'North Celestial Pole Cluster',
    'NGC1333': 'Embryo Nebula',
    'NGC1746': 'Cluster of Clusters',
    'NGC1788': 'Cosmic Bat Nebula',
    'NGC1960': 'Pinwheel Cluster',
    'NGC2168': 'Shoe-Buckle Cluster',
    'NGC2170': 'Angel Nebula',
    'NGC2682': 'King Cobra Cluster',
    'NGC5457': 'Pinwheel Galaxy',
    'NGC6334': "Cat's Paw Nebula",
    'NGC6355':  'The Snake Nebula',
    'NGC 6979': "Fleming's Triangular Nebula",
    'NGC7078': 'Great Pegasus Cluster',
    'NGC7092': 'Pyramid Cluster',
    'NGC7380': 'Wizard Nebula',
    'NGC7822': "Question Mark Nebula",
    'IC0410': "Tadpole Nebula",
    'IC0423': 'Tear Drop Nebula',
    'IC0446': 'Coyote Cloud',
    'IC0447': "Dreyer's Nebula",
    'IC1318': 'Gamma Cygni Nebula',
    'IC1396': "Elephant's Trunk Nebula",
    'IC1613': 'Cetus Dwarf Galaxy',
    'IC4665': 'Summer Beehive Cluster',
    'IC4756': 'Graff’s Cluster',
}

MESSIER_CORRECTIONS = {"NGC1432": 45}

# === LOOKUP: HERSCHEL 400 ===
HERSCHEL_400 = {
    "NGC0040": 1, "NGC0129": 2, "NGC0136": 3, "NGC0157": 4, "NGC0185": 5, "NGC0205": 6, "NGC0225": 7, "NGC0246": 8, "NGC0247": 9, "NGC0253": 10, "NGC0278": 11, "NGC0288": 12, "NGC0381": 13, "NGC0404": 14, "NGC0436": 15, "NGC0457": 16, "NGC0488": 17, "NGC0524": 18, "NGC0559": 19, "NGC0584": 20, "NGC0596": 21, "NGC0598": 22, "NGC0613": 23, "NGC0615": 24, "NGC0637": 25, "NGC0651": 26, "NGC0654": 27, "NGC0659": 28, "NGC0663": 29, "NGC0720": 30, "NGC0752": 31, "NGC0772": 32, "NGC779": 33, "NGC869": 34, "NGC884": 35, "NGC891": 36, "NGC0908": 37, "NGC0936": 38, "NGC1022": 39, "NGC1023": 40, "NGC1027": 41, "NGC1052": 42, "NGC1055": 43, "NGC1084": 44, "NGC1245": 45, "NGC1342": 46, "NGC1407": 47, "NGC1444": 48, "NGC1501": 49, "NGC1502": 50, "NGC1513": 51, "NGC1528": 52, "NGC1535": 53, "NGC1545": 54, "NGC1647": 55, "NGC1664": 56, "NGC1788": 57, "NGC1817": 58, "NGC1857": 59, "NGC1907": 60, "NGC1931": 61, "NGC1961": 62, "NGC1964": 63, "NGC1980": 64, "NGC1999": 65, "NGC2022": 66, "NGC2024": 67, "NGC2126": 68, "NGC2129": 69, "NGC2158": 70, "NGC2169": 71, "NGC2185": 72, "NGC2186": 73, "NGC2194": 74, "NGC2204": 75, "NGC2215": 76, "NGC2232": 77, "NGC2244": 78, "NGC2251": 79, "NGC2264": 80, "NGC2266": 81, "NGC2281": 82, "NGC2286": 83, "NGC2301": 84, "NGC2304": 85, "NGC2311": 86, "NGC2324": 87, "NGC2335": 88, "NGC2343": 89, "NGC2353": 90, "NGC2354": 91, "NGC2355": 92, "NGC2360": 93, "NGC2362": 94, "NGC2371": 95, "NGC2372": 96, "NGC2392": 97, "NGC2395": 98, "NGC2403": 99, "NGC2419": 100, "NGC2420": 101, "NGC2421": 102, "NGC2422": 103, "NGC2423": 104, "NGC2438": 105, "NGC2440": 106, "NGC2479": 107, "NGC2482": 108, "NGC2489": 109, "NGC2506": 110, "NGC2509": 111, "NGC2527": 112, "NGC2539": 113, "NGC2548": 114, "NGC2567": 115, "NGC2571": 116, "NGC2613": 117, "NGC2627": 118, "NGC2655": 119, "NGC2681": 120, "NGC2683": 121, "NGC2742": 122, "NGC2768": 123, "NGC2775": 124, "NGC2782": 125, "NGC2787": 126, "NGC2811": 127, "NGC2841": 128, "NGC2859": 129, "NGC2903": 130, "NGC2950": 131, "NGC2964": 132, "NGC2974": 133, "NGC2976": 134, "NGC2985": 135, "NGC3034": 136, "NGC3077": 137, "NGC3079": 138, "NGC3115": 139, "NGC3147": 140, "NGC3166": 141, "NGC3169": 142, "NGC3184": 143, "NGC3190": 144, "NGC3193": 145, "NGC3198": 146, "NGC3226": 147, "NGC3227": 148, "NGC3242": 149, "NGC3245": 150, "NGC3277": 151, "NGC3294": 152, "NGC3310": 153, "NGC3344": 154, "NGC3377": 155, "NGC3379": 156, "NGC3384": 157, "NGC3395": 158, "NGC3412": 159, "NGC3414": 160, "NGC3432": 161, "NGC3486": 162, "NGC3489": 163, "NGC3504": 164, "NGC3521": 165, "NGC3556": 166, "NGC3593": 167, "NGC3607": 168, "NGC3608": 169, "NGC3610": 170, "NGC3613": 171, "NGC3619": 172, "NGC3621": 173, "NGC3626": 174, "NGC3628": 175, "NGC3631": 176, "NGC3640": 177, "NGC3655": 178, "NGC3665": 179, "NGC3675": 180, "NGC3686": 181, "NGC3726": 182, "NGC3729": 183, "NGC3810": 184, "NGC3813": 185, "NGC3877": 186, "NGC3893": 187, "NGC3898": 188, "NGC3900": 189, "NGC3912": 190, "NGC3938": 191, "NGC3941": 192, "NGC3945": 193, "NGC3949": 194, "NGC3953": 195, "NGC3962": 196, "NGC3982": 197, "NGC3992": 198, "NGC3998": 199, "NGC4026": 200, "NGC4027": 201, "NGC4030": 202, "NGC4036": 203, "NGC4038": 204, "NGC4041": 205, "NGC4051": 206, "NGC4085": 207, "NGC4088": 208, "NGC4102": 209, "NGC4111": 210, "NGC4143": 211, "NGC4147": 212, "NGC4150": 213, "NGC4151": 214, "NGC4179": 215, "NGC4203": 216, "NGC4214": 217, "NGC4216": 218, "NGC4245": 219, "NGC4251": 220, "NGC4258": 221, "NGC4261": 222, "NGC4273": 223, "NGC4274": 224, "NGC4278": 225, "NGC4281": 226, "NGC4293": 227, "NGC4303": 228, "NGC4314": 229, "NGC4346": 230, "NGC4350": 231, "NGC4361": 232, "NGC4365": 233, "NGC4371": 234, "NGC4394": 235, "NGC4414": 236, "NGC4419": 237, "NGC4429": 238, "NGC4435": 239, "NGC4438": 240, "NGC4442": 241, "NGC4448": 242, "NGC4449": 243, "NGC4450": 244, "NGC4459": 245, "NGC4473": 246, "NGC4477": 247, "NGC4478": 248, "NGC4485": 249, "NGC4490": 250, "NGC4494": 251, "NGC4526": 252, "NGC4527": 253, "NGC4535": 254, "NGC4536": 255, "NGC4546": 256, "NGC4548": 257, "NGC4550": 258, "NGC4559": 259, "NGC4565": 260, "NGC4570": 261, "NGC4594": 262, "NGC4596": 263, "NGC4618": 264, "NGC4631": 265, "NGC4636": 266, "NGC4643": 267, "NGC4654": 268, "NGC4656": 269, "NGC4660": 270, "NGC4665": 271, "NGC4666": 272, "NGC4689": 273, "NGC4697": 274, "NGC4698": 275, "NGC4699": 276, "NGC4725": 277, "NGC4753": 278, "NGC4754": 279, "NGC4762": 280, "NGC4781": 281, "NGC4800": 282, "NGC4845": 283, "NGC4856": 284, "NGC4866": 285, "NGC4900": 286, "NGC4958": 287, "NGC4995": 288, "NGC5005": 289, "NGC5033": 290, "NGC5054": 291, "NGC5195": 292, "NGC5248": 293, "NGC5273": 294, "NGC5322": 295, "NGC5363": 296, "NGC5364": 297, "NGC5466": 298, "NGC5473": 299, "NGC5474": 300, "NGC5557": 301, "NGC5566": 302, "NGC5576": 303, "NGC5631": 304, "NGC5634": 305, "NGC5676": 306, "NGC5689": 307, "NGC5694": 308, "NGC5746": 309, "NGC5846": 310, "NGC5866": 311, "NGC5897": 312, "NGC5907": 313, "NGC5982": 314, "NGC6118": 315, "NGC6144": 316, "NGC6171": 317, "NGC6207": 318, "NGC6217": 319, "NGC6229": 320, "NGC6235": 321, "NGC6284": 322, "NGC6287": 323, "NGC6293": 324, "NGC6304": 325, "NGC6316": 326, "NGC6342": 327, "NGC6355": 328, "NGC6356": 329, "NGC6369": 330, "NGC6401": 331, "NGC6426": 332, "NGC6440": 333, "NGC6445": 334, "NGC6451": 335, "NGC6514": 336, "NGC6517": 337, "NGC6520": 338, "NGC6522": 339, "NGC6528": 340, "NGC6540": 341, "NGC6543": 342, "NGC6544": 343, "NGC6553": 344, "NGC6568": 345, "NGC6569": 346, "NGC6583": 347, "NGC6624": 348, "NGC6629": 349, "NGC6633": 350, "NGC6638": 351, "NGC6642": 352, "NGC6645": 353, "NGC6664": 354, "NGC6712": 355, "NGC6755": 356, "NGC6756": 357, "NGC6781": 358, "NGC6802": 359, "NGC6818": 360, "NGC6823": 361, "NGC6826": 362, "NGC6830": 363, "NGC6834": 364, "NGC6866": 365, "NGC6882": 366, "NGC6885": 367, "NGC6905": 368, "NGC6910": 369, "NGC6934": 370, "NGC6939": 371, "NGC6940": 372, "NGC6946": 373, "NGC7000": 374, "NGC7006": 375, "NGC7008": 376, "NGC7009": 377, "NGC7044": 378, "NGC7062": 379, "NGC7086": 380, "NGC7128": 381, "NGC7142": 382, "NGC7160": 383, "NGC7209": 384, "NGC7217": 385, "NGC7243": 386, "NGC7296": 387, "NGC7331": 388, "NGC7380": 389, "NGC7448": 390, "NGC7479": 391, "NGC7510": 392, "NGC7606": 393, "NGC7662": 394, "NGC7686": 395, "NGC7723": 396, "NGC7727": 397, "NGC7789": 398, "NGC7790": 399, "NGC7814": 400
    }

# === LOOKUP: CALDWELL ===
CALDWELL = {
    # IC Caldwell
    "IC0342": 5, "IC0405": 31, "IC1613": 51, "IC2391": 85, "IC2602": 102, "IC2944": 100, "IC5146": 19,
    # NGC Caldwell
    "NGC0040": 2, "NGC0055": 72, "NGC0104": 106, "NGC0147": 17, "NGC0185": 18, "NGC0188": 1,
    "NGC0246": 56, "NGC0247": 62, "NGC0253": 65, "NGC0300": 70, "NGC0362": 104, "NGC0457": 13,
    "NGC0559": 8, "NGC0663": 10, "NGC0752": 28, "NGC0891": 23, "NGC1097": 67, "NGC1261": 87,
    "NGC1275": 24, "NGC1851": 73, "NGC2070": 103, "NGC2237": 49, "NGC2239": 50, "NGC2261": 46,
    "NGC2360": 58, "NGC2362": 64, "NGC2392": 39, "NGC2403": 7, "NGC2419": 25, "NGC2477": 71,
    "NGC2506": 54, "NGC2516": 96, "NGC2775": 48, "NGC2867": 90, "NGC3115": 53, "NGC3132": 74,
    "NGC3195": 109, "NGC3201": 79, "NGC3242": 59, "NGC3372": 92, "NGC3532": 91, "NGC3626": 40,
    "NGC3766": 97, "NGC4038": 60, "NGC4039": 61, "NGC4236": 3, "NGC4244": 26, "NGC4372": 108,
    "NGC4449": 21, "NGC4559": 36, "NGC4565": 38, "NGC4609": 98, "NGC4631": 32, "NGC4697": 52,
    "NGC4755": 94, "NGC4833": 105, "NGC4884": 35, "NGC4945": 83, "NGC5005": 29, "NGC5128": 77,
    "NGC5139": 80, "NGC5248": 45, "NGC5286": 84, "NGC5694": 66, "NGC5823": 88, "NGC6025": 95,
    "NGC6087": 89, "NGC6101": 107, "NGC6124": 75, "NGC6193": 82, "NGC6231": 76, "NGC6302": 69,
    "NGC6352": 81, "NGC6397": 86, "NGC6541": 78, "NGC6543": 6, "NGC6729": 68, "NGC6744": 101,
    "NGC6752": 93, "NGC6822": 57, "NGC6826": 15, "NGC6882": 37, "NGC6888": 27, "NGC6934": 47,
    "NGC6946": 12, "NGC6960": 34, "NGC6992": 33, "NGC7000": 20, "NGC7006": 42, "NGC7009": 55,
    "NGC7023": 4, "NGC7243": 16, "NGC7293": 63, "NGC7331": 30, "NGC7479": 44, "NGC7635": 11,
    "NGC7662": 22, "NGC7814": 43,
    # Addendum Caldwell (bez NGC/IC prefixu)
    "C9": 9, "C14": 14, "C41": 41, "C99": 99,
}

def safe_float(value) -> float | None:
    """Bezpieczna konwersja na float. Zwraca None dla pustych/błędnych."""
    if value is None:
        return None
    s = str(value).strip()
    if not s or s.lower() in ("nan", "none", ""):
        return None
    try:
        return float(s)
    except ValueError:
        return None

def sexa_to_deg(sexa: str, is_ra: bool) -> float | None:
    """
    Konwersja HH:MM:SS.SS (RA) lub +/-DD:MM:SS.SS (Dec) na stopnie.
    Używana wyłącznie dla NGC/IC (OpenNGC podaje współrzędne w tej notacji).
    Dla katalogów VizieR (stopnie) używaj convert_ra_dec().
    """
    if not sexa or sexa.strip() == "":
        return None
    s = sexa.strip().replace(" ", "")
    sign = 1
    if s and s[0] in "+-":
        if s[0] == "-":
            sign = -1
        s = s[1:]
    parts = s.split(":")
    if len(parts) != 3:
        return None
    try:
        h_d = float(parts[0])
        m   = float(parts[1])
        sec = float(parts[2])
    except ValueError:
        return None
    val = h_d + m / 60.0 + sec / 3600.0
    if is_ra:
        val *= 15.0     # RA: godziny → stopnie
    else:
        val *= sign     # Dec: zachowaj znak
    return val

def get_magnitude(vmag, bmag, jmag, hmag, kmag, surfbr, obj_type: str) -> float | None:
    """
    Zwraca magnitude według priorytetu astronomicznego:
      1. V-Mag (zawsze priorytet — wizualne = najbardziej porównywalne)
      2. Galaktyki ('G', 'GPair', 'GTrpl', 'GGroup'): SurfBr jeśli V-Mag puste
      3. Inne: najniższa wartość z B/J/H/K (najniższa = najjaśniejszy obiekt)
    Zwraca None jeśli żadne magnitudo niedostępne.
    """
    v = safe_float(vmag)
    if v is not None:
        return v

    GALAXY_TYPES = {"G", "GPAIR", "GTRPL", "GGROUP"}
    if obj_type.upper() in GALAXY_TYPES:
        sb = safe_float(surfbr)
        if sb is not None:
            return sb

    mags = [safe_float(x) for x in [bmag, jmag, hmag, kmag]]
    valid = [m for m in mags if m is not None]
    return min(valid) if valid else None

def convert_ra_dec(
    df: pd.DataFrame,
    ra_col: str,
    dec_col: str,
    unit_type: str = "sexagesimal",
) -> pd.DataFrame:
    """
    Konwertuje RA/Dec z kolumn DataFrame na stopnie dziesiętne.
    Używana dla katalogów VizieR (_RAJ2000, _DEJ2000 już w stopniach).
    Dla NGC/IC używaj sexa_to_deg() wiersz po wierszu.

    unit_type: 'sexagesimal' (godziny:minuty:sekundy) lub 'deg' (stopnie)
    """
    if df.empty or ra_col not in df.columns or dec_col not in df.columns:
        return df
    df = df.dropna(subset=[ra_col, dec_col]).copy()
    try:
        unit = (u.hourangle, u.deg) if unit_type == "sexagesimal" else u.deg
        c = SkyCoord(
            ra=df[ra_col].astype(str).values,
            dec=df[dec_col].astype(str).values,
            unit=unit,
            frame="icrs",
        )
        df["ra"]  = np.round(c.ra.deg, 5)
        df["dec"] = np.round(c.dec.deg, 5)
    except Exception as e:
        print(f"  [!] Błąd konwersji RA/DEC: {e}")
    return df

def process_common_names(text_series: pd.Series) -> str:
    """
    Scala nazwy zwyczajowe z kilku obiektów w klastrze:
      1. Rozbija stringi po przecinkach.
      2. Usuwa 'The ' z początku (case-insensitive).
      3. Sortuje od najdłuższej (priorytet dłuższych fraz).
      4. Odrzuca nazwę jeśli ZAWIERA SIĘ w innej, już zaakceptowanej nazwie.
         Np. 'Orion Nebula' odpada jeśli jest 'Great Orion Nebula'.
         Ale 'Triangulum Galaxy' i 'Triangulum Pinwheel' zostają obie.
    """
    unique_candidates: set[str] = set()
    for val in text_series.dropna().astype(str):
        for p in val.split(","):
            p = p.strip()
            if not p or p.lower() in ("nan", "none", ""):
                continue
            if p.lower().startswith("the "):
                p = p[4:].strip()
            if p:
                unique_candidates.add(p)

    sorted_candidates = sorted(unique_candidates, key=len, reverse=True)
    kept: list[str] = []
    for candidate in sorted_candidates:
        if not any(candidate.lower() in existing.lower() for existing in kept):
            kept.append(candidate)

    return ", ".join(sorted(kept))
    
def build_famous_columns(df: pd.DataFrame) -> pd.DataFrame:
    """
    Dodaje kolumny messier_nr, herschel_nr, caldwell_nr (int, 0 = brak)
    na podstawie kolumny 'id' (pełna nazwa obiektu np. 'NGC7000').
    Dla NGC: uzupełnia też z kolumny 'M' (bezpośredni numer Messiera
    z OpenNGC) oraz z MESSIER_CORRECTIONS (korekty ręczne).
    """
    # Herschel 400: ID obiektu → numer H
    # HERSCHEL_400 używa formatu "NGC0040" (z wiodącymi zerami)
    df["herschel_nr"] = (
        df["id"]
        .str.strip()
        .map(HERSCHEL_400)
        .fillna(0)
        .astype(int)
    )

    # Caldwell: ID obiektu → numer C
    df["caldwell_nr"] = (
        df["id"]
        .str.strip()
        .map(CALDWELL)
        .fillna(0)
        .astype(int)
    )

    # Messier: najpierw z kolumny 'M' (jeśli istnieje — tylko NGC/IC)
    if "M" in df.columns:
        df["messier_nr"] = (
            pd.to_numeric(df["M"], errors="coerce")
            .fillna(0)
            .astype(int)
        )
        # Korekty ręczne (np. NGC1432 = M45, często brak w OpenNGC)
        for ngc_id, m_num in MESSIER_CORRECTIONS.items():
            mask = df["id"].str.strip() == ngc_id
            df.loc[mask, "messier_nr"] = m_num
    else:
        # Katalogi VizieR nie mają kolumny M — startują z 0
        # Smart Merge uzupełni jeśli wchłoną obiekt NGC z messier_nr > 0
        df["messier_nr"] = 0

    return df

# === 1. POBIERANIE DANYCH ===
def fetch_data() -> dict[str, pd.DataFrame]:
    print("[INFO] Pobieranie katalogów")
    
    # NGC/IC z lokalnego pliku 
    print(f"       > NGC/IC z {PATHS.ngc_csv}")
    try:
        # Używamy surowego pliku NGC.csv, 
        df_ngc = pd.read_csv(PATHS.ngc_csv, sep=";")
    except FileNotFoundError:
        print_step(f"[!] Brak pliku {PATHS.ngc_csv}!")
        df_ngc = pd.DataFrame()
        
    v_std = Vizier(row_limit=-1, columns=['**', '_RAJ2000', '_DEJ2000'])
    
    def get_v(cat_id: str, name: str) -> pd.DataFrame:
        print(f"       > {name} ({cat_id})")
        try:
            cats = v_std.get_catalogs(cat_id)
            return cats[0].to_pandas() if cats else pd.DataFrame()
        except Exception as e:
            print(f"        [!] Błąd połączenia: {e}")
            return pd.DataFrame()

    data = {
        "ngc": df_ngc,
        "sh2": get_v('VII/20/catalog', 'Sharpless'),
        "barn": get_v('VII/220A', 'Barnard'),
        "rcw": get_v('VII/216', 'RCW'),
        "pgc": get_v('VII/119', 'PGC'),
        "ldn": get_v('VII/7A', 'Lynds Dark'),
        "lbn": get_v('VII/9', 'Lynds Bright'),
        "ced": get_v('VII/231', 'Cederblad'),
    }
    
    print(f"\n       PODSUMOWANIE POBIERANIA:")
    total_rows = 0
    for k, v in data.items():
        cnt = len(v)
        print(f"       {k.upper():8}: {fmt(cnt):>8} wierszy")
        total_rows += cnt
    print(" " * 7 + "=" * 26)
    print(f"       Razem   : {fmt(total_rows):>8} wierszy\n")
    
    return data

# === 2. PEŁNA NORMALIZACJA ===
def normalize_all(raw: dict[str, pd.DataFrame]) -> pd.DataFrame:
    print("[INFO] Normalizacja katalogów (unifikacja kolumn)")
    results = []

    # 1. NGC/IC 
    if not raw['ngc'].empty:
        df = raw['ngc'].copy()
        print("       > Normalizacja NGC/IC")
        
        # Filtrowanie śmieci przed normalizacją 
        invalid_types = {"Dup", "NonEx"}
        df = df[~df["Type"].isin(invalid_types)].copy()
        
        # RA/Dec z notacji seksagesymalnej
        df["ra"] = df["RA"].apply(lambda x: sexa_to_deg(x, is_ra=True))
        df["dec"] = df["Dec"].apply(lambda x: sexa_to_deg(x, is_ra=False))
        
        # Magnitude 
        mags = []
        for _, row in df.iterrows():
            m = get_magnitude(
                row.get("V-Mag"), row.get("B-Mag"), row.get("J-Mag"),
                row.get("H-Mag"), row.get("K-Mag"), row.get("SurfBr"),
                str(row.get("Type", ""))
            )
            mags.append(m)
        df["mag"] = mags
        
        df["id"] = df["Name"]
        df["size"] = df["MajAx"].apply(safe_float)
        df["type"] = df["Type"]
        df["common_names"] = df["Common names"]
        def _add_ngc_extra_common_names(row):
            """Dodaje NGC_EXTRA_COMMON_NAMES do istniejących common names"""
            orig = row['Common names'] or ''
            extra = NGC_EXTRA_COMMON_NAMES.get(row['Name'].strip(), '')
            if extra:
                return f"{orig}, {extra}".strip(', ')
            return orig
        
        df['common_names'] = df.apply(_add_ngc_extra_common_names, axis=1)
        # Kolumny famous (Messier, Herschel, Caldwell)
        df = build_famous_columns(df)
        
        df['catalog'] = 'ngc'
        df['catalog_priority'] = CATALOG_PRIORITY['ngc']
        results.append(df)
        
        # 2. SHARPLESS
        df = raw['sh2'].copy()
        if not df.empty:
            print("       > Normalizacja Sharpless")
            df = convert_ra_dec(df, '_RAJ2000', '_DEJ2000', 'deg')
            df['id'] = df['Sh2'].apply(lambda x: f"Sh2-{int(x)}" if pd.notnull(x) else "")
            df['size'] = pd.to_numeric(df['Diam'], errors='coerce')
            df['mag'], df['type'], df['extra_info'] = np.nan, 'HII', ""
            df['common_names'] = df['Sh2'].apply(
                lambda x: SH2_COMMON_NAMES.get(f"Sh2-{int(x)}", "") if pd.notnull(x) else ""
            )
            df['catalog'] = 'sh2'
            df['catalog_priority'] = CATALOG_PRIORITY['sh2']
            results.append(df)
        
        # 3. BARNARD
        df = raw['barn'].copy()
        if not df.empty:
            print("       > Normalizacja Barnard")
            df = convert_ra_dec(df, '_RAJ2000', '_DEJ2000', 'deg')
            df['id'] = df['Barn'].apply(lambda x: f"B{str(x).strip()}")
            df['size'] = pd.to_numeric(df['Diam'], errors='coerce')
            df['mag'], df['type'], df['extra_info'] = np.nan, 'DN', ""
            df['common_names'] = ""
            df['catalog'] = 'barn'
            df['catalog_priority'] = CATALOG_PRIORITY['barn']
            results.append(df)
        
        # 4. RCW
        df = raw['rcw'].copy()
        if not df.empty:
            print("       > Normalizacja RCW")
            df = convert_ra_dec(df, '_RAJ2000', '_DEJ2000', 'deg')
            df['id'] = df['RCW'].apply(lambda x: f"RCW{str(x).strip()}")
            df['size'] = pd.to_numeric(df['MajAxis'], errors='coerce')
            df['mag'], df['type'], df['extra_info'] = np.nan, 'HII', ""
            df['common_names'] = ""
            df['catalog'] = 'rcw'
            df['catalog_priority'] = CATALOG_PRIORITY['rcw']
            results.append(df)
        
        # 5. CEDERBLAD
        df = raw['ced'].copy()
        if not df.empty:
            print("       > Normalizacja Cederblad")
            df = convert_ra_dec(df, '_RAJ2000', '_DEJ2000', 'deg')
            ids = []
            for _, row in df.iterrows():
                num = str(row['Ced']).strip()
                sub = str(row.get('m_Ced', "")).strip()
                sub = sub if pd.notnull(row.get('m_Ced')) and sub != 'nan' else ""
                ids.append(f"Ced{num}{sub}")
            df['id'] = ids
            
            d1 = pd.to_numeric(df['Dim1'], errors='coerce')
            d2 = pd.to_numeric(df['Dim2'], errors='coerce')
            df['size'] = np.where((d2.notna()) & (d2 > 0), (d1 + d2) / 2, d1)
            df['mag'] = pd.to_numeric(df['vmag'], errors='coerce')
            df['type'], df['extra_info'] = 'NB', ""
            df['common_names'] = ""
            df['catalog'] = 'ced'
            df['catalog_priority'] = CATALOG_PRIORITY['ced']
            results.append(df)
        
        # 6. LBN
        df = raw['lbn'].copy()
        if not df.empty:
            print("       > Normalizacja LBN")
            df = convert_ra_dec(df, '_RAJ2000', '_DEJ2000', 'deg')
            df['id'] = df['Seq'].apply(lambda x: f"LBN{x}")
            df['size'] = pd.to_numeric(df['Diam1'], errors='coerce')
            df['mag'], df['type'], df['extra_info'] = np.nan, 'NB', ""
            df['common_names'] = ""
            df['catalog'] = 'lbn'
            df['catalog_priority'] = CATALOG_PRIORITY['lbn']
            results.append(df)
        
        # 7. LDN
        df = raw['ldn'].copy()
        if not df.empty:
            print("       > Normalizacja LDN")
            df = convert_ra_dec(df, '_RAJ2000', '_DEJ2000', 'deg')
            df['id'] = df['LDN'].apply(lambda x: f"LDN{x}")
            df['size'] = np.sqrt(pd.to_numeric(df['Area'], errors='coerce')) * 60
            df['mag'], df['type'], df['extra_info'] = np.nan, 'DN', ""
            df['common_names'] = ""
            df['catalog'] = 'ldn'
            df['catalog_priority'] = CATALOG_PRIORITY['ldn']
            results.append(df)
        
        # 8. PGC 
        df = raw['pgc'].copy()
        if not df.empty:
            print("       > Normalizacja PGC")
            df = convert_ra_dec(df, '_RAJ2000', '_DEJ2000', 'deg')
            df['id'] = df['PGC'].apply(lambda x: f"PGC{int(x)}" if pd.notnull(x) else "")
            # UWAGA: MajAxis w PGC VII/119 jest W ARCMIN (wg VizieR ReadMe)
            df['size'] = pd.to_numeric(df['MajAxis'], errors='coerce')
            df['mag'] = pd.to_numeric(df['Btot'], errors='coerce')
            df['type'], df['extra_info'] = 'G', ""
            df['common_names'] = ""
            df['catalog'] = 'pgc'
            df['catalog_priority'] = CATALOG_PRIORITY['pgc']
            results.append(df)
            
        # SKLEJANIE
        final_cols = [
            'id', 'ra', 'dec', 'size', 'mag', 'type', 'extra_info', 'common_names',
            'catalog', 'catalog_priority', 'messier_nr', 'caldwell_nr', 'herschel_nr'
        ]
        
        output = []
        for d in results:
            for c in final_cols:
                if c not in d.columns:
                    # Jeśli kolumny int sławy nie ma, wypełnij 0
                    if c.endswith('_nr'):
                        d[c] = 0
                    else:
                        d[c] = np.nan
            output.append(d[final_cols])
            
        unified = pd.concat(output, ignore_index=True)
        
        # Podstawowa konwersja na typ numeryczny dla pewności
        for c in ['messier_nr', 'caldwell_nr', 'herschel_nr']:
            unified[c] = pd.to_numeric(unified[c], errors='coerce').fillna(0).astype(int)
            
        return unified

#  === 3. MERGE ===
def show_merge_statistics(df: pd.DataFrame) -> float:
    """
    Liczy statystyki par dla tolerancji od 1' do 20'.
    Na tej podstawie automatycznie wybiera tolerancję (naturalny skok).
    
    Drukuje tabelkę, żeby decyzja była audytowalna.
    Zwraca wybraną tolerancję w stopniach.
    """
    print_step("Analiza tolerancji dla Smart Merge")
    
    coords = SkyCoord(
        ra=df['ra'].values * u.deg,
        dec=df['dec'].values * u.deg,
    )
    
    # Progi do przetestowania (w arcminach)
    tolerances_arcmin = [1, 2, 3, 5, 8, 12, 20]
    
    results = []
    for tol in tolerances_arcmin:
        idx1, idx2, _, _ = search_around_sky(
            coords, coords, (tol / 60.0) * u.deg
        )
        # Pary cross-katalogowe (różne katalogi, różne ID)
        mask_cross = (idx1 != idx2) & (
            df['catalog'].values[idx1] != df['catalog'].values[idx2]
        )
        n_cross = mask_cross.sum() // 2  # //2 bo pary symetryczne
        
        # Buduj graf i policz klastry
        g = nx.Graph()
        g.add_nodes_from(range(len(df)))
        for i1, i2 in zip(idx1[mask_cross], idx2[mask_cross]):
            g.add_edge(int(i1), int(i2))
        clusters = list(nx.connected_components(g))
        n_multi = sum(1 for c in clusters if len(c) > 1)
        
        # Klastry z "podejrzanym" mieszaniem typów:
        # np. emisja (HII/NB/PN/SNR) + ciemna mgławica (DN)
        EMISSION = {"HII", "NB", "PN", "SNR"}
        DARK     = {"DN"}
        n_suspicious = 0
        for cluster in clusters:
            if len(cluster) < 2:
                continue
            types_in = set(df['type'].values[list(cluster)])
            if EMISSION & types_in and DARK & types_in:
                n_suspicious += 1
        
        results.append({
            'tol': tol,
            'pairs': n_cross,
            'clusters': n_multi,
            'suspicious': n_suspicious,
        })
    
    # Tabelka
    print(f"\n       {'Tol':>5} | {'Pary X-kat':>10} | {'Klastry':>8} | {'Podejrzane':>10}")
    print("       " + "-" * 42)
    for r in results:
        print(
            f"       {r['tol']:>4}' | "
            f"{fmt(r['pairs']):>10} | "
            f"{fmt(r['clusters']):>8} | "
            f"{fmt(r['suspicious']):>10}"
        )
    
    # Automatyczny wybór: pierwszy skok > 50% wzrostu liczby podejrzanych
    # względem poprzedniej wartości — cofamy się o jeden krok
    chosen = tolerances_arcmin[1]  # fallback: 2'
    for i in range(1, len(results)):
        prev = results[i - 1]['suspicious']
        curr = results[i]['suspicious']
        if prev > 0 and curr / prev > 1.5:
            chosen = tolerances_arcmin[i - 1]
            break
    else:
        # Jeśli nie ma wyraźnego skoku, bierz 3' (sprawdzona wartość z 1_)
        chosen = 3
    
    print(f"\n       [AUTO] Wybrana tolerancja: {chosen}' (arcmin)")
    print(f"       [AUTO] Tolerancja w stopniach: {chosen / 60.0:.5f}°")
    return chosen / 60.0

def smart_merge(df: pd.DataFrame, tolerance_deg: float) -> pd.DataFrame:
    """
    Entity resolution przez graf sąsiedztwa + connected components.
    Lider klastra wybierany przez catalog_priority (wyższy = ważniejszy),
    potem rozmiar (większy = lepszy).
    """
    df = df.dropna(subset=['ra', 'dec']).reset_index(drop=True)
    print_step(
        f"Smart Merge (tolerancja {tolerance_deg * 60:.2f} arcmin, "
        f"{fmt(len(df))} obiektów)"
    )
    
    coords = SkyCoord(ra=df['ra'].values * u.deg, dec=df['dec'].values * u.deg)
    idx1, idx2, _, _ = search_around_sky(coords, coords, tolerance_deg * u.deg)
    
    g = nx.Graph()
    g.add_nodes_from(range(len(df)))
    for i1, i2 in zip(idx1, idx2):
        if i1 != i2:
            g.add_edge(i1, i2)
    
    clusters = list(nx.connected_components(g))
    n_single  = sum(1 for c in clusters if len(c) == 1)
    n_multi   = sum(1 for c in clusters if len(c) > 1)
    print(f"       Klastry jednoobiektowe : {fmt(n_single)}")
    print(f"       Klastry wieloobiektowe : {fmt(n_multi)}")
    print(f"       Łącznie klastrów       : {fmt(len(clusters))}")

    merged_rows = []
    
    for cluster in tqdm(clusters,
                    desc="       Merge klastrów",
                    unit="klaster",
                    colour="green",
                    ncols=119):
        subset = df.iloc[list(cluster)]
        
        # CATALOG_PRIORITY: wyższy = ważniejszy → ascending=[False, False]
        subset_sorted = subset.sort_values(
            by=['catalog_priority', 'size'],
            ascending=[False, False],
            na_position='last',
        )
        master = subset_sorted.iloc[0].copy()
        
        # Współrzędne — od lidera
        master['ra']  = subset_sorted.iloc[0]['ra']
        master['dec'] = subset_sorted.iloc[0]['dec']
        
        # Rozmiar — max z klastra z limitem
        sizes = subset['size'].dropna()
        if not sizes.empty:
            reasonable = sizes[sizes <= MAX_OBJECT_SIZE_ARCMIN]
            master['size'] = reasonable.max() if not reasonable.empty \
                             else MAX_OBJECT_SIZE_ARCMIN
        
        # Mag — tylko od lidera (nie bierzemy od sąsiadów — to byłoby błędne)
        # master['mag'] pozostaje niezmieniony (już jest od lidera)
        
        # extra_info: TYLKO IDs wchłoniętych (bez M/C/H — mają własne kolumny)
        def clean_ids(series: pd.Series) -> set:
            unique = set()
            for val in series.dropna().astype(str):
                for p in val.split(','):
                    p = p.strip()
                    if p and p.lower() != 'nan':
                        unique.add(p)
            return unique
        
        master_id = str(master['id']).strip()
        all_ids   = clean_ids(subset['id'])
        all_ids.discard(master_id)
        existing  = clean_ids(subset['extra_info'])
        master['extra_info'] = ", ".join(sorted(all_ids | existing))
        
        # common_names — scalanie z logiką usuwania substring
        master['common_names'] = process_common_names(subset['common_names'])
        
        # Flagi sławy — bierzemy MAX z klastra, nie tylko od lidera
        # Jeśli Sh2-49 wchłonie NGC6611 (M16), to M16 ma messier_nr=16
        # Klaster powinien odziedziczyć to po wchłoniętym NGC, nie je stracić
        for col in ['messier_nr', 'caldwell_nr', 'herschel_nr']:
            max_val = subset[col].max()
            master[col] = int(max_val) if pd.notnull(max_val) else 0
        
        merged_rows.append(master.to_dict())
    
    result = pd.DataFrame(merged_rows)
    print(f"\n       Po merge: {fmt(len(result))} obiektów.")
    return result

def merge_by_common_names(df: pd.DataFrame) -> pd.DataFrame:
    """
    Fuzja obiektów o tej samej nazwie zwyczajowej (common_names).

    Obiekty z różnych katalogów opisujące ten sam fizyczny obiekt
    (np. NGC1499 i Sh2-220 jako "California Nebula") zostają złączone
    w jeden rekord przed klastrowaniem przestrzennym.

    Zasady fuzji:
    - Lider:       najwyższy catalog_priority (NGC > IC > Sh2 > ...)
    - size:          max ze wszystkich wierszy w grupie
    - mag:         min (najjaśniejsza wartość, ignorując NaN)
    - extra_info:  unia wszystkich ID ze wszystkich wierszy grupy
    - common_names: unia nazw przez process_common_names()
    - messier/caldwell/herschel_nr: max (OR logiczne)
    - Współrzędne i typ: z lidera
    """
    if df.empty:
        return df

    df = df.copy().reset_index(drop=True)

    def normalize_name(name: str) -> str:
        return re.sub(r'\s+', ' ', name.strip().lower())

    def get_names_set(val) -> set:
        if pd.isna(val) or str(val).strip().lower() in ('nan', 'none', ''):
            return set()
        return {
            normalize_name(n) for n in str(val).split(',')
            if len(normalize_name(n)) >= 5  # Ochrona przed krótkimi generycznymi nazwami
        }

    def clean_ids(series: pd.Series) -> set:
        unique = set()
        for val in series.dropna().astype(str):
            for p in val.split(','):
                p = p.strip()
                if p and p.lower() not in ('nan', 'none', ''):
                    unique.add(p)
        return unique

    # Buduj graf: dwa wiersze → krawędź, jeśli mają wspólną nazwę zwyczajową
    g = nx.Graph()
    g.add_nodes_from(range(len(df)))

    name_to_indices: dict = {}
    for i, row in df.iterrows():
        for name in get_names_set(row.get('common_names')):
            name_to_indices.setdefault(name, []).append(i)

    for name, indices in name_to_indices.items():
        if len(indices) > 1:
            for j in range(len(indices)):
                for k in range(j + 1, len(indices)):
                    g.add_edge(indices[j], indices[k])

    clusters = list(nx.connected_components(g))
    n_multi = sum(1 for c in clusters if len(c) > 1)
    print_step(f"Fuzja common_names: {n_multi} grup do scalenia")

    merged_rows = []

    for cluster in clusters:
        subset = df.iloc[sorted(cluster)]

        if len(subset) == 1:
            merged_rows.append(subset.iloc[0].to_dict())
            continue

        # Lider: najwyższy catalog_priority, przy remisie największy size
        subset_sorted = subset.sort_values(
            by=['catalog_priority', 'size'],
            ascending=[False, False],
            na_position='last'
        )
        master = subset_sorted.iloc[0].copy()

        # extra_info: unia wszystkich ID (bez ID lidera)
        all_ids = clean_ids(subset['id'])
        master_id = str(master.get('id', '')).strip()
        all_ids.discard(master_id)
        existing_extra = clean_ids(subset['extra_info'])
        master['extra_info'] = ','.join(sorted(all_ids | existing_extra))

        # size: max
        sizes = pd.to_numeric(subset['size'], errors='coerce').dropna()
        if not sizes.empty:
            master['size'] = float(sizes.max())

        # mag: min (najjaśniejsza), ignorując NaN
        mags = pd.to_numeric(subset['mag'], errors='coerce').dropna()
        if not mags.empty:
            master['mag'] = float(mags.min())

        # Sława: OR logiczne (max)
        for col in ['messier_nr', 'caldwell_nr', 'herschel_nr']:
            if col in subset.columns:
                vals = pd.to_numeric(subset[col], errors='coerce').fillna(0)
                master[col] = int(vals.max())

        # Nazwy zwyczajowe: unia przez istniejącą funkcję
        master['common_names'] = process_common_names(subset['common_names'])

        merged_rows.append(master.to_dict())

    result = pd.DataFrame(merged_rows).reset_index(drop=True)
    n_before, n_after = len(df), len(result)
    print(f"       {n_before} → {n_after} obiektów (scalono {n_before - n_after} duplikatów nazewniczych)")
    return result

# === 4. UZUPEŁNIANIE DANYCH ===
def impute_size_and_mag(df: pd.DataFrame) -> pd.DataFrame:
    """
    Imputacja brakujących wartości Size i Mag na podstawie mediany
    wg typu obiektu.
    
    Zasady (omówione w analizie):
      - size: imputujemy dla wszystkich typów (geometria)
      - mag:  imputujemy TYLKO dla typów emisyjnych (HII, NB, PN, SNR, G itd.)
              NIE imputujemy dla DN (ciemne mgławice — mag nie ma sensu fizycznego)
    
    Dodawane kolumny flag:
      - size_imputed (bool)
      - mag_imputed  (bool)
    """
    df = df.copy()
    df['size_imputed'] = False
    df['mag_imputed']  = False
    
    # Typy, dla których mag NIE ma sensu fizycznego
    NO_MAG_TYPES = {"DN"}
    
    # Statystyki do logowania
    size_filled = 0
    mag_filled  = 0
    
    # Wylicz mediany per typ
    size_medians = df.groupby('type')['size'].median()
    mag_medians  = df.groupby('type')['mag'].median()
    
    print_step("Imputacja brakujących Size i Mag (mediana wg typu)")
    print(f"\n       {'Typ':>8} | {'Mediana size':>12} | {'Mediana mag':>11}")
    print("       " + "-" * 36)
    for t in sorted(df['type'].dropna().unique()):
        s = size_medians.get(t, float('nan'))
        m = mag_medians.get(t, float('nan'))
        s_str = f"{s:.1f}'" if not np.isnan(s) else "   brak"
        m_str = f"{m:.2f} " if not np.isnan(m) else "   brak"
        print(f"       {t:>8} | {s_str:>12} | {m_str:>11}")
    
    for idx, row in df.iterrows():
        obj_type = str(row.get('type', '')).strip().upper()
        
        # --- SIZE ---
        if pd.isna(row['size']) or row['size'] <= 0:
            median_s = size_medians.get(row['type'])
            if median_s is not None and not np.isnan(median_s):
                df.at[idx, 'size'] = round(median_s, 2)
                df.at[idx, 'size_imputed'] = True
                size_filled += 1
        
        # --- MAG ---
        if obj_type in NO_MAG_TYPES:
            # Ciemne mgławice — nie imputujemy, zostawiamy NaN świadomie
            continue
        if pd.isna(row['mag']):
            median_m = mag_medians.get(row['type'])
            if median_m is not None and not np.isnan(median_m):
                df.at[idx, 'mag'] = round(median_m, 2)
                df.at[idx, 'mag_imputed'] = True
                mag_filled += 1
    
    print(f"\n       Uzupełniono size : {fmt(size_filled)} obiektów")
    print(f"       Uzupełniono mag  : {fmt(mag_filled)} obiektów")
    print(f"       Bez imputacji mag (DN) : pozostają z NaN")
    return df

# === 5. FILTR KOŃCOWY ===
def apply_final_filters(df: pd.DataFrame, min_size: float, max_mag: float) -> pd.DataFrame:
    """
    Ostatni filtr przed zapisem gotowego katalogu.
    Odrzuca obiekty zbyt małe (rozmiar < min_size) lub zbyt ciemne (mag > max_mag).
    
    UWAGA DLA MAGNITUDO:
    Dla ciemnych mgławic (DN) mag jest NaN (zgodnie z decyzją o braku imputacji).
    Warunek (mag <= max_mag) odrzuciłby wszystkie DN!
    Dlatego warunek to: (mag is NaN) LUB (mag <= max_mag).
    """
    print_step(f"Filtr docelowy (Size >= {min_size}', Mag <= {max_mag})")
    
    # 1. Filtr rozmiaru: zachowaj NaN (choć imputacja powinna je wyeliminować) lub >= min_size
    mask_size = df['size'].isna() | (df['size'] >= min_size)
    df_filtered = df[mask_size].copy()
    dropped_size = len(df) - len(df_filtered)
    
    # 2. Filtr jasności: zachowaj NaN (niezbędne dla DN) lub <= max_mag
    mask_mag = df_filtered['mag'].isna() | (df_filtered['mag'] <= max_mag)
    df_final = df_filtered[mask_mag].copy()
    dropped_mag = len(df_filtered) - len(df_final)
    
    print(f"       Odrzucono ze względu na mały rozmiar : {fmt(dropped_size)}")
    print(f"       Odrzucono ze względu na słabą jasność: {fmt(dropped_mag)}")
    
    return df_final

# === 6. GŁÓWNY PUNKT WEJŚCIA ===
def main():
   
    # 1. Pobranie z VizieR / lokalnego NGC.csv
    print_green(f"\n{'=' * 34} KROK 1: Pobieranie bazy danych (ENTITY RESOLUTION) {'=' * 33}")
    raw = fetch_data()
    
    # 2. Normalizacja jednostek, nazw, dodanie kolumn flag (M/H/C)
    print_green(f"{'=' * 49} KROK 2: Normalizacja {'=' * 48}")
    full = normalize_all(raw)
    
    # 3. Filtr śmieci (przed merge)
    print_green(f"{'=' * 49} KROK 3: Filtrowanie {'=' * 49}")
    mask_trash = (
        (full['mag'].notna() & (full['mag'] > TRASH_MAG)) |
        (full['size'].notna() & (full['size'] < TRASH_SIZE))
    )
    df_clean = full[~mask_trash].copy()
    print(
        f"       Odrzucono {fmt(mask_trash.sum())} obiektów "
        f"(Mag > {TRASH_MAG} lub Size < {TRASH_SIZE}')."
    )
    print(f"       Do analizy merge trafia: {fmt(len(df_clean))} obiektów.")
    
    # 4. Automatyczna ocena tolerancji i Smart Merge
    best_tol_deg = show_merge_statistics(df_clean)
    final_merged = smart_merge(df_clean, best_tol_deg)
    final_merged = merge_by_common_names(final_merged)
    # 5. Imputacja brakujących Size i Mag (tylko na połączonych obiektach)
    print_green(f"{'=' * 42} KROK 4: Imputacja danych mag/size {'=' * 42}")
    final_imputed = impute_size_and_mag(final_merged)
    
    # 6. Filtr końcowy (zastępuje: input("Minimalny rozmiar..."))
    print_green(f"{'=' * 49} KROK 5: Filtrowanie {'=' * 49}")
    final_filtered = apply_final_filters(
        final_imputed, 
        min_size=MIN_SIZE_ARCMIN, 
        max_mag=MAX_MAG
    )
    
    # 7. Zapis do CSV
    # Porządek kolumn w pliku wyjściowym
    print_green(f"{'=' * 49} KROK 2: Zapis danych {'=' * 48}")
    output_cols = [
        'id', 'ra', 'dec', 'size', 'size_imputed', 'mag', 'mag_imputed', 
        'type', 'catalog', 'catalog_priority', 
        'messier_nr', 'caldwell_nr', 'herschel_nr', 
        'extra_info', 'common_names'
    ]
    
    # Upewniamy się, że nie zapisujemy niepotrzebnych kolumn
    final_filtered = final_filtered[output_cols]
    
    # Sortujemy przed zapisem: najpierw Messier, potem Caldwell, Herschel, a reszta po priorytecie
    # To sprawi, że w pliku najciekawsze obiekty będą na początku
    final_filtered = final_filtered.sort_values(
        by=['messier_nr', 'caldwell_nr', 'herschel_nr', 'catalog_priority', 'size'],
        ascending=[False, False, False, False, False]
    )
    
    print_green(f"[INFO] ZAPIS DO PLIKU: {PATHS.catalog_full}")
    final_filtered.to_csv(PATHS.catalog_full, index=False)
    print(f"       Sukces! Zapisano {fmt(len(final_filtered))} obiektów.")
    
    answer = input('\n\n[INFO] Uruchomić pełną analizę katalogu? [y/n] (domyślnie y): ').strip().lower()
    if answer in ["y", ""]:
        analyze_catalog()
    else:
        print_step("Pominięto analizę.")

if __name__ == "__main__":
    main()

