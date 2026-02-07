#!/usr/bin/env python3

import csv
from pathlib import Path

BASE_DIR = Path(__file__).resolve().parent
INPUT_FILE = BASE_DIR / "OpenNGC" / "NGC.csv"
OUTPUT_FILE = BASE_DIR / "updated_ngc.csv"

# ==================== LOOKUPS z lookups.py ====================

# Lista Herschel 400 po numerach NGC (poprawiona wersja)
herschel_set = {"NGC0040": 1, "NGC0129": 2, "NGC0136": 3, "NGC0157": 4, "NGC0185": 5, "NGC0205": 6, "NGC0225": 7, "NGC0246": 8, "NGC0247": 9, "NGC0253": 10, "NGC0278": 11, "NGC0288": 12, "NGC0381": 13, "NGC0404": 14, "NGC0436": 15, "NGC0457": 16, "NGC0488": 17, "NGC0524": 18, "NGC0559": 19, "NGC0584": 20, "NGC0596": 21, "NGC0598": 22, "NGC0613": 23, "NGC0615": 24, "NGC0637": 25, "NGC0651": 26, "NGC0654": 27, "NGC0659": 28, "NGC0663": 29, "NGC0720": 30, "NGC0752": 31, "NGC0772": 32, "NGC779": 33, "NGC869": 34, "NGC884": 35, "NGC891": 36, "NGC0908": 37, "NGC0936": 38, "NGC1022": 39, "NGC1023": 40, "NGC1027": 41, "NGC1052": 42, "NGC1055": 43, "NGC1084": 44, "NGC1245": 45, "NGC1342": 46, "NGC1407": 47, "NGC1444": 48, "NGC1501": 49, "NGC1502": 50, "NGC1513": 51, "NGC1528": 52, "NGC1535": 53, "NGC1545": 54, "NGC1647": 55, "NGC1664": 56, "NGC1788": 57, "NGC1817": 58, "NGC1857": 59, "NGC1907": 60, "NGC1931": 61, "NGC1961": 62, "NGC1964": 63, "NGC1980": 64, "NGC1999": 65, "NGC2022": 66, "NGC2024": 67, "NGC2126": 68, "NGC2129": 69, "NGC2158": 70, "NGC2169": 71, "NGC2185": 72, "NGC2186": 73, "NGC2194": 74, "NGC2204": 75, "NGC2215": 76, "NGC2232": 77, "NGC2244": 78, "NGC2251": 79, "NGC2264": 80, "NGC2266": 81, "NGC2281": 82, "NGC2286": 83, "NGC2301": 84, "NGC2304": 85, "NGC2311": 86, "NGC2324": 87, "NGC2335": 88, "NGC2343": 89, "NGC2353": 90, "NGC2354": 91, "NGC2355": 92, "NGC2360": 93, "NGC2362": 94, "NGC2371": 95, "NGC2372": 96, "NGC2392": 97, "NGC2395": 98, "NGC2403": 99, "NGC2419": 100, "NGC2420": 101, "NGC2421": 102, "NGC2422": 103, "NGC2423": 104, "NGC2438": 105, "NGC2440": 106, "NGC2479": 107, "NGC2482": 108, "NGC2489": 109, "NGC2506": 110, "NGC2509": 111, "NGC2527": 112, "NGC2539": 113, "NGC2548": 114, "NGC2567": 115, "NGC2571": 116, "NGC2613": 117, "NGC2627": 118, "NGC2655": 119, "NGC2681": 120, "NGC2683": 121, "NGC2742": 122, "NGC2768": 123, "NGC2775": 124, "NGC2782": 125, "NGC2787": 126, "NGC2811": 127, "NGC2841": 128, "NGC2859": 129, "NGC2903": 130, "NGC2950": 131, "NGC2964": 132, "NGC2974": 133, "NGC2976": 134, "NGC2985": 135, "NGC3034": 136, "NGC3077": 137, "NGC3079": 138, "NGC3115": 139, "NGC3147": 140, "NGC3166": 141, "NGC3169": 142, "NGC3184": 143, "NGC3190": 144, "NGC3193": 145, "NGC3198": 146, "NGC3226": 147, "NGC3227": 148, "NGC3242": 149, "NGC3245": 150, "NGC3277": 151, "NGC3294": 152, "NGC3310": 153, "NGC3344": 154, "NGC3377": 155, "NGC3379": 156, "NGC3384": 157, "NGC3395": 158, "NGC3412": 159, "NGC3414": 160, "NGC3432": 161, "NGC3486": 162, "NGC3489": 163, "NGC3504": 164, "NGC3521": 165, "NGC3556": 166, "NGC3593": 167, "NGC3607": 168, "NGC3608": 169, "NGC3610": 170, "NGC3613": 171, "NGC3619": 172, "NGC3621": 173, "NGC3626": 174, "NGC3628": 175, "NGC3631": 176, "NGC3640": 177, "NGC3655": 178, "NGC3665": 179, "NGC3675": 180, "NGC3686": 181, "NGC3726": 182, "NGC3729": 183, "NGC3810": 184, "NGC3813": 185, "NGC3877": 186, "NGC3893": 187, "NGC3898": 188, "NGC3900": 189, "NGC3912": 190, "NGC3938": 191, "NGC3941": 192, "NGC3945": 193, "NGC3949": 194, "NGC3953": 195, "NGC3962": 196, "NGC3982": 197, "NGC3992": 198, "NGC3998": 199, "NGC4026": 200, "NGC4027": 201, "NGC4030": 202, "NGC4036": 203, "NGC4038": 204, "NGC4041": 205, "NGC4051": 206, "NGC4085": 207, "NGC4088": 208, "NGC4102": 209, "NGC4111": 210, "NGC4143": 211, "NGC4147": 212, "NGC4150": 213, "NGC4151": 214, "NGC4179": 215, "NGC4203": 216, "NGC4214": 217, "NGC4216": 218, "NGC4245": 219, "NGC4251": 220, "NGC4258": 221, "NGC4261": 222, "NGC4273": 223, "NGC4274": 224, "NGC4278": 225, "NGC4281": 226, "NGC4293": 227, "NGC4303": 228, "NGC4314": 229, "NGC4346": 230, "NGC4350": 231, "NGC4361": 232, "NGC4365": 233, "NGC4371": 234, "NGC4394": 235, "NGC4414": 236, "NGC4419": 237, "NGC4429": 238, "NGC4435": 239, "NGC4438": 240, "NGC4442": 241, "NGC4448": 242, "NGC4449": 243, "NGC4450": 244, "NGC4459": 245, "NGC4473": 246, "NGC4477": 247, "NGC4478": 248, "NGC4485": 249, "NGC4490": 250, "NGC4494": 251, "NGC4526": 252, "NGC4527": 253, "NGC4535": 254, "NGC4536": 255, "NGC4546": 256, "NGC4548": 257, "NGC4550": 258, "NGC4559": 259, "NGC4565": 260, "NGC4570": 261, "NGC4594": 262, "NGC4596": 263, "NGC4618": 264, "NGC4631": 265, "NGC4636": 266, "NGC4643": 267, "NGC4654": 268, "NGC4656": 269, "NGC4660": 270, "NGC4665": 271, "NGC4666": 272, "NGC4689": 273, "NGC4697": 274, "NGC4698": 275, "NGC4699": 276, "NGC4725": 277, "NGC4753": 278, "NGC4754": 279, "NGC4762": 280, "NGC4781": 281, "NGC4800": 282, "NGC4845": 283, "NGC4856": 284, "NGC4866": 285, "NGC4900": 286, "NGC4958": 287, "NGC4995": 288, "NGC5005": 289, "NGC5033": 290, "NGC5054": 291, "NGC5195": 292, "NGC5248": 293, "NGC5273": 294, "NGC5322": 295, "NGC5363": 296, "NGC5364": 297, "NGC5466": 298, "NGC5473": 299, "NGC5474": 300, "NGC5557": 301, "NGC5566": 302, "NGC5576": 303, "NGC5631": 304, "NGC5634": 305, "NGC5676": 306, "NGC5689": 307, "NGC5694": 308, "NGC5746": 309, "NGC5846": 310, "NGC5866": 311, "NGC5897": 312, "NGC5907": 313, "NGC5982": 314, "NGC6118": 315, "NGC6144": 316, "NGC6171": 317, "NGC6207": 318, "NGC6217": 319, "NGC6229": 320, "NGC6235": 321, "NGC6284": 322, "NGC6287": 323, "NGC6293": 324, "NGC6304": 325, "NGC6316": 326, "NGC6342": 327, "NGC6355": 328, "NGC6356": 329, "NGC6369": 330, "NGC6401": 331, "NGC6426": 332, "NGC6440": 333, "NGC6445": 334, "NGC6451": 335, "NGC6514": 336, "NGC6517": 337, "NGC6520": 338, "NGC6522": 339, "NGC6528": 340, "NGC6540": 341, "NGC6543": 342, "NGC6544": 343, "NGC6553": 344, "NGC6568": 345, "NGC6569": 346, "NGC6583": 347, "NGC6624": 348, "NGC6629": 349, "NGC6633": 350, "NGC6638": 351, "NGC6642": 352, "NGC6645": 353, "NGC6664": 354, "NGC6712": 355, "NGC6755": 356, "NGC6756": 357, "NGC6781": 358, "NGC6802": 359, "NGC6818": 360, "NGC6823": 361, "NGC6826": 362, "NGC6830": 363, "NGC6834": 364, "NGC6866": 365, "NGC6882": 366, "NGC6885": 367, "NGC6905": 368, "NGC6910": 369, "NGC6934": 370, "NGC6939": 371, "NGC6940": 372, "NGC6946": 373, "NGC7000": 374, "NGC7006": 375, "NGC7008": 376, "NGC7009": 377, "NGC7044": 378, "NGC7062": 379, "NGC7086": 380, "NGC7128": 381, "NGC7142": 382, "NGC7160": 383, "NGC7209": 384, "NGC7217": 385, "NGC7243": 386, "NGC7296": 387, "NGC7331": 388, "NGC7380": 389, "NGC7448": 390, "NGC7479": 391, "NGC7510": 392, "NGC7606": 393, "NGC7662": 394, "NGC7686": 395, "NGC7723": 396, "NGC7727": 397, "NGC7789": 398, "NGC7790": 399, "NGC7814": 400}


# Caldwell: mapowanie pełnej nazwy NGC/IC na numer C
caldwell_dict = {
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

# ==================== FUNKCJE POMOCNICZE ====================

def sexa_to_deg(sexa, is_ra):
    """Konwersja HH:MM:SS.SS (RA) lub +/-DD:MM:SS.SS (Dec) na stopnie."""
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
        m = float(parts[1])
        sec = float(parts[2])
    except ValueError:
        return None

    val = h_d + m / 60.0 + sec / 3600.0

    if is_ra:
        val *= 15.0  # RA w godzinach -> stopnie
    else:
        val *= sign  # Dec ze znakiem

    return val


def safe_float(value):
    """Bezpieczna konwersja na float."""
    if not value or value.strip() == "":
        return None
    try:
        return float(value)
    except ValueError:
        return None


def get_magnitude(vmag, bmag, jmag, hmag, kmag, surfbr, obj_type):
    """
    Zwraca magnitude według priorytetu:
    1. V-Mag (zawsze priorytet)
    2. Dla galaktyk (Type='G'): SurfBr jeśli V-Mag puste
    3. Dla innych: najniższa wartość z B/J/H/K-Mag
    """
    v = safe_float(vmag)
    if v is not None:
        return v

    # Galaktyki - użyj SurfBr
    if obj_type == "G":
        sb = safe_float(surfbr)
        if sb is not None:
            return sb

    # Inne obiekty - wybierz najniższą magnitudę
    mags = []
    for mag_val in [bmag, jmag, hmag, kmag]:
        m = safe_float(mag_val)
        if m is not None:
            mags.append(m)

    if mags:
        return min(mags)  # Najniższa wartość = najjaśniejszy

    return None


def build_extra_info(name, m_id, ngc_id, ic_id):
    """
    Buduje extra_info z indeksów M, H, C
    """
    extras = []

    # Messier
    if m_id and m_id.strip():
        extras.append(f"M{int(m_id.strip())}")

    # Herschel 400
    # Sprawdź czy NGC jest w Herschel
    full_name = ""
    if ngc_id and ngc_id.strip():
        full_name = f"NGC{ngc_id.strip().zfill(4)}"
    else:
        full_name = name.strip()
    
    if full_name in herschel_set:
           h_num = herschel_set[full_name]
           extras.append(f"H{int(h_num)}")
    

    # Caldwell
    # Trzeba sprawdzić pełną nazwę NGC/IC
    full_name = ""
    if ngc_id and ngc_id.strip():
        full_name = f"NGC{ngc_id.strip().zfill(4)}"
    elif ic_id and ic_id.strip():
        full_name = f"IC{ic_id.strip().zfill(4)}"
    else:
        full_name = name.strip()

    if full_name in caldwell_dict:
        c_num = caldwell_dict[full_name]
        extras.append(f"C{int(c_num)}")

    return ",".join(extras)


# ==================== GŁÓWNA LOGIKA ====================

def process_catalog():
    """Przetwarza NGC.csv i tworzy unified catalog."""

    if not INPUT_FILE.exists():
        print(f"BŁĄD: Brak pliku {INPUT_FILE}")
        return

    output_rows = []
    skipped = 0

    with INPUT_FILE.open("r", encoding="utf-8") as f:
        reader = csv.DictReader(f, delimiter=";")

        for row in reader:
            name = row.get("Name", "").strip()
            if not name:
                skipped += 1
                continue

            obj_type = row.get("Type", "").strip()
            ra_sexa = row.get("RA", "").strip()
            dec_sexa = row.get("Dec", "").strip()

            # Konwersja RA/Dec
            ra_deg = sexa_to_deg(ra_sexa, is_ra=True)
            dec_deg = sexa_to_deg(dec_sexa, is_ra=False)

            # Magnitude
            vmag = row.get("V-Mag", "")
            bmag = row.get("B-Mag", "")
            jmag = row.get("J-Mag", "")
            hmag = row.get("H-Mag", "")
            kmag = row.get("K-Mag", "")
            surfbr = row.get("SurfBr", "")

            mag = get_magnitude(vmag, bmag, jmag, hmag, kmag, surfbr, obj_type)

            # Size (MajAx)
            size = safe_float(row.get("MajAx", ""))

            # Extra info
            m_id = row.get("M", "")
            ngc_id = row.get("NGC", "")
            ic_id = row.get("IC", "")
            extra_info = build_extra_info(name, m_id, ngc_id, ic_id)

            # Common names
            common_names = row.get("Common names", "").strip()

            # Dodaj do wynikowej listy
            output_rows.append({
                "id": name,
                "extra_info": extra_info,
                "type": obj_type,
                "ra": round(ra_deg, 5) if ra_deg is not None else "",
                "dec": round(dec_deg, 5) if dec_deg is not None else "",
                "mag": round(mag, 2) if mag is not None else "",
                "size": round(size, 1) if size is not None else "",
                "base_score": "",  # Puste pole
                "common_names": common_names,
            })

    # Zapisz do pliku
    if output_rows:
        with OUTPUT_FILE.open("w", encoding="utf-8", newline="") as f:
            fieldnames = ["id", "extra_info", "type", "ra", "dec", "mag", "size", "base_score", "common_names"]
            writer = csv.DictWriter(f, fieldnames=fieldnames)
            writer.writeheader()
            writer.writerows(output_rows)

        print(f"✓ Zapisano {len(output_rows)} obiektów do {OUTPUT_FILE}")
        if skipped > 0:
            print(f"  Pominięto {skipped} wierszy bez nazwy")
    else:
        print("BŁĄD: Brak danych do zapisania")


if __name__ == "__main__":
    process_catalog()
