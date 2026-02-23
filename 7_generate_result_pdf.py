#!/usr/bin/env python3
# =============================================================================
# 7_generate_result_pdf.py  –  Astrophotography Planner  |  Etap 7: Tworzenie Plannera w PDF
# =============================================================================
"""
Etap po wygenerowaniu dwóch składowych plików PDF.
1) wczytuje vis_data.json,
3) dodanie strony tytułowej i informacji, 
2) tworzenie pdf z Astrophotography Planner na dany rok.
"""

import json
from pypdf import PdfWriter, PdfReader
from reportlab.pdfgen import canvas
from reportlab.pdfbase import pdfmetrics
from reportlab.pdfbase.ttfonts import TTFont
from reportlab.pdfbase.pdfmetrics import registerFontFamily
from reportlab.lib.pagesizes import A4
from reportlab.lib.units import mm

import os
from datetime import date
import matplotlib
mpl_font_dir = os.path.join(matplotlib.get_data_path(), 'fonts', 'ttf')

# Import shared.py
from shared import (PATHS, 
                    fmt, print_step, print_green,
                    load_vis_data, save_vis_data, 
                    UserConfig, CameraConfig,
                    MONTH_NAMES_PRINT, TYPE_NAMES, CATALOG_PRIORITY, CAT_ORDER,
                    get_engine_raw_hash, get_engine_final_hash, 
                    load_observing_data, extract_famous_labels,)

vis_path = PATHS.vis_data
with open(vis_path, 'r') as f:
    dane = json.load(f)
user_config = UserConfig.from_vis_data(dane)

year_val = user_config.year
miasto = user_config.location.name
lat_lon = f"Lat.: {user_config.location.lat:.2f}°N, Lon.: {user_config.location.lon:.2f}°E"
miasto_wyd = miasto.split(',', 1)[0].strip()
dzisiaj = date.today()
dzien = dzisiaj.day
miesiac_pl = MONTH_NAMES_PRINT[dzisiaj.month]
rok = dzisiaj.year
miejsce_wydania = f"{miasto_wyd}\n{dzien} {miesiac_pl} {rok} roku"
cam = user_config.camera
focal_length = cam.lens_focal_length
fov_width, fov_height, _ = cam.calculate_fov()  # Dokładny wzór z shared [file:1]
sensor_info = f"{cam.sensor_width:.1f}×{cam.sensor_height:.1f} mm ({cam.sensor_cols}×{cam.sensor_rows} px, {cam.sensor_pitch:.2f} µm)"
min_alt = user_config.min_altitude
sun_limit = abs(user_config.sun_limit)
min_hours = user_config.min_hours
bortle = f"{user_config.bortle_range[0]}–{user_config.bortle_range[1]}"
fov_info = f"{fov_width:.2f}° × {fov_height:.2f}°"

# Policz liczbę obiektów - jeśli masz w JSON, użyj; tu przykład
objects_list = dane.get('objects', [])
total_objects = len([obj for obj in objects_list if obj.get('selected')])

# Ustawienia plików
tytul = "Astrophotography"
tytul1 = "Planner"
tytul2 = f"{dane['year']}"
output_title = "tytulowa.pdf"
output_info = "info_page.pdf"
pdf1 = f"Astrophotography_Planner_{year_val}_1.pdf"
pdf2 = f"Astrophotography_Planner_{year_val}_2.pdf"
output = f"Astrophotography_Planner_{year_val}_{miasto_wyd}.pdf"

# Rejestracja czcionek
# Rejestracja
pdfmetrics.registerFont(TTFont("DejaVu", os.path.join(mpl_font_dir, "DejaVuSans.ttf")))
pdfmetrics.registerFont(TTFont("DejaVu-Bold", os.path.join(mpl_font_dir, "DejaVuSans-Bold.ttf")))
pdfmetrics.registerFont(TTFont("DejaVu-Italic", os.path.join(mpl_font_dir, "DejaVuSans-Oblique.ttf")))
pdfmetrics.registerFont(TTFont("DejaVu-BoldItalic", os.path.join(mpl_font_dir, "DejaVuSans-BoldOblique.ttf")))

# Spięcie w rodzinę
registerFontFamily("DejaVu",
    normal="DejaVu",
    bold="DejaVu-Bold",
    italic="DejaVu-Italic",
    boldItalic="DejaVu-BoldItalic"
)

# ========== STRONA TYTUŁOWA ==========
c = canvas.Canvas(output_title, pagesize=A4)
width, height = A4
c.setTitle(tytul)

# Tytuł główny
c.setFont('DejaVu-Bold', 48)
c.setFillColorRGB(0.07, 0.07, 0.25)
c.drawCentredString(width / 2, height / 2 + 280, tytul)
c.drawCentredString(width / 2, height / 2 + 230, tytul1)
c.setFont('DejaVu-Bold', 34)
c.drawCentredString(width / 2, height / 2 + 150, tytul2)

# Miejsce wydania
c.setFillColor('black')
c.setFont('DejaVu', 10)
y_pos = 60
lines = miejsce_wydania.split('\n')
line_height = 12
for i, line in enumerate(lines):
    c.drawCentredString(width / 2, y_pos - i * line_height, line)

# Parametry - sekcja na dole strony
y_start = 390  # Pozycja startowa dla parametrów
x_left = 56.69
line_spacing = 16

c.setFillColorRGB(0.07, 0.07, 0.25)
y_current = y_start - line_spacing - 4

params_text = [
   ("Lokalizacja:", f"{miasto}"),
    ("", f"{lat_lon}"),  # cała linia italic (puste pole 1)
    ("Bortle:", f"{bortle}"),
    "",  # pusta linia
    ("Ogniskowa:", f"{focal_length:.0f} mm"),
    ("Sensor:", f"{sensor_info}"),
    ("Pole widzenia (FOV):", f"{fov_info}"),
    ("Filtry wąskopasmowe:", f"{'Tak' if user_config.has_narrowband else 'Nie'}"),
    "",
    ("Minimalna wysokość nad horyzontem:", f"{min_alt:.1f}°"),
    ("Zmierzch liczony jako Słońce poniżej", f"{sun_limit:.1f}°"),
    ("Minimalne okno obserwacji:", f"{min_hours:.1f}h"),
    ("Priorytet katalogów Messier, Caldwell, Herschel:", f"{'Tak' if user_config.prefer_famous else 'Nie'}"),
    "",
    ("Liczba obiektów w atlasie:", f"{total_objects}"),
]

for item in params_text:
    if item == "":  # pusta linia
        y_current -= line_spacing * 0.5
    else:
        label, value = item
        x_pos = x_left + 10
        
        if label:  # jeśli jest etykieta
            c.setFont('DejaVu', 10)
            c.drawString(x_pos, y_current, label)
            # oblicz szerokość etykiety + odstęp
            x_pos += pdfmetrics.stringWidth(label + " ", 'DejaVu', 10)
        
        # wartość italic
        c.setFont('DejaVu-BoldItalic', 10)
        c.drawString(x_pos, y_current, value)
        
        y_current -= line_spacing

c.save()

# ========== STRONA INFORMACYJNA (druga strona) ==========
c2 = canvas.Canvas(output_info, pagesize=A4)
c2.setTitle("Informacje")

# Nagłówek
c2.setFont('DejaVu-Bold', 16)
c2.setFillColorRGB(0.07, 0.07, 0.25)
c2.drawCentredString(width / 2, height - 40.03, "Informacje o atlasie")

# Tekst wyjaśniający
c2.setFont('DejaVu', 10)
c2.setFillColor('black')

info_text = """
Strony z wykresami miesięcznymi przedstawiają obiekty przypisane do wariantów na 
dany miesiąc. Prezentowana jest noc nowiu. Dlatego nie każdy obiekt na tych wykre-
sach spełnia warunek minimalnego okna obserwacji. Należy odszukać go na stronie 
jemu poświęconej i tam sprawdzić, która noc w danym miesiącu pozwala wykorzystać
okno obserwacji.

Dane typu: RA, Dec, size, mag są wynikiem automatycznego agregowania wielu obiek- 
tów z kadru. Dlatego mogą się nieco różnić od tego, co znajdziesz w innych atlasach. 
Różnica jednak nie powinna być znacząca.

Nazwa obiektu i nazwa zwyczajowa to również efekt łączenia różnych rekordów z ka-
talogów astronomicznych. Dlatego ich wzajemna adekwatność może czasem być 
umowna. Ale obiekt opisany przez jedną lub drugą nazwę z pewnością jest w kadrze.

Indeksy to obiekty, których wpisy zostały połączone w tym kadrze. Ich liczba jest 
zmienna, a tu została ograniczona do pierwszych 24. Ich kolejność odzwierciedla 
popularność katalogu. Stąd najpierw są NGC, IC, Sh2, a później pozostałe. 



Wszelkie sugestie zmian i poprawek mile widziane.
Podziel się swoją opinią: morus@dominikanie.pl. 
"""

# Rysowanie tekstu z podziałem na linie
y_text = height - 120
left_margin = 56.69
right_margin = width - 56.69
max_width = right_margin - left_margin

for line in info_text.strip().split('\n'):
    line = line.strip()
    if line:
        # Podział długich linii
        words = line.split()
        current_line = ""
        for word in words:
            test_line = current_line + " " + word if current_line else word
            if pdfmetrics.stringWidth(test_line, "DejaVu", 11) < max_width:
                current_line = test_line
            else:
                c2.drawString(left_margin, y_text, current_line)
                y_text -= 16
                current_line = word
        if current_line:
            c2.drawString(left_margin, y_text, current_line)
            y_text -= 16
    else:
        y_text -= 8  # Pusta linia = mniejszy odstęp

c2.save()

# ========== MERGOWANIE ==========
merger = PdfWriter()
merger.append(PdfReader(output_title, "rb"))  # Strona tytułowa
merger.append(PdfReader(output_info, "rb"))   # Strona informacyjna
merger.append(PdfReader(pdf1, "rb"))
merger.append(PdfReader(pdf2, "rb"))
merger.add_blank_page(width=595.28, height=841.89)
with open(output, "wb") as f:
    merger.write(f)

# Czyszczenie
os.remove(output_title)
os.remove(output_info)
print(f"[INFO] Pliki połączono z tytułem i stroną informacyjną w {output}.")
