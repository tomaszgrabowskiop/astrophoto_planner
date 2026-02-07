import json

from pypdf import PdfWriter, PdfReader
from reportlab.pdfgen import canvas
from reportlab.lib.pagesizes import A4
from reportlab.pdfbase import pdfmetrics
from reportlab.pdfbase.ttfonts import TTFont
#from reportlab.lib.colors import blue  # Dla niebieskiego
from reportlab.pdfbase.pdfmetrics import stringWidth

import os
from reportlab.pdfbase.pdfmetrics import registerFontFamily
from datetime import date

MONTH_NAMES = {
    1: "stycznia",
    2: "lutego",
    3: "marca",
    4: "kwietnia",
    5: "maja",
    6: "czerwca",
    7: "lipca",
    8: "sierpnia",
    9: "września",
    10: "października",
    11: "listopada",
    12: "grudnia",
}

with open('vis_data.json', 'r') as f:  
    dane = json.load(f)

miasto = dane['location']['name'] 
lat_lon = f"Lat.: {dane['location']['lat']:.2f}°N, Lon.: {dane['location']['lon']:.2f}°E"
dzisiaj = date.today()
dzien = dzisiaj.day
miesiac_pl = MONTH_NAMES[dzisiaj.month]
rok = dzisiaj.year
miejsce_wydania= f"{miasto}\n{lat_lon}\n{dzien} {miesiac_pl} {rok} roku"

# Ustawienia tytułu
tytul = "Astrophotography"
tytul1 = "Planner"
tytul2  = f"{dane['year']}"
output_title = "tytulowa.pdf"
pdf1 = "Astrophotography_Planner_2026_1.pdf"
pdf2 = "Astrophotography_Planner_2026_2.pdf"
output = "Astrophotography_Planner_2026.pdf"

# Rejestracja polskiej czcionki TTF (ścieżka macOS; dostosuj jeśli brak)
font_path = "/System/Library/Fonts/Helvetica.ttc"
pdfmetrics.registerFont(TTFont("HelveticaPL", font_path))
pdfmetrics.registerFont(TTFont("HelveticaPL-Bold", font_path))
pdfmetrics.registerFont(TTFont("HelveticaPL-Italic", font_path))

# Tworzenie strony tytułowej
c = canvas.Canvas(output_title, pagesize=A4)
width, height = A4
c.setTitle(tytul)

# Tytuł: niebieski, polska czcionka bold
c.setFont('HelveticaPL-Bold', 48)
c.setFillColorRGB(0.07, 0.07, 0.25)  # Niebieski kolor tekstu
c.drawCentredString(width / 2, height / 2 + 120, tytul)  # Pierwsza linia wyżej
c.drawCentredString(width / 2, height / 2 + 70, tytul1)  # Druga linia niżej
c.setFont('HelveticaPL-Bold', 34)
c.drawCentredString(width / 2, height / 2 + 20, tytul2)  # Trzecia linia niżej

c.setFillColor('black')  # Reset do czarnego
c.setFont('HelveticaPL-Italic', 12)

# Data: czarna, normalna czcionka
y_pos = 108

lines = miejsce_wydania.split('\n')
line_height = 14

for i, line in enumerate(lines):
    c.drawCentredString(width / 2, y_pos - i * line_height, line)

c.save()

# Mergowanie (jak wcześniej)...

# Mergowanie
merger = PdfWriter()
merger.append(PdfReader(output_title, "rb"))  # Tytuł na początku
merger.add_blank_page(width=595.28, height=841.89)  # A4 w punktach [web:72]
merger.append(PdfReader(pdf1, "rb"))
merger.append(PdfReader(pdf2, "rb"))
merger.add_blank_page(width=595.28, height=841.89)  # A4 w punktach [web:72]
with open(output, "wb") as f:
    merger.write(f)

# Czyszczenie
os.remove(output_title)
print(f"Pliki połączono z tytułem w {output}")
