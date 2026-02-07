from pypdf import PdfWriter, PdfReader
from reportlab.pdfgen import canvas
from reportlab.lib.pagesizes import A4
from reportlab.pdfbase import pdfmetrics
from reportlab.pdfbase.ttfonts import TTFont
from reportlab.lib.colors import blue  # Dla niebieskiego
import os
from reportlab.pdfbase.pdfmetrics import registerFontFamily

# Ustawienia tytułu
tytul = "Astrophotography"
tytul1 = "Planner"
tytul2  = "2026"
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
c.drawCentredString(width / 2, height / 2 + 20, tytul2)  # Druga linia niżej

c.setFillColor('black')  # Reset do czarnego

# Data: czarna, normalna czcionka
c.setFont('HelveticaPL-Italic', 12)
c.drawCentredString(width / 2, 50, "Kraków, 6 lutego 2026 roku")

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
