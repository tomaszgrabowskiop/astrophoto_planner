# Astrophotography Planner & Atlas Generator

Zestaw skryptów w języku Python służący do generowania spersonalizowanego rocznego planera astronomicznego oraz atlasu obiektów głębokiego nieba (DSO). Dla niezaawansowanych amatorów nieba. Atomatycznie wybiera z bazy 90k+ obiektów te, którymi warto się zająć w bieżącym roku. Obiekty są "rozkładane" na przestrzeni roku tak, żeby stworzyć najbardziej zachłanny plan astrofotograficzny.

Możesz wywołać `run.sh`, który założy środowisko i będzie uruchamiał kolejne skrypty. Pozwoli Ci też pominąć wybrane kroki. Przydatne, jeśli dokonasz przeliczeń, które są najbardziej czasochłonne i nie będziesz chciał ich powtarzać.

Pełna dokumentacja znajduje się w pliku: ![doc/AstroPhotography Planner – Dokumentacja techniczna.md]
**Powodzenia!**

---

## Szczegóły działania skryptów

System pobiera dane z katalogów astronomicznych, **filtruje** je pod kątem lokalizacji obserwatora i posiadanego sprzętu (teleskop/kamera), oblicza precyzyjną widoczność na dany rok, a następnie generuje profesjonalny plik PDF zawierający: 

1. Przegląd roczny (kiedy obserwować dany obiekt).

![Widok propozycji miesięcznych](doc/Planner_roczny.png)

2. Szczegółowe strony dla każdego obiektu (wykresy wysokości, kadry FOV, mapy kontekstowe). 

![Strona obiektu](doc/Strona_obiektu.png)

---

## 🚀 Możliwości

- **Agregacja** danych: łączy katalogi NGC/IC, Sharpless (Sh2), RCW, Barnard, LBN, LDN, Cederblad i PGC. 
- **Inteligentne** filtrowanie: wybiera obiekty na podstawie szerokości geograficznej, minimalnej wysokości nad horyzontem, jasności (Mag), rozmiaru oraz skali Bortle. 
- Symulacja FOV: generuje symulacje kadru (Field of View) dla kamery i teleskopu przy użyciu biblioteki `starplot`. 
- Obliczenia astronomiczne: wylicza widoczność w ciągu roku w zależności od podanego progu nocy, widoczność bez Księżyca, wysokość górowania w przykładowej nocy na tle zmierzchu i schodu (cywilnych, żeglarskich i astronomicznych). 
- Format PDF: generuje gotowy do druku atlas w formacie A4. 

---

## 🛠️ Wymagania

Projekt wymaga Pythona 3.10+ oraz następujących bibliotek: 

```bash
pip install pandas numpy astropy astroplan matplotlib reportlab pypdf tqdm astroquery starplot networkx
```

**Uwaga:** Biblioteka `starplot` może wymagać dodatkowej konfiguracji (pobrania danych gwiazd). 

---

## 📂 Struktura plików i dane wejściowe

Aby rozpocząć, upewnij się, że posiadasz plik źródłowy dla katalogu NGC (używany w kroku 0): 

- `OpenNGC/NGC.csv` – plik CSV z danymi OpenNGC (wymagany przez skrypt `0_opracuj_katalog_ngc.py`). 

---

## ⚙️ Instrukcja użycia (krok po kroku)

Skrypty są ponumerowane, aby ułatwić zachowanie odpowiedniej kolejności wykonywania operacji. 

### Krok 0: Przygotowanie bazy NGC

Uruchom: 

```bash
python 0_opracuj_katalog_ngc.py
```

- Parsuje surowy plik CSV z OpenNGC. 
- Tworzy plik `updated_ngc.csv`. 

### Krok 1: Pobieranie i unifikacja katalogów

Uruchom: 

```bash
python 1_generuj_katalog_astro.py
```

- Pobiera dane z serwisu VizieR (Sharpless, Barnard, LDN, itp.). 
- Łączy je z bazą NGC. 
- Wykonuje „Smart Merge” (łączenie duplikatów i obiektów blisko siebie). 
- Tworzy plik `katalog_astro_full.csv`. 
- Opcjonalnie: uruchom `analiza_katalog.py`, aby sprawdzić statystyki bazy. 

### Krok 2: Konfiguracja i selekcja obiektów

Uruchom: 

```bash
python 2_ograniczenie_katalogu.py
```

- Interaktywny skrypt: pyta o lokalizację, parametry teleskopu/kamery, filtry (Ha/OIII) oraz minimalną wysokość obiektu, próg wysokości słońca, czas trwania okna obserwacyjnego. 
- Filtruje bazę pod kątem parametrów użytkownika i określonego minimalnego rozmiaru i jasnosci dla obiektów. 
- Tworzy plik konfiguracyjny `vis_data.json` z kandydatami do atlasu. 

### Krok 3: Silnik obliczeniowy (Engine)

Uruchom: 

```bash
python 3_wyliczenia.py
```

- Wykonuje ciężkie obliczenia astronomiczne dzięki AstroPy (równolegle na wielu rdzeniach CPU). 
- Wylicza dokładną widoczność minuta po minucie dla całego roku. 
- Zapisuje wyniki do `observing_data.pkl`.
- Możliwe wywołanie bez powtarzania obliczeń (sam sprawdza, czy zmieniły się parametry).
- Podjęćie obliczeń tylko w zakresie maski (parametry: wysokość nad horyzontem, długość okna obserwacyhnego, określenie zmierzchu/świtu).
- Pełne obliczenia dla obiektów, których nie było wcześniej i przy zmianie lokalizacji.  

### Krok 4: Plan roczny i wybór wariantów

Uruchom: 

```bash
python 4_plan_roczny.py
```

- Analizuje dane z kroku 3. 
- Przydziela obiekty do miesięcy (Warianty A, B, C), aby zbalansować sesje obserwacyjne. 
- Generuje **Część** 1 PDF: `Astrophotography_Planner_2026_1.pdf` (wykresy roczne). 
- Aktualizuje `vis_data.json` o flagę `selected`. 

### Krok 5: Generowanie map nieba

Uruchom: 

```bash
python 5_fov_and_maps.py
```

- Korzysta z biblioteki `starplot`. 
- Generuje pliki PNG w katalogu `starplots/`: 
  - Kadry optyczne (symulacja kamery). 
  - Mapy kontekstowe (szersze pole widzenia). 

### Krok 6: Generowanie stron obiektów

Uruchom: 

```bash
python 6_drukuj_strony_obiektow.py
```

- Składa szczegółowe strony dla każdego wybranego obiektu. 
- Zawiera wykresy wysokości w nocy, wykres rocznej widoczności, wykres liczby godzin z/bez księżyca oraz wygenerowane mapy. 
- Tworzy **Część** 2 PDF: `Astrophotography_Planner_2026_2.pdf`. 

### Krok 7: Finalizacja

Uruchom: 

```bash
python 7_polacz_pliki_pdf.py
```

- Generuje stronę tytułową. 
- Łączy część 1 i część 2 w jeden plik. 
- Wynik końcowy: `Astrophotography_Planner_2026.pdf`. 

---

## 📝 Uwagi dodatkowe

- Czcionki: skrypt `7_polacz_pliki_pdf.py` jest skonfigurowany pod system macOS (`/System/Library/Fonts/Helvetica.ttc`); na Windows lub Linux należy edytować ścieżkę do czcionek. 
- Wydajność: krok 3 i 5 wykorzystują wielowątkowość (`multiprocessing`). Mimo to obliczenia AstroPy i generowanie map może zająć sporo czasu w zależności od liczby obiektów i wydajności komputera.
- TimeZone: użytkownik może wybrać TimeZone. Wpływa na obliczenia, wykresy miesięczne i przykładowej nocy również na wykresie widoczności w skali roku jest zaznaczone przesunięcie godzinowe.  
- Lokalizacja: domyślnie ustawiony jest rok 2026 i lokalizacja w Polsce; można to zmienić w trakcie działania skryptu nr 2 lub edytując stałe w plikach. 

---

## 📄 Licencja

Projekt do użytku własnego. Korzysta z danych OpenNGC oraz serwisów VizieR. 
