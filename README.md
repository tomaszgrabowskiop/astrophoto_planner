# Astrophotography Planner & Atlas Generator

Zestaw zaawansowanych skryptów w języku Python służący do generowania spersonalizowanego rocznego planera astronomicznego oraz atlasu obiektów głębokiego nieba (DSO). [file:1]

System pobiera dane z katalogów astronomicznych, filtruje je pod kątem lokalizacji obserwatora i posiadanego sprzętu (teleskop/kamera), oblicza precyzyjną widoczność na dany rok, a następnie generuje profesjonalny plik PDF zawierający: [file:1]

1. Przegląd roczny (kiedy obserwować dany obiekt). [file:1]
2. Szczegółowe strony dla każdego obiektu (wykresy wysokości, kadry FOV, mapy kontekstowe). [file:1]

---

## 🚀 Możliwości

- **Agregacja** danych: łączy katalogi NGC/IC, Sharpless (Sh2), RCW, Barnard, LBN, LDN, Cederblad i PGC. [file:1]
- **Inteligentne** filtrowanie: wybiera obiekty na podstawie szerokości geograficznej, minimalnej wysokości nad horyzontem, jasności (Mag), rozmiaru oraz skali Bortle. [file:1]
- Symulacja FOV: generuje symulacje kadru (Field of View) dla kamery i teleskopu przy użyciu biblioteki `starplot`. [file:1]
- Obliczenia astronomiczne: wylicza okna obserwacyjne (godziny bez Księżyca, wysokość górowania). [file:1]
- Format PDF: generuje gotowy do druku atlas w formacie A4. [file:1]

---

## 🛠️ Wymagania

Projekt wymaga Pythona 3.10+ oraz następujących bibliotek: [file:1]

```bash
pip install pandas numpy astropy astroplan matplotlib reportlab pypdf tqdm astroquery starplot networkx
```

**Uwaga:** Biblioteka `starplot` może wymagać dodatkowej konfiguracji (pobrania danych gwiazd). [file:1]

---

## 📂 Struktura plików i dane wejściowe

Aby rozpocząć, upewnij się, że posiadasz plik źródłowy dla katalogu NGC (używany w kroku 0): [file:1]

- `OpenNGC/NGC.csv` – plik CSV z danymi OpenNGC (wymagany przez skrypt `0_opracuj_katalog_ngc.py`). [file:1]

---

## ⚙️ Instrukcja użycia (krok po kroku)

Skrypty są ponumerowane, aby ułatwić zachowanie odpowiedniej kolejności wykonywania operacji. [file:1]

### Krok 0: Przygotowanie bazy NGC

Uruchom: [file:1]

```bash
python 0_opracuj_katalog_ngc.py
```

- Parsuje surowy plik CSV z OpenNGC. [file:1]
- Tworzy plik `updated_ngc.csv`. [file:1]

### Krok 1: Pobieranie i unifikacja katalogów

Uruchom: [file:1]

```bash
python 1_generuj_katalog_astro.py
```

- Pobiera dane z serwisu VizieR (Sharpless, Barnard, LDN, itp.). [file:1]
- Łączy je z bazą NGC. [file:1]
- Wykonuje „Smart Merge” (łączenie duplikatów i obiektów blisko siebie). [file:1]
- Tworzy plik `katalog_astro_full.csv`. [file:1]
- Opcjonalnie: uruchom `analiza_katalog.py`, aby sprawdzić statystyki bazy. [file:1]

### Krok 2: Konfiguracja i selekcja obiektów

Uruchom: [file:1]

```bash
python 2_ograniczenie_katalogu.py
```

- Interaktywny skrypt: pyta o lokalizację, parametry teleskopu/kamery, filtry (Ha/OIII) oraz minimalną wysokość obiektu. [file:1]
- Filtruje bazę pod kątem używanego sprzętu. [file:1]
- Tworzy plik konfiguracyjny `vis_data.json` z kandydatami do atlasu. [file:1]

### Krok 3: Silnik obliczeniowy (Engine)

Uruchom: [file:1]

```bash
python 3_wyliczenia.py
```

- Wykonuje ciężkie obliczenia astronomiczne (równolegle na wielu rdzeniach CPU). [file:1]
- Wylicza dokładną widoczność minuta po minucie dla całego roku. [file:1]
- Zapisuje wyniki do `observing_data.pkl`. [file:1]

### Krok 4: Plan roczny i wybór wariantów

Uruchom: [file:1]

```bash
python 4_plan_roczny.py
```

- Analizuje dane z kroku 3. [file:1]
- Przydziela obiekty do miesięcy (Warianty A, B, C), aby zbalansować sesje obserwacyjne. [file:1]
- Generuje **Część** 1 PDF: `Astrophotography_Planner_2026_1.pdf` (wykresy roczne). [file:1]
- Aktualizuje `vis_data.json` o flagę `selected`. [file:1]

### Krok 5: Generowanie map nieba

Uruchom: [file:1]

```bash
python 5_fov_and_maps.py
```

- Korzysta z biblioteki `starplot`. [file:1]
- Generuje pliki PNG w katalogu `starplots/`: [file:1]
  - Kadry optyczne (symulacja kamery). [file:1]
  - Mapy kontekstowe (szersze pole widzenia). [file:1]

### Krok 6: Generowanie stron obiektów

Uruchom: [file:1]

```bash
python 6_drukuj_strony_obiektów.py
```

- Składa szczegółowe strony dla każdego wybranego obiektu. [file:1]
- Zawiera wykresy wysokości w noc nowiu, wykres roczny, statystyki godzinowe oraz wygenerowane mapy. [file:1]
- Tworzy **Część** 2 PDF: `Astrophotography_Planner_2026_2.pdf`. [file:1]

### Krok 7: Finalizacja

Uruchom: [file:1]

```bash
python 7_połącz_pliki_pdf.py
```

- Generuje stronę tytułową. [file:1]
- Łączy część 1 i część 2 w jeden plik. [file:1]
- Wynik końcowy: `Astrophotography_Planner_2026.pdf`. [file:1]

---

## 📝 Uwagi dodatkowe

- Czcionki: skrypt `7_połącz_pliki_pdf.py` jest skonfigurowany pod system macOS (`/System/Library/Fonts/Helvetica.ttc`); na Windows lub Linux należy edytować ścieżkę do czcionek. [file:1]
- Wydajność: krok 3 i 5 wykorzystują wielowątkowość (`multiprocessing`); generowanie map może zająć kilka minut w zależności od liczby obiektów. [file:1]
- Lokalizacja: domyślnie ustawiony jest rok 2026 i lokalizacja w Polsce; można to zmienić w trakcie działania skryptu nr 2 lub edytując stałe w plikach. [file:1]

---

## 📄 Licencja

Projekt do użytku własnego. Korzysta z danych OpenNGC oraz serwisów VizieR. [file:1]
