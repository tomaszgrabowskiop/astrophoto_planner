# 🔭 Astrophotography Planner & Atlas Generator

Kompleksowy zestaw narzędzi w Pythonie służący do generowania spersonalizowanego, rocznego planera i atlasu astrofotograficznego. System automatycznie pobiera dane o obiektach, filtruje je pod kątem Twojego sprzętu i lokalizacji, oblicza widoczność na przestrzeni roku, a następnie generuje profesjonalny dokument PDF gotowy do druku.

## 🚀 Możliwości

*   **Agregacja danych:** Pobieranie i łączenie katalogów (NGC/IC, Sharpless, Barnard, RCW, PGC, LBN, LDN, Cederblad) z bazy VizieR i OpenNGC.
*   **Smart Merge:** Inteligentne łączenie dublujących się obiektów (np. Mgławica Kalifornia jako NGC 1499 i Sh2-220) oraz klastrowanie obiektów mieszczących się w jednym kadrze (FOV).
*   **Punktacja (Scoring):** Ocenianie obiektów na podstawie ich rodzaju, jasności, rozmiaru oraz "sławy" (Messier, Caldwell, Herschel 400).
*   **Obliczenia astronomiczne:** Precyzyjne wyliczanie wysokości nad horyzontem i godzin obserwacyjnych dla każdej nocy w roku z uwzględnieniem faz księżyca i zanieczyszczenia światłem (Bortle).
*   **Planowanie roczne:** Wykorzystanie algorytmu optymalizacji (Hungarian Algorithm) do przydzielenia najlepszych obiektów do optymalnych miesięcy obserwacyjnych.
*   **Generowanie map:** Tworzenie symulacji pola widzenia (FOV) oraz map kontekstowych (star hopping) przy użyciu biblioteki `starplot`.
*   **Output:** Finalny plik PDF zawierający harmonogram roczny, szczegółowe karty obiektów oraz mapy.

## 🛠️ Wymagania i Instalacja

Projekt wymaga Pythona 3.10+ oraz szeregu bibliotek astronomicznych i graficznych.

1.  **Sklonuj repozytorium:**
    ```bash
    git clone https://github.com/twoje-konto/astro-planner.git
    cd astro-planner
    ```

2.  **Zainstaluj zależności:**
    Zaleca się użycie wirtualnego środowiska (venv).
    ```bash
    pip install -r requirements.txt
    ```
    *Główne biblioteki to: `astropy`, `astroplan`, `starplot`, `pandas`, `matplotlib`, `reportlab`, `networkx`.*
    
    **UWAGA: pandas musi być poniżej wersji 3.0!**

3.  **Przygotowanie danych:**
    Upewnij się, że w folderze `data/` znajduje się plik `NGC.csv` (baza OpenNGC), jeśli skrypt go nie pobierze automatycznie.

##Jak używać?

Proces składa się z 7 kroków, które należy uruchamiać sekwencyjnie. Każdy skrypt korzysta z danych wygenerowanych przez poprzedni.

### Krok 1: Budowa katalogu
```bash
python 1_build_catalog.py
```
*   Pobiera dane z VizieR i łączy je w jeden spójny plik CSV (`katalog_astro_full.csv`).
*   Dokonuje normalizacji nazw, współrzędnych i "Entity Resolution" (łączenie tych samych obiektów z różnych katalogów).
*   *Interakcja:* Pyta o progi filtrowania (jasność mag, rozmiar kątowy).

### Krok 2: Konfiguracja i Scoring
```bash
python 2_plan_and_score.py
```
*   **Kluczowy etap konfiguracji użytkownika.**
*   Pyta o: Rok, Lokalizację (współrzędne), Parametry kamery/teleskopu (do obliczenia FOV), Limity horyzontu i Zanieczyszczenie nieba (Bortle).
*   Grupuje obiekty w kadry pasujące do Twojego sensora.
*   Zapisuje wynik w `vis_data.json`.

### Krok 3: Obliczenia astronomiczne (Compute Engine)
```bash
python 3_compute.py
```
*   Najbardziej czasochłonny etap. Wykorzystuje wielordzeniowość (multiprocessing).
*   Liczy wysokość obiektów dla każdej nocy w roku.
*   Cache'uje wyniki w `observing_data.pkl`, aby przy kolejnych uruchomieniach liczyć tylko zmiany.

### Krok 4: Selekcja i Plan roczny
```bash
python 4_select_objects.py
```
*   Dzieli obiekty na warianty (A, B, C) i przypisuje je do najlepszych miesięcy.
*   Generuje pierwszą część PDF: wykresy zbiorcze na każdy miesiąc.
*   Modyfikuje `vis_data.json` dodając flagę `selected`.

### Krok 5: Generowanie map nieba
```bash
python 5_generate_fov_and_ctx.py
```
*   Generuje pliki PNG dla każdego wybranego obiektu.
*   **FOV:** Symulacja kadru Twojej kamery.
*   **Context:** Mapa szerszego pola (szukacz/star hopping).
*   Pliki lądują w folderze `data/starplots/`.

### Krok 6: Tworzenie stron obiektów
```bash
python 6_generate_objects_pages.py
```
*   Generuje szczegółowe strony PDF dla każdego obiektu (wykres wysokości w noc nowiu, wykres roczny, mapy, dane techniczne).

### Krok 7: Finalizacja (Scalanie PDF)
```bash
python 7_generate_result_pdf.py
```
*   Dodaje stronę tytułową i informacyjną.
*   Łączy wszystkie wygenerowane wcześniej PDF-y w jeden kompletny plik: `Astrophotography_Planner_ROK_MIASTO.pdf`.

---

## 📂 Struktura plików

*   `shared.py` – Plik konfiguracyjny współdzielony przez wszystkie skrypty. Zawiera ścieżki (`PATHS`), klasy konfiguracyjne (`UserConfig`, `CameraConfig`) oraz stałe astronomiczne.
*   `analyse_catalog.py` – Narzędzie pomocnicze do statystycznej analizy pliku `katalog_astro_full.csv` (sprawdza kompletność danych, rozkład jasności itp.).
*   `data/` – Katalog roboczy (tworzony automatycznie), w którym przechowywane są pliki tymczasowe, cache obliczeń, grafiki i wyniki.

## ⚙️ Konfiguracja zaawansowana

Większość parametrów podaje się interaktywnie w kroku 2. Jednak stałe systemowe można edytować w pliku `shared.py` oraz na początku poszczególnych skryptów (np. wagi punktacji w `2_plan_and_score.py`).

*   **CATALOG_PRIORITY:** W `shared.py` określa, który katalog jest "ważniejszy" przy łączeniu nazw (np. NGC > IC > Sh2).
*   **SCORING:** W `2_plan_and_score.py` znajdują się tabele punktacji (np. bonusy za obiekty Messiera, punkty za typ obiektu vs filtry wąskopasmowe).

## 📄 Licencja

Projekt przeznaczony do użytku prywatnego i edukacyjnego. Nie jestem programistą ani astronomem. Zgłaszaj poprawki lub błędy, jeśli je zauważysz. 

Generowane mapy korzystają z danych z VizieR. [OpenNGC](https://github.com/mattiaverga/OpenNGC/tree/master) to osobny projekt częściowo wykorzystywany w tym repozytorium.