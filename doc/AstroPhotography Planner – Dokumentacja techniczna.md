# AstroPhotography Planner – Dokumentacja techniczna

## 1. Przegląd systemu

AstroPhotography Planner to wieloetapowy pipeline w Pythonie służący do generowania rocznego planu astrofotografii oraz drukowalnego atlasu PDF na podstawie wielu katalogów DSO, parametrów lokalizacji i konfiguracji sprzętu.
System łączy dane z katalogów NGC/IC (OpenNGC) oraz VizieR (Sharpless, Barnard, RCW, LDN, LBN, Cederblad, PGC), przetwarza je do wspólnego formatu, wylicza widoczność roczną obiektów, optymalizuje ich przydział do miesięcy i generuje mapy nieba oraz strony atlasu.

### 1.1. Etapy pipeline’u

Pipeline składa się z następujących kroków, realizowanych kolejno przez osobne skrypty:

    0_opracuj_katalog_ngc.py       → preprocessing NGC/IC (OpenNGC → updated_ngc.csv)
    1_generuj_katalog_astro.py     → pobranie katalogów, normalizacja, Smart Merge → katalog_astro_full.csv
    2_ograniczenie_katalogu.py     → konfiguracja, scoring, selekcja → vis_data.json
    3_wyliczenia.py                → silnik widoczności (RAW/FINAL cache) → observing_data.pkl
    4_plan_roczny.py               → optymalny plan roczny + monthly_overview.pdf + selected w vis_data.json
    5_fov_and_maps.py              → generacja kadrów FOV i map kontekstowych → starplots/*.png
    6_drukuj_strony_obiektow.py    → strony obiektów do PDF (warianty) → Astrophotography_Planner_2026_*.pdf
    7_polacz_pliki_pdf.py          → strona tytułowa + merge → Astrophotography_Planner_2026.pdf

---

## 2. Przepływ danych i formaty plików

### 2.1. Pliki wejściowe

- `OpenNGC/NGC.csv` – oryginalny katalog OpenNGC wykorzystywany jako źródło NGC/IC.
- Katalogi VizieR (pobierane dynamicznie): Sharpless, Barnard, RCW, LDN, LBN, Cederblad, PGC.

### 2.2. Pliki pośrednie

- `updated_ngc.csv` – katalog NGC/IC po preprocessingu, z polami `id`, `extra_info`, `type`, `ra`, `dec`, `mag`, `size`, `base_score`, `common_names`.
- `katalog_astro_full.csv` – zunifikowany katalog wszystkich obiektów po normalizacji i Smart Merge, z kolumnami m.in. `id`, `ra`, `dec`, `size`, `mag`, `type`, `extra_info`, `common_names`, `priority`, `catalog`.
- `vis_data.json` – główny plik konfiguracyjno‑katalogowy zawierający sekcje:
  - `location` – nazwa miejsca, szerokość/długość geograficzna, wysokość, strefa czasowa,
  - `year` – rok planowania,
  - `parameters` – m.in. `minalt`, `sunlimit`, informacje o filtrach, klasie Bortle i kamerze,
  - `objects` – lista obiektów z polami katalogowymi oraz polem `score` i opcjonalnym `selected`,
  - `parameters_hash` – hash parametrów użytych w silniku widoczności.
- `observing_data_raw.pkl` – cache surowych danych wysokości obiektu i czasów crossingów Słońca/obiektu, niezależny od progów jakości.
- `observing_data.pkl` – cache FINAL z wyliczonymi polami `q_hours`, `m_hours`, `qual_segments`, `sun_pts` i `tz_offset` dla każdego dnia roku.
- `observing_data_final.hash` – hash parametrów (`minalt`, `sunlimit`, lokalizacja, rok) użytych przy generowaniu FINAL, służy do wykrywania potrzeby ponownego przeliczenia.
- `observing_data_raw.hash` – hash parametrów RAW (lokalizacja + rok), służy do walidacji cache RAW
- `ENGINE_STATE_PKL` – zapis stanu silnika (poprzednia lokalizacja, rok i parametry), wykorzystywany do logów diagnostycznych.

### 2.3. Pliki wyjściowe

- `monthly_overview.pdf` – PDF z miesięcznym przeglądem widoczności podczas nocy nowiu, z wykresami dla wariantów A/B/C.
- `starplots/{id}.png` – kadry FOV dla wybranych obiektów na podstawie konfiguracji kamery.
- `starplots/context_{id}.png` – szerokie mapy kontekstowe dla tych samych obiektów.
- `Astrophotography_Planner_2026_1.pdf`, `Astrophotography_Planner_2026_2.pdf` – części atlasu z indywidualnymi stronami obiektów.
- `Astrophotography_Planner_2026.pdf` – finalny atlas po dodaniu strony tytułowej i scaleniu wszystkich części.

---

## 3. Opis modułów

## 3. Opis modułów

### 3.1. `0_opracuj_katalog_ngc.py`

**Cel:** preprocessing katalogu OpenNGC do formatu wewnętrznego `updated_ngc.csv`.

Moduł wczytuje plik `OpenNGC/NGC.csv` i iteruje po wierszach, budując listę rekordów z wybranymi polami.  
Dla każdego obiektu wyprowadza identyfikator `id` (NGC/IC), typ obiektu (`type`), współrzędne `ra`/`dec` w stopniach, jasność `mag`, rozmiar kątowy `size` w arcmin, dodatkowe identyfikatory w polu `extra_info` oraz nazwy zwyczajowe w polu `common_names`.  
Pole `extra_info` uzupełniane jest na podstawie predefiniowanych zbiorów Messier, Caldwell i Herschel, co umożliwia późniejsze bonusy w systemie punktacji.  
Na końcu moduł zapisuje dane do `updated_ngc.csv` z nagłówkiem i raportuje liczbę przetworzonych oraz pominiętych wierszy.

---

### 3.2. `1_generuj_katalog_astro.py`

**Cel:** pobranie katalogów z VizieR, ich normalizacja oraz konsolidacja w jeden katalog przy użyciu algorytmu Smart Merge.

#### 3.2.1. Pobieranie danych

Funkcja `fetch_data()` wczytuje lokalny `updated_ngc.csv` oraz pobiera z VizieR katalogi Sharpless, Barnard, RCW, PGC, LDN, LBN i Cederblad przy użyciu `Vizier`.  
Dla każdego katalogu wyliczana i wypisywana jest liczba wierszy oraz suma wszystkich rekordów.  
Rezultat zapisywany jest do `katalog_astro_full.csv` po zakończeniu normalizacji i Smart Merge.

#### 3.2.2. Normalizacja kolumn

Funkcja `normalize_all(raw)` unifikuje wszystkie katalogi do wspólnego zestawu kolumn (`id`, `ra`, `dec`, `size`, `mag`, `type`, `extra_info`, `common_names`, `priority`, `catalog`).  
Do konwersji RA/DEC z różnych formatów (sekstagesymalny/stopnie) używany jest `SkyCoord` z odpowiednimi jednostkami, a wyniki są zaokrąglane do pięciu miejsc po przecinku.  
Dla każdego źródła zastosowane są specyficzne reguły:

- Sharpless: generowane identyfikatory `Sh2-XX`, rozmiar z `Diam`, typ `HII`, nazwy zwyczajowe z mapy `SH2_COMMON_NAMES`.
- Barnard: `id=BXX`, rozmiar z `Diam`, typ `DN`.
- RCW: `id=RCWXX`, rozmiar z `MajAxis`, typ `HII`.
- Cederblad: złożone `id` z `Ced` i `m_Ced`, rozmiar z średniej `Dim1/Dim2`, `mag` z `vmag`, typ `NB`.
- LBN: `id=LBNSeq`, rozmiar z `Diam1`, typ `NB`.
- LDN: `id=LDNXX`, rozmiar z `sqrt(Area)*60`, typ `DN`.
- PGC: `id=PGCXXXX`, rozmiar z `MajAxis/60`, `mag` z `Btot`, typ `Gx`.

#### 3.2.3. Smart Merge

Funkcja `smart_merge(df, tolerance_deg)` konsoliduje obiekty w odległości mniejszej niż zadana tolerancja kątowa, traktując je jako potencjalne duplikaty.  
Na bazie `search_around_sky` budowany jest graf sąsiedztwa, w którym wierzchołki to obiekty, a krawędzie łączą obiekty występujące w odległości mniejszej niż `tolerance_deg`.  
Dla każdej składowej spójnej wybierany jest lider na podstawie priorytetu katalogu (NGC/IC > Sh2 > RCW > Barnard > LDN/LBN > Cederblad > PGC) oraz rozmiaru kątowego, a jego współrzędne stają się reprezentatywne dla całej grupy.  
Rozmiar wynikowy przyjmowany jest jako maksimum spośród sensownych rozmiarów w klastrze, ograniczonych limitem 300 arcmin, co zapobiega powstawaniu obiektów o absurdalnie dużych wymiarach.  
Identyfikatory i pola `extra_info` są łączone w unikalne listy, a `common_names` są czyszczone i deduplikowane przez funkcję `process_common_names`, która usuwa krótsze nazwy zawarte w dłuższych.  
Wynikowe rekordy zapisywane są do `katalog_astro_full.csv`.


---

### 3.3. `2_ograniczenie_katalogu.py`

**Cel:** interaktywna konfiguracja lokalizacji, parametrów obserwacyjnych i sprzętu oraz selekcja i punktacja obiektów do pliku `vis_data.json`.

Moduł definiuje klasę `Config` z domyślnymi ścieżkami plików oraz parametrami kamery (ogniskowa, rozmiar matrycy, pitch, rozdzielczość), które mogą być nadpisane danymi od użytkownika.  
Parametry domyślne obejmują m.in. nazwę katalogu roboczego, ścieżkę do `katalog_astro_full.csv` oraz wartości startowe dla kamery (np. ogniskowa 300 mm, typowa matryca z danymi pitch i rozdzielczości).  

Użytkownik w trybie konsolowym podaje:  
- nazwę lokalizacji,  
- szerokość i długość geograficzną,  
- wysokość n.p.m.,  
- rok planowania (domyslnie 2026),  
- minimalną wysokość obiektu (`minalt` w stopniach),  
- limit wysokości Słońca (`sunlimit` w stopniach, odpowiadający typowi zmierzchu: civil, nautical, astronomical),  
- informację o filtrach narrowband (np. obecność filtrów Ha, OIII, SII, Hb),  
- klasę Bortle miejsca obserwacji.  

Na podstawie współrzędnych moduł automatycznie wylicza i ustala nazwę strefy czasowej za pomocą `TimezoneFinder`, a następnie wypełnia w `vis_data.json` strukturę `location` z polami `name`, `lat`, `lon`, `height` oraz `time_zone`.  

#### 3.3.1. System punktacji

System punktacji bazuje na kilku tabelach: `SCORING_FILTER`, `SCORING_BORTLE`, `SCORE_FAMOUS` i `SCORE_DATA_QUALITY`, które są zdefiniowane w postaci stałe słowników wewnątrz modułu.  

- `SCORING_FILTER` przydziela różne punkty dla obiektów w zależności od typu i używanych filtrów. Na przykład mgławice NB/HII uzyskują wyższe punkty w przypadku filtrów narrowband, a galaktyki w przypadku filtrów typu triband (np. L, R, G, B lub „L‑RGB” z moderowanym wpływem L).  
- `SCORING_BORTLE` uwzględnia ciemność nieba: ciemne mgławice otrzymują wysoki bonus dla Bortle ≤ 5, a punkty są stopniowo zmniejszane do zera dla Bortle > 6–7, co odzwierciedla utratę kontrastu na jasnym niebie.  
- `SCORE_FAMOUS` dodaje punkty za przynależność do katalogów Messier, Caldwell i Herschel 400, rozpoznawaną na podstawie tokenów w polu `extra_info` (np. `"messier"`, `"caldwell"`, `"herschel"`).  
- `SCORE_DATA_QUALITY` rozróżnia pomiary bezpośrednie od szacunków: obiekty z danymi „mierzone” (`measured`) otrzymują wyższy bonus punktowy, natomiast te o danych „estymowane” (`estimated`) dostają niższe punkty lub nawet kary.  

Dodatkowo wyliczana jest jasność powierzchniowa obiektu na podstawie przybliżonej powierzchni koła o średnicy `size_arcmin` i jasności `mag`.  
Wzór ma postać typu `mu ≈ mag + 2.5 * log10(π * (size_arcmin/2)**2)`, a wynik jest wykorzystywany do filtrowania „zbyt cichych” powierzchniowo mgławic, które praktycznie nie są widoczne w danym warunkach obserwacyjnych.  

#### 3.3.2. Wyjście `vis_data.json`

Na wyjściu powstaje słownik `vis_data` z sekcjami:  

- `location` – opis miejsca obserwacji (nazwa, współrzędne, wysokość, strefa czasowa),  
- `year` – rok planowania (np. 2026),  
- `parameters` – lista parametrów obserwacyjnych (`minalt`, `sunlimit`, `bortle`, informacje o filtrach, parametry kamery, ewentualnie flagi konfiguracyjne),  
- `objects` – lista wszystkich obiektów z `katalog_astro_full.csv` po przefiltrowaniu i punktacji.  

Struktura `objects` w `vis_data.json` zawiera dla każdego obiektu:  

- pola katalogowe (np. `id`, `ra`, `dec`, `size`, `mag`, `type`, `extra_info`, `common_names`, `catalog`),  
- `score` – suma punktów z `SCORING_FILTER`, `SCORING_BORTLE`, `SCORE_FAMOUS`, `SCORE_DATA_QUALITY` oraz ewentualnych modyfikacji z jasności powierzchniowej,  
- `surface_brightness` – obliczona jasność powierzchniowa,  
- `selected` – opcjonalnie pole `null` albo `undefined` na tym etapie, a następnie wypełniane w `4_plan_roczny.py`.  

Funkcja `add_parameters_hash_to_output` dokleja do słownika `vis_data` pole `parameters_hash`, które jest wartością MD5 z konkatenacji `minalt`, `sunlimit`, rok, `lat`, `lon` oraz `bortle`.

Ten hash jest później wykorzystywany w `3_wyliczenia.py` do szybkiego wykrywania, czy zmiany parametrów jakości obserwacji wymagają ponownego wyliczenia widoczności (RAW/FINAL).  
Cała struktura jest zapisywana do `vis_data.json` z nagłówkiem `{"generated_at": ...}` oraz komentarzem w komentarzach JSON wskazującym na etap kroku `2_ograniczenie_katalogu.py`.

---

### 3.4. `3_wyliczenia.py`

**Cel:** obliczanie rocznej widoczności obiektów z wykorzystaniem cache RAW/FINAL oraz logiki incrementalnej.

Moduł definiuje parametry globalne czasu (`H_START`, `H_END`, `N_SAMPLES`) oraz liczby próbek do detekcji crossingów (`CROSSING_SAMPLES`), a także ścieżki do plików cache.
Dataclass `RawObjectData` przechowuje surowe wysokości obiektu (`o_alt_all`) oraz listy czasów crossingów Słońca i obiektu względem horyzontu.

#### 3.4.1. Etap RAW

Moduł został zoptymalizowany poprzez wydzielenie obliczeń crossingów Słońca do etapu globalnego, wykonywanego raz przed przetwarzaniem obiektów. Funkcja compute_raw_data teraz najpierw oblicza precomputed_sun_pts – listę czasów wschodu/zachodu Słońca dla wszystkich 365 dni – a następnie przekazuje je do workerów poprzez partial. Dzięki temu crossingi Słońca nie są liczone redundantnie dla każdego obiektu osobno, co znacząco przyspiesza etap RAW.

* Funkcja `precomputed_sun_pts: List[List[datetime]]` – pre-obliczone crossingi Słońca (0°) dla każdego dnia, przekazywane z głównego wątku zamiast obliczania w pętli.
* Funkcja `compute_raw_data(json_path, object_limit, max_workers)` wczytuje `vis_data.json`, określa lokalizację, rok i listę obiektów ograniczoną do `object_limit`.
* Dla każdego dnia roku tworzony jest czas południa (`days_noon`), a następnie siatka czasów nocnych przez dodanie wektora przesunięć godzinowych `t_night_offsets_hours`.
* Dla każdego obiektu obliczana jest macierz wysokości `o_alt_all` w układzie `AltAz`, a także czasy crossingów 0° dla Słońca i obiektu przy użyciu funkcji `get_crossings`.
* Jeśli `RAW_DATA_PKL` istnieje, moduł wczytuje istniejące dane i liczy tylko brakujące obiekty, po czym nadpisuje plik z pełnym zbiorem RAW.

#### 3.4.1.1
Wprowadzono osobny mechanizm hashowania dla danych RAW niezależny od parametrów jakości (`minalt`, `sunlimit`). Funkcja `get_raw_params_hash` generuje hash MD5 z lokalizacji (lat, lon) i roku, zapisując go w pliku `observing_data_raw.hash`. Cache RAW jest invalidowany tylko gdy zmienią się parametry geometryczne (lokalizacja lub rok), a nie przy każdej zmianie progów obserwacyjnych. Dzięki temu ponowne przeliczenie FINAL z innym `minalt` lub `sunlimit` nie wymaga przetwarzania RAW od nowa.

#### 3.4.2. Etap FINAL

* Funkcja `process_raw_to_final` przekształca dane RAW jednego obiektu w listę rekordów dziennych zawierających m.in. `q_hours`, `m_hours`, `qual_segments` oraz zmodyfikowane `sun_pts` względem progu `sunlimit`.
Maska jakości tworzona jest jako warunek spełnienia minimalnej wysokości obiektu oraz odpowiedniego zanurzenia Słońca pod horyzontem.
* Funkcja `reprocess_to_final` wylicza finalne dane dla wszystkich obiektów, generując wcześniej globalne siatki wysokości Słońca i Księżyca dla całego roku, co minimalizuje liczbę transformacji układu odniesienia.
Po zakończeniu przetwarzania zapisuje `observing_data.pkl` oraz aktualny hash parametrów do `observing_data_final.hash`.

Obie funkcje zostały zrefaktoryzowane do pracy z `ProcessPoolExecutor`. Funkcja `process_raw_to_final` została wydzielona jako worker function i uruchamiana jest równolegle dla wszystkich obiektów. Wspólne obliczenia wysokości Słońca i Księżyca (`sun_alt_all`, `moon_alt_all`) są wykonywane raz w głównym procesie i przekazywane do workerów poprzez `partial`, co minimalizuje overhead transformacji układu odniesienia.

#### 3.4.3. Logika incrementalna

* Funkcja `should_reprocess` porównuje aktualny hash parametrów z zapisanym w `observing_data_final.hash`, decydując czy wymagane jest pełne przeliczenie FINAL.
* `reprocess_missing_final` potrafi dogenerować brakujące obiekty do istniejącego FINAL cache, scalając nowe wyniki z dotychczasowymi i zwracając tylko podzbiór dla wybranych `target_ids`.
* Funkcja `run_engine_from_vis_json` steruje całością procesu, obsługując inteligentne cache'owanie na dwóch poziomach: RAW (zależny od geometrii) i FINAL (zależny od progów jakości). Przed rozpoczęciem obliczeń raportuje poprzedni i bieżący stan silnika (lokalizacja, rok, parametry). Pyta użytkownika o limit obiektów do przeliczenia. Obsługuje flagę --force-all, która usuwa wszystkie cache'e i wymusza pełne przeliczenie.

---

### 3.5. `4_plan_roczny.py`

**Cel:** obliczenie miesięcznych parametrów widoczności, optymalny przydział obiektów do miesięcy i wariantów oraz zapis wyboru w `vis_data.json` i zgenerowanie PDF przeglądowego (`monthly_overview.pdf`).

Moduł `4_plan_roczny.py` jest odpowiedzialny za przejście od surowych danych widoczności (`observing_data.pkl`) do konkretnego rocznego planu, wyrażonego jako przydział obiektów do miesięcy i wariantów (A/B/C).  
Na wyjściu moduł aktualizuje plik `vis_data.json`, dodając w każdym obiekcie pole `selected` z informacją o wariant, miesiącu i dacie przypisania, oraz generuje wykresy widoczności w formacie `monthly_overview.pdf`.

Dataclass `MonthlyAssignment` opisuje pojedynczy wariant (np. wariant A, B, C) w postaci struktury składającej się z:
- `name` – nazwa wariantu (np. `"A"`, `"B"`, `"C"`),
- `month_to_ids` – słownik, w którym kluczem jest numer miesiąca (1–12), a wartością lista `id` obiektów przypisanych do tego wariantu i miesiąca.  

Struktura ta jest uogólnieniem zadanego problemu optymalizacji: system określa, które obiekty i w jakich miesiącach mają być przydzielone do którego wariantu, z uwzględnieniem ograniczeń pojemnościowych.

Funkcje pomocnicze `load_vis_data` i `load_observing_data` wczytują odpowiednio:
- `vis_data.json` – zawierający obiekty, ich `score`, parametry lokalizacji, rok planowania oraz bieżące obliczenia z `2_ograniczenie_katalogu.py`,
- `observing_data.pkl` – z pre‑obliczonymi danymi widoczności dla każdego obiektu (pola `q_hours`, `m_hours`, `qual_segments`, `sun_pts` itp. dla każdego dnia roku).  

Wszystkie kolejne obliczenia są wykonywane w przestrzeni obiektów o `score > 0`, filtrując te zbyt słabie punktowane lub z niską wizualną wartością w danym setupie obserwacyjnym.

#### 3.5.1. Metryki widoczności

* Funkcja `compute_monthly_best_q_hours(vis_data, observing_data, threshold=None)` tworzy tablicę / DataFrame z kolumnami `id`, `month`, `best_q_hours`.  
Dla każdego obiektu i miesiąca wyliczana jest najlepsza wartość `q_hours` spośród wszystkich nocy w danym miesiącu, spełniających warunki jakości obserwacji (minimalna wysokość obiektu, odpowiedni typ zmierzchu).  
Pole `best_q_hours` reprezentuje więc „najlepszą możliwą długość okna obserwacyjnego” w miesiącu dla danego obiektu.

* Dodatkowo funkcja `compute_yearly_annual_vis(vis_data, observing_data, q_hours_threshold=0.5)` zlicza liczbę nocy w całym roku, w których `q_hours` dla obiektu przekracza zadany próg `q_hours_threshold` (domyślnie 0.5 godziny).  
Wynik jest zapisywany w kolumnie `annual_vis` DataFrame, reprezentując „liczbę użytecznych nocy” w roku dla danego obiektu – im większa wartość, tym częściej obiekt jest w dobrych warunkach obserwacyjnych.

* Obie metryki (`best_q_hours` oraz `annual_vis`) są używane w kolejnym kroku jako podstawowe dane do określenia jakości potencjalnego przypisania obiektu do miesiąca/wariantu.  

Warto zaznaczyć, że `annual_vis` może być używana jako dodatkowy filtr – obiekty o bardzo niskiej liczbie użytecznych nocy mogą być zredukowane w liczbie lub w ogóle wykluczone z planu rocznego.

#### 3.5.2. Algorytm przydziału (Hungarian)

Funkcja `build_monthly_variants` wykorzystuje logikę hybrydową dwupoziomową opartą na _Hungarian Algorithm_ do optymalnego przydziału obiektów:

**Podział na dwie grupy**:

1. ELITA (Score > Mediana): obiekty o najwyższych punktach, priorytetyzowane przez Score^3 × Quality_Ratio. Dostają offset 1 000 000 000, co gwarantuje, że zawsze są wybierane jako pierwsze.
1. RESZTA (Score ≤ Mediana): obiekty o niższym score, priorytetyzowane przez prestiż katalogu (NGC/IC > Sh2 > RCW > LBN > Cederblad > PGC > Barnard > LDN). Score jest ignorowany. Waga katalogowa × 1000 + jakość widoczności × 100.

* Funkcja pomocnicza `get_catalog_weight` zwraca wagę numeryczną (10-90) na podstawie przynależności do katalogu, skanując zarówno `id` jak i `common_names` obiektu.

**Stałe warianty**: Liczba wariantów jest teraz sztywno ustalona na 3 (A, B, C), parametr block_size został usunięty z logiki (zachowany tylko dla kompatybilności interfejsu).

Na etapie wstępnej konfiguracji użytkownik lub konfiguracja (np. z pliku lub stałych w kodzie) definiuje:
`monthly_capacities` – słownik `month -> int`, mówiący o maksymalnej liczbie obiektów przypisanych w danym miesiącu (np. 9 obiektów na miesiąc),
`variants` – lista nazw wariantów (np. `["A", "B", "C"]`), co razem definiuje łączną liczbę slotów dla miesiąca: `total_monthly_slots = len(variants) * monthly_capacity`.

Całkowita liczba slotów w roku to `n_slots = sum(monthly_capacities.values()) * len(variants)`.

Funkcja:

- buduje listę wszystkich „slotów” w roku, reprezentowanych jako trójki `(id_variant, month, slot_index)`,
- tworzy macierz kosztów o wymiarach `n_objects × n_slots`,

Macierz inicjalizowana jest dużą wartością `INVALID_COST` dla wszystkich par `(obiekt, slot)`, zgodnie z konwencją `scipy.optimize.linear_sum_assignment` – węzły niedozwolone mają bardzo wysoki koszt, aby nie być wybierane.

Dla każdej pary `(obiekt, slot)`:

- jeśli obiekt w danym `month` ma `best_q_hours > min_avg_q_hours` (domyślnie jakiś graniczny próg, np. 1.0 godzina), to slot jest uznawany za „dopuszczalny”,
- wtedy:
  - obliczany jest `quality_ratio = best_q_hours / max_possible_hours` dla tego miesiąca (maksymalna możliwa jakość jako 100%),
  - obliczana jest waga `weighted_score = score^3 * quality_ratio`,
  - koszt przypisania jest ustawiany na `-weighted_score`, tak aby minimalizacja sumy kosztów odpowiadała maksymalizacji sumy `weighted_score`.

W praktyce:

- obiekty z wysokim `score` są bardzo mocno nagradzane (wzrost trzecią potęgą),
- obiekty z dobrą jakością widoczności w danym miesiącu (duży `quality_ratio`) otrzymują wyższe `weighted_score`,
- niezalecani obiekty (niższe `score`, niskie `best_q_hours` lub poniżej `min_avg_q_hours`) albo nie są w ogóle przypisane, albo trafiają dopiero „na ostatnią ławkę” przy wypełnianiu planu.

Funkcja `linear_sum_assignment` z modułu `scipy.optimize` zwraca optymalne przypisanie minimalizujące sumę kosztów, co równe jest maksymalizacji globalnej sumy `weighted_score` przy zachowaniu ograniczeń pojemnościowych na sloty miesięczne.  
Wynik jest przetworzony do listy `MonthlyAssignment` po jednym dla każdego wariantu, zawierającej `id` obiektów przypisanych do konkretnych miesięcy.

Dodatkowo, moduł pozwala na konfigurację ograniczeń typu:
- minimalna liczba obiektów w roku na dany typ (np. minimalna liczba mgławic, galaktyk, gromad),
- maksymalna liczba obiektów tego samego typu w danym miesiącu,
co jest realizowane przez dodatkowe warstwy filtrowania i korygujące iteracje po wyliczeniu pierwotnego rozwiązania.

#### 3.5.3. Raport i zapis wyników

Po wyznaczeniu optymalnego przydziału moduł generuje zestaw statystyk opisujących uzyskany plan roczny.

Statystyki obejmują m.in.:

- medianę `score` z przypisanych obiektów – wskazuje na ogólne poziome „wartości obserwacyjnej” planu,
- udział obiektów z górnej połowy rankingów `score` (np. top 50% obiektów z bazy) – pokazuje, jak bardzo plan koncentruje się na obiektach wysokiego priorytetu,
- udział obiektów spoza topu `score` – wskazuje, czy w planie zostają włączone „mniej ważne” cele, które jednak w danym roku mają wyjątkowo dobre warunki widoczności,
- średnią jakość okna obserwacyjnego (średnie `q_hours` z przypisanych slotów) – pokazuje, jak „wygodne” są warunki obserwacyjne dla danego planu.

Funkcja `save_selected_to_vis_data(vis_data, monthly_assignments, year)`:

- przejmuje strukturę `vis_data` z wcześniejszego wczytania,
- dla wszystkich obiektów zeruje lub czyszcza pola `selected` (ustawiając `null` albo usuwając),
- dla każdego obiektu przypisanego w którymś z wariantów i miesięcy:
  - tworzy strukturę `{"variant": "A", "month": 3, "assignment_date": "2026-03-01"}` (data jest przykładowa, zwykle wskazuje pierwszą noc w danym miesiącu),
  - zapisuje ją pod kluczem `selected` w danym obiekcie w `vis_data["objects"]`.
Na koniec funckja wywołuje `append_famous_labels_to_ids(vis_data)`, która dodaje do nazwy obiektu w formacie `"(M, C, H)"` informację o przynależności do Messier, Caldwell i Herschel (jeśli obiekt jest w jednym lub większej liczbie z tych katalogów), co jest później wykorzystywane w nagłówkach stron atlasu.

Moduł generuje również raport w formie `monthly_overview.pdf`.

Dla każdego miesiąca i każdego wariantu (A/B/C) w PDF robi się wykres:

- oś X – godziny nocy (lokalne),
- oś Y – wysokość obiektu, Słońca i Księżyca w stopniach,
- zaznaczone okna jakości obserwacji (z tłem kodującym fazę zmierzchu: civil, nautical, astronomical),
- na wykresie zaznaczone są punkty centralne nocy nowiu, w których obiekt osiąga `best_q_hours` w danym miesiącu.  

Wszystkie wykresy są ułożone na stronach A4 w formacie siatki (np. po 2–3 miesiące na stronie), a na końcu PDF znajduje się podsumowanie liczby obiektów per wariant oraz podstawowe statystyki planu rocznego (mediana, udział top, średnia jakość itp.).

W efekcie `4_plan_roczny.py` dostarcza zarówno strukturalny plan roczny (coded w `vis_data.json`), jak i wygodny, wizualny przegląd PDF, który jest później wykorzystywany w krokach `5_fov_and_maps.py` i `6_drukuj_strony_obiektow.py` do generowania stron atlasu.

* Funkcja `save_selected_to_vis_data` wywołuje `append_famous_labels_to_ids`, która dla wybranych obiektów dopisuje etykiety M/C/H do pola `name` (ale nie `id`), np. 'NGC 7000 (M33, C20)'

Raport końcowy zawiera:

1. Medianę Score i podział Elita vs Reszta
1. Listę odrzuconych obiektów z Top (score > mediana) z przyczynami
1. Listę obiektów użytych z puli poniżej mediany, posortowanych wg prestiżu katalogu
1. Średnią jakość okna obserwacyjnego (% względem najlepszej możliwej nocy)
1. Ocenę planu (WYBITNA / BARDZO DOBRA / DOBRA / KOMPROMISOWA)
1. Miesięczne obciążenie kalendarza ze statusem (ELITA / DOBRE / WYPEŁNIACZE)
1. Listę kompromisów – obiektów Top przypisanych do gorszych miesięcy"

### 3.5.4
Moduł generuje dodatkową stronę ze spisem wybranych obiektów na początku PDF. Funkcja `generate_summary_page` tworzy układ dwukolumnowy z tabelami zawierającymi nazwę obiektu, miesiąc przypisania i wariant. Obiekty są sortowane alfabetycznie. Wykorzystuje `matplotlib.table` z customowym formatowaniem.

---

### 3.6. `5_fov_and_maps.py`

**Cel:** generacja kadrów FOV i map kontekstowych dla obiektów oznaczonych jako `selected` w `vis_data.json` z użyciem biblioteki starplot.

`preload_open_ngc` tworzy przykładową mapę w projekcji Mercatora, ładując do pamięci katalog `OPEN_NGC`, gwiazdozbiory i gwiazdy do ~11 mag, co działa jako warm‑up backendu starplot.
`build_objects_from_vis_data` filtruje sekcję `objects` pliku `vis_data.json`, pozostawiając tylko obiekty z niepustym polem `selected` i zwraca je jako DataFrame.
`get_camera_params` pobiera konfigurację kamery z `vis_data["parameters"]["camera"]`, wypisuje ją w formie czytelnego komunikatu i zwraca jako słownik przekazywany do workerów.
`create_camera_object` tworzy obiekt `Camera` na podstawie słownika, z polami ogniskowej, pitch, rozmiaru matrycy i rozdzielczości.

Funkcje `_worker_render_fov` i `_worker_render_context` tworzą odpowiednio silniki `OpticMapEngine` i `ContextMapEngine` ze zdefiniowanymi stylami i zapisują wygenerowane obrazy do katalogu `starplots`.
`generate_fov_pngs` oraz `generate_context_pngs` budują listę zadań dla obiektów, dla których brakuje odpowiednich plików PNG, a następnie renderują je równolegle przy użyciu `ProcessPoolExecutor`, raportując postęp przez `tqdm`.
Funkcja `main` spina wszystko: preload danych, wczytanie `vis_data`, wybór obiektów i wywołanie generatorów FOV oraz map kontekstowych.

---

### 3.7. `6_drukuj_strony_obiektow.py`

**Cel:** wygenerowanie stron atlasu w formacie A4 z informacjami o obiekcie, wykresem nocy nowiu, kadrem FOV i rocznym wykresem widoczności.

Moduł `6_drukuj_strony_obiektow.py` odpowiada za konwersję obiektów z `vis_data.json` (w tym ich pola `selected` z `4_plan_roczny.py`) na wizualne strony A4 w formacie PDF, zgodne z konwencją atlasu astrofotograficznego.  
Każda strona odpowiada jednemu obiektowi i składa się z czterech głównych bloków: nagłówka z metadanymi, wykresu nocy nowiu, kadru FOV oraz rocznego wykresu widoczności.  
Moduł używa `matplotlib` do tworzenia wykresów oraz `reportlab` do układu strony i generowania finalnego PDF, a także wykorzystuje pliki `starplots/{id}.png` z `5_fov_and_maps.py` jako kadrów FOV.

#### 3.7.1. Ustawienia strony i layout

Moduł definiuje parametry strony A4 w centymetrach i calach, zgodne z typową konwencją PDF (np. 21 cm × 29.7 cm, 8.27 in × 11.69 in).  
Określane są również:
- uniwersalne `margins` (np. 2 cm od każdej krawędzi),
- `page_width` i `page_height` w calach,
- `plot_area` – obszar roboczy wewnątrz marginesów, który jest wykorzystywany do układu `matplotlib`.

Cała strona jest dzielona na cztery sekcje:
- nagłówek na górze strony,
- główny wykres nocy nowiu (w środku),
- kadr FOV z prawej lub lewej strony,
- roczny wykres widoczności na dole (lub w dalszej części strony).

Obszar dla każdej sekcji jest wyrażony jako `Rect` w współrzędnych „figury” Matplotlib, co pozwala na dokładne pozycjonowanie wszystkich elementów.

#### 3.7.2. Konfiguracja z `vis_data.json`

Funkcja `apply_config_from_vis_data(vis_data, config)` wczytuje z `vis_data.json` wartości konfiguracyjne i nadpisuje domyślne parametry modułu, takie jak:
- lokalizacja obserwacji (`location["name"]`, `location["lat"]`, `location["lon"]`),
- rok planowania (`year`),
- minimalna wysokość dla wykresów (`parameters["minalt"]`),
- parametry kamery (`parameters["camera"]` – ogniskowa, rozmiar matrycy, pitch, rozdzielczość, które są używane w `compute_diagonal_fov_arcmin`).

Wszystkie te wartości są używane do spójności między planem rocznym a stronami atlasu: np. wykres roczny używa tej samej wartości `minalt` i `sunlimit`, co moduł `3_wyliczenia.py`.

#### 3.7.3. Obliczanie FOV kamery

Funkcja `compute_diagonal_fov_arcmin(camera_params)` pobiera parametry kamery (np. `focal_length_mm`, `sensor_width_mm`, `sensor_height_mm`) i oblicza przekątną pola widzenia w arcmin i stopniach.  
Algorytmy:

- obliczony jest rozmiar przekątnej matrycy:  
  `diag_sensor = sqrt(sensor_width**2 + sensor_height**2)`,
- długość ogniskowej jest używana do wyliczenia kąta widzenia:  
  `fov_deg = 2 * arctan(diag_sensor / (2 * focal_length)) * 180/pi`,
- przekątna FOV w arcmin:  
  `fov_arcmin = fov_deg * 60`.

Wynik (`fov_arcmin`, `fov_deg`) jest później używany:
- jako tekst w tytule kadru FOV (`"FOV: {:.1f}′ (diagonal)"` lub podobnie),
- do kalibracji skali i legendy w kadrze.

#### 3.7.4. Formatowanie indeksów i nagłówka

Funkcja `format_indeksy(common_names, extra_info, max_length=120)` służy do przycinania długich ciągów indeksów i aliasów obiektu, aby zmieściły się w jednej linii nagłówka na stronie A4.  
Przykładowe operacje:
- łączenie pól `common_names` i `extra_info` w jeden tekst (np. `"Sh2-155, IC 1805, Heart Nebula"`),
- przycinanie do `max_length` znaków, z dodaniem `…` na końcu, jeśli tekst jest zbyt długi,
- usuwanie powtórzonych lub zduplikowanych nazw.

Wyjściowy tekst jest używany w nagłówku strony jako linia indeksów pod głównym `id` i nazwą obiektu.

#### 3.7.5. Budowa strony obiektu

Główna logika jest zawarta w funkcji `draw_object_page(pdf_canvas, obj, vis_data, observing_data, starplot_dir, year=2026)`.

Dla pojedynczego obiektu `obj` z `vis_data["objects"]` funkcja tworzy stronę A4 składającą się z:

1. **Nagłówek**  
   - `id` obiektu (np. `NGC 2237, Sh2-27, M 42` – z użyciem `append_famous_labels_to_ids` z `4_plan_roczny.py`),
   - nazwa obiektu (pierwszy element z `common_names`),
   - typ obiektu (np. `HII`, `DN`, `Gx`),
   - współrzędne `RA` i `Dec` w formacie `hh:mm:ss` / `dd:mm:ss`,
   - rozmiar kątowy `size` (np. `80′`),
   - jasność `mag` (np. `4.5`),
   - skrócone `indeksy` z `format_indeksy`,
   - informacja o `variant` i `month` z `selected` (np. `Wariant A, Marzec 2026`).

2. **Wykres nocy nowiu**  
   Funkcja tworzy `axes` w obszarze środkowym strony i rysuje:
   - krzywą wysokości obiektu w ciągu nocy (godziny lokalne na osi X, wysokość w stopniach na osi Y),
   - krzywą wysokości Słońca (z oznaczeniem granic zmierzchu: `minalt` oraz `sunlimit`),
   - półcześnie `qual_segments`, z tłem kodującym fazy nocy (np. ciemne tło dla nocy astronomicznej),
   - oznaczenie momentu nowiu (np. `New Moon`),
   - siatkę godzinową i linie graniczne (np. `minalt`, `sunlimit`).

   Wartości `q_hours` dla tej konkretnej nocy nowiu są pobierane z `observing_data` i zaznaczone na wykresie jako „okno jakości obserwacji”.

3. **Kadr FOV**  
   - funkcja wczytuje obraz `starplots/{id}.png` z katalogu `starplot_dir`,
   - skaluje go do odpowiedniego rozmiaru wewnątrz `plot_area` strony (np. prostokąt 6-in × 6-in),
   - dodaje tytuł nad kadrem zawierający informację o przekątnej FOV (np. `FOV: 120.5′ (diagonal)`),
   - ewentualnie dodaje etykiety kierunków `N/E/S/W` i skalę kątową w arcminach.

   Kadr ten pokazuje okolice obiektu z widoku kamery/obserwatora, z uwzględnieniem rzeczywistych gwiazd i obiektów wokół obiektu głównego.

4. **Roczny wykres widoczności**  
   Funkcja generuje `axes` z osiami:
   - X: dni roku (np. 1–365, z oznaczeniem miesięcy),
   - Y: godziny (0–24) z przesunięciem czasowym lokalnym (`tz_offset` z `observing_data`).

   Na wykresie:
   - zaznaczone jest tło reprezentujące nocne interwały (np. między `sunlimit` a `minalt`),
   - nałożone są `qual_segments` z `observing_data` dla tego obiektu – segmenty godzin o `q_hours > 0`,
   - oznaczone są noce nowiu (np. po 28 dniach),
   - ewentualnie zaznaczone są punkty `selected` z `vis_data` (np. okrągła ramka w dniu przypisania do wariantu/miesiąca).

   W rezultacie powstaje mapka „dzień–godzina”, z której widać, w jakie okna w roku dany obiekt ma dobre warunki obserwacyjne.

#### 3.7.6. Ostrzeżenia i czyszczenie wyjścia

Ograniczenia `Astropy` dotyczące modeli IERS i długiego zakresu czasu (np. dla 2026 i później) są obsługiwane przez globalne ustawienie `warnings.filterwarnings("ignore", module="astropy.iers")` na początku skryptu.  
To zapobiega „zatłaczeniu” logów komunikatami typu `Downloading remote IERS` i podobnymi, szczególnie przy generowaniu dużych partii stron obiektów.

Po zakończeniu `draw_object_page` moduł:
- zapisuje aktualną `figura` Matplotlib jako `Canvas` w `pdf_canvas`,
- dodaje strony `Astrophotography_Planner_2026_1.pdf` (część z rocznymi wykresami widoczności) oraz `Astrophotography_Planner_2026_2.pdf` (część z bardziej szczegółowymi strukturami, zależnie od konfiguracji),
- zamyka `plt.close("all")` po zakończeniu, aby uniknąć przecieku pamięci przy wielu stronach.

W efekcie `6_drukuj_strony_obiektow.py` dostarcza zestaw stron A4, które są potem scalane w `7_polacz_pliki_pdf.py` z tytułem i pustymi stronami do finalnego atlasu `Astrophotography_Planner_2026.pdf`.


---

### 3.8. `7_polacz_pliki_pdf.py`

**Cel:** wygenerowanie strony tytułowej atlasu i połączenie częściowych PDF w jeden końcowy plik.

Moduł wczytuje `vis_data.json`, aby pobrać nazwę lokalizacji, współrzędne i rok, następnie na tej podstawie buduje tekst „miejsce wydania” z aktualną datą.
Rejestrowane są czcionki Helvetica (z systemowego pliku TTF) w wariantach normal, bold i italic, co pozwala na poprawne wyświetlanie polskich znaków.
Przy użyciu ReportLab generowana jest strona tytułowa z napisem „Astrophotography Planner {rok}” oraz informacją o miejscu i dacie, zapisana jako `tytulowa.pdf`.
Następnie, przy użyciu `PdfWriter`, moduł składa finalny atlas: dodaje stronę tytułową, pustą stronę A4, dwie części atlasu (`Astrophotography_Planner_2026_1.pdf` i `_2.pdf`) oraz końcową pustą stronę, zapisując wynik do `Astrophotography_Planner_2026.pdf`.
Na koniec tymczasowy plik tytułowy jest usuwany, a w logu wypisywana jest informacja o pomyślnym scaleniu plików.

---

## 4. Konfiguracja i uruchamianie

### 4.1. Główne punkty konfiguracji

- `2_ograniczenie_katalogu.py` – ustawienia lokalizacji, roku, progów `minalt` i `sunlimit`, klasy Bortle oraz parametrów kamery.
- `3_wyliczenia.py` – domyślny limit liczby obiektów do przeliczenia, liczba workerów `max_workers`, flaga `force_all` do wymuszenia pełnego przeliczenia.
- `4_plan_roczny.py` – liczba wariantów, pojemność miesięczna (`per_month_capacity`), próg `min_avg_q_hours` używany przy budowie macierzy kosztów.

### 4.2. Typowa sekwencja uruchomienia

1. `python 0_opracuj_katalog_ngc.py` – jednorazowo, po pobraniu OpenNGC.
2. `python 1_generuj_katalog_astro.py` – pobranie katalogów z VizieR i Smart Merge.
3. `python 2_ograniczenie_katalogu.py` – interaktywna konfiguracja, scoring i zapis `vis_data.json`.
4. `python 3_wyliczenia.py [--force-all]` – budowa cache RAW/FINAL widoczności z opcjonalnym wymuszeniem pełnego przeliczenia.
5. `python 4_plan_roczny.py` – optymalny plan roczny i aktualizacja `vis_data.json` o `selected`.
6. `python 5_fov_and_maps.py` – generacja FOV i map kontekstowych dla wybranych obiektów.
7. `python 6_drukuj_strony_obiektow.py` – wydruk stron atlasu dla wariantów.
8. `python 7_polacz_pliki_pdf.py` – wygenerowanie finalnego atlasu PDF.

---

## 5. Możliwe kierunki rozwoju

- Wydzielenie konfiguracji do wspólnego pliku `config.yaml` lub `config.toml` i rezygnacja z części interaktywnych zapytań `input()` na rzecz odczytu konfiguracji.
- Dodanie testów jednostkowych dla kluczowych funkcji takich jak `smart_merge`, `analyze_famous_status`, `mask_to_segments` oraz budowa macierzy kosztów w `build_monthly_variants`.
- Uspójnienie logowania i wprowadzenie poziomów logów (INFO/WARN/ERROR) zamiast wyłącznie `print`, co ułatwi analizę pracy systemu.
