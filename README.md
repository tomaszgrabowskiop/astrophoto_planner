# Astrophotography Planner & Atlas Generator

Zestaw skryptów w Python służący do generowania spersonalizowanego rocznego planera astronomicznego oraz atlasu obiektów głębokiego nieba (DSO). Dla niezaawansowanych amatorów nieba. Atomatycznie wybiera z bazy 90k+ obiektów te, których parametry określa użytkownik. Obiekty są "rozkładane" na przestrzeni roku tak, żeby stworzyć najbardziej zachłanny plan astrofotograficzny.

[Przykładowy planer](doc/AstroPhotography_Planner_2026_Poznań_compressed.pdf) - tu w wersji skompresowanej - nie jest tak ładny, jak oryginał, ale daje obraz całości.

Możesz wywołać `run.sh`, który założy środowisko i będzie uruchamiał kolejne skrypty. Pozwoli Ci też pominąć wybrane kroki. Nie jest najwygodniejszą formą używania skryptów, ale moze się przydać.

Pełna dokumentacja znajduje się w pliku:
[Dokumentacja techniczna w formacie md](doc/Dokumentacja_techniczna.md)

**Powodzenia!**

---

## Szczegóły działania skryptów

System **pobiera dane** z katalogów astronomicznych, **filtruje** je pod kątem lokalizacji obserwatora i posiadanego sprzętu (teleskop/kamera), **oblicza** precyzyjną widoczność na dany rok, **dystrybuuje** po miesiącach, a następnie generuje profesjonalny plik PDF zawierający: 

1. Przegląd roczny (kiedy obserwować dany obiekt).

![Widok propozycji miesięcznych](doc/Planner_roczny.png)

2. Szczegółowe strony dla każdego obiektu (wykres wysokości, kadr FOV, widoczność w ciągu roku, noce z/bez Księżyca, mapa kontekstowa). 

![Strona obiektu](doc/Strona_obiektu.png)

---

## 🚀 Możliwości

- **Agregacja** danych: łączy katalogi NGC/IC, Sharpless (Sh2), RCW, Barnard, LBN, LDN, Cederblad i PGC. 
- **Inteligentne filtrowanie**: wybiera obiekty na podstawie roku, szerokości geograficznej, długości okna obserwacyjnego, minimalnej wysokości obiektu nad horyzontem, określenia zmierzchu, jasności (Mag), rozmiaru (Size) oraz skali Bortle i ewentualności używania filtrów wąskopasmowych.
- **Symulacja FOV**: generuje symulacje kadru (Field of View) dla kamery i teleskopu przy użyciu biblioteki `starplot`.
- **Mapa kontekstowa** dająca szersze pole widzenia obiektu na niebie. Mapy wykorzystują różną projekcję w zależności od wysokosci obiektu.
- **Obliczenia astronomiczne**: wylicza widoczność w ciągu roku w zależności od podanego progu zmierzchu i długości okna obserwacyjnego; widoczność w skali roku, godziny z/bez Księżyca, wysokość górowania w przykładowej nocy na tle zmierzchu i wschodu (cywilnych, żeglarskich i astronomicznych). 
- Format PDF: generuje gotowy do druku atlas w formacie A4. 

---

## 🛠️ Wymagania

Projekt wymaga Pythona 3.10+ oraz następujących bibliotek: 

```bash
pip install pandas numpy astropy astroplan matplotlib reportlab pypdf tqdm astroquery starplot networkx
```

**Uwaga:** Biblioteka `starplot` może wymagać dodatkowej konfiguracji (pobrania danych gwiazd). Skrypt pobiera je automatycznie, ale może się okazać, że przy wyborze specyficznych parametrów, będzie trzeba "dociągnąć" coś jeszcze. 

---

## 📂 Struktura plików i dane wejściowe

Aby rozpocząć, upewnij się, że posiadasz plik źródłowy dla katalogu NGC. 

- `OpenNGC/NGC.csv` – plik CSV z poszerzonymi o dodatkowe nazwy zwyczajowe danymi OpenNGC (wymagany przez skrypt `0_opracuj_katalog_ngc.py`).
- Jeżeli chcesz dodać kolejne nazwy zwyczajowe, możesz wyedytować `uzupelnij_openngc` i uruchomić. Powstanie nowy plik `NGC_updated.csv`, którym możesz zastąpić `OpenNGC/NGC.csv`

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
- Wykonuje „Smart Merge” - łączenie duplikatów i obiektów blisko siebie według oczekiwania użytkownika (zakres od 1 do 60 arcmin).
- Usuwa "szum" - według kryteriów użytkownika. Obiekty mniejsze niż Size i ciemniejsze niż Mag. (Domyślnie: size mniejszy niż 5' i magnitudo ciemniejsze niż 18)
- Tworzy plik `katalog_astro_full.csv`. 
- Opcjonalnie: uruchom `analiza_katalog.py`, aby sprawdzić statystyki bazy. 

### Krok 2: Konfiguracja i selekcja obiektów

Uruchom: 

```bash
python 2_ograniczenie_katalogu.py
```

- Interaktywny skrypt: pyta o rok, lokalizację, Bortle, parametry teleskopu/kamery, filtry narrowband oraz minimalną wysokość obiektu, próg wysokości słońca, czas trwania okna obserwacyjnego. Pytania o filtry i Bortle służą selekcji obiektów, które warto fotografować bądź na ciemnym niebie, bądź jedynie z filtrami narrowband.
- Filtruje bazę według parametrów użytkownika i określonego minimalnego rozmiaru i jasnosci dla obiektów. 
- Tworzy plik `vis_data.json` z kandydatami do atlasu. 

### Krok 3: Silnik obliczeniowy (Engine)

Uruchom: 

```bash
python 3_wyliczenia.py
```

- Wykonuje ciężkie obliczenia astronomiczne dzięki AstroPy (równolegle na wielu rdzeniach CPU). To zajmuje czas!
- Wylicza dokładną widoczność minuta po minucie dla całego roku. 
- Zapisuje wyniki do `observing_data.pkl`.
- Możliwe przyrostowe **wywołanie bez powtarzania obliczeń** (sam sprawdza, czy zmieniły się parametry). 
	- **UWAGA**: Za każdym razem kiedy zostanie uruchomiony `2_ograniczenie_katalogu.py`, skrypt `3_wyliczenia` sprawdza nowe parametry. 
	- Jeżeli lokalizacja lub rok zostały zmienione, wówczas wszystkie wybrane obiekty są przeliczane ponownie "od zera". 
	- Jeżeli zmieniły się parametry okna obserwacyjnego, próg zmierzchu, wówczas obliczenia ograniczają się do zmiany maski (trwaja krócej).
	- Jeżeli w json pojawiły się nowe pozycje, wyliczenia są ograniczone tylko do tych obiektów (o ile nie zmieniła się lokalizacja lub rok).
- Podjęćie obliczeń tylko w zakresie maski (parametry: wysokość nad horyzontem, długość okna obserwacyhnego, określenie zmierzchu/świtu).
- Pełne obliczenia dla obiektów, których nie było wcześniej i przy zmianie lokalizacji.  

### Krok 4: Plan roczny i wybór wariantów

Uruchom: 

```bash
python 4_plan_roczny.py
```

- Analizuje dane z kroku 3. 
- Przydziela obiekty do miesięcy (Warianty A, B, C), aby zbalansować sesje obserwacyjne. Użytkownik określa ile obiektów chce przypisać do każdego z wariantów.
- Generuje **Część** 1 PDF: `Astrophotography_Planner_ROK_1.pdf` (wykresy miesieczne, spis obiektów). 
- Aktualizuje `vis_data.json` o flagę `selected`. 

### Krok 5: Generowanie map nieba

Uruchom: 

```bash
python 5_fov_and_maps.py
```

- Korzysta z biblioteki `starplot`. 
- Generuje pliki PNG w katalogu `starplots/`: 
  - FOV czyli kadry optyczne (symulacja kamery). 
  - Mapy kontekstowe (szersze pole widzenia).

Pliki są w generowane w wysokiej rozdzielczości. Można to zmienić wewnątrz skryptu.

### Krok 6: Generowanie stron obiektów

Uruchom: 

```bash
python 6_drukuj_strony_obiektow.py
```

- Składa szczegółowe strony dla każdego wybranego obiektu. 
- Zawierają: 
	- wykresy wysokości w najlepszej nocy (miesiąc według wariantu, w którym obiekt się pojawia) plus informacja, która noc w roku daje najdłuższe okno obserwacyjne, 
	- FOV, 
	- wykres rocznej widoczności z zaznaczeniem nocy według podanego wcześniej progu (w `2_ograniczenie_katalogu.py`), uwzględnia zmianę czasu z/na letni/zimowy,
	- wykres liczby godzin z/bez księżyca w ciągu roku,
	- mapę kontekstową. 
- Tworzy **Część** 2 PDF: `Astrophotography_Planner_ROK_2.pdf`. 

### Krok 7: Finalizacja

Uruchom: 

```bash
python 7_polacz_pliki_pdf.py
```

- Generuje stronę tytułową. 
- Łączy część 1 i część 2 w jeden plik. 
- Wynik końcowy: `Astrophotography_Planner_ROK_MIASTO.pdf`. 

---

## 📝 Uwagi dodatkowe

- Wydajność: krok 3 i 5 wykorzystują wielowątkowość (`multiprocessing`). Mimo to obliczenia AstroPy i generowanie map może zająć sporo czasu w zależności od liczby obiektów i wydajności komputera.
- TimeZone jest przypisywane na podstawie lokalizacji. Wpływa na obliczenia oraz  na wykres widoczności w skali roku.  

---

## 📄 Licencja

Projekt do użytku własnego. Korzysta z danych OpenNGC oraz serwisów VizieR.

[OpenNGC](https://github.com/mattiaverga/OpenNGC/tree/master) to osobny projekt częściowo wykorzystywany w tym repozytorium.

---

Jeżeli masz propozycje zmian, widzisz błędy, napisz: <morus@dominikanie.pl>. 
