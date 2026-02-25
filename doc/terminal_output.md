
**% python3 1_build_catalog.py**
```
=======================================================================================================================
============================================ PARAMETRY FILTROWANIA KATALOGU =========================================== 
=======================================================================================================================
  (naciśnij Enter, aby zaakceptować wartość domyślną)
  Minimalny rozmiar obiektu  [arcmin] [domyślnie: 5.0]: 
  MInimalna jasność (granica słaba) [mag] [domyślnie: 17.0]: 18
  Maksymalny dopuszczalny rozmiar obiektu (cap) [arcmin] [domyślnie: 300.0]: 

================================== KROK 1: Pobieranie bazy danych (ENTITY RESOLUTION) =================================
[INFO] Pobieranie katalogów
       > NGC/IC z data/NGC.csv
       > Sharpless (VII/20/catalog)
       > Barnard (VII/220A)
       > RCW (VII/216)
       > PGC (VII/119)
       > Lynds Dark (VII/7A)
       > Lynds Bright (VII/9)
       > Cederblad (VII/231)

       PODSUMOWANIE POBIERANIA:
       NGC     :   13 969 wierszy
       SH2     :      313 wierszy
       BARN    :      349 wierszy
       RCW     :      181 wierszy
       PGC     :   77 141 wierszy
       LDN     :    1 791 wierszy
       LBN     :    1 125 wierszy
       CED     :      330 wierszy
       ==========================
       Razem   :   95 199 wierszy

================================================= KROK 2: Normalizacja ================================================
[INFO] Normalizacja katalogów (unifikacja kolumn)
       > Normalizacja NGC/IC
       > Normalizacja Sharpless
       > Normalizacja Barnard
       > Normalizacja RCW
       > Normalizacja Cederblad
       > Normalizacja LBN
       > Normalizacja LDN
       > Normalizacja PGC
================================================= KROK 3: Filtrowanie =================================================
       Odrzucono 10 871 obiektów (Mag > 24.0 lub Size < 0.5').
       Do analizy merge trafia: 79 723 obiektów.

[INFO] Analiza tolerancji dla Smart Merge

         Tol | Pary X-kat |  Klastry | Podejrzane
       ------------------------------------------
          1' |      8 851 |    7 917 |         21
          2' |     10 442 |    8 234 |         31
          3' |     11 826 |    8 138 |         41
          5' |     14 983 |    7 866 |         63
          8' |     20 267 |    7 545 |         88
         12' |     28 237 |    7 073 |        138
         20' |     46 996 |    6 270 |        195

       [AUTO] Wybrana tolerancja: 3' (arcmin)
       [AUTO] Tolerancja w stopniach: 0.05000°

[INFO] Smart Merge (tolerancja 3.00 arcmin, 79 723 obiektów)
       Klastry jednoobiektowe : 52 542
       Klastry wieloobiektowe : 11 584
       Łącznie klastrów       : 64 126
       Merge klastrów: 100%|███████████████████████████████████████████████| 64126/64126 [02:02<00:00, 523.34klaster/s]

       Po merge: 64 126 obiektów.

[INFO] Fuzja common_names: 12 grup do scalenia
       64126 → 64114 obiektów (scalono 12 duplikatów nazewniczych)
========================================== KROK 4: Imputacja danych mag/size ==========================================

[INFO] Imputacja brakujących Size i Mag (mediana wg typu)

            Typ | Mediana size | Mediana mag
       ------------------------------------
              * |         1.1' |      11.50 
             ** |         0.9' |      11.17 
           *Ass |         4.8' |      10.43 
           Cl+N |         6.9' |       9.88 
             DN |        13.7' |        brak
            EmN |         3.2' |      12.21 
              G |         1.0' |      15.50 
            GCl |         3.8' |      10.36 
         GGroup |         1.1' |      23.82 
          GPair |         1.2' |      23.25 
          GTrpl |         1.0' |      22.82 
            HII |        12.0' |      11.00 
             NB |        35.0' |       6.66 
            Neb |        15.1' |       8.97 
           Nova |         brak |       7.00 
            OCl |         4.6' |       9.10 
          Other |         1.4' |       8.90 
             PN |         0.9' |      11.25 
            RfN |         5.0' |       9.00 
            SNR |        25.1' |       8.40 

       Uzupełniono size : 17 462 obiektów
       Uzupełniono mag  : 7 298 obiektów
       Bez imputacji mag (DN) : pozostają z NaN
================================================= KROK 5: Filtrowanie =================================================

[INFO] Filtr docelowy (Size >= 5.0', Mag <= 18.0)
       Odrzucono ze względu na mały rozmiar : 60 128
       Odrzucono ze względu na słabą jasność: 4
================================================= KROK 2: Zapis danych ================================================
[INFO] ZAPIS DO PLIKU: data/katalog_astro_full.csv
       Sukces! Zapisano 3 982 obiektów.
```

**[INFO] Uruchomić pełną analizę katalogu? [y/n] (domyślnie y):**
```
═══════════════════════════════════════════════════════════════════════════════════════════════════════════════════════
       ANALIZA KATALOGU ASTRONOMICZNEGO  ›  data/katalog_astro_full.csv
═══════════════════════════════════════════════════════════════════════════════════════════════════════════════════════

[01] OGÓLNE STATYSTYKI
     ──────────────────────────────────────────────────────────────────────────────────────────────────────────────────
     Liczba obiektów : 3,982
     Liczba kolumn   : 15

[02] JASNOŚĆ CAŁKOWITA – kolumna 'mag'
     ──────────────────────────────────────────────────────────────────────────────────────────────────────────────────
     Obecne   : 2,270 (57.0%)
     Brakujące: 1,712 (43.0%)

     Rozkład jasności:
       mag < 5        :    85 (  3.7%)  – bardzo jasne
       5  ≤ mag < 8   : 1,203 ( 53.0%)  – jasne
       8  ≤ mag < 12  :   932 ( 41.1%)  – średnie
       mag ≥ 12       :    50 (  2.2%)  – słabe

[03] ROZMIAR KĄTOWY – kolumna 'size' [arcmin]
     ──────────────────────────────────────────────────────────────────────────────────────────────────────────────────
     Obecne   : 3,980 (99.9%)
     Brakujące: 2 (0.1%)

     Rozkład rozmiaru:
       < 5'             :     0 (  0.0%)  – bardzo małe [crop]
       5' – 20'         : 2,084 ( 52.4%)  – małe/średnie
       20' – 60'        : 1,062 ( 26.7%)  – dobre
       60' – 180'  (>1°):   664 ( 16.7%)  – duże
       ≥ 180'           :   170 (  4.3%)  – bardzo duże

[08] KOMPLETNOŚĆ DANYCH
     ──────────────────────────────────────────────────────────────────────────────────────────────────────────────────
     Pełne  (ra + dec + type + mag + size)           2,268 / 3,982  ( 57.0%)
     Pozycja + typ + rozmiar  (ra + dec + type + size)  3,980 / 3,982  ( 99.9%)
     Tylko pozycja  (ra + dec + type)                3,982 / 3,982  (100.0%)
     Mag ✓  ale brak size                                2 / 3,982  (  0.1%)
     Size ✓  ale brak mag                            1,712 / 3,982  ( 43.0%)

[05] KATALOG ŹRÓDŁOWY – kolumna 'catalog'
     ──────────────────────────────────────────────────────────────────────────────────────────────────────────────────
     ldn             1,423 ( 35.7%)
     lbn               931 ( 23.4%)
     ngc               806 ( 20.2%)
     barn              289 (  7.3%)
     sh2               213 (  5.3%)
     ced               162 (  4.1%)
     rcw               117 (  2.9%)
     pgc                41 (  1.0%)

[06] PRZYNALEŻNOŚĆ DO KATALOGÓW – Messier / Caldwell / Herschel
     ──────────────────────────────────────────────────────────────────────────────────────────────────────────────────
     Messier  (M1–M110)  :   99  (  2.5%)
     Caldwell (C1–C109)  :   79  (  2.0%)
     Herschel (H1–H400)  :  186  (  4.7%)

[07] POPULARNE NAZWY – kolumna 'common_names'
     ──────────────────────────────────────────────────────────────────────────────────────────────────────────────────
     Z nazwą popularną: 132 (3.3%)

     Przykłady (pierwsze 10):
       NGC4594          Sombrero Galaxy                           [G]
       NGC5457          Pinwheel Galaxy                           [G]
       NGC4254          Coma Pinwheel, Virgo Cluster Pinwheel     [G]
       NGC4486          Virgo Galaxy                              [G]
       NGC5236          Southern Pinwheel Galaxy                  [G]
       NGC3034          Cigar Galaxy                              [G]
       NGC3031          Bode's Galaxy                             [G]
       NGC2682          King Cobra Cluster                        [OCl]
       NGC4826          Black Eye Galaxy, Evil Eye Galaxy         [G]
       NGC5055          Sunflower Galaxy                          [G]

[10] STATYSTYKI PER KLASA OBIEKTU
     ──────────────────────────────────────────────────────────────────────────────────────────────────────────────────
  Klasa                                     N      %   mag śr  mag med   size śr   size med
  ─────────────────────────────────────────────────────────────────────────────────────────
  Galaktyki (G / GPair / GGroup)          398   10.0%     10.1     10.3     10.9'       6.8'
  Mgławice Planetarne (PN)                  3    0.1%      8.0      7.4      9.8'       6.7'
  Mgławice Emisyjne (HII / EmN)           363    9.1%     10.9     11.0     51.6'      20.0'
  Mgławice Refleksyjne (NB / RfN)        1113   28.0%      6.6      6.7     62.9'      38.0'
  Ciemne Mgławice (DN / DrkN)            1712   43.0%        –        –     34.0'      16.9'
  Pozostałości Supernowej (SNR)             8    0.2%      8.5      8.4     63.2'      37.6'

═══════════════════════════════════════════════════════════════════════════════════════════════════════════════════════
```
**% python3 2_plan_and_score.py**
```
=======================================================================================================================
KROK 1: WYBÓR ROKU I LOKALIZACJI
=======================================================================================================================
Enter current year lub wpisz wybrany pomiędzy 2000 a 2100 (Enter = 2026): 
Rok: 2026

[INFO] Wybierz opcję lub wpisz własne miejsce:
1. Poznań, Polska
2. Kraków, Polska
3. Wpisz lokalizację
[1/2/3, domyślnie 1]: 
Wybrano Poznań, Polska (52.4095°N, 16.9319°E)

=======================================================================================================================
KROK 2: PARAMETRY WIDOCZNOŚCI I FOV
=======================================================================================================================
A. Minimalna wysokość obiektu nad horyzontem [domyślnie 25°]: 22.5

B. Minimalna liczba godzin, którą obiekt jest widoczny w nocy powyżej progu wysokości [domyślnie 3]: 2.5

C. Ciemność nieba – kąt słońca pod horyzontem:
   • zmierzch cywilny (-6°)
   • zmierzch żeglarski (-12°)
   • zmierzch astronomiczny (-18°)
   Możesz wpisać dowolną sensowną wartość (liczba stopni pod horyzontem) [domyślnie 12]: 

D. Określenie FOV
   Użyć domyślnego setupu RedCat61 + ASI2600MC Pro? [t/n, domyślnie t]: 

[INFO] Szacowany FOV: 4.49° x 3.00°

E. Minimalny rozmiar obiektu jako % krótszego boku FOV [domyślnie 10]: 
[INFO] Minimalny rozmiar obiektu (arcmin): 18.0'.

F. Skala Bortle – określenie zanieczyszczenia światłem:
   1) Bortle 1–3 (wieś, ciemne niebo, pomijalne LP)
   2) Bortle 4–5 (przedmieścia, umiarkowane LP)
   3) Bortle 6–7 (miasto, silne LP)
   4) Bortle 8–9 (centrum miasta, ekstremalne LP)
   Wybierz Twój stopień zanieczyszczenia nieba światłem [domyślnie 4]: 

G. Czy zamierzasz korzystać z filtrów narrowband (H/O/S)? [t/n, domyślnie n]: 

H. Czy w wyborze obiektów premiować katalogi Messier/Caldwell/Herschel? [t/n, domyślnie t]: 

=======================================================================================================================
KROK 3. ŁADOWANIE I FILTROWANIE
=======================================================================================================================

[INFO] Filtr geometryczny (Alt >= 22.5°, Visible >= 2.5h)

[INFO] Po filtrze min. wysokości (22.5°) pozostało 2 811 z 3 982 obiektów.

[INFO] Generowanie SMART GRID dla roku 2026, limit słońca: -12.0°
       SmartGrid: 100%|██████████████████████████████████████████████████████████████| 365/365 [00:17<00:00, 20.87day/s]
       ✓ Wygenerowano 12 967 punktów pomiarowych dla 365 nocy (~3242h obserwacji).

[INFO] Obliczanie macierzowe dla 2 811 obiektów. (To może potrwać kilkanaście sekund.)
       ✓ Obliczono widoczność per noc i per miesiąc.
       ✓ Spośród 2 811 obiektów, 2 759 przynajmniej w jedną noc są widoczne przez >= 2.5h.

=======================================================================================================================
KROK 4. PUNKTACJA (SCORING) OBIEKTÓW
=======================================================================================================================

[INFO] Zastosowano ostre cięcie dla obiektw z wynikiem < 10.0 pkt.
[INFO] Pozostało 2 759 obiektów z 2 759.

[INFO] Soft-cut: lista < 5000 obiektów — pomijam cięcie percentylowe.

=======================================================================================================================
KROK 5. ŁĄCZENIE OBIEKTÓW W KADRY (CLUSTERING)
=======================================================================================================================

[INFO] Optymalizacja kadrów (Clustering z dynamicznym FOV)
       • FOV sprzętu: 4.49° x 3.00°
       • Promień klastrowania ustalony na: 2.40°
       ✓ Zredukowano kadrów z 2 759 do 618.

[INFO] Ostateczna selekcja punktowa zoptymalizowanych kadrów...

[INFO] Zastosowano ostre cięcie dla obiektw z wynikiem < 50.0 pkt.
       Pozostało 522 obiektów z 618.
[INFO] Wartości wyliczone przez cięcie adaptacyjne: final_score >= p30 (59.0).
       Pozostało 371 obiektów z 522.

[INFO] Hard Filter Size >= 18.0'
       Pozostało 214 finalnych kadrów.

=======================================================================================================================
KROK 6. ZAPIS DANYCH DO VIS_DATA.JSON
=======================================================================================================================

=======================================================================================================================
 Nr ID           Extra                     Common name               Typ           RA     Dec     Mag     Size  Score
-----------------------------------------------------------------------------------------------------------------------
  1. NGC1976      NGC1973, NGC1975, NGC1977 Great Orion Nebula, Lower Cl+N        83.8    -5.4    ●4.0  ●240.0'  110.0
  2. NGC6705      B104, B107, B110, B111, B Amas de l'Ecu de Sobieski OCl        282.8    -6.3    ●5.8  ●120.0'  110.0
  3. NGC7654      NGC7538, NGC7635, Sh2-157 Bubble Nebula             OCl        351.2    61.6    ●6.9  ●224.5'  110.0
  4. NGC1912      NGC1907, NGC1931, NGC1960 Fly Nebula, Pinwheel Clus OCl         82.2    35.9    ●6.4  ●300.0'  110.0
  5. NGC2323      NGC2306, NGC2309, NGC2335 Seagull Nebula head       OCl        105.7    -8.4    ●5.9   ●90.0'  110.0
  6. NGC7092      NGC7082, Sh2-124, LBN416, Pyramid Cluster           OCl        323.0    48.4    ●4.6  ●147.0'  110.0
  7. NGC7023      Sh2-136, LBN468, LBN475,  Ghost Nebula, Iris Nebula Neb        315.4    68.2    ●7.2  ●150.0'  105.0
  8. IC0405       NGC1893, IC0410, Sh2-229, Flaming Star Nebula, Tadp Neb         79.1    34.4   ●10.0  ●240.0'  105.0
  9. NGC6694      NGC6712, B314, LDN433.0,                            OCl        281.3    -9.4    ●8.9  ●103.9'  105.0
 10. NGC7000      NGC6989, NGC6997, NGC7039 N American & Pelican Nebu HII        314.8    44.5    ●4.0  ●240.0'  105.0
 11. NGC6888      IC4996, Sh2-104, Sh2-105, Crescent Nebula           HII        303.0    38.4    ●7.4  ●300.0'  105.0
 12. NGC0559      NGC0654, Sh2-187, LBN630,                           OCl         22.4    63.3    ●9.5   ●31.9'  105.0
 13. NGC2239      NGC2236, NGC2237, NGC2238 Rosette A, Rosette B, Ros Cl+N        98.0     4.9    ●4.8  ●264.3'  105.0
 14. IC5146       NGC7209, Sh2-125, LBN424, Cocoon Nebula             Cl+N       328.4    47.3    ●7.2   ●26.8'  103.5
 15. NGC2099      IC0439, Ced60, B34, LDN15                           OCl         88.1    32.6    ●5.6   ●35.0'  102.0
 16. NGC2682      NGC2664, NGC2678          King Cobra Cluster        OCl        132.8    11.8    ●6.9   ●33.0'  101.0
 17. NGC6838      Sh2-84, LBN131, LDN745.0,                           GCl        298.4    18.8    ●6.1   ●61.2'  101.0
 18. NGC6171      LBN12, LBN16, LDN145.0, L                           GCl        248.1   -13.1    ●8.8  ●150.0'  100.0
 19. NGC2264      NGC2251, NGC2259, Sh2-273 Christmas Tree Cluster, F Cl+N       100.2     9.9    ●3.9  ●250.0'  100.0
 20. NGC1528      NGC1491, NGC1513, Sh2-206                           OCl         63.8    51.2    ●6.4   ●50.0'  100.0
 21. NGC5457      NGC5447, NGC5450, NGC5453 Pinwheel Galaxy           G          210.8    54.3    ●7.9   ●28.5'  100.0
 22. NGC7380      Sh2-139, Sh2-142, Sh2-151 Wizard Nebula             Cl+N       341.8    58.1    ●7.2  ●180.0'  100.0
 23. NGC2343      NGC2353, Sh2-293, Sh2-294 Seagull Nebula wings      OCl        107.0   -10.6    ●6.7  ●200.0'  100.0
 24. NGC6910      IC1311, Sh2-108, Sh2-109, Sadr Region               OCl        305.8    40.8    ●7.4  ●300.0'  100.0
 25. NGC7044      Sh2-119, LBN340, LBN360,                            OCl        318.3    42.5   ●12.0  ●171.1'  100.0
 26. NGC6664      NGC6649, Sh2-56, Sh2-57,                            OCl        279.1    -8.2    ●7.8  ●173.9'  100.0
 27. NGC0129      IC0010, Sh2-173, Sh2-177,                           OCl          7.5    60.2    ●6.5   ●40.0'  100.0
 28. NGC1027      NGC0896, IC1795, IC1805,  Heart Nebula, Maffei 2, S OCl         40.6    61.6    ●6.7  ●207.8'  100.0
 29. NGC7142      NGC7129, IC5132, IC5133,                            OCl        326.3    65.8    ●9.3  ●120.0'  100.0
 30. NGC0381      NGC0366, IC0059, IC0063,                            OCl         17.1    61.6    ●9.3  ●199.0'  100.0
 31. NGC6866      LBN232, LBN237, LBN238, L                           OCl        301.0    44.2    ●7.6  ●300.0'  100.0
 32. NGC6823      Sh2-86, Sh2-87, Sh2-88, L                           Cl+N       295.8    23.3    ●7.1   ●43.6'   99.5
 33. NGC1039      NGC1003, PGC10052                                   OCl         40.5    42.7    ●5.2   ●22.5'   98.5
 34. NGC2194      NGC2169, Sh2-267, Sh2-268                           OCl         93.4    12.8    ●8.5   ●60.0'   98.5
 35. NGC2168      NGC2158, LDN1564.0        Shoe-Buckle Cluster       OCl         92.3    24.3    ●5.1   ●24.0'   98.5
 36. NGC0457      NGC0436, Sh2-188, LBN633, Owl Cluster               OCl         19.9    58.3    ●6.4   ●63.2'   97.5
 37. NGC7078      LBN152                    Great Pegasus Cluster     GCl        322.5    12.2    ●6.3  ●190.0'   97.0
 38. NGC6853      LDN787.0, LDN797.0, LDN80 Dumbbell Nebula           PN         299.9    22.7    ●7.4  ●293.9'   96.5
 39. NGC3031      NGC2976, NGC3034, NGC3077 Bode's Galaxy, Cigar Gala G          148.9    69.1    ●6.9   ●24.9'   96.5
 40. NGC1952      Sh2-243, Sh2-244, LBN830, Crab Nebula               SNR         83.6    22.0    ●8.4   ●69.7'   95.5
 41. NGC1342      LBN716, LBN718, LBN719, L                           OCl         52.9    37.4    ●6.7  ●210.0'   95.5
 42. NGC2632                                Beehive, Praesepe Cluster OCl        130.1    19.7    ●3.1  ●108.6'   95.0
 43. NGC1432      NGC1435, IC0336, IC0349,  Barnard's Merope Nebula,  HII         56.5    24.4   ○11.0  ●300.0'   95.0
 44. NGC2548                                                          OCl        123.4    -5.8    ●5.8   ●28.2'   95.0
 45. NGC2301      Sh2-284, LBN983, LBN984,  Great Bird Cluster        OCl        102.9     0.5    ●6.0   ●80.0'   95.0
 46. NGC6755      Sh2-72, RCW179, LBN104, B                           OCl        287.0     4.3    ●7.5   ●70.2'   94.0
 47. NGC0224      NGC0205, NGC0221, PGC2429 Andromeda Galaxy          G           10.7    41.3    ●3.4  ●189.1'   92.5
 48. NGC0188      LBN617                    North Celestial Pole Clus OCl         11.9    85.3    ●8.1  ●300.0'   92.0
 49. NGC6882      LDN809.0, LDN814.0, LDN81                           OCl        303.0    26.5   ●14.1  ●189.7'   91.5
 50. NGC0752                                                          OCl         29.4    37.8    ●5.7   ●39.0'   90.0
 51. NGC1664      LBN755, B25, LDN1461.0, L                           OCl         72.8    43.7    ●7.6  ●129.2'   90.0
 52. NGC6960      NGC6974, NGC6979, IC1340, Filamentary Nebula, Veil  SNR        311.5    30.6    ●7.0  ●210.0'   90.0
 53. NGC6939      NGC6946, PGC65001, B150,  Fireworks Galaxy          OCl        307.9    60.7    ●7.8   ●60.0'   90.0
 54. NGC7789      LBN562, LDN1253.0, LDN125                           OCl        359.4    56.7    ●6.7  ●110.0'   88.0
 55. NGC6940      Sh2-102, LDN846.0                                   OCl        308.6    28.3    ●6.3  ●103.9'   87.5
 56. NGC6992      NGC6995, Sh2-103, LBN191  Eastern Veil, Network Neb SNR        314.1    31.7    ●7.0  ●210.0'   87.0
 57. NGC1647      LDN1558.0                                           OCl         71.5    19.1    ●6.4   ●50.7'   85.5
 58. NGC0598      PGC5818                   Triangulum Galaxy, Triang G           23.5    30.7    ●5.8   ●68.7'   85.0
 59. Sh2-134      NGC7261, Sh2-135, LBN474,                           HII        332.9    59.4   ○11.0  ●190.0'   85.0
 60. NGC2174      NGC2175, IC2159, Sh2-247, Monkey Head Nebula        Neb         92.3    20.7    ●6.8   ●40.0'   85.0
 61. IC0348       LBN749, LBN758, Ced18a, C omi Per Cloud             Cl+N        56.1    32.2    ○9.9  ●150.0'   85.0
 62. IC0434       NGC1990, NGC2023, NGC2024 Alnilam, Flame Nebula, Or HII         85.3    -2.5   ●11.0  ●146.7'   85.0
 63. Sh2-265      Sh2-263, LBN866, LBN867,                            HII         79.7     7.4   ○11.0   ●70.0'   85.0
 64. Sh2-113      Sh2-114, LBN296, LBN299,                            HII        320.2    38.1   ○11.0  ●230.0'   85.0
 65. IC1396       Sh2-131, LBN451, LBN452,  Elephant's Trunk Nebula   Cl+N       324.7    57.5    ○9.9  ●170.0'   85.0
 66. Sh2-171      NGC7762, NGC7822, LBN580, Question Mark Nebula      HII          1.2    67.2   ○11.0  ●180.0'   85.0
 67. IC2087       IC2088, LBN793, LBN812, L                           Neb         70.0    25.7    ●6.3  ●300.0'   85.0
 68. NGC2112      NGC2064, NGC2067, NGC2068                           OCl         88.4     0.4    ●9.1  ●140.0'   85.0
 69. NGC1333      LBN734, LBN740, LBN741, L Embryo Nebula             Cl+N        52.2    31.4   ●10.9   ●77.3'   85.0
 70. IC0447       IC0446, LBN887, LBN895, L Coyote Cloud, Dreyer's Ne HII         97.8     9.9    ●7.7  ●175.0'   85.0
 71. Sh2-249      IC0443, IC0444, Sh2-248,  Gem A, Jellyfish Nebula   HII         95.2    23.1   ○11.0   ●80.0'   85.0
 72. NGC2374      NGC2359, NGC2361, NGC2396 Thor's Helmet             OCl        111.0   -13.3    ●8.0   ●35.0'   85.0
 73. NGC6871      IC1310, Sh2-101, LBN162,  Tulip Nebula              OCl        301.5    35.8    ●5.2  ●270.0'   85.0
 74. NGC6604      NGC6625, NGC6631, Sh2-53, Cauda                     OCl        274.5   -12.2    ●6.5  ●259.5'   85.0
 75. NGC1624      Sh2-211, Sh2-212, LBN722,                           Cl+N        70.2    50.5   ●11.8   ●96.7'   82.5
 76. LBN857       Ced47, Ced51, B30, B31, B                           NB          80.4    11.4    ○6.7   ●70.0'   82.0
 77. Sh2-202      LBN677, LBN681, LBN682, L                           HII         49.7    59.6   ○11.0  ●240.0'   81.5
 78. Sh2-82       LBN128, LBN129, Ced168, L                           HII        292.6    18.3   ○11.0  ●261.5'   80.5
 79. NGC6819      LBN186, LBN196, LBN200, L Foxhead Cluster           OCl        295.3    40.2    ●7.3   ●70.0'   80.0
 80. NGC2183      LBN990, LBN999, Ced65, Ce                           HII         92.7    -6.2   ●15.2  ●106.3'   80.0
 81. IC1613       PGC3844                   Cetus Dwarf Galaxy        G           16.2     2.1    ●9.5   ●18.3'   80.0
 82. NGC2403      NGC2404, PGC21396                                   G          114.2    65.6    ●8.4   ●23.4'   80.0
 83. NGC4236      PGC39346                                            G          184.2    69.5    ●9.8   ●23.5'   80.0
 84. IC0426       IC0423, LBN918, Ced52, Ce Tear Drop Nebula          Neb         84.1    -0.3    ○9.0  ●300.0'   79.5
 85. LBN404       LBN407, Ced181, B155, B15                           NB         321.0    44.7    ○6.7  ●120.0'   78.5
 86. IC1470       Sh2-154, Sh2-156, LBN521,                           HII        346.3    60.2   ●11.5   ●75.7'   78.5
 87. NGC1909      LBN959, LBN968, Ced41a, C Witch Head Nebula         RfN         76.2    -7.3    ○9.0  ●180.0'   78.0
 88. NGC0772      LBN688, PGC7525                                     G           29.8    19.0   ●10.3   ●60.0'   77.5
 89. Sh2-123      LBN405, LBN414, Ced183b,                            HII        325.6    44.5   ○11.0  ●270.0'   77.5
 90. LBN526       LBN532, Ced194, Ced196, L                           NB         338.4    66.8    ○6.7   ●50.0'   77.0
 91. NGC1499      Sh2-220, Ced26, LDN1449.0 California Nebula         Neb         60.8    36.4    ●5.0  ●300.0'   76.5
 92. NGC6760      NGC6749, B138, LDN622.0,                            GCl        287.8     1.0    ●9.8  ●180.0'   76.0
 93. LBN420       Ced183a, Ced183c, Ced183d                           NB         331.3    42.7    ○6.7  ●190.0'   76.0
 94. NGC0744      LBN637, LBN638, LBN640                              OCl         29.6    55.5    ●7.9   ●35.0'   76.0
 95. Ced186       Ced200, Ced202, LDN1091.0                           NB         330.3    54.6   ●12.0  ●199.0'   75.0
 96. Sh2-118      LBN366, Ced176f                                     HII        324.3    40.2   ○11.0  ●300.0'   74.0
 97. NGC2126      LBN764, LBN765                                      OCl         90.6    49.9    ○9.1   ●50.0'   74.0
 98. NGC7063      LBN355, Ced180                                      OCl        321.1    36.5    ●7.0   ●90.0'   74.0
 99. NGC1662      NGC1663, LDN1561.0                                  OCl         72.1    10.9    ●6.4   ●31.4'   73.5
100. NGC6811      LBN225                                              OCl        294.3    46.4    ●6.8  ●210.0'   72.0
101. IC4756       LDN630.0, LDN633.0, LDN63 Graff’s Cluster           OCl        279.7     5.5    ●4.6   ●32.9'   71.5
102. NGC1746      LDN1542.0, LDN1544.0      Cluster of Clusters       OCl         76.0    23.8    ●6.1   ●47.5'   71.0
103. IC4665                                 Summer Beehive Cluster    OCl        266.6     5.6    ●4.2   ●24.6'   70.0
104. LBN788       LBN799, LBN800, LBN814, L                           NB          66.8    26.1    ○6.7  ●300.0'   70.0
105. LBN278       LBN298, LBN300, LBN301, L                           NB         304.7    43.2    ○6.7  ●300.0'   70.0
106. Sh2-155      Sh2-141, LBN509, LBN510,  Cave Nebula               HII        344.2    62.6   ○11.0  ●300.0'   70.0
107. Sh2-107      LBN230, LBN231, LBN235, L                           HII        310.7    36.3   ○11.0  ●180.0'   70.0
108. NGC0281      IC1590, Sh2-184, LBN599,  Pac Man Nebula            HII         13.2    56.6   ○11.0  ●180.0'   70.0
109. Sh2-115      LBN344, LBN352, LBN357, L                           HII        308.6    46.9   ○11.0  ●240.0'   70.0
110. NGC6991      IC5076, LBN361, LBN388, L                           OCl        313.7    47.5    ○9.1  ●240.0'   70.0
111. LBN114       LBN115, LBN119, LBN120, L                           NB         307.9    -2.5    ○6.7  ●140.0'   70.0
112. Sh2-69       Sh2-66, RCW176, RCW177, L                           HII        281.1    -0.3   ○11.0   ●77.8'   70.0
113. NGC6793      LBN130, LBN133, LBN134, C                           OCl        290.8    22.1    ○9.1  ●169.7'   70.0
114. Sh2-64       RCW174, LBN90, LBN95, LDN                           HII        277.9    -1.9   ○11.0  ●247.4'   70.0
115. Sh2-98       Sh2-99, LBN154, LBN155, L                           HII        299.7    31.4   ○11.0  ●300.0'   70.0
116. Sh2-92       Sh2-89, Sh2-90, LBN144, L                           HII        296.7    28.2   ○11.0   ●94.9'   70.0
117. Sh2-145      Sh2-140, Sh2-150, LBN500,                           HII        336.4    64.3   ○11.0   ●94.9'   70.0
118. IC5068       LBN265, LBN275, LBN277, L                           HII        312.6    42.5   ○11.0  ●238.3'   70.0
119. IC1287       NGC6639, Sh2-55, LBN73, L                           RfN        277.9   -10.8    ●6.1   ●60.0'   70.0
120. IC0360       LBN775, LBN777, Ced27, Ce                           Neb         62.3    26.1    ○9.0  ●180.0'   69.0
121. Sh2-126      LBN428, LBN429, LBN430, L                           HII        338.4    38.6   ○11.0  ●160.0'   69.0
122. Sh2-261      Sh2-254, LBN858, LBN862,  Lower's Nebula            HII         92.2    15.8   ○11.0  ●162.2'   69.0
123. LBN199       LBN207, LBN213, LBN217, L                           NB         300.2    39.5    ○6.7  ●160.0'   68.5
124. Sh2-170      Sh2-168, Sh2-169, LBN568,                           HII          0.4    64.6   ○11.0   ●30.0'   68.0
125. Sh2-112      LBN315, LBN317, LBN330, L                           HII        308.5    45.7   ○11.0  ●180.0'   67.5
126. Sh2-67       Sh2-65, RCW175, LBN91, LB                           HII        282.4    -2.4   ○11.0   ●91.0'   67.5
127. LBN198       LBN210, LBN211, LBN218, L                           NB         307.2    35.7    ○6.7  ●240.0'   67.0
128. Sh2-276      LBN942, LBN956, LBN957, L Barnard's Loop            HII         81.9    -4.0   ○11.0  ●300.0'   67.0
129. Sh2-280      Sh2-282, LBN970, LBN971,                            HII         98.6     2.5   ○11.0   ●40.0'   67.0
130. Sh2-221      NGC1798, Sh2-217, LBN745,                           HII         75.4    46.3   ○11.0  ●120.0'   67.0
131. Sh2-227      Sh2-225, Sh2-228, LBN778,                           HII         80.0    39.0   ○11.0   ●30.0'   67.0
132. Sh2-129      LBN445, LBN446, LBN449, L Flying Bat Nebula         HII        317.9    60.0   ○11.0  ●150.0'   67.0
133. LBN250       LBN256, LBN266, LBN269, L                           NB         313.7    36.2    ○6.7  ●130.0'   67.0
134. LBN10        LBN11, LBN15, LBN18, LBN1                           NB         238.4    -4.7    ○6.7  ●140.0'   66.5
135. NGC7640      LBN499, LBN503, LBN507, P                           G          350.5    40.8   ●11.0   ●60.0'   66.5
136. IC0359A      LBN782, LBN785, Ced30, Ce                           RfN         64.7    28.3    ○9.0   ●51.1'   66.5
137. Sh2-132      LBN471, LBN473, Ced192, C                           HII        334.7    56.1   ○11.0   ●90.0'   66.0
138. Sh2-278      LBN907, LBN915, LBN919, L                           HII         80.0    -5.7   ○11.0  ●200.0'   65.5
139. Sh2-68       LBN93, LBN96, LBN97, LBN9                           HII        276.3     0.9   ○11.0   ●68.9'   65.5
140. LBN792       B23, B24, B26, B27, B221,                           NB          72.0    29.8    ○6.7   ●61.5'   65.5
141. RCW4         LBN1026, LBN1028, LBN1029                           HII        111.0    -8.6   ○11.0  ●140.0'   65.0
142. Sh2-251      LBN842, LBN846, LBN847, L                           HII         68.2     5.9   ○11.0  ●270.0'   65.0
143. NGC0100      LBN570, LBN574, PGC1525                             G            6.0    16.5   ●13.2   ●80.0'   64.5
144. IC2067       NGC1579, Sh2-222, Ced35,                            Neb         67.7    35.4   ●11.5  ○126.4'   64.5
145. IC1276       NGC6539, LDN408.0, LDN418                           GCl        272.7    -7.2   ○10.4   ●44.3'   64.0
146. LBN422       B357, B359, LDN1033.0, LD                           NB         312.6    56.8    ○6.7   ●55.9'   64.0
147. IC0448       LBN929, LBN947, Ced79, Ce                           HII         98.2     7.4   ○11.0  ●300.0'   63.5
148. Sh2-232      Sh2-235, LBN789, LBN805,                            HII         85.6    36.2   ○11.0  ●128.7'   63.5
149. LBN878       LBN879, Ced59, B35, B36,                            NB          86.2     9.2    ○6.7  ●120.0'   63.5
150. Sh2-27       LBN21, LBN22, LBN32, LBN3                           HII        249.3   -10.6   ○11.0  ●300.0'   63.5
151. LBN283       LBN289, LBN295, LBN304, L                           NB         313.7    38.2    ○6.7  ●169.7'   63.5
152. Sh2-216      LBN742, LBN744, B15, B16,                           HII         71.2    46.8   ○11.0   ●80.0'   63.0
153. LBN883       LBN884, LBN888, LBN890, L                           NB          79.4     4.2    ○6.7  ●120.0'   63.0
154. LBN723       LBN727, LBN729, LBN732, L                           NB         121.8    61.4    ○6.7   ●55.0'   63.0
155. Sh2-183      Sh2-181, LBN614, LBN618,                            HII         13.5    65.7   ○11.0   ●98.2'   63.0
156. LBN559       LBN560, LBN561, LBN563, L                           NB         359.1    49.6    ○6.7  ●140.0'   63.0
157. LBN442       LBN448, LBN450, LBN458, L                           NB         339.6    41.1    ○6.7  ●290.0'   63.0
158. Sh2-111      LBN282, LBN287, LBN306, C                           HII        325.5    30.1   ○11.0   ●90.0'   63.0
159. Sh2-239      LBN817, LBN819, LBN821, L                           HII         67.8    18.1   ○11.0  ●150.0'   62.5
160. Sh2-210      LBN712, B12, B13, LDN1402                           HII         67.7    52.6   ○11.0   ●26.9'   62.5
161. LBN528       LBN535, LBN541, LBN546, L                           NB         331.3    70.9    ○6.7  ●155.0'   62.5
162. NGC1337      LBN870, PGC12916                                    G           52.0    -8.4   ●12.9   ●65.0'   62.5
163. NGC1788      LBN910, LBN923, Ced40, LD Cosmic Bat Nebula         RfN         76.7    -3.3    ●5.8  ●120.0'   62.0
164. LBN603       LBN608, LBN609, LBN612, L                           NB          10.7    52.3    ○6.7   ●85.0'   62.0
165. Sh2-75       Sh2-76, RCW181, LBN110, L                           HII        284.8     7.1   ○11.0  ●268.3'   61.5
166. Sh2-94       Sh2-96, LBN148, LBN150, L                           HII        292.0    31.5   ○11.0   ●50.0'   61.5
167. Sh2-200      Sh2-201, LBN674, LBN676,                            HII         47.7    62.8   ○11.0   ●20.0'   61.5
168. LBN545       LBN550, LBN552, LBN555, L                           NB         318.6    77.0    ○6.7  ●145.0'   61.5
169. LBN828       LBN831, LBN832, LBN834                              NB          63.4    10.2    ○6.7  ●140.0'   61.0
170. LBN464       LBN465, LBN469, LBN481                              NB         292.9    69.9    ○6.7   ●90.0'   61.0
171. NGC6847      Sh2-97, LBN151, Ced172                              Cl+N       299.2    30.2    ○9.9   ●35.0'   61.0
172. Sh2-223      Sh2-224, LBN768, LBN769                             HII         79.3    42.2   ○11.0   ●70.0'   61.0
173. LBN893       LBN896, LBN897, LBN900                              NB          68.1    -4.9    ○6.7   ●60.0'   61.0
174. Sh2-23       Sh2-24, LBN13, LBN17                                HII        243.4    -8.4   ○11.0   ●50.0'   61.0
175. LBN690       LBN693, LBN699, LBN700                              NB          83.8    66.5    ○6.7   ●90.0'   61.0
176. LBN165       LBN173, LBN176, LBN184                              NB         296.0    36.6    ○6.7  ●195.0'   61.0
177. LBN454       LBN460, LBN462, LBN463                              NB         341.8    41.8    ○6.7  ●180.0'   61.0
178. LBN1006      LBN1007, LBN1009, LBN1010                           NB          83.3   -12.3    ○6.7  ●200.0'   61.0
179. LBN982       LBN985, LBN988, LBN989                              NB          88.4    -5.5    ○6.7  ●180.0'   61.0
180. NGC7708      LBN573, LBN575, LDN1259.0                           OCl        353.8    72.8    ○9.1   ●40.0'   60.5
181. Sh2-218      LBN750, LBN751, LDN1460.0                           HII         85.3    52.2   ○11.0  ●160.0'   60.5
182. PGC13826                                                         G           56.7    68.1    ●9.1   ●20.9'   60.0
183. PGC54074                                                         G          227.2    67.2   ●12.8   ●30.1'   60.0
184. PGC60095                                                         G          260.1    57.9   ●11.7   ●35.8'   60.0
185. LBN292       LBN305, LBN309, B348                                NB         306.2    43.2    ○6.7   ●90.0'   60.0
186. Sh2-110      NGC7037, LBN254                                     HII        320.2    32.5   ○11.0  ●100.0'   60.0
187. LBN531       LBN538, Ced201, B175                                NB         333.3    70.3    ○6.7   ●60.0'   60.0
188. LBN818       Ced50, Ced56, B226                                  NB          84.8    30.7    ○6.7   ●35.0'   60.0
189. NGC6828      B338, B339, LDN668.0, LDN                           OCl        297.6     7.9    ○9.1   ●60.0'   60.0
190. NGC2143      LBN908, LDN1611.0, LDN161                           OCl         90.8     5.8    ○9.1  ●240.0'   59.5
191. Sh2-91       LBN147, LDN811.0, LDN812.                           HII        293.9    29.6   ○11.0  ●134.2'   59.5
192. Sh2-241      LBN824, LBN825, LDN1557.0                           HII         91.0    30.2   ○11.0   ●61.2'   59.5
193. LBN634       LBN635, LBN636, LDN1333.0                           NB          36.3    75.4    ○6.7   ●23.0'   59.5
194. IC1369       B361, LDN959.0, LDN963.0,                           OCl        318.0    47.8    ○9.1   ●39.3'   59.5
195. LBN628       LBN631, LBN632, LDN1320.0                           NB          66.6    86.0    ○6.7  ●105.0'   59.5
196. Sh2-205      LBN696, LBN701, LDN1391.0                           HII         59.0    53.2   ○11.0  ●120.0'   59.5
197. Sh2-245      LBN835, LBN836            Fishhook Nebula           HII         60.6     4.1   ○11.0  ●300.0'   59.0
198. Sh2-264      LBN865, Ced54             Angelfish Nebula          HII         83.8     9.9   ○11.0  ●270.0'   59.0
199. IC0341       Ced19a, Ced19b                                      Neb         55.2    22.0    ○9.0  ●134.9'   59.0
200. LBN743       LBN753, LDN1453.0, LDN145                           NB          44.0    20.2    ○6.7  ●140.0'   59.0
=======================================================================================================================
                                                                           ... i 14 więcej kadrów (limit wyświetlania).
=======================================================================================================================

[INFO] Zapisano 214 zoptymalizowanych kadrów do data/vis_data.json
```
**% python3 3_compute.py**
```
[INFO] Brak poprzedniego stanu silnika (PATHS.engine_state).
[INFO] W vis_data.json jest 214 obiektów.
       Ile obiektów przeliczyć? [domyślnie 108, Enter = wszystkie]: 
[INFO] Do przeliczenia: 214 obiektów.
=======================================================================================================================
ENGINE: Smart Cache System
=======================================================================================================================
[INFO] Lokalizacja: Poznań, Polska (52.41°, 16.93°)
[INFO] Rok: 2026
[INFO] Parametry: minimalna wysokość obiektu: 22.5°, wysokość słońca: -12.0°
[INFO] Cache RAW:   data/observing_data_raw.pkl
[INFO] Cache FINAL: data/observing_data.pkl
=======================================================================================================================
[CACHE] Brak hash lub RAW cache - obliczam od nowa.
[RAW]   Brak ważnego RAW cache - przygotowuję przeliczenie wszystkich obiektów.
[RAW]   Brakujących obiektów w RAW cache: 214.
[COMPUTE] Obliczam dane Słońca dla 365 dni.
          Sun Calc: 100%|████████████████████████████████████████████████████████████| 365/365 [00:43<00:00,  8.31it/s]
[COMPUTE] Obliczam RAW dla 214 brakujących obiektów...
          Raw compute: 100%|████████████████████████████████████████| 214/214 [34:30<00:00,  9.68s/obj, ostatni=LBN126]
[SAVE]  Raw data zapisane do data/observing_data_raw.pkl (214 obiektów).
[SAVE]  RAW hash zapisany: ac206b38...
[COMPUTE] Parametry/lokalizacja/rok się zmieniły lub brak FINAL cache – pełne przeliczenie FINAL.
          Przeliczanie FINAL (Parallel): konwersja 214 obiektów.
          Reprocess: 100%|██████████████████████████████████████████████████████████| 214/214 [00:05<00:00, 35.71obj/s]
[COMPUTE] Przetworzono 214 obiektów.
[SAVE] Final data zapisane do: data/observing_data.pkl
[SAVE] Hash zapisany do: data/observing_data_final.hash
[INFO] Raw objects: 214. Final objects: 214.
=======================================================================================================================
[INFO] Silnik zakończył pracę.
=======================================================================================================================
```
**% python3 4_select_objects.py**
```
=======================================================================================================================
       OPTYMALIZACJA MIESIĘCZNA I GENEROWANIE PLANERA
=======================================================================================================================
[INFO] Rok obserwacji: 2026
[INFO] Lokalizacja: Poznań, Polska (52.41°N, 16.93°E)
[INFO] Strefa czasowa: Europe/Warsaw
[INFO] Minimalna wysokość nad horyzontem: 22.5°
[INFO] Wymagane okno widoczności: 2.5h
[INFO] Limit ciemności nieba (Słońce): -12.0°
=======================================================================================================================
[INFO] Miesięczna macierz najlepszych nocy: 2568 rekordów.
[INFO] Macierz widoczności rocznej: 214 zliczonych obiektów.
[INFO] Całkowita pula obiektów kandydujących z JSON: 214

[USER] Podaj liczbę obiektów w pojedynczym wariancie na miesiąc [Enter = 3]: 

[INFO] Uruchamiam algorytm optymalizacji globalnej (Hungarian Algorithm)...
[INFO] Obiektów do rozplanowania: 213
[INFO] Dostępnych slotów: 108 (3 warianty x 12 miesięcy x 3)
[INFO] Mediana Score: 70.00 (Granica logiki Elita vs Reszta)

=======================================================================================================================
 RAPORT KOŃCOWY PO PRZYDZIALE DO WARIANTÓW (OPTIMAL + PRESTIGE LOGIC)
=======================================================================================================================
[INFO] Top (score > mediana) – 102 obiektów:
       Przypisanych: 101/102 (99.0%)
       Odrzuconych:  1/102 (1.0%)
[INFO] Odrzucone obiekty z Top (score > mediana):
       • NGC2374  Score:  85.0 | Najlepsze: 01 (2.4h), 02 (2.4h), 03 (2.4h)
[INFO] Wszystkie sloty zostały wypełnione!

[INFO] Obiekty użyte z puli poniżej mediany score (Sortowane wg Prestiżu):
       Łącznie: 7 obiektów z tej puli.
        1. IC4665       (Waga:  90) Score: 70.0 | Najlepsze: 05 (4.2h), 07 (3.9h), 04 (3.8h)
        2. NGC7640      (Waga:  90) Score: 66.5 | Najlepsze: 10 (10.7h), 11 (10.0h), 09 (9.9h)
        3. IC1276       (Waga:  90) Score: 64.0 | Najlepsze: 07 (3.2h), 06 (3.0h), 05 (2.7h)
        4. NGC6847      (Waga:  90) Score: 61.0 | Najlepsze: 08 (7.2h), 09 (7.2h), 10 (6.5h)
        5. NGC6828      (Waga:  90) Score: 60.0 | Najlepsze: 08 (5.8h), 09 (5.2h), 07 (5.1h)
        6. IC1369       (Waga:  90) Score: 59.5 | Najlepsze: 09 (9.8h), 10 (9.6h), 11 (8.7h)
        7. Sh2-27       (Waga:  70) Score: 63.5 | Najlepsze: 05 (3.6h), 06 (3.5h), 04 (2.8h)

[INFO] Łącznie przypisanych obiektów: 108/214 (50.5%)
=======================================================================================================================
 PODSUMOWANIE JAKOŚCI PLANU (213 obiektów / 108 slotów)
=======================================================================================================================
[INFO] Średnia jakość okna obserwacyjnego: 86.5% (BARDZO DOBRA).
       Średnia jakość okna obserwacyjnego względem najlepszej możliwej w roku: 87%.

[INFO] Obciążenie kalendarza (Średni Score obiektów w miesiącu):
       Miesiąc    Śr. Score  Liczba   Status
       ---------------------------------------------
       Sty          94.4      9 szt.   🔥 ELITA (Top Obiekty)
       Lut          91.1      9 szt.   🔥 ELITA (Top Obiekty)
       Mar          95.9      9 szt.   🔥 ELITA (Top Obiekty)
       Kwi          82.0      9 szt.   🔥 ELITA (Top Obiekty)
       Maj          77.9      9 szt.   ✨ DOBRE (Solidne)
       Cze          69.9      9 szt.   ✨ DOBRE (Solidne)
       Lip          87.7      9 szt.   🔥 ELITA (Top Obiekty)
       Sie          92.8      9 szt.   🔥 ELITA (Top Obiekty)
       Wrz          98.7      9 szt.   🔥 ELITA (Top Obiekty)
       Paź          88.3      9 szt.   🔥 ELITA (Top Obiekty)
       Lis          87.9      9 szt.   🔥 ELITA (Top Obiekty)
       Gru         100.2      9 szt.   🔥 ELITA (Top Obiekty)

[WARN] Kompromisy (Top Obiekty przesunięte do gorszych miesięcy):
       • NGC0188  (Score 92): Miesiąc 04 (8.3h) -> Zamiast 12 (13.4h). Strata: -38%
       • NGC6939  (Score 90): Miesiąc 04 (8.3h) -> Zamiast 12 (13.4h). Strata: -38%
       • Sh2-134  (Score 85): Miesiąc 05 (6.0h) -> Zamiast 11 (12.3h). Strata: -51%
       • IC1396   (Score 85): Miesiąc 05 (6.0h) -> Zamiast 10 (11.4h). Strata: -47%
       • Sh2-171  (Score 85): Miesiąc 04 (8.3h) -> Zamiast 12 (13.4h). Strata: -38%
       • Sh2-202  (Score 82): Miesiąc 04 (6.9h) -> Zamiast 12 (13.4h). Strata: -48%
       • NGC6819  (Score 80): Miesiąc 05 (5.2h) -> Zamiast 09 (7.8h). Strata: -34%
       • NGC2403  (Score 80): Miesiąc 04 (8.3h) -> Zamiast 12 (13.4h). Strata: -38%
       • NGC4236  (Score 80): Miesiąc 04 (8.3h) -> Zamiast 12 (13.4h). Strata: -38%
=======================================================================================================================


[INFO] Zakończono tworzenie wariantów (max. 9 DSO w każdym miesiącu w sumie dla wszystkich wariantów 108 obiektów).
[INFO] Zapisano flagi 'selected' do 108 obiektów w data/vis_data.json
[INFO] Statystyki: 108/214 obiektów wybranych (warianty A B C)

[INFO] Rozkład obiektów w utworzonych wariantach:
       Wariant A: łącznie przypisano 36 obiekty.
       Wariant B: łącznie przypisano 36 obiekty.
       Wariant C: łącznie przypisano 36 obiekty.
Generowanie PDF.
=======================================================================================================================

[INFO] Planer roczny zapisany jako: Astrophotography_Planner_2026_1.pdf
```
**% python3 5_generate_fov_and_ctx.py**
```
=======================================================================================================================
GENEROWANIE GRAFIK: FOV & MAPY KONTEKSTOWE
=======================================================================================================================

[INFO] Pobieram katalogi dla StarPlot.

[INFO] Konfiguracja kamery:
       ogniskowa:           300.0 mm,
       rozmiar sensora:     23.5 x 15.7 mm,
       rozdzielczość:       6248 x 4176 px,
       wielkość piksela:    3.76 µm,
       pole widzenia:       4.49° x 3.00° (przekątna 5.39°).

[INFO] Wybrano 108 obiektów (wszystkie z 'selected' w vis_data.json).

[INFO] Generowanie kadrów FOV: 108 do zrobienia, 0 pominięto (już istnieją). 
       Engine FOV: 100%|████████████████████████████████████████████| 108/108 [04:02<00:00,  2.24s/obj, ostatni=IC1369]
[INFO] Mamy łącznie 108 kadrów FOV.
[INFO] Generowanie map kontekstowych: 108 do zrobienia, 0 pominięto (już istnieją).
       Engine CTX: 100%|███████████████████████████████████████████| 108/108 [10:09<00:00,  5.64s/obj, ostatni=NGC6847]
[INFO] Mamy łącznie 108 map kontekstowych.

[INFO] Gotowe. Wygenerowane pliki znajdują się w folderze /Users/wdrodze/Public/astrophoto_planner/data/starplots.
```
**% python3 6_generate_objects_pages.py**
```
[INFO] Inicjalizacja danych.
       Generowanie CMAP: Sh2-82: 100%|███████████████████████████████████████████| 108/108 [10:59<00:00,  6.11s/obiekt]
[INFO] Proces zakończony. Plik wynikowy: Astrophotography_Planner_2026_2.pdf.
```
**% python3 7_generate_result_pdf.py**
```
[INFO] Pliki połączono z tytułem i stroną informacyjną w Astrophotography_Planner_2026_Poznań.pdf.
```
