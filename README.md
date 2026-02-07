# Astrophotography Planner & Atlas Generator

Zestaw zaawansowanych skryptów w języku Python służący do generowania spersonalizowanego rocznego planera astronomicznego oraz atlasu obiektów głębokiego nieba (DSO).

System pobiera dane z katalogów astronomicznych, filtruje je pod kątem lokalizacji obserwatora i posiadanego sprzętu (teleskop/kamera), oblicza precyzyjną widoczność na dany rok, a następnie generuje profesjonalny plik PDF zawierający:
1. Przegląd roczny (kiedy obserwować dany obiekt).
2. Szczegółowe strony dla każdego obiektu (wykresy wysokości, kadry FOV, mapy kontekstowe).

## 🚀 Możliwości

*   **Agregacja danych:** Łączy katalogi NGC/IC, Sharpless (Sh2), RCW, Barnard, LBN, LDN, Cederblad i PGC.
*   **Inteligentne filtrowanie:** Wybiera obiekty na podstawie Twojej szerokości geograficznej, minimalnej wysokości nad horyzontem, jasności (Mag), rozmiaru oraz zanieczyszczenia nieba (skala Bortle).
*   **Symulacja FOV:** Generuje symulacje kadru (Field of View) dla Twojej kamery i teleskopu przy użyciu biblioteki `starplot`.
*   **Obliczenia astronomiczne:** Precyzyjne wyliczanie okien obserwacyjnych (godziny bez Księżyca, wysokość górowania).
*   **Format PDF:** Generuje gotowy do druku atlas w formacie A4.

## 🛠️ Wymagania

Projekt wymaga Pythona 3.10+ oraz następujących bibliotek:

```bash
pip install pandas numpy astropy astroplan matplotlib reportlab pypdf tqdm astroquery starplot networkx