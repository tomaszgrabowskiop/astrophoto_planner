#!/bin/bash

# Katalog projektu = katalog skryptu
PROJ_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
cd "$PROJ_DIR" || exit 1

VENV_NAME="dso_env"
REQUIREMENTS="requirements.txt"

# Kolejność uruchamiania
PY_FILES=(
  "0_opracuj_katalog_ngc.py"
  "1_generuj_katalog_astro.py"
  "2_ograniczenie_katalogu.py"
  "3_wyliczenia.py"
  "4_plan_roczny.py"
  "5_fov_and_maps.py"
  "6_drukuj_strony_obiektow.py"
  "7_polacz_pliki_pdf.py"
)

STEP_DESC=(
  "Krok 0: Przygotowanie katalogu NGC — parsuje OpenNGC/NGC.csv. Tworzy: updated_ngc.csv."
  "Krok 1: Pobieranie i unifikacja katalogów — VizieR + smart merge. Tworzy: katalog_astro_full.csv."
  "Krok 2: Konfiguracja i selekcja obiektów — interaktywny filtr lokalizacja/limity alt, mag, size... Tworzy: vis_data.json."
  "Krok 3: Silnik obliczeniowy AstroPy — widoczność obiektów przez cały rok. Tworzy: observing_data.pkl. [UWAGA! Zajmuje sporo czasu]"
  "Krok 4: Plan roczny — przydział obiektów do miesięcy. Tworzy: Astrophotography_Planner_2026_1.pdf."
  "Krok 5: StarPlot: generowanie widoków FOV i map kontekstowych w PNG."
  "Krok 6: Generowanie stron obiektów. Tworzy: Astrophotography_Planner_2026_2.pdf."
  "Krok 7: Finalizacja — Scalenie PDF do Astrophotography_Planner_2026.pdf."
)

echo "🚀 Setup & run yearly_plan w $PROJ_DIR"

# 1. Tworzenie venv (jeśli brak)
if [ ! -d "$VENV_NAME" ]; then
  echo "📦 Tworzę środowisko $VENV_NAME..."
  python3 -m venv "$VENV_NAME" || { echo "❌ Nie udało się utworzyć venv"; exit 1; }
fi

# 2. Aktywacja venv
# shellcheck source=/dev/null
source "$VENV_NAME/bin/activate"

# 3. Instalacja pakietów
if [ -f "$REQUIREMENTS" ]; then
  echo "📥 Instaluję pakiety z $REQUIREMENTS..."
  pip install --upgrade pip
  pip install -r "$REQUIREMENTS"
else
  echo "⚠️ Brak $REQUIREMENTS – pomijam instalację pakietów."
fi

echo
echo "Dostępne kroki:"
for i in "${!PY_FILES[@]}"; do
  printf "  %d) %s - %s\n" "$i" "${PY_FILES[$i]}" "${STEP_DESC[$i]}"
done
echo

# Wybór kroku startowego
while true; do
  read -p "Od którego kroku mam zacząć? [0-7]: " START
  if [[ "$START" =~ ^[0-7]$ ]]; then
    break
  else
    echo "Podaj liczbę z zakresu 0–7."
  fi
done

CURRENT=$START
TOTAL=${#PY_FILES[@]}

while [ $CURRENT -lt $TOTAL ]; do
  py_file="${PY_FILES[$CURRENT]}"
  desc="${STEP_DESC[$CURRENT]}"

  echo
  echo "----------------------------------------"
  echo "$desc"
  echo "Plik: $py_file"
  echo "Krok: $CURRENT z $((TOTAL-1))"
  echo "----------------------------------------"

  if [ -f "$py_file" ]; then
    echo "⚡ Uruchamiam $py_file..."
    python "$py_file"
    STATUS=$?
    if [ $STATUS -ne 0 ]; then
      echo "❌ Błąd w $py_file (kod $STATUS) – zatrzymuję pipeline."
      deactivate
      exit $STATUS
    fi
  else
    echo "⚠️ Plik nie istnieje: $py_file"
  fi

    # Po ostatnim kroku po prostu wychodzimy
  if [ $CURRENT -ge $((TOTAL-1)) ]; then
    echo "✅ Osiągnięto ostatni krok ($CURRENT)."
    break
  fi

  # Co zrobić dalej:
  while true; do
    read -p "Co dalej? [n]astępny, [p]owtórz, [0-7 wybierz krok], [z]akończ: " next
    case "$next" in
      n|N)
        CURRENT=$((CURRENT+1))
        break
        ;;
      p|P)
        # powtórz ten sam krok
        # (CURRENT pozostaje bez zmian)
        break
        ;;
      z|Z)
        echo "🛑 Przerywam działanie sekwencji."
        deactivate
        exit 0
        ;;
      [0-7])
        # przejście do wybranego numeru kroku
        if [ "$next" -ge 0 ] && [ "$next" -lt "$TOTAL" ]; then
          CURRENT=$next
          break
        else
          echo "Dozwolone numery kroków: 0–$((TOTAL-1))."
        fi
        ;;
      *)
        echo "Wpisz: n (następny), p (powtórz), 0–7 (wybierz krok), z (zakończ)."
        ;;
    esac
  done

done

deactivate
echo "🎉 Zakończono działanie sekwencji kroków."
