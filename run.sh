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
  "6_drukuj_strony_obiektów.py"
  "7_połącz_pliki_pdf.py"
)

echo "🚀 Setup & run yearly_plan w $PROJ_DIR"

# 1. Tworzenie venv (jeśli brak)
if [ ! -d "$VENV_NAME" ]; then
  echo "📦 Tworzę środowisko $VENV_NAME..."
  python3 -m venv "$VENV_NAME" || { echo "❌ Nie udało się utworzyć venv"; exit 1; }
fi

# 2. Aktywacja venv
source "$VENV_NAME/bin/activate"

# 3. Instalacja pakietów
if [ -f "$REQUIREMENTS" ]; then
  echo "📥 Instaluję pakiety z $REQUIREMENTS..."
  pip install --upgrade pip
  pip install -r "$REQUIREMENTS"
else
  echo "⚠️ Brak $REQUIREMENTS – pomijam instalację pakietów."
fi

# 4. Uruchamianie skryptów po kolei
for py_file in "${PY_FILES[@]}"; do
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
done

echo "✅ Wszystkie skrypty zakończone pomyślnie."
deactivate