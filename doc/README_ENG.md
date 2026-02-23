
# Astrophotography Planner

**Astrophotography Planner** is a comprehensive Python tool designed to generate a personalized, year-long imaging schedule and atlas. It tailors targets specifically to your geographical location, time, camera equipment, and local sky conditions (Bortle scale).

![Sample page](doc/strona_a)

The tool automates the entire process: from aggregating astronomical catalogs and calculating complex visibility windows (considering moon phases and twilight), to generating field-of-view (FOV) simulations and compiling a print-ready PDF book. [Here is an example, compressed however.](doc/AP_compressed.pdf)

![Month view](doc/planer_miesieczny.png)
---

## 📂 Project Structure & Workflow

The project is organized into sequential steps. You should run the scripts in the order numbered 1 through 7.

### ⚙️ Configuration (`shared.py`)
**Central Configuration Hub**
This file acts as the repository for all settings.

*   **UserConfig:** Stores user preferences (Location, Year).
*   **CameraConfig:** Defines equipment parameters (Focal length, Sensor size).
*   **PATHS:** Manages input/output directories.
*   **Constants:** Astronomical constants, catalog priorities, and color definitions.
*   *Note:* Changes here affect the entire pipeline.

---

### 1️⃣ Data Aggregation (`1_build_catalog.py`)
**Build and Aggregate Astronomical Catalog**
Prepares the master raw database (`katalog_astro_full.csv`).

*   **Sources:** Merges the base list (OpenNGC) with supplementary data from VizieR (Sharpless, Barnard, LDN, LBN, RCW, etc.).
*   **Normalization:** Unifies units, coordinates, and naming conventions.
*   **Entity Resolution:** Inteligent deduplication (e.g., recognizing that NGC 1499 and Sh2-220 are the same object).
*   **Filtering:** Interactively filters out objects that are too small or faint based on user preference.
*   *Requires:* Active internet connection.

### 2️⃣ Scoring & Setup (`2_plan_and_score.py`)
**Session Configuration and Scoring System**
Functions as both the configuration interface and the object assessment engine.

*   **Interactive Setup:** Prompts for Year, Location (Lat/Lon), Equipment (FOV calculation), and constraints (Min Altitude, Light Pollution/Bortle).
*   **Scoring Algorithm:** Assigns points to every object based on:
    *   Object Type vs. Filters (e.g., boosting emission nebulae for narrowband).
    *   Surface brightness and size.
    *   Inclusion in "Best of" lists (Messier, Caldwell, Herschel 400).
*   *Output:* `vis_data.json` containing scored objects.

### 3️⃣ Computations (`3_compute.py`)
**Compute Engine (The Heavy Lifter)**
Performs precise astronomical calculations using **AstroPy**. This is the most time-consuming step.

*   **Visibility Analysis:** Calculates altitude for every object, for every night of the year.
*   **Window Calculation:** Accounts for:
    *   Twilight (Astronomical/Nautical).
    *   Moon Phase & Position (excludes times when the Moon interferes).
    *   Hardware limits (horizon obstructions).
*   **Imaging Hours:** Sums up total quality imaging hours per year.

*   *Performance:* Uses Multiprocessing and caches results (`observing_data.pkl`) to speed up subsequent runs.

### 4️⃣ Scheduling (`4_select_objects.py`)
**Object Selection and Schedule Optimization**
Solves the resource allocation problem (assigning targets to specific nights/months).

*   **Distribution:** Splits objects into priority groups based on Step 2 scores.
*   **Allocation Algorithm:**
    *   Assigns objects to months with the best visibility.
    *   Balances the load (Slots A, B, C) to ensure a diverse plan.
    *   Prioritizes rare targets (short visibility windows) over circumpolar ones.

*Output:* Generates a report with planning statistics and a pdf file with yearly schedulle.

### 5️⃣ Visualization (`5_generate_fov_and_ctx.py`)
**Sky Map Generator**
Creates visual aids for all selected objects using the `starplot` library.

*   **FOV (Field of View):** Simulates your specific camera frame, showing how the target fits and marking other deep-sky objects in the frame.
*   **Context Maps:** Wide-field charts to assist with star hopping and general orientation.

### 6️⃣ Page Rendering (`6_generate_objects_pages.py`)
**Object Page Renderer**
Compiles detailed PDF datasheets for every planned target.

*   **Content:**
    *   Names, Type, Magnitude, Size and other objects in the frame.
    *   **Graphs:** Annual visibility, Altitude during the best night in selected month, Moon-free hours.
    *   **Maps:** FOV and Context images generated in Step 5.
    *   **Metadata:** Best night of the year, window duration, framing notes.

### 7️⃣ Final Assembly (`7_generate_result_pdf.py`)
**PDF Merger**
The final step that unifies all components into a professional document.

*   **Components:**
    *   Title pages and Legend.
    *   Monthly/Annual Calendars.
    *   Alphabetically ordered collection of Object Pages (from Step 6).
*   *Result:* `Astrophotography_Planner_[YEAR]_[CITY].pdf` ready for printing.

---

## 🚀 How to Run

1.  Install dependencies (see `requirements.txt`).
2.  Run the scripts sequentially:

```bash
python 1_build_catalog.py
python 2_plan_and_score.py
python 3_compute.py
python 4_select_objects.py
python 5_generate_fov_and_ctx.py
python 6_generate_objects_pages.py
python 7_generate_result_pdf.py
```

## 📝 Requirements

*   Python 3.x
*   AstroPy
*   Pandas
*   Starplot
*   Matplotlib
*   (See `requirements.txt` for full list)
