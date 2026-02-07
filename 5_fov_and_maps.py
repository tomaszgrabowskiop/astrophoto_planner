#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
Etap po zrobieniu wykresów dla wybranych nocy, wtedy tylko te obiekty są rysowane. 
1) wczytuje vis_data.json,
2) generuje FOV oraz context map PNG dla obiektów, które mają uzupełnione pole 'selected'
   – zgodnie z silnikami OpticMapEngine/ContextMapEngine
"""
import multiprocessing
import multiprocessing as mp
from multiprocessing import Pool

from tqdm import tqdm
from concurrent.futures import ProcessPoolExecutor, as_completed

import json
import math
import gc
from pathlib import Path
from typing import Dict, Any, List

import pandas as pd

from starplot import OpticPlot, MapPlot, DSO, LambertAzEqArea, Mercator, Miller, DsoType, _
from starplot.styles import PlotStyle, extensions, LineStyle
from starplot.models import Camera
from starplot.callables import color_by_bv
from starplot.data.catalogs import OPEN_NGC
import matplotlib.pyplot as plt

# --- KONFIGURACJA PLIKÓW ---

VIS_DATA_PATH = "vis_data.json"  
STARPLOTS_DIR = Path("starplots")
STARPLOTS_DIR.mkdir(exist_ok=True)


# --- Wczytanie vis_data.json i budowa listy obiektów ---
def preload_open_ngc():
    """
    Mały preload: ściąga katalogi:
    - gwiazdozbiory
    - OPEN_NGC (DSO)
    - gwiazdy bigsky do ~11 mag
    """
    style = PlotStyle().extend(extensions.MAP)

    # Prosta projekcja Mercatora wokół RA=0, Dec=0
    proj = Mercator(center_ra=0.0, center_dec=0.0)

    plot = MapPlot(
        projection=proj,
        ra_min=0, ra_max=10,
        dec_min=-5, dec_max=5,
        style=style,
        resolution=50,  # mało, żeby było szybko
    )
    plot.dsos(catalog=OPEN_NGC)
    plot.constellations()               
    plot.stars(where=[_.magnitude < 11])

    
def _worker_render_fov(args):
    """
    Worker rozpakowuje argumenty, tworzy obiekty i renderuje.
    args: (obj, cam_params, out_dir)
    """
    obj, cam_params, out_dir = args
    
    try:
        camera = create_camera_object(cam_params)
        style = make_optic_style()
        engine = OpticMapEngine(camera, style, out_dir)
        return engine.render_field(obj)
    except Exception as e:
        return None


def _worker_render_context(args):
    """
    args: (obj, out_dir)
    """
    obj, out_dir = args
    
    try:
        style = make_context_style()
        engine = ContextMapEngine(style, out_dir)
        return engine.render_context(obj)
    except Exception as e:
        return None

def load_vis_data(path: str = VIS_DATA_PATH) -> Dict[str, Any]:
    p = Path(path)
    if not p.exists():
        raise FileNotFoundError(f"Brak pliku {path}")
    with p.open("r", encoding="utf-8") as f:
        return json.load(f)

def build_objects_from_vis_data(
    data: Dict[str, Any],
) -> pd.DataFrame:
    """
    Buduje DataFrame z sekcji 'objects' w vis_data.json,
    filtrując tylko obiekty, które mają uzupełnione pole 'selected'.
    """
    objs = []
    for obj in data.get("objects", []):
        if obj.get("selected"):
            objs.append(obj)

    df = pd.DataFrame(objs)
    return df

def _display_name(row: pd.Series) -> str:
    cn = row.get("commonname")
    if pd.notna(cn) and str(cn).strip():
        return str(cn)
    nm = row.get("name")
    if pd.notna(nm) and str(nm).strip():
        return str(nm)
    return str(row.get("id"))


# --- KROK 2: KONFIGURACJA KAMERY ---
def get_camera_params(path: str = VIS_DATA_PATH) -> dict:
    """
    Zwraca surowe parametry kamery (słownik), bezpieczne do przesłania do workera.
    Wyświetla informacje o wczytanej konfiguracji.
    """
    p = Path(path)
    if not p.exists():
        raise FileNotFoundError(f"Brak pliku {path}")

    with p.open("r", encoding="utf-8") as f:
        data = json.load(f)
    
    camcfg = data.get("parameters", {}).get("camera", {})

    print(
        "Konfiguracja kamery z vis_data.json: "
        f"f={camcfg['lens_focal_length']} mm, "
        f"sensor {camcfg['sensor_width']} x {camcfg['sensor_height']} mm, "
        f"{camcfg['sensor_cols']} x {camcfg['sensor_rows']} px, "
        f"pitch {camcfg['sensor_pitch']} µm"
    )

    return camcfg

def create_camera_object(camcfg: dict) -> Camera:
    """
    Tworzy obiekt Camera na podstawie słownika.
    Ta funkcja będzie wywoływana WEWNĄTRZ procesu workera.
    """
    return Camera(
        lens_focal_length=camcfg["lens_focal_length"],
        sensor_pitch=camcfg["sensor_pitch"],
        sensor_rows=camcfg["sensor_rows"],
        sensor_cols=camcfg["sensor_cols"],
        sensor_width=camcfg["sensor_width"],
        sensor_height=camcfg["sensor_height"],
    )

def ensure_starplots_dir() -> None:
    STARPLOTS_DIR.mkdir(exist_ok=True)


# --- STYL I HELPERY ---
def make_optic_style() -> PlotStyle:
    """
    Styl dla kadrów optycznych – przejęty z bak_3_drukuj_atlas.py (make_optic_style).
    """
    style = PlotStyle().extend(
        extensions.BLUE_NIGHT,
        extensions.OPTIC,
    )

    style.dso_galaxy.label.font_size = 28
    style.dso_galaxy.label.font_color = "33691E"
    style.dso_galaxy.label.font_weight = 900

    style.dso_nebula.label.font_size = 28
    style.dso_nebula.label.font_color = "FF6F00"
    style.dso_nebula.label.font_weight = 900

    style.dso_open_cluster.label.font_size = 28

    return style

def make_context_style() -> PlotStyle:
    """
    Styl dla MapPlot – zgodnie z przykładem map Cassiopei z bak_3_drukuj_atlas.py.
    """
    style = PlotStyle().extend(
        extensions.BLUE_MEDIUM,
        extensions.MAP,
    )
    return style

def generate_circle_coords_deg(ra_center_deg, dec_center_deg, size_arcmin, steps=100):
    """
    Generuje okrąg w stopniach z korekcją cos(dec) dla RA.
    Wzięte z bak_3_drukuj_atlas.py.
    """
    if size_arcmin <= 0:
        return [], []

    radius_deg = (size_arcmin / 2.0) / 60.0
    ra_points = []
    dec_points = []

    dec_rad = math.radians(dec_center_deg)
    cos_dec = math.cos(dec_rad)
    if abs(cos_dec) < 0.001:
        cos_dec = 0.001

    for i in range(steps + 1):
        angle = 2 * math.pi * i / steps
        d_dec = radius_deg * math.cos(angle)
        d_ra = (radius_deg * math.sin(angle)) / cos_dec

        new_dec = dec_center_deg + d_dec
        new_dec = max(-90.0, min(90.0, new_dec))

        ra_points.append(ra_center_deg + d_ra)
        dec_points.append(new_dec)

    return ra_points, dec_points


# --- SILNIK OPTICPLOT (FOV) ---
class OpticMapEngine:
    """
    Silnik do generowania kadrów kamery (OpticPlot).
    Brak funkcji pomocniczych typu obrysy/krzyżyki – czysty kadr optyczny.

    Obiekty wejściowe: dict z polami name, ra, dec (w stopniach).
    """

    def __init__(self, camera: Camera, style: PlotStyle, out_dir: Path):
        self.camera = camera
        self.style = style
        self.out_dir = out_dir

    @staticmethod
    def _dso_label(dso: DSO) -> str:
        """
        Preferuje common_names, potem NGC/IC, potem name.
        """
        if dso.common_names:
            return dso.common_names[0]
        elif dso.ngc:
            return f"NGC{dso.ngc}"
        elif dso.ic:
            return f"IC{dso.ic}"
        else:
            return dso.name

    def render_field(self, obj: Dict[str, Any]) -> Path | None:
        """
        - używa OpticPlot,
        - FOV wynika z kamery i ogniskowej,
        - brak niestandardowych obrysów (są w mapie kontekstowej).
        """
        name = obj.get("id") or obj.get("name")
        name = str(name)
        ra_center_deg = float(obj["ra"])  # RA w stopniach
        dec_center = float(obj["dec"])    # Dec w stopniach

        try:
            p = OpticPlot(
                ra=ra_center_deg,
                dec=dec_center,
                optic=self.camera,
                style=self.style,
                resolution=4600,
                autoscale=True,
                raise_on_below_horizon=False,
            )

            p.stars(
                where=[_.magnitude < 14],
                color_fn=color_by_bv,
                bayer_labels=True,
            )

            p.dsos(
                where_labels=[True],
                label_fn=self._dso_label,
                catalog=OPEN_NGC,
            )

            p.open_clusters(
                label_fn=self._dso_label,
                where_true_size=[True],
                catalog=OPEN_NGC,
            )

            p.nebula(
                label_fn=self._dso_label,
                where_true_size=[True],
                catalog=OPEN_NGC,
            )

            p.galaxies(
                label_fn=self._dso_label,
                where_true_size=[True],
                catalog=OPEN_NGC,
            )

            p.legend(
                style__alignment="left",
                style__location="upper left",
                style__padding_x=50,
                style__padding_y=50,
                style__num_columns=1,
            )

            out_png = self.out_dir / f"{name}.png"
            p.export(str(out_png), padding=0, transparent=True)
            return out_png

        except Exception as e:
            print(f"Błąd przy obiekcie {name}: {e}")
            return None

        finally:
            plt.close("all")
            if p is not None:
                del p
            gc.collect()


# --- SILNIK MAPPLOT (MAPA KONTEKSTOWA) ---
class ContextMapEngine:
    """
    Silnik do generowania map kontekstowych (MapPlot 20x20°).
    Tu można korzystać z funkcji pomocniczych (np. obrysów).
    """

    def __init__(self, style: PlotStyle, out_dir: Path):
        self.style = style
        self.out_dir = out_dir

    def _choose_projection(self, ra_center_deg: float, dec_center: float):
        """
        Wybór projekcji w zależności od |Dec|.
        """
        abs_dec = abs(dec_center)
        if abs_dec <= 60:
            return Mercator(center_ra=ra_center_deg, center_dec=dec_center)
        elif abs_dec <= 75:
            return Miller(center_ra=ra_center_deg, center_dec=dec_center)
        else:
            return LambertAzEqArea(center_ra=ra_center_deg, center_dec=dec_center)

    def render_context(self, obj: Dict[str, Any]) -> Path | None:
        """
        Druga mapa: MapPlot 20x20 stopni.

        Obrysy (okrąg + krzyżyk) dla dużych / niestandardowych obiektów.
        """
        name = obj.get("id") or obj.get("name")
        name = str(name)
        ra_center_deg = float(obj["ra"])   # RA w stopniach
        dec_center = float(obj["dec"])     # Dec w stopniach

        # 20x20 stopni – prosty box wokół środka
        ra_min = ra_center_deg - 10.0
        ra_max = ra_center_deg + 10.0
        dec_min = dec_center - 10.0
        dec_max = dec_center + 10.0

        # Przycięcie Dec do zakresu -90..+90
        dec_min = max(-90.0, dec_min)
        dec_max = min(90.0, dec_max)

        proj = self._choose_projection(ra_center_deg, dec_center)

        try:
            p = MapPlot(
                projection=proj,
                ra_min=ra_min,
                ra_max=ra_max,
                dec_min=dec_min,
                dec_max=dec_max,
                style=self.style,
                resolution=4600,
                scale=1.2,
            )

            # Konstelacje
            p.constellations()

            # Gwiazdy – umiarkowany limit
            p.stars(
                where=[
                    _.magnitude < 12,
                ],
                bayer_labels=True,
                flamsteed_labels=True,
            )

            # DSO – filtr typów i jasności jak w bak_3
            p.dsos(
                where=[
                    _.type.isin([
                        DsoType.OPEN_CLUSTER.value,
                        DsoType.GALAXY.value,
                        DsoType.GALAXY_PAIR.value,
                        DsoType.GALAXY_TRIPLET.value,
                        DsoType.GROUP_OF_GALAXIES.value,
                        DsoType.GLOBULAR_CLUSTER.value,
                        DsoType.NEBULA.value,
                        DsoType.EMISSION_NEBULA.value,
                        DsoType.REFLECTION_NEBULA.value,
                        DsoType.PLANETARY_NEBULA.value,
                        DsoType.HII_IONIZED_REGION.value,
                        DsoType.STAR_CLUSTER_NEBULA.value,
                        DsoType.DARK_NEBULA.value,]
                    ),
                    (_.magnitude.isnull()) | (_.magnitude < 16),
                ],
                where_true_size=[_.size > 1],
            )

            p.constellation_labels()

            p.gridlines(
                dec_locations=[d for d in range(-90, 91, 2)],
                ra_locations=[d for d in range(int(ra_min), int(ra_max), 2)],
                ra_formatter_fn=lambda ra: f"{round(ra * 15)}\u00B0",
            )

            # ===== OBRYS NIESTANDARDOWY =====
            try:
                size_arcmin = float(obj.get("size", 0) or 0)
            except Exception:
                size_arcmin = 0.0

            if 0 < size_arcmin < 3000:
                c_ra_deg, c_dec_deg = generate_circle_coords_deg(
                    ra_center_deg, dec_center, size_arcmin
                )

                if c_ra_deg:
                    outline_style = LineStyle(
                        color="#e74c3c", width=2, style="dashed", zorder=10
                    )
                    outline_coords = list(zip(c_ra_deg, c_dec_deg))
                    p.line(style=outline_style, coordinates=outline_coords, label=name)

                    cross_size_deg = 0.15
                    cos_dec = math.cos(math.radians(dec_center))
                    if cos_dec < 0.01:
                        cos_dec = 0.01

                    d_ra_cross = (cross_size_deg / 2.0) / cos_dec
                    d_dec_cross = cross_size_deg / 2.0

                    cross_style = LineStyle(
                        color="#e74c3c", width=1, style="solid", zorder=10
                    )

                    h_cross = [
                        (ra_center_deg - d_ra_cross, dec_center),
                        (ra_center_deg + d_ra_cross, dec_center),
                    ]

                    v_start = max(-90.0, dec_center - d_dec_cross)
                    v_end = min(90.0, dec_center + d_dec_cross)
                    v_cross = [(ra_center_deg, v_start), (ra_center_deg, v_end)]

                    p.line(style=cross_style, coordinates=h_cross)
                    p.line(style=cross_style, coordinates=v_cross)

            out_png = self.out_dir / f"context_{name}.png"
            p.export(str(out_png), padding=0.5, transparent=False)
            return out_png

        except Exception as e:
            print(f"Błąd przy mapie kontekstowej dla {name}: {e}")
            return None

        finally:
            plt.close("all")
            if p is not None:
                del p
            gc.collect()


# --- GENERATORY FOV/CTX DLA LISTY OBIEKTÓW ---
def generate_fov_pngs(objs: List[Dict[str, Any]]) -> List[Path | None]:
    ensure_starplots_dir()
    
    cam_params = get_camera_params(VIS_DATA_PATH)
    tasks = [(obj, cam_params, STARPLOTS_DIR) for obj in objs]

    print(f"\nGenerowanie {len(objs)} kadrów FOV (Równolegle)...")
    all_results: Dict[str, Path] = {}

    with ProcessPoolExecutor(max_workers=mp.cpu_count()) as executor:
        future_to_objid = {
            executor.submit(_worker_render_fov, args): obj["id"]
            for args, obj in zip(tasks, objs)
        }

        with tqdm(total=len(objs), desc="Engine FOV", unit="obj") as pbar:
            for future in as_completed(future_to_objid):
                obj_id = future_to_objid[future]
                try:
                    out_path = future.result()
                except Exception as e:
                    print(f"Błąd przy obiekcie {obj_id}: {e}")
                    out_path = None

                if out_path is not None:
                    all_results[obj_id] = out_path

                pbar.set_postfix_str(f"ostatni={obj_id}")
                pbar.update(1)

    valid_paths = list(all_results.values())
    print(f"Wygenerowano {len(valid_paths)} kadrów FOV.")
    return valid_paths


def generate_context_pngs(objs: List[Dict[str, Any]]) -> List[Path | None]:
    ensure_starplots_dir()
    
    tasks = [(obj, STARPLOTS_DIR) for obj in objs]
    print(f"Generowanie {len(objs)} map kontekstowych (Równolegle)...")

    all_results: Dict[str, Path] = {}

    with ProcessPoolExecutor(max_workers=mp.cpu_count()) as executor:
        future_to_objid = {
            executor.submit(_worker_render_context, args): obj["id"]
            for args, obj in zip(tasks, objs)
        }

        with tqdm(total=len(objs), desc="Engine CTX", unit="obj") as pbar:
            for future in as_completed(future_to_objid):
                obj_id = future_to_objid[future]
                try:
                    out_path = future.result()
                except Exception as e:
                    print(f"Błąd przy mapie kontekstowej dla {obj_id}: {e}")
                    out_path = None

                if out_path is not None:
                    all_results[obj_id] = out_path

                pbar.set_postfix_str(f"ostatni={obj_id}")
                pbar.update(1)

    valid_paths = list(all_results.values())
    print(f"Wygenerowano {len(valid_paths)} map kontekstowych.\n")
    return valid_paths


# --- MAIN ---
def main() -> None:
    print("Pobieram katalogi dla StarPlot.")
    preload_open_ngc()
    vis_data = load_vis_data()
    df_objects = build_objects_from_vis_data(vis_data)

    if df_objects.empty:
        print("Brak obiektów z uzupełnionym polem 'selected' w vis_data.json.")
        return

    selected_objects = df_objects.to_dict(orient="records")
    print(f"Wybrano {len(selected_objects)} obiektów (wszystkie z 'selected' w vis_data.json).")
    
    generate_fov_pngs(selected_objects)
    generate_context_pngs(selected_objects)

    print("Gotowe: FOV i context PNG wygenerowane w katalogu 'starplots/'.")


if __name__ == "__main__":
    main()