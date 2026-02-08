import napari
import numpy as np
import pandas as pd

from magicgui import magic_factory
from shapely.geometry import Point, Polygon
from napari.layers import Shapes, Points

from spatialpipe.neighborhoods import compute_region_expression

from matplotlib.figure import Figure
from matplotlib.backends.backend_qt5agg import FigureCanvasQTAgg as FigureCanvas
from qtpy.QtWidgets import QWidget, QVBoxLayout


# =========================================================
# Helpers
# =========================================================

def _current_viewer():
    try:
        return napari.current_viewer()
    except Exception:
        return None


def _get_gene_choices(_=None):
    viewer = _current_viewer()
    if viewer is None:
        return []

    if "tissue" not in viewer.layers:
        return []

    return viewer.layers["tissue"].metadata.get("genes", [])


def _collect_polygons(viewer):
    polygons, names = [], []

    for layer in viewer.layers:
        if isinstance(layer, Shapes):
            for i, shape in enumerate(layer.data):
                if len(shape) >= 3:
                    polygons.append(Polygon(shape))
                    names.append(f"{layer.name}_{i+1}")

    return polygons, names


def _points_inside_polygons(points, polygons):
    mask = np.zeros(len(points), dtype=bool)
    for poly in polygons:
        mask |= np.array(
            [poly.contains(Point(x, y)) for y, x in points],
            dtype=bool,
        )
    return mask


# =========================================================
# NAPARI PLUGIN WIDGET (MUST BE MagicFactory)
# =========================================================

@magic_factory(
    genes={"choices": _get_gene_choices, "label": "Select genes"},
    call_button="Compute",
)
def region_expression_widget(
    viewer: napari.Viewer,
    genes: list[str],
):
    # ---------- safety ----------
    if viewer is None or not genes:
        return

    if "tissue" not in viewer.layers:
        print("⚠️ Missing 'tissue' layer")
        return

    tissue: Points = viewer.layers["tissue"]

    if "sdata" not in tissue.metadata:
        print("⚠️ tissue.metadata['sdata'] missing")
        return

    polygons, names = _collect_polygons(viewer)
    if not polygons:
        print("⚠️ Draw at least one polygon")
        return

    sdata = tissue.metadata["sdata"]

    # ---------- compute ----------
    records = []
    for poly, name in zip(polygons, names):
        expr = compute_region_expression(sdata, poly, genes)
        for g, v in expr.items():
            records.append(
                {"region": name, "gene": g, "expression": float(v)}
            )

    df = pd.DataFrame(records)
    if df.empty:
        print("⚠️ No expression found")
        return

    # ---------- build Qt plot ----------
    fig = Figure(figsize=(6, 4))
    canvas = FigureCanvas(fig)
    ax = fig.add_subplot(111)

    for region in df.region.unique():
        sub = df[df.region == region]
        ax.bar(sub.gene, sub.expression, label=region, alpha=0.75)

    ax.set_ylabel("Expression")
    ax.set_title("Spatial Region Expression")
    ax.legend()
    fig.tight_layout()

    plot_widget = QWidget()
    layout = QVBoxLayout(plot_widget)
    layout.addWidget(canvas)

    viewer.window.add_dock_widget(
        plot_widget,
        area="right",
        name="Region Expression Plot",
    )

    # ---------- highlight ----------
    mask = _points_inside_polygons(tissue.data, polygons)

    if "Selected points" in viewer.layers:
        viewer.layers.remove("Selected points")

    viewer.add_points(
        tissue.data[mask],
        name="Selected points",
        size=6,
        face_color="red",
        opacity=0.9,
    )

    df.to_csv("region_expression.csv", index=False)
    print("✅ Saved region_expression.csv")
