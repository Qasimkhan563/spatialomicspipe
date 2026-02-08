import napari
import scanpy as sc
from spatialpipe.io import load_visium_to_spatialdata
from spatialpipe.neighborhoods import compute_region_expression
from shapely.geometry import Polygon
import numpy as np


class SpatialPipeWidget:
    def __init__(self, viewer: napari.Viewer):
        self.viewer = viewer
        self.sdata = None

    def load_visium(self):
        adata = sc.datasets.visium_sge()
        self.sdata = load_visium_to_spatialdata(adata)

        coords = adata.obsm["spatial"]
        self.viewer.add_points(coords, size=3, name="spots")

    def compute_region_stats(self):
        shapes = self.viewer.layers["Shapes"].data
        if len(shapes) == 0:
            print("Draw a polygon first!")
            return

        poly = Polygon(shapes[-1])
        genes = ["Reln", "Slc17a7", "Gad1"]

        expr = compute_region_expression(self.sdata, poly, genes)
        print("Region expression:")
        print(expr)


def napari_experimental_provide_dock_widget():
    return SpatialPipeWidget
