import napari
import squidpy as sq
from spatialpipe.io import load_visium_to_spatialdata

adata = sq.datasets.visium(sample_id="V1_Mouse_Brain_Sagittal_Posterior")
sdata = load_visium_to_spatialdata(adata)

viewer = napari.Viewer()

points = viewer.add_points(
    adata.obsm["spatial"],
    size=5,
    name="Spots",
)

# attach genes + SpatialData to layer metadata
points.metadata["genes"] = adata.var_names.tolist()
points.metadata["sdata"] = sdata

viewer.add_shapes(name="Shapes")

napari.run()
