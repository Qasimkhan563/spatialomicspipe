import squidpy as sq
from spatialpipe.io import load_visium_to_spatialdata


def test_load_visium():
    """
    Test that Visium data can be converted into a SpatialData object
    containing a 'spots' table.
    """
    adata = sq.datasets.visium()
    sdata = load_visium_to_spatialdata(adata)

    assert "spots" in sdata.tables
