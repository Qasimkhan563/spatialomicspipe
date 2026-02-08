import squidpy as sq
import dask.array as da
from spatialpipe.io import load_visium_to_spatialdata
from spatialpipe.lazy import make_spatialdata_lazy


def test_lazy_backend():
    """
    Test that a SpatialData object can be converted to a lazy
    Dask-backed representation.
    """
    adata = sq.datasets.visium()
    sdata = load_visium_to_spatialdata(adata)

    lazy = make_spatialdata_lazy(sdata)

    X = lazy.tables["spots"].X
    assert isinstance(X, da.Array)
