import zarr
import dask.array as da
import numpy as np
from spatialdata import SpatialData
from anndata import AnnData
from pathlib import Path
from typing import Tuple


def to_lazy_anndata(
    adata: AnnData,
    chunks: Tuple[int, int] = (1000, 1000),
) -> AnnData:
    """
    Convert AnnData.X to a Dask-backed array.

    Parameters
    ----------
    adata : AnnData
        Input AnnData object.
    chunks : tuple[int, int]
        Chunk size for Dask array.

    Returns
    -------
    AnnData
        AnnData with X stored as a Dask array.
    """
    if not isinstance(adata.X, da.Array):
        # Handle sparse or dense matrices
        if hasattr(adata.X, "toarray"):
            array = adata.X.toarray()
        else:
            array = np.asarray(adata.X)

        adata.X = da.from_array(array, chunks=chunks)

    return adata


def make_spatialdata_lazy(
    sdata: SpatialData,
    table_key: str = "spots",
    chunks: Tuple[int, int] = (1000, 1000),
) -> SpatialData:
    """
    Convert a SpatialData table to a lazy (Dask-backed) representation.

    Parameters
    ----------
    sdata : SpatialData
        Input SpatialData object.
    table_key : str
        Table key inside SpatialData (default: "spots").
    chunks : tuple[int, int]
        Chunk size for Dask array.

    Returns
    -------
    SpatialData
        New SpatialData object with lazy-backed AnnData table.
    """
    adata = sdata.tables[table_key].copy()
    adata = to_lazy_anndata(adata, chunks=chunks)

    lazy_sdata = SpatialData(
        tables={table_key: adata},
    )

    return lazy_sdata


def save_lazy_spatialdata(
    sdata: SpatialData,
    path: str | Path,
    overwrite: bool = True,
):
    """
    Save SpatialData to Zarr with chunking preserved.

    Parameters
    ----------
    sdata : SpatialData
        SpatialData object to save.
    path : str or Path
        Output directory.
    overwrite : bool
        Whether to overwrite existing store.
    """
    path = Path(path)
    store = zarr.DirectoryStore(path)
    root = zarr.group(store=store, overwrite=overwrite)

    sdata.write(root)
