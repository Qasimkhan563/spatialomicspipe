import scanpy as sc
from spatialdata import SpatialData
from spatialdata.models import TableModel
from anndata import AnnData
import zarr
from pathlib import Path
from typing import Union


def load_visium_to_spatialdata(data: Union[str, Path, AnnData]) -> SpatialData:
    """
    Convert Visium spatial transcriptomics data (path or AnnData)
    into a SpatialData object with region and instance annotations.

    Parameters
    ----------
    data : str, Path, or AnnData
        Path to Visium dataset or preloaded AnnData object.

    Returns
    -------
    SpatialData
        SpatialData object containing a 'spots' table.
    """
    if isinstance(data, AnnData):
        adata = data
    elif isinstance(data, (str, Path)):
        adata = sc.read_visium(data)
    else:
        raise TypeError("data must be a path (str/Path) or an AnnData object")

    # Ensure unique identifiers
    adata.obs_names_make_unique()
    adata.var_names_make_unique()

    # Create instance IDs
    adata.obs["spot_id"] = adata.obs_names.astype(str)

    # Map in_tissue (0/1) -> semantic region labels
    adata.obs["region"] = adata.obs["in_tissue"].map(
        {1: "tissue", 0: "background"}
    ).astype(str)

    # Only use regions that actually exist
    regions = sorted(adata.obs["region"].dropna().unique().tolist())

    table = TableModel.parse(
        adata,
        region=regions,
        region_key="region",
        instance_key="spot_id",
    )

    sdata = SpatialData(tables={"spots": table})
    return sdata


def save_spatialdata_to_omezarr(
    sdata: SpatialData,
    store_path: Union[str, Path],
    overwrite: bool = True,
):
    """
    Save a SpatialData object to an OME-Zarr compatible store.

    Parameters
    ----------
    sdata : SpatialData
        SpatialData object to save.
    store_path : str or Path
        Output directory for the Zarr store.
    overwrite : bool
        Whether to overwrite an existing store.
    """
    store_path = Path(store_path)
    store = zarr.DirectoryStore(store_path)
    root = zarr.group(store=store, overwrite=overwrite)
    sdata.write(root)
