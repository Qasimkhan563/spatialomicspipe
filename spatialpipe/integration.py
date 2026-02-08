from spatialdata import SpatialData
from anndata import AnnData
import anndata as ad
from typing import List


def merge_spatialdata(
    sdata_list: List[SpatialData],
    table_key: str = "spots",
) -> SpatialData:
    """
    Merge multiple SpatialData objects into a single SpatialData.

    Parameters
    ----------
    sdata_list : list[SpatialData]
        List of SpatialData objects to merge.
    table_key : str
        Table key to merge (default: "spots").

    Returns
    -------
    SpatialData
        New SpatialData object with merged table and a `sample_id` column.

    Notes
    -----
    Each input SpatialData is assumed to contain the same table_key.
    A `sample_id` column is added to track the origin of each observation.
    """
    if len(sdata_list) == 0:
        raise ValueError("sdata_list must contain at least one SpatialData object")

    tables: list[AnnData] = []

    for i, sdata in enumerate(sdata_list):
        if table_key not in sdata.tables:
            raise KeyError(f"Table key '{table_key}' not found in SpatialData object")

        adata: AnnData = sdata.tables[table_key].copy()
        adata.obs["sample_id"] = f"sample_{i}"
        tables.append(adata)

    merged = ad.concat(tables, join="outer", axis=0)

    merged_sdata = SpatialData(
        tables={table_key: merged},
    )

    return merged_sdata
