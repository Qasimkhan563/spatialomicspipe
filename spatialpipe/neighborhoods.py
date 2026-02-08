import geopandas as gpd
from shapely.geometry import Point, Polygon
import pandas as pd
from spatialdata import SpatialData
from anndata import AnnData
import numpy as np


def spots_to_geodataframe(
    sdata: SpatialData,
    table_key: str = "spots",
) -> gpd.GeoDataFrame:
    """
    Convert a SpatialData table into a GeoDataFrame using spatial coordinates.

    Parameters
    ----------
    sdata : SpatialData
        Input SpatialData object.
    table_key : str
        Key of the table inside SpatialData (default: "spots").

    Returns
    -------
    GeoDataFrame
        GeoDataFrame with spot metadata and Point geometries.
    """
    adata: AnnData = sdata.tables[table_key]

    coords = adata.obsm["spatial"]

    gdf = gpd.GeoDataFrame(
        adata.obs.copy(),
        geometry=[Point(float(x), float(y)) for x, y in coords],
        crs="EPSG:4326",  # arbitrary planar CRS for demo purposes
    )

    return gdf


def compute_region_expression(
    sdata: SpatialData,
    polygon: Polygon,
    genes: list[str],
    table_key: str = "spots",
    agg: str = "mean",
) -> pd.Series:
    """
    Aggregate gene expression inside a polygon region.

    Parameters
    ----------
    sdata : SpatialData
        Input SpatialData object.
    polygon : shapely.geometry.Polygon
        Region of interest.
    genes : list[str]
        Genes to aggregate.
    table_key : str
        Table key inside SpatialData.
    agg : str
        Aggregation method: "mean" or "sum".

    Returns
    -------
    pd.Series
        Aggregated expression values indexed by gene.
    """
    adata: AnnData = sdata.tables[table_key]
    gdf = spots_to_geodataframe(sdata, table_key)

    mask = gdf.geometry.within(polygon)

    # Handle empty regions robustly
    if int(mask.sum()) == 0:
        return pd.Series([0.0] * len(genes), index=genes, dtype=float)

    subset = adata[mask, genes]

    if agg == "mean":
        values = subset.X.mean(axis=0)
    elif agg == "sum":
        values = subset.X.sum(axis=0)
    else:
        raise ValueError("agg must be 'mean' or 'sum'")

    # Handle sparse or dense matrices
    if hasattr(values, "A1"):
        values = values.A1
    else:
        values = np.asarray(values).ravel()

    return pd.Series(values, index=genes, dtype=float)


def build_neighborhood_graph(
    sdata: SpatialData,
    radius: float,
    table_key: str = "spots",
) -> pd.DataFrame:
    """
    Build a radius-based spatial neighborhood graph.

    Parameters
    ----------
    sdata : SpatialData
        Input SpatialData object.
    radius : float
        Distance threshold for neighborhood edges.
    table_key : str
        Table key inside SpatialData.

    Returns
    -------
    pd.DataFrame
        Edge list with columns:
        - spot_i
        - spot_j
        - distance
    """
    adata: AnnData = sdata.tables[table_key]
    coords = adata.obsm["spatial"]

    edges = []
    n = coords.shape[0]

    for i in range(n):
        for j in range(i + 1, n):
            dist = float(np.linalg.norm(coords[i] - coords[j]))
            if dist <= radius:
                edges.append((adata.obs_names[i], adata.obs_names[j], dist))

    return pd.DataFrame(edges, columns=["spot_i", "spot_j", "distance"])
