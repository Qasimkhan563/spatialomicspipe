import squidpy as sq
import numpy as np
from shapely.geometry import Polygon
from spatialpipe.io import load_visium_to_spatialdata
from spatialpipe.neighborhoods import compute_region_expression


def test_region_expression():
    """
    Test that gene expression can be aggregated inside a polygon
    region and returns non-negative values.
    """
    adata = sq.datasets.visium()
    sdata = load_visium_to_spatialdata(adata)

    coords = adata.obsm["spatial"]

    # Build polygon that definitely contains some spots
    minx, miny = coords.min(axis=0)
    maxx, maxy = coords.max(axis=0)

    poly = Polygon(
        [
            (minx, miny),
            (maxx, miny),
            (maxx, maxy),
            (minx, maxy),
        ]
    )

    genes = adata.var_names[:5].tolist()
    expr = compute_region_expression(sdata, poly, genes)

    assert len(expr) == 5
    assert (expr.values >= 0).all()
