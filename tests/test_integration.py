import squidpy as sq
from spatialpipe.io import load_visium_to_spatialdata
from spatialpipe.integration import merge_spatialdata


def test_merge_spatialdata():
    """
    Test that multiple SpatialData objects can be merged and
    that sample identity is preserved via a `sample_id` column.
    """
    adata1 = sq.datasets.visium()
    adata2 = sq.datasets.visium()

    sdata1 = load_visium_to_spatialdata(adata1)
    sdata2 = load_visium_to_spatialdata(adata2)

    merged = merge_spatialdata([sdata1, sdata2])

    assert "spots" in merged.tables
    assert "sample_id" in merged.tables["spots"].obs

    assert merged.tables["spots"].n_obs == (
        sdata1.tables["spots"].n_obs + sdata2.tables["spots"].n_obs
    )
