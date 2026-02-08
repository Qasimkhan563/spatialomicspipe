from pathlib import Path
import yaml
from typing import Union
from spatialpipe.io import load_visium_to_spatialdata, save_spatialdata_to_omezarr
from spatialdata import SpatialData


def run_ingest_pipeline(config_path: Union[str, Path]) -> SpatialData:
    """
    Run an ingestion pipeline from a YAML configuration file.

    Parameters
    ----------
    config_path : str or Path
        Path to YAML configuration file with input/output settings.

    Returns
    -------
    SpatialData
        Ingested SpatialData object.
    """
    config_path = Path(config_path)

    if not config_path.exists():
        raise FileNotFoundError(f"Config file not found: {config_path}")

    with open(config_path, "r") as f:
        config = yaml.safe_load(f)

    try:
        input_path = Path(config["input"]["path"])
        output_path = Path(config["output"]["omezarr"])
    except KeyError as e:
        raise KeyError("Config must contain 'input.path' and 'output.omezarr'") from e

    sdata = load_visium_to_spatialdata(input_path)
    save_spatialdata_to_omezarr(sdata, output_path)

    return sdata
