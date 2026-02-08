# SpatialPipe

**SpatialPipe** is a scalable, NGFF-native framework for spatial omics analysis built on
the `SpatialData` ecosystem.

It provides:
- Conversion of spatial transcriptomics data to `SpatialData`
- Polygon-based spatial neighborhood queries
- Multi-sample integration
- Lazy, Dask-backed computation
- Zarr/OME-NGFF compatible storage

Designed as a lightweight pipeline framework for reproducible spatial omics workflows.

---

## Features

✔ SpatialData ingestion (Visium)  
✔ Region-based gene expression queries  
✔ Spatial neighborhood graphs  
✔ Multi-sample merging  
✔ Lazy Dask-backed backend  
✔ Zarr persistence  
✔ Tested API  

---

## Installation

```bash
git clone https://github.com/yourname/spatialdata-pipelines
cd spatialdata-pipelines
python -m poetry install
