# SpatialOmicsPipe

> **Status:** 🚧 Active development (research / prototype stage)

SpatialOmicsPipe is an **NGFF-native spatial omics analysis framework** built on top of **SpatialData**, **OME-Zarr**, and the geospatial Python ecosystem. It focuses on **interactive, region-based spatial analysis** with first‑class support for **napari** as a visual analytics frontend.

This repository represents a **working, tested development version** of SpatialOmicsPipe. Core functionality is implemented and validated locally, while packaging, UX polish, and CI/CD are still evolving.

---

## 🔬 Project Motivation

Modern spatial transcriptomics produces *spatially indexed, multi‑modal data* that is difficult to explore with traditional single‑cell pipelines alone.

SpatialOmicsPipe was created to:

- Treat spatial omics data as **geometric objects** (points, regions, neighborhoods)
- Leverage **SpatialData / NGFF** as the canonical data model
- Enable **interactive region selection** and analysis inside napari
- Bridge **scanpy / squidpy analytics** with **geospatial workflows**

---

## ✨ Key Features

### Core Library (`SpatialPipe`)

- ✅ NGFF‑native data handling via **SpatialData**
- ✅ Region‑based expression computation
- ✅ Polygon‑driven spatial queries
- ✅ Compatible with Visium / spot‑based data

### napari Plugin (`napari_SpatialPipe`)

- ✅ Auto‑discovered napari plugin (manifest‑based)
- ✅ Interactive polygon drawing (Shapes layer)
- ✅ Region‑wise gene expression computation
- ✅ Visual highlighting of selected spatial points
- ✅ CSV export of region expression results

### Developer‑Focused

- ✅ Poetry‑managed dependencies
- ✅ Python ≥ 3.11
- ✅ Modular package layout
- ✅ Ready for future CI / Dockerization

---

## 📁 Repository Structure

```text
spatialdata-pipelines/
│
├── SpatialPipe/                 # Core analysis library
│   ├── neighborhoods/           # Region / neighborhood computations
│   ├── io/                      # Data loading helpers
│   ├── cli.py                   # CLI entry point
│   └── __init__.py
│
├── napari_SpatialPipe/           # napari plugin package
│   ├── widget.py                # Interactive napari widget
│   ├── napari.yaml              # napari plugin manifest
│   └── __init__.py
│
├── tests/                        # Unit & integration tests
│
├── pyproject.toml                # Project + dependency specification
├── poetry.lock                  # Locked dependency graph
├── README.md                    # (this file)
└── .gitignore
```

---

## ⚙️ Installation (Development Setup)

> **Recommended:** Local development install using Poetry

### 1️⃣ Clone the repository

```bash
git clone https://github.com/Qasimkhan563/spatialomicspipe.git
cd spatialomicspipe
```

### 2️⃣ Install dependencies

```bash
poetry install
```

> ⚠️ This project targets **Python 3.11–3.13**. Older versions are not supported.

### 3️⃣ Activate the environment

```bash
poetry shell
```

---

## ▶️ Running napari + SpatialOmicsPipe

### Launch napari

```bash
napari
```

SpatialOmicsPipe should automatically appear under:

> **Plugins → SpatialOmicsPipe → Spatial Region Expression**

---

## 🧪 Example Workflow (Tested)

The following workflow has been **executed and validated locally**:

### 1️⃣ Load Visium demo data

```python
import scanpy as sc
from SpatialPipe.io import load_visium_to_spatialdata

adata = sc.datasets.visium_sge()
sdata = load_visium_to_spatialdata(adata)
```

### 2️⃣ Add spatial points to napari

```python
import napari
viewer = napari.current_viewer()

points = sdata.tables["spots"]
viewer.add_points(
    points.obsm["spatial"],
    name="tissue",
    size=6,
    face_color="white",
)

# attach SpatialData to layer metadata
viewer.layers["tissue"].metadata["sdata"] = sdata
viewer.layers["tissue"].metadata["genes"] = list(points.var_names)
```

### 3️⃣ Draw regions

- Select **Shapes** layer
- Draw one or more **polygons** over the tissue

### 4️⃣ Run region expression

- Open the **Spatial Region Expression** widget
- Select one or more genes
- Click **Compute**

**Results:**
- Bar plot of expression per region
- Highlighted spatial points inside regions
- `region_expression.csv` exported

---

## 📊 Outputs

- 📈 Region‑wise gene expression plots
- 🧬 CSV export of expression values
- 🔴 Highlighted spatial points per region

---

## 🧠 Design Decisions

- **SpatialData first** — no ad‑hoc coordinate handling
- **napari Shapes** used as the canonical region definition
- **No hidden global state** — all data flows through layer metadata
- **Fail‑safe widgets** — errors print warnings instead of crashing napari

---

## ⚠️ Known Limitations (Current)

- ❌ UX still experimental (research‑grade)
- ❌ Widget API may change
- ❌ Large datasets not yet performance‑optimized

These are **planned improvements**, not architectural blockers.

---

## 🛣️ Roadmap (Planned)

- [ ] Dockerized reproducible environment
- [ ] GitHub Actions CI
- [ ] Multi‑region comparison UI
- [ ] Neighborhood / buffer‑based analysis
- [ ] Performance optimizations for large tissues
- [ ] Documentation website

---

## 🧪 Development Status

> **This repository represents a validated development prototype.**

- Core algorithms implemented
- Plugin integration tested locally
- Active iteration ongoing

The project is **not abandoned**, **not broken**, and **not theoretical** — it is in an development phase.

---

## 👤 Author

**Muhammad Qasim**  
Geospatial Data Professional | EO, GIS & Applied Analytics

📧 khanjiqasim@gmail.com

---

## 📜 License

MIT License

---

## 🤝 Contributions

Issues, discussions, and suggestions are welcome.

