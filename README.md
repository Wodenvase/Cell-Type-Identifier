# 🧬 scRNA-seq Cell Type Identifier 🔬

A multi-notebook scRNA-seq analysis pipeline that takes raw single-cell gene expression data and automatically identifies cell types using dimensionality reduction, clustering, and marker-gene analysis.  
This repository is structured as a set of small, focused Jupyter notebooks so each stage of the workflow is easy to read, run, and modify.

---

## 📂 Project structure
```
.
├── data/                                  # raw & intermediate .h5ad files (git-lfs or download externally)
├── figures/                               # saved figures and plots
├── notebooks/
│   ├── 01_Data_Acquisition_and_Loading.ipynb
│   ├── 02_Data_Preprocessing_and_QC.ipynb
│   ├── 03_Dimensionality_Reduction.ipynb
│   ├── 04_Clustering.ipynb
│   └── 05_Visualization_and_Interpretation.ipynb
├── environment.yml                         # optional conda env file
├── requirements.txt                        # pip installable requirements
└── README.md
```
---

## 🧾 Overview of notebooks

- **01_Data_Acquisition_and_Loading.ipynb**  
  Download/load the dataset into an `AnnData` object (`.h5ad`). Inspect metadata, counts, and data shape. Example: downloading 10x "3k PBMCs" and creating `pbmc3k_raw.h5ad`.

- **02_Data_Preprocessing_and_QC.ipynb**  
  Quality control and filtering: remove low-quality cells & lowly-expressed genes, compute mitochondrial gene fraction, normalize counts, log-transform, and identify highly variable genes. Saves `pbmc3k_processed.h5ad`.

- **03_Dimensionality_Reduction.ipynb**  
  Perform PCA for linear reduction, compute neighborhood graph, and embed with UMAP for visualization. Save UMAP coordinates into the `AnnData` object.

- **04_Clustering.ipynb**  
  Build the neighborhood graph (if not already), run Leiden clustering across chosen resolutions and attach cluster labels to `adata.obs`. Explore cluster stability.

- **05_Visualization_and_Interpretation.ipynb**  
  Identify cluster marker genes, produce cluster heatmaps / dotplots, and manually or semi-automatically annotate clusters (B cell, T cell, NK, Monocyte, etc.). Save final annotated object `pbmc3k_final_annotated.h5ad` and high-resolution figures to `figures/`.

---

## 📊 Dataset
This pipeline uses the 10x Genomics "3k PBMCs from a Healthy Donor" (approx. 2,700 cells). The notebooks show how to load this dataset via `scanpy.datasets` or download the official 10x files and convert to `AnnData`.

---

## ⚙️ Requirements & Installation

**Using pip (recommended in a virtualenv):**
```bash
python -m venv .venv
source .venv/bin/activate        # macOS / Linux
# .venv\Scripts\activate         # Windows (PowerShell/CMD)

pip install --upgrade pip
pip install -r requirements.txt
```
requirements.txt
```
scanpy>=1.9
anndata
python-igraph
leidenalg
pandas
numpy
matplotlib
seaborn
scikit-learn
jupyterlab
```
Optional (conda) environment.yml:
```
name: scrnaseq
channels:
  - conda-forge
dependencies:
  - python=3.10
  - scanpy
  - anndata
  - python-igraph
  - leidenalg
  - pandas
  - numpy
  - matplotlib
  - seaborn
  - scikit-learn
  - jupyterlab
```

⸻

🚀 How to run
	1.	Clone the repo and create the environment (see above).
	2.	Start Jupyter Lab or Notebook:

jupyter lab
# or
jupyter notebook

	3.	Open and execute notebooks in order: 01_ → 02_ → 03_ → 04_ → 05_.

Non-interactive / batch execution

If you want to run notebooks headlessly and save outputs, use papermill or nbconvert. Example with papermill:

pip install papermill
papermill notebooks/01_Data_Acquisition_and_Loading.ipynb notebooks/out/01_out.ipynb


⸻

🔁 Reproducibility & tips
	•	Fix random seeds where relevant (PCA/UMAP/leiden) to help reproducibility.
	•	Keep .h5ad intermediate files in data/ to avoid re-running expensive steps.
	•	If memory is limited, subsample cells for exploratory analysis or use scanpy.pp.subsample.
	•	If UMAP produces different layouts between runs, increase min_dist/n_neighbors or use sc.tl.pca(…).obsm['X_pca'] as initialization with a fixed seed.

⸻

✅ Expected outputs
	•	data/pbmc3k_raw.h5ad — raw loaded dataset
	•	data/pbmc3k_processed.h5ad — filtered & normalized data (post-QC)
	•	data/pbmc3k_final_annotated.h5ad — final AnnData with .obs['cell_type'] annotations
	•	figures/ — UMAP plots, marker heatmaps, dotplots, cluster diagnostics

⸻

🧭 Quick example: minimal code snippet (Scanpy)
```
import scanpy as sc

# load example (or replace with data path)
adata = sc.datasets.pbmc3k_processed()   # or sc.read_10x_h5(path_to_10x_h5)
sc.pp.filter_cells(adata, min_genes=200)
sc.pp.filter_genes(adata, min_cells=3)
adata.var['mt'] = adata.var_names.str.startswith('MT-')
sc.pp.calculate_qc_metrics(adata, qc_vars=['mt'], inplace=True)
sc.pp.normalize_total(adata, target_sum=1e4)
sc.pp.log1p(adata)
sc.pp.highly_variable_genes(adata, n_top_genes=2000)
sc.tl.pca(adata)
sc.pp.neighbors(adata)
sc.tl.umap(adata)
sc.tl.leiden(adata, key_added='leiden')
sc.tl.rank_genes_groups(adata, 'leiden', method='t-test')
```

⸻

🛠️ Troubleshooting
	•	Import errors (leiden/igraph): Make sure python-igraph and leidenalg are installed from conda-forge or compatible wheels. If using pip, install python-igraph first then leidenalg.
	•	Slow UMAP / memory issues: Try reducing n_pcs or subsampling cells. Use sc.pp.highly_variable_genes to reduce gene count.
	•	Cluster over/under-splitting: Try different resolution values in sc.tl.leiden(adata, resolution=0.5).

⸻

🧾 Citation & Further reading

If you use this pipeline in work or publications, cite Scanpy and the primary algorithms (UMAP, Leiden). See the notebooks for references and links to original papers.

⸻

🧑‍💻 Extending this project
	•	Add automated cluster annotation using reference datasets (e.g., SingleR, scmap, or Azimuth).
	•	Integrate batch-correction methods (Harmony, BBKNN) for multi-sample inputs.
	•	Convert parts of the pipeline into reusable Python scripts or a CLI.
	•	Create a lightweight Streamlit / Dash / Panel app to interactively explore clusters & markers.

⸻

📬 Contact

Questions, suggestions or contributions welcome — open an issue or send a PR.

⸻
