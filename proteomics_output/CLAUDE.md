# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## Project Overview

PTX3 proteomics analysis pipeline — compares PTX3-treated vs control human cells (60h, 3 vs 3 replicates) using label-free quantitative mass spectrometry data. Produces publication-quality figures, pathway enrichment analysis, and a comprehensive text report.

## Running the Scripts

```bash
# Step 1: Run main analysis (generates individual plots, enrichment CSVs, DEP tables, report)
python proteomics_analysis.py

# Step 2: Run panel/bubble generator (reads enrichment CSVs from step 1, generates composite figures)
python make_panels.py
```

**Execution order matters**: `proteomics_analysis.py` must run first — it produces the `enrichment_*.csv` files that `make_panels.py` reads. Both scripts require network access (Enrichr API calls in step 1).

## Architecture

- **`proteomics_analysis.py`** — Main pipeline. Reads raw CSV input, performs differential expression classification (based on pre-computed `Sig` column), generates 6 individual publication figures (volcano, heatmap, bar chart, correlation, distributions, scatter), queries the Enrichr REST API for KEGG/GO enrichment, saves enrichment CSVs + DEP tables + full report.
- **`make_panels.py`** — Post-processing. Loads the enrichment CSVs produced by the main script, generates composite figures: a 5-panel summary (`00_MAIN_PANEL`), GO BP bubble plots, KEGG combined bar chart. Also overwrites `ANALYSIS_REPORT.txt` with an expanded version including biological interpretation.

## Key Conventions

- Input data path is hardcoded: `D:\桌面\2026.3.12 pr-ptx3 60h 蛋白质谱3vs3.csv`
- Output directory is hardcoded: `D:\桌面\proteomics_output`
- All plots are saved as both PDF (vector) and PNG (300 DPI)
- Output files are numbered `00_` through `14_` for ordering
- Sample columns: `CTRL-1/2/3` and `PTX3-1/2/3`; derived columns: `log2FC`, `P.value`, `Sig` (up/down), `mean c`, `mean ptx3`
- Enrichment uses the Enrichr API (`maayanlab.cloud/Enrichr`) with libraries: `KEGG_2021_Human`, `GO_Biological_Process_2021`, `GO_Molecular_Function_2021`
- Matplotlib backend is set to `Agg` (non-interactive) for headless rendering
- Color scheme: Up-regulated = `#D73027` (red), Down-regulated = `#4575B4` (blue), NS = `#BBBBBB`/`#CCCCCC` (gray)

## Dependencies

Python packages: `numpy`, `pandas`, `matplotlib`, `seaborn`, `scipy`. No requirements.txt exists — install manually via pip.
