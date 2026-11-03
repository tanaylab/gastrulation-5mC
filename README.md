# DNMT embryo paper: figure code

Analysis accompanying: **A weak-interaction model defines the cell-autonomous function of DNA methylation in gastrulation**.
Saifeng Cheng, Roni Stok Ranen, Marthe Marie Behrends, Yoav Mayshar, Ofir Raz, Shuailin Du, Aviezer Lifshitz, Ying Zhu, Evgheni Casimov, Tal Sokolov, Akhiad Bercovich, Raz Ben-Yair, Oren Ben-Kiki, Netta Reines, Ayelet-Hashahar Orenbuch, Jinfang Wang, Yonatan Stelzer, and Amos Tanay.

## Layout

- `*.ipynb`: Fig 1-5, Extended Data and Supplementary table notebooks; run from the repo root
- `E9.5_multiome.Fig6/`: Fig 6 (`Fig6.ipynb`, E9.5 multiome of Dnmt-DKO chimera head parts) and its functions (`Fig6_func.r`, `e95_func.r`)
- `prep/`: preprocessing notebooks that build the inputs of the figures (`E9.5_*` for Fig 6)
- `scripts/`: R functions the notebooks source; notebooks define no functions themselves (except the python notebook `prep/embryo_profiles`)
- `paths.R`: locations of the data and external resources; edit this first

How paths work:
- Each R notebook starts with `source("paths.R")` (`source("../paths.R")` in `prep/` and `E9.5_multiome.Fig6/`). `paths.R` finds the repo root itself.
- Root notebooks read `./data/...` from the repo root (link `data`, see Data); notebooks in `prep/` and `E9.5_multiome.Fig6/` change into the original analysis directory first (`embdnmt_dir`, or `embdnmt_dir/revision_e95` for E9.5).
- Everything written goes through `out_path()` (in `paths.R`) into `git_test/` (gitignored).
- Notebooks read cached intermediates (e.g. `data/output/mod.RDS`, `chisq_screen.RDS`, `mod_e95.RDS`) from the original directory. To use the copies regenerated under `git_test/` by `prep/` or the figure notebooks, set `use_regenerated_inputs <- TRUE` in `paths.R`. The reads go through `in_path()`, which takes the `git_test/` copy when it exists and the original otherwise.

Figures:
- Each figure notebook starts with `save_figure = F`. Figures are then only drawn in the notebook; set it to `T` to also write every figure file.
- All figures go through `show_fig()` (in `paths.R`), which draws the figure and, when `save_figure` is `T`, saves it: `show_fig(file, plot)` for ggplot/grid objects, `show_fig(file, draw = {...})` for base graphics and for functions that save through their own `filename` argument (pheatmap).
- Saved but not drawn inline (`show_fig(..., show = FALSE)`): the Fig 4 cluster profiles (~380 small plots), to keep the notebook small, and the stand-alone colour scales (Fig 4, 5a, 5c, 6e).
- In Jupyter, each inline figure has the size in inches of its saved file, so it has the same proportions and layout. Its resolution keeps the longest side at most 300 pixels (`fig_inline_px`, `fig_inline_dpi` in `paths.R`).
- The stored notebooks were run top to bottom with `save_figure = F`; their outputs are the figures only.

## Figure → notebook

- **Fig 1**: `Fig1.ipynb` (`scripts/Fig1_func.r`; Fig 1i: `scripts/single_cell_func.r`)
- **Fig 2**: `Fig2.ipynb` (`scripts/rna_projection_func.r`; Fig 2e: `scripts/single_cell_func.r`)
- **Fig 3**: `Fig3.ipynb` (`scripts/lineage_peaks_func.r`, `scripts/ATAC_preprocessing.r`)
- **Fig 4** (spatial): `Fig4.ipynb` (`scripts/ATAC_preprocessing.r`, `scripts/ATAC_spatial.r`)
- **Fig 4F, 5D, S7H** (genome browser tracks): `Fig4.ipynb`, last section
- **Fig 5A**: `Fig5.ipynb` (`scripts/lineage_peaks_func.r`, `scripts/motif_screen_func.r`)
- **Fig 5B** (TEs): `Fig5.ipynb` (`scripts/ATAC_preprocessing.r`, `scripts/Fig5_func.r`)
- **Fig 5C, ED 8c, ED 8d**: `Fig5.ipynb` (`scripts/Fig5_func.r`, `scripts/lineage_peaks_func.r`, `scripts/rna_projection_func.r`; ED 8c: `scripts/single_cell_func.r`)
- **Fig 6, ED 10c-h**: `E9.5_multiome.Fig6/Fig6.ipynb` (`E9.5_multiome.Fig6/Fig6_func.r`, `E9.5_multiome.Fig6/e95_func.r`, `scripts/lineage_peaks_func.r`, `scripts/ATAC_preprocessing.r`)

### Extended Data and Supplementary tables

- **ED 2c, 2d, 3b, 3c, 4b, 5a-c, 6a** (single-cell RNA): `ExtData_RNA.ipynb` (`scripts/ExtData_RNA_func.r`; ED 5a: `scripts/single_cell_func.r`)
- **ED 4a, 6b** (single-cell ATAC): `ExtData_ATAC.ipynb` (`scripts/ATAC_preprocessing.r`)
- **Supplementary Tables 1-3**: `Supp_tables.ipynb`
- Other Extended Data panels are in the figure notebooks (ED 2e, 3d, 3f: Fig 1; ED 6d: Fig 3; ED 7: Fig 4; ED 8: Fig 5; ED 10c-h: Fig 6).
- Not in this repo: Fig 3b, ED 6c and ED 9 (FOXA2 and published ChIP-seq).

### Preprocessing (upstream of the figures)

- **ATAC import** from the multiome pipeline: `prep/atac_tracks.ipynb` (`scripts/ATAC_import.r`)
- **ATAC model** (`data/output/mod.RDS`): `prep/atac_model.ipynb` (`scripts/ATAC_preprocessing.r`)
- **GC per embryo**: `prep/embryo_profiles.ipynb` (python)
- **Transcriptional time per embryo** (for Fig S2D): `prep/transcriptional_time.ipynb` (`scripts/transcriptional_time_func.r`); runs in `dnmt_final_dir`
- **KO vs host pseudo-bulk** (`data/output/mko.RDS`, `mhost.RDS`) and **chi-square screen** (`data/output/chisq_screen.RDS`): `prep/ko_vs_host.ipynb` (`scripts/single_cell_func.r`; figure notebooks load both with `prep_single_cell()`)
- **Made outside this code** (E9.5 inputs in `revision_e95/data/`; metacell building is not part of this repo):
  - `dnmt_dko_brain_cells.h5ad`, `dnmt_dko_brain_combined_cells.h5ad`: metacells (`metacells` python package) of the Dnmt-DKO and host fractions, then concatenated.
  - `dnmt_dko_brain_metacells_gotg95_projection.h5ad` (and `..._projection_weights.csv`): projection of the combined metacells on the E9.5 atlas with `mc.pl.projection_pipeline`.
  - `annot_2909.csv`: the final metacell cell-state annotation, curated manually in MCView.
  - `data/output/cre_peaks_for_energ_passed_threshold.RDS` (input of `prep/E9.5_motif_energies.ipynb`): the 500 bp CREs whose accessibility varies across E9.5 states by > 3 (max - min of the CRE x state matrix).
- **E9.5 ATAC tracks and peaks**: `prep/E9.5_atac_tracks.ipynb`, `prep/E9.5_marginal_tracks.ipynb`, `prep/E9.5_atac_peaks.ipynb` (`scripts/ATAC_preprocessing.r`, `E9.5_multiome.Fig6/e95_func.r`)
- **E9.5 motif energies of the CREs** (`energ_cre_*.RDS`, for Fig 6): `prep/E9.5_motif_energies.ipynb`

## Data

The data are not in the repo (`data/` ~88 GB, `revision_e95/data/` ~11 GB). Set `embdnmt_dir` in `paths.R` to the original analysis directory, which holds `data/` and `revision_e95/data/`. The root notebooks read `./data/` from the repo root, so link it: `data` -> `embdnmt_dir/data`. Notebooks in `prep/` and `E9.5_multiome.Fig6/` change into `embdnmt_dir` (or `embdnmt_dir/revision_e95`) themselves.

Main inputs:
- `data/MC_models/*.h5ad`: metacell models per genotype and the WT atlas
- `data/output/*.RDS`: cached ATAC model (`mod.RDS`), motif and trend tables, screens
- `data/ref/`, `data/summary/`, `data/download/`: annotations, TSS clusters, published gene lists (Weber et al.)
- `revision_e95/data/`: E9.5 DKO brain cells, metacells, annotations

Raw and processed sequencing data: GEO (accession in the paper).

### External resources

Set in `paths.R`, except a few single-use inputs still referenced by absolute path:
- misha genome database `mm10` (Tanay lab trackdb)
- `dnmt_final_dir`: parent project with the embryo metadata, the metacell database and the transcriptional times for Fig S2D
- WT gastrulation atlas `spatwt` (iteration 10.5.5)
- `prego`, `iceqream`, `mcATAC` (tanaylab R packages) and motif models from the motif_reg project
- `embdnmt_atac.r` / `embdnmt_atac_mod.r` from the original project's shared scripts (`prep/atac_tracks.ipynb`)

## Known issues

- **Jittered points**: the chimerism plot (ED 3b), ED 2d and Fig 2b, d now use jitter with seed 1, so they are reproducible. The published panels used unseeded jitter, so their point positions differ slightly; the values do not.

## Environment

R notebooks ran on Jupyter IRkernel with R 4.3–4.4. `renv.lock` pins the R 4.4.1 environment the notebooks were last run with (Bioconductor 3.19; tanaylab packages by GitHub commit). Restore it with `renv::restore()`. Main packages:
- `misha`, `misha.ext`, `tgstat`, `tgutil`, `tglkmeans`, `metacell`, `anndata`, `prego`, `iceqream`, `mcATAC`
- `data.table`, `tidyverse`, `ggplot2`, `ggrepel`, `ggpubr`, `patchwork`, `cowplot`, `pheatmap`, `ComplexHeatmap`
- `GenomicRanges`, `Rsamtools`, `Matrix`

`mcATAC` is not on CRAN: `remotes::install_github("tanaylab/mcATAC")`.

The R `anndata` package needs a Python with `anndata` installed; point `RETICULATE_PYTHON` to it.

Running notebooks headless (e.g. with `nbclient`): R's `png()` defaults to X11 and fails without a display. Add `options(bitmapType = "cairo")` to an R profile.

The python notebook (`prep/embryo_profiles`) uses `metacells`, `anndata`, `scanpy`, `numpy`, `pandas`, `scipy`. It originally ran on an `ipyflow` kernel (Python 3.9.5) that no longer exists, so its package versions are not known.
