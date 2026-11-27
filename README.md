Do T cells change with age, or just their mix?

Single-cell analysis of EBV-specific CD8 T cells from 7 young and 7 older donors.
[Your name], [course / institute], [month year]

Data
- GEO accession **GSE275633** (BEAM-T single-cell RNA, TCR and surface protein data)
- 14 donors (7 young, 7 older), 3 experiments (Exp01 to Exp03), HLA-A*02:01-positive
- Source paper: Sturmlechner, Jain et al., Nat Commun 2025, doi:10.1038/s41467-025-61627-y
- Raw data are not in this repository. Download them from GEO into `data/`.

Questions and main results
1. Mix or change? Older donors have fewer stem-like/naive cells (about 31% to 10%).
   Inside a cell type, only effector memory shows significant gene changes (10 genes higher
   in older donors, BCL2-led survival score). This is exploratory (6 young vs 7 older donors).
2. Clone spread: no detectable age difference in how TCR clones spread across cell
   states, or in clonality.
3. LMP2: not analysed in depth. Only 912 cells, 8 of 14 donors with 20 or more cells,
   and 2 to 278 cells per donor.

How to run
Run the scripts in `scripts/` in numbered order from the project folder (R, Seurat).
Scripts read and write objects in `results/`, which are not uploaded because of their size,
so they are recreated by running the pipeline.

| Scripts | Purpose |
|---|---|
| 01-04 | Setup, cleaning of Exp01 to Exp03 |
| 05-08 | Merge, QC filtering, normalization, PCA checks |
| 09-12 | Antigen metadata check, doublets, clustering |
| 13-15 | Cluster composition, extra markers, cell state labels |
|16 | Harmony Check |
| 17 | Donor-level cell-state composition |
| 18-20, 23 | Q1: pseudobulk, survival score, robustness |
| 21, 22, 24 | Q2: TCR clones, spread across states, clonality |
| 25 | Q3: LMP2 compared with other antigens |


## Limits
14 donors in 3 batches, small cell states, gene-set choice, and very uneven LMP2 cell
numbers. Statistics are done at donor level, not per cell.
