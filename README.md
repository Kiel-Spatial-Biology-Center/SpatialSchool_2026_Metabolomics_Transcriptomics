# Spatial Multi-Omics: MALDI Metabolomics + Xenium Transcriptomics (Mouse Colon)

A five-script pipeline integrating MALDI imaging mass spectrometry
(metabolomics) with 10x Xenium spatial transcriptomics on adjacent mouse
colon swiss roll sections, using Cardinal, Seurat, and SpaMTP.

## Running order

| `001_spatial_metabolomics_prep.R` | Raw MALDI processing, annotation, and class-sized subsetting (two parts in one file) | Cluster (needs real RAM) |
| `002_spatial_transcriptomics_prep.R` | Xenium load, QC, clustering, marker genes | Cluster or laptop |
| `003b_manual_alignment.R` | Aligns MALDI pixel coordinates into Xenium coordinate space | Cluster or laptop, **before 003** |
| `003_integration_spamtp.R` | Maps metabolomics onto Xenium cells, clusters each modality, cell-type annotation, per-cell-type metabolite signatures, gene-metabolite correlation | Cluster or laptop, **after 003b** |


`003b` must be run once before `003` — it produces
`data/processed/sm_transformed_manual_new.rds`, which `003` loads directly.

## Why the pipeline looks like this

- **`AlignSpatialOmics()` doesn't work on Xenium data.** It's written for
  Visium's image class (`@scale.factors` slot); Seurat's `FOV` class (used
  for Xenium) doesn't have that slot at all — a structural incompatibility,
  not a parameter to fix. `003b` does the same job (landmark-based affine
  registration) manually instead — a standard approach in the field, not a
  downgrade; even SpaMTP's own published Xenium+MALDI case study used an
  external tool for this step rather than `AlignSpatialOmics()`.
- **Both slides are adjacent, not identical, sections.** Even with correct
  alignment, you're comparing two different physical slices of tissue —
  expect real biological/positional noise, not pixel-perfect correspondence.
  Worth saying explicitly to students.

##  Limitations

- **Large Xenium objects (~230k+ cells) can overflow R's protection stack**
  (`protect(): protection stack overflow`) on operations that touch the full
  cell-boundary polygon data, even with plenty of RAM available. Fixed in two
  ways used throughout this pipeline: (1) launching R with
  `R --max-ppsize=500000` (or setting `R_MAX_PPSIZE=500000` in `.Renviron`)
  when possible, and (2) stripping `@images` before a `subset()` call and
  reattaching afterward, which avoids the expensive operation entirely rather
  than just giving it more stack room.

## Exercises

Each of `001`, `002`, and `003` has 1–3 blocks marked `## >>> EXERCISE N <<<`
directly in the script, at the point in the analysis they're most relevant.
They're meant to be done in place — re-running a small piece of the
preceding code with a changed parameter or a different gene/cluster, then
comparing the output — not separate assignments.

| Script | Exercise | What it explores |
|---|---|---|
| 001 | 1 | Effect of a stricter annotation-specificity cutoff on the final feature panel |
| 001 | 2 | Which lipid classes dominate this DAN-negative MALDI dataset, and why |
| 002 | 1 | Sensitivity of cell retention to QC thresholds on a small custom gene panel |
| 002 | 2 | Spatial coherence of a marker gene from a different cluster |
| 003 | 1 | Effect of clustering resolution on transcriptomics cluster count/quality |
| 003 | 2 | Critiquing and revising a specific cell-type annotation call |
| 003 | 3 | Repeating the per-cell-type metabolite signature analysis for a new cell type |

## Data

Scripts assume `data/metabolomics/` (raw MALDI `.imzML`) and
`data/xenium/` (Xenium `outs` folder, the classic flat-file bundle — **not**
the Xenium Explorer `.zarr.zip` bundle, which `LoadXenium()` cannot read)
exist under the working directory set near the top of each script
(`setwd(...)` — update this path for your own environment).
