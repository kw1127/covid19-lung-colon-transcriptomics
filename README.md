# Transcriptional profiling of fatal COVID-19 cases using lung and colon samples

A tissue-adjusted differential expression analysis of post-mortem lung and
colon from fatal COVID-19 cases and controls, exploring which parts of the
host response are shared across organs and which are organ-specific.

## What I did

Starting from bulk RNA-seq of paired lung and colon samples (n = 38), I:

- ran QC and model diagnostics (PCA, sample-distance clustering, library-size
  and dispersion checks) to justify a tissue-adjusted design
- fitted a pooled, tissue-adjusted DESeq2 model to isolate the disease effect
  common to both organs
- re-analysed each tissue separately and used a tissue-by-disease interaction
  test to classify genes as shared, lung-specific, colon-specific, or discordant
- interpreted the results with GO over-representation analysis, GSEA, and
  enrichment-map networks (clusterProfiler)

## What I found

A concordant pro-fibrotic, innate-inflammatory core (SPP1, SERPINE2, PLOD2,
S100A8) sits on top of substantially different organ-specific programmes:

- **Lung** — fibrillar-collagen deposition (COL1A1/3A1) and an antimicrobial
  immune response
- **Colon** — TGF-β and vascular remodelling, with a marked loss of gut-defence
  (DEFA6, REG3A) and adaptive-immune transcripts

A theme common to both organs was a broad down-regulation of immune-cell
identity, consistent with the immunosuppression seen in severe disease.

The full write-up, all ten figures, and a detailed limitations section are in
the **PDF report** in this repository.

## Repository contents

- `Script.R` — the full analysis pipeline
- `Transcriptomic_analysis_of_fatal_COVID_19_cases_in_lung_and_colon_samples.pdf` — report
- `figures/` — figures 1–10

## Data

Preprocessed counts and metadata from the EMBL-EBI Gene Expression Atlas,
from Wu et al. (2020), *PNAS* (DOI: 10.1073/pnas.2018030117).
