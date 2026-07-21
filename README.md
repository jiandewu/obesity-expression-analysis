# TNBC Obesity Expression Analysis

R workflows for studying BMI-associated gene-expression differences using
TNBC samples and normal reference samples from GEO datasets GSE76124 and
GSE102088.

## Repository layout

- `analysis/bmi_analysis_preprocessed.R`: analysis starting from prepared
  expression tables.
- `analysis/bmi_analysis_from_geo.R`: preparation and analysis starting from
  GEO series matrices and probe-to-gene mappings.
- `environment.yml`: Conda environment for the R dependencies.

Both workflows use edgeR and limma-voom, compare lean samples with obese and
overweight samples, and produce diagnostic, MDS, volcano, and heatmap plots.

## Data policy

Downloaded expression matrices, phenotype attributes, sample-level metadata,
generated results, and manuscripts are not tracked in this repository.

Place independently downloaded inputs under `data/`:

```text
data/
|-- GSE76124/
|   |-- 00_genes_ID.txt
|   |-- 01_GSE76124_series_matrix_165_bmisorted.txt
|   `-- 03_data_exp_genes_165_bmisorted.txt
`-- GSE102088/
    |-- 00_GSE102088_series_matrix.txt
    |-- 01_genes_ID.txt
    `-- 03_filtered.csv
```

The preprocessed workflow needs the two `03_*` files. The from-GEO workflow
needs the `00_*` and `01_*` files.

## Environment

```bash
conda env create -f environment.yml
conda activate tnbc-obesity-expression
```

## Configuration

Set a different private input directory without editing the scripts:

```bash
export OBESITY_DATA_DIR=/path/to/input
```

If unset, the scripts use `data/` relative to the repository root.

## Run

```bash
Rscript analysis/bmi_analysis_preprocessed.R
Rscript analysis/bmi_analysis_from_geo.R
```

## Analysis groups

BMI values are converted into three categories:

- `Lean`: BMI up to 24.99
- `OverWeight`: BMI from 25.00 through 29.99
- `Obese`: BMI 30.00 or greater

The fitted model includes BMI category and evaluates `Lean - Obese` and
`Lean - OverWeight` contrasts.

## Reproducibility notes

- Source data are not bundled; download them through an appropriate GEO
  workflow and retain the filenames expected above.
- The from-GEO script relies on fixed metadata row positions in the source
  series matrices. Verify those positions if the downloaded representation
  changes.
- The scripts retain the original normalization, filtering, statistical
  thresholds, and plotting choices except where model inconsistencies were
  corrected during repository extraction.

## Privacy

No sample-level phenotype files, expression matrices, generated results,
manuscripts, or private paths are included in the repository or its retained
history.
