# Integrated Analysis of DNA Methylation and Transcriptomic Data

R code for the data overview and initial preparation of DNA methylation microarray and RNA-seq data, extracted and adapted from **Data_Overview.html** by Yining Wang.

The current workflow covers data import, sample annotation, raw methylation preprocessing, and RNA count-table preparation. Statistical integration, differential methylation, differential expression, and enrichment analysis are not implemented in this report-derived script.

## Datasets

| Modality | GEO accession | Samples reported in the HTML |
| --- | --- | --- |
| DNA methylation microarray | [GSE189148](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE189148) | 88 |
| RNA-seq | [GSE189149](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE189149) | 48 |

These totals describe the original report. The script prints the dimensions of the data actually loaded. It does not assume that the two modalities contain identical individuals or paired samples.

## Repository contents

```text
Data_Overview.R       Executable R workflow adapted from the report
Data_Overview.html    Original rendered report
README.md            Setup, execution, outputs, and interpretation notes
```

Raw inputs and generated results are not included in the repository.

## Requirements

Install R and the following Bioconductor packages from an R console:

```r
install.packages("BiocManager")
BiocManager::install(c("GEOquery", "Biobase", "minfi", "S4Vectors"))
```

Use package versions compatible with your R/Bioconductor installation. The source HTML does not provide a complete software environment or lockfile. Internet access is required to download GEO metadata. Loading the full methylation dataset requires substantial memory; exact memory requirements have not been measured.

## Prepare the data

Download the methylation IDAT archive and the RNA-seq count table from the corresponding GEO records. Extract the IDAT files before running the script. The default layout is:

```text
data/
  GSE189148_RAW/
    <sample IDAT files>
  GSE189149_counts_table.txt
```

The count table must be tab-delimited, with a header row, gene identifiers in the first column, and numeric sample counts in subsequent columns. IDAT sample names must begin with the original ten-character GSM identifiers, as assumed in the source report.

## Run

From the repository root:

```bash
Rscript Data_Overview.R
```

Alternatively, open the project folder in RStudio and run:

```r
source("Data_Overview.R")
```

To use another input location, configure paths before sourcing:

```r
Sys.setenv(
  METHYLATION_IDAT_DIR = "D:/RStudiodata/Methylation/GSE189148_RAW",
  RNA_COUNTS_FILE = "D:/RStudiodata/RNA/GSE189149_counts_table.txt",
  DATA_OVERVIEW_OUTPUT_DIR = "results",
  GEO_CACHE_DIR = "data/geo_cache"
)
source("Data_Overview.R")
```

The script checks package availability, input paths, required metadata fields, methylation sample identifiers, and numeric count values. It stops when methylation samples cannot be matched safely to annotation.

## Workflow

1. Reproduce the study overview table from the HTML.
2. Download GSE189148 metadata and extract group, stimulation status, individual, sex, and allergy annotation.
3. Read IDAT intensity data with `minfi::read.metharray.exp()`, align metadata by sample identifier, and apply `minfi::preprocessRaw()`.
4. Define helpers for M values and Beta values using the original offsets:
   - `calculate_m_values(mSet, offset = 1)`: `log2((M + 1) / (U + 1))`.
   - `calculate_beta_values(mSet, offset = 100)`: `M / (M + U + 100)`.
5. Read the RNA count table, retain the first occurrence of each gene identifier, round numeric counts, and assign gene row names.
6. Download GSE189149 metadata and extract organism and allergy annotation.
7. Export prepared tables and R session information.

Here, `M` and `U` denote methylated and unmethylated intensities. The helper functions are defined without automatically allocating the full derived matrices. To calculate them interactively after sourcing, call `calculate_m_values(mSet)` or `calculate_beta_values(mSet)`.

## Outputs

By default, outputs are written to `results/`:

| File | Contents |
| --- | --- |
| `study_overview.csv` | Dataset accessions and sample totals from the original report |
| `methylation_sample_annotation.csv` | Methylation annotation in IDAT sample order; first column contains sample IDs |
| `rna_counts_rounded.tsv` | Gene identifiers and rounded RNA counts, retaining the first duplicate gene row |
| `rna_sample_annotation.csv` | RNA annotation with GEO sample IDs in the first column |
| `sessionInfo.txt` | R version, platform, and package session information |

GEO metadata files are cached under `data/geo_cache/`. Re-running the script overwrites the named output tables. Large methylation objects remain in memory when the script is sourced; they are not automatically saved to disk.

## Interpretation and reproducibility notes

- The `NA` group value is the literal string used by the source grouping rule, not an R missing value. Group assignment preserves the original substring match on the `indiidual:ch1` metadata field.
- Metadata field spellings and positional allergy fields are inherited from the report. Required fields are checked explicitly; inspect the GEO metadata if the schema changes.
- `preprocessRaw()` provides raw methylated/unmethylated intensities. This workflow does not perform a full quality-control, normalization, probe-filtering, or batch-correction pipeline.
- Keeping the first duplicate gene row and rounding values reproduces the report's choices. Their suitability should be assessed before downstream modeling, especially if the input represents estimated counts.
- RNA annotation is exported in GEO order. If count-column names differ from GEO IDs, the script warns; resolve sample mapping explicitly before using it in a statistical design.
- Raw inputs were not supplied with the HTML. The complete workflow has not been executed against the original data in this conversion environment.

## Changes from the rendered report

The script excludes rendered console output, adds explicit package namespaces and configurable paths, validates inputs, aligns methylation annotation by identifier, and exports reusable tables. The original M/Beta formulas and RNA preprocessing choices are retained. The original HTML remains available for reference.
