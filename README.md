# Integrated Analysis of DNA Methylation and Transcriptomic Data

Reproducible R workflows for DNA methylation, RNA-seq, and gene-level result integration in a study of food allergy.....

**Scope:** compare activation responses in food-allergic (FA) and non-food-allergic (NA) participants, analyze RNA expression differences, and identify overlapping genes. This is result-level integration, not a matched multi-omics correlation or a causal model.

## Study data

| Modality | GEO accession | Samples in the original overview |
| --- | --- | --- |
| DNA methylation | [GSE189148](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE189148) | 88 |
| RNA-seq | [GSE189149](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE189149) | 48 |

The original RNA follow-up notebook selected 36 samples using positional group/sex vectors. This repository does **not** guess which 36 samples those were: supply an explicit sample sheet and document the cohort selection. RNA and methylation samples are not assumed to be paired across assays.

Only source code, configuration defaults, documentation, and synthetic test code are versioned. Input data, caches, generated plots, tables, saved R objects, gene lists, and session history are excluded.

## Repository layout

```text
Data_Overview.R         Optional initial data inspection and annotation
all_for_DNA.R           Methylation QC, normalization, paired models, enrichment
all_for_rnaseq.R        RNA import, filtering, TMM, edgeR, annotation, enrichment
integrate.R            Gene overlap for FA activation and interaction contrasts
config.R               Shared settings and thresholds
R/helpers.R            Validation, identifier matching, annotation, effect helpers
install_packages.R     Explicit dependency installation
tests/test_workflow.R  Offline synthetic regression and workflow checks
```

The three follow-up R scripts replace the executable content of the supplied R Markdown files; they do not require knitting, LaTeX, or hidden objects from an earlier R session.

## Install

Use a compatible R/Bioconductor installation:

```bash
Rscript install_packages.R
```

This installs GEOquery, minfi, ENmix, limma, edgeR, tximport, missMethyl, human annotation packages, and supporting dependencies. Analysis scripts never install packages automatically. The default array configuration is EPIC with hg19 annotation; verify this against your input platform. Other platforms require matching minfi annotation/manifest packages and a supported missMethyl array setting.

No exact environment lockfile was supplied. Each completed stage records its configuration and session information. GEO metadata and KEGG enrichment may require network access.

## Prepare inputs

Run commands from the repository root. Keep private/local inputs under `data/`:

```text
data/
  GSE189148_RAW/                 Extracted methylation IDAT files
  GSE189149_counts_table.txt     Optional gene-level RNA count matrix
  rna_samples.tsv               Curated RNA sample metadata
  tx2gene.tsv                   Required for Kallisto import
  kallisto_quant/...            Kallisto abundance TSVs
```

### RNA sample sheet

Create a tab-delimited file with these columns:

| Column | Meaning |
| --- | --- |
| `sample_id` | Unique sample identifier; must match count-table columns in counts mode |
| `group` | Literal `NA` or `FA`; `NA` is a group label, not a missing value |
| `sex` | Verified sex annotation with consistent category labels |
| `individual` | Verified participant identifier |
| `quant_file` | Path to that sample's Kallisto abundance TSV; required only in Kallisto mode |

File paths are relative to the repository root unless absolute paths are supplied. Do not copy a group/sex vector based on directory order. Sample sheets may contain additional covariates referenced by the configured design.

The default RNA design is `~ sex + group`, with `NA` as the reference. Repeated RNA individuals are rejected; choose an appropriate independent cohort or implement a model for within-person dependence. The script also rejects confounded designs and missing covariates.

### Transcript-to-gene mapping

Kallisto mode requires a TSV with `TXNAME` and `GENEID` columns. Use the annotation release that matches the transcriptome used to build the Kallisto index. Transcript bars and trailing numeric versions are stripped consistently; ambiguous transcript-to-gene mappings are rejected.

### Alternative gene-count input

Set `config$rna_input <- "counts"` to use a TSV whose first column contains unique Ensembl gene IDs and remaining columns contain nonnegative counts with sample IDs as column headers. This mode assumes the input is appropriate for count-based modeling. TPM/FPKM values are not valid substitutes. For transcript-derived estimates with sample-dependent gene-length effects, use Kallisto mode or prepare an appropriate offset-aware workflow rather than treating uncorrected estimates as ordinary counts.

Duplicate gene identifiers are rejected in the follow-up RNA analysis rather than silently discarding rows. The earlier overview script retains its original first-duplicate/rounding behavior for inspection; its rounded export is not the default modeling input.

### Methylation metadata

By default, GSE189148 metadata is downloaded from GEO, retaining the original `indiidual:ch1`, stimulation, and sex fields. Group assignment preserves the original individual-name substring rule; verify it against the study records. Override with `config$methylation_samples` pointing to a TSV containing `sample_id`, `group`, `status`, `individual`, and optionally `sex`. Status must be `U` (unstimulated) or `S` (stimulated).

IDAT basenames must start with the ten-character GSM identifiers used in the source data. All identifiers are checked and annotation is explicitly reordered. After QC, each included participant must have exactly one U/S pair within one group. Incomplete pairs stop the analysis for an explicit cohort decision; they are not silently treated as independent samples.

## Configuration

Create `config.local.R` (ignored by Git) containing only overrides, for example:

```r
config$idat_dir <- "data/GSE189148_RAW"
config$rna_samples <- "data/rna_samples.tsv"
config$rna_input <- "kallisto"
config$tx2gene_file <- "data/tx2gene.tsv"
config$results_dir <- "results"
config$fdr <- 0.05
config$min_delta_beta <- 0.05
config$min_abs_rna_logfc <- 0.5
```

Use `PROJECT_CONFIG` to select another override file. Supported path environment variables are `METHYLATION_IDAT_DIR`, `RNA_COUNTS_FILE`, `GEO_CACHE_DIR`, and `ANALYSIS_RESULTS_DIR`. Local overrides take precedence over defaults.

Optional methylation settings include a trusted normalized MethylSet RDS, a one-CpG-ID-per-line probe exclusion list, sex-chromosome filtering, and ENmix core count. The normalized-RDS route requires matching metadata for exactly its included samples and bypasses raw QC/normalization; the object must already have undergone appropriate preprocessing. Rebuild saved objects when inputs or preprocessing change.

The initial `Data_Overview.R` remains an independent inspection step using its documented environment variables, including `DATA_OVERVIEW_OUTPUT_DIR`; it does not load the follow-up configuration or save prerequisite analysis objects.

## Execution order

```bash
# Optional original overview; requires both original modalities.
Rscript Data_Overview.R

# Follow-up analysis stages; each can run in a fresh R session.
Rscript all_for_DNA.R
Rscript all_for_rnaseq.R
Rscript integrate.R
```

The overview is not a required computational dependency: follow-up stages read configured inputs directly. Integration reads the full tested result tables generated by the two analysis stages. Re-running a stage overwrites its named outputs; use separate results directories for different settings. Set `config$run_enrichment <- FALSE` for a run without GO/KEGG requests.

## Analysis logic

### Methylation

1. Read IDATs and align GEO or user-supplied metadata by sample identifier.
2. Generate raw QC plots and run ENmix QC.
3. Remove ENmix quality failures and nonfinite raw probes; apply ENmix OOB background correction, RELIC dye correction, and quantile2 normalization.
4. Check complete U/S pairs, finite normalized values, and design rank. Optional user-provided probe exclusions and sex-chromosome filtering are applied before testing.
5. Model M values using participant fixed effects, activation, and an FA-by-activation term. Participant effects absorb time-invariant subject characteristics; the group main effect is not separately identifiable in this paired design.
6. Test FA activation, NA activation, and the direct interaction `(FA_S - FA_U) - (NA_S - NA_U)`.
7. Calculate Beta differences and attach them to sorted model results **by CpG ID**.
8. Apply BH-adjusted P < 0.05 and absolute Beta difference > 0.05 by default. For the interaction, the effect threshold uses the difference of Beta changes.
9. Expand CpG-to-gene annotations without treating island annotations as gene-aligned lists; retain the original probe identifiers.
10. Run directional GO/KEGG enrichment for FA and interaction CpGs using missMethyl, with **all tested CpGs** as the background and probe-number/multi-gene bias adjustment.

The FA-significant minus NA-significant gene list remains an explicitly exploratory output. A significant result in one group and a nonsignificant result in another does not establish a between-group difference; use the direct interaction test.

QC should still be reviewed by the analyst. Cross-reactive/SNP probe exclusion lists, batch effects, and platform-specific filters require study-specific decisions and are not invented by this conversion.

### RNA-seq

1. Validate the curated sample sheet and construct an identifiable design.
2. Import Kallisto estimates with tximport `lengthScaledTPM`, avoiding an additional length offset, or read an appropriate gene-count matrix.
3. Filter low expression using the actual design, apply TMM normalization, estimate dispersion, and fit a robust edgeR quasi-likelihood model.
4. Test `groupFA` by default; positive logFC means higher expression in FA after adjustment for sex.
5. Map Ensembl IDs to gene symbols and Entrez IDs using the installed human annotation package. Version suffixes are removed for lookup; ambiguous mappings use the package's first mapping and should be reviewed.
6. Select genes using the same FDR and absolute logFC thresholds in both directions.
7. Run GO/KEGG enrichment on unique tested Entrez genes. Where multiple tested Ensembl IDs map to one Entrez ID, select the most highly expressed representative before selection; mean log-expression is the aligned detectability covariate. Correct enrichment P values by BH within each collection/direction.

The old RNA FDR thresholds of 0.4/0.45 and asymmetric effect cutoffs were replaced by configurable, symmetric defaults (FDR 0.05; absolute logFC > 0.5). Corrected results may differ substantially from the original exported tables.

### Integration

For both FA activation and the response interaction, construct the shared universe from **all tested, annotated genes** present in both modalities. Reapply the configured thresholds, deduplicate gene sets, export overlap counts, and join detailed results explicitly by gene symbol.

CpG-to-gene and Ensembl-to-symbol mappings may be many-to-many; detail tables preserve those rows. Gene summaries average the selected CpG effects and flag mixed directions. A gene with at least one selected CpG is not a gene-level FDR-controlled methylation test.

Methylation measures within-person activation or a difference in activation responses; RNA measures a between-group expression contrast by default. Overlap is exploratory biological prioritization, not evidence of matched molecular correlation, methylation-driven expression, or promoter-specific repression. The optional Fisher overlap P value is descriptive and unadjusted; dependent genes and annotation can invalidate simple independence assumptions.

## Local outputs

| Directory | Main outputs |
| --- | --- |
| `results/` | Optional overview tables |
| `results/methylation/` | Full FA/NA/interaction CpG tables, annotated gene tables, selected gene summaries, exploratory FA-only list, QC/diagnostic PDFs, enrichment tables, retained sample/design tables |
| `results/rnaseq/` | Full edgeR table, sample/design tables, QC and top-gene plots, directional enrichment tables |
| `results/integration/` | Shared universes, overlap gene lists/details, gene summaries, counts plots, descriptive overlap summaries |

Each follow-up stage records `run_config.R` and `sessionInfo.txt`. Empty selected sets produce header-only overlap tables and skip enrichment where appropriate; they are not evidence of a successful replication of previous biological claims.

## Validation and limitations

Run offline checks from the project root:

```bash
Rscript --vanilla tests/test_workflow.R
```

Tests exercise metadata reordering, missing/incomplete samples, paired interaction direction, CpG effect alignment, multi-gene annotation, common-universe overlap, and synthetic RNA, normalized-methylation, and integration workflows. These checks use synthetic inputs and temporary outputs, not study data. The suite passed locally with R 4.4.1 and the installed Bioconductor packages; raw IDAT preprocessing and remote enrichment are outside this offline test.

The original archive did not include raw assay inputs, the curated RNA sample sheet, the Kallisto reference mapping, or all intermediate methylation objects. Full real-data preprocessing, remote enrichment, and biological conclusions have not been reproduced. Review QC and cohort definitions before interpreting results.

## Method references

- [ENmix documentation](https://bioconductor.org/packages/release/bioc/html/ENmix.html)
- [tximport workflow](https://bioconductor.org/packages/release/bioc/vignettes/tximport/inst/doc/tximport.html)
- [limma documentation](https://bioconductor.org/packages/release/bioc/html/limma.html)
- [missMethyl workflow](https://bioconductor.org/packages/release/bioc/vignettes/missMethyl/inst/doc/missMethyl.html)
