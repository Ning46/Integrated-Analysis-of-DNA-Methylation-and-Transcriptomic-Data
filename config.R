# Run scripts from the repository root. Override settings in config.local.R
# (ignored by Git), or set PROJECT_CONFIG to another R configuration file.
config <- list(
  idat_dir = Sys.getenv("METHYLATION_IDAT_DIR", "data/GSE189148_RAW"),
  geo_cache = Sys.getenv("GEO_CACHE_DIR", "data/geo_cache"),
  results_dir = Sys.getenv("ANALYSIS_RESULTS_DIR", "results"),
  methylation_array = "EPIC",
  methylation_samples = "", # optional TSV: sample_id, group, status, individual, sex
  methylation_normalized_rds = "", # optional trusted normalized MethylSet
  excluded_cpg_file = "", # optional one-probe-ID-per-line exclusion list
  remove_sex_chromosomes = FALSE,
  n_cores = 2L,
  rna_input = "kallisto", # alternatively: "counts"
  rna_samples = "data/rna_samples.tsv",
  tx2gene_file = "data/tx2gene.tsv",
  rna_counts_file = Sys.getenv("RNA_COUNTS_FILE", "data/GSE189149_counts_table.txt"),
  rna_design = "~ sex + group",
  rna_coefficient = "groupFA", # positive logFC = FA higher than NA
  fdr = 0.05,
  min_delta_beta = 0.05,
  min_abs_rna_logfc = 0.5,
  run_enrichment = TRUE
)
local_config <- Sys.getenv("PROJECT_CONFIG", "config.local.R")
if (file.exists(local_config)) sys.source(local_config, envir = environment())
stopifnot(is.list(config), config$fdr > 0, config$fdr < 1,
          config$min_delta_beta >= 0, config$min_abs_rna_logfc >= 0)
