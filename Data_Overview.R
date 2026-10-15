# Data overview: DNA methylation and transcriptomic data
# Adapted from Data_Overview.html by Yining Wang.
# Run from the project root with: Rscript Data_Overview.R

required_packages <- c("GEOquery", "Biobase", "minfi", "S4Vectors")
missing_packages <- required_packages[!vapply(
  required_packages, requireNamespace, logical(1), quietly = TRUE
)]
if (length(missing_packages)) {
  stop("Missing Bioconductor packages: ", paste(missing_packages, collapse = ", "),
       ". See README.md for installation instructions.", call. = FALSE)
}

# Set these environment variables to override the default project paths.
idat_path <- Sys.getenv("METHYLATION_IDAT_DIR", "data/GSE189148_RAW")
counts_path <- Sys.getenv("RNA_COUNTS_FILE", "data/GSE189149_counts_table.txt")
output_dir <- Sys.getenv("DATA_OVERVIEW_OUTPUT_DIR", "results")
geo_cache <- Sys.getenv("GEO_CACHE_DIR", "data/geo_cache")

if (!dir.exists(idat_path)) stop("IDAT directory not found: ", idat_path)
if (!length(list.files(idat_path, pattern = "\\.idat(\\.gz)?$",
                       recursive = TRUE, ignore.case = TRUE))) {
  stop("No IDAT files found in: ", idat_path)
}
if (!file.exists(counts_path)) stop("RNA count table not found: ", counts_path)
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(geo_cache, recursive = TRUE, showWarnings = FALSE)

require_columns <- function(x, columns, context) {
  missing <- setdiff(columns, colnames(x))
  if (length(missing)) {
    stop(context, " is missing columns: ", paste(missing, collapse = ", "))
  }
}

# 1. Study overview -----------------------------------------------------------
# Sample totals below are reported in the original HTML, not inferred here.
table_general_info <- rbind(
  c("GSE189148", "GSE189149"), c(88, 48), c(88, 48)
)
colnames(table_general_info) <- c("Microarray Data", "RNAseq Data")
rownames(table_general_info) <- c(
  "Accession ID", "Number of Samples (original)", "Number of Samples (analysed)"
)
print(table_general_info)
write.csv(table_general_info, file.path(output_dir, "study_overview.csv"))

# 2. Methylation sample annotation -------------------------------------------
geo <- GEOquery::getGEO(
  GEO = "GSE189148", destdir = geo_cache, GSEMatrix = TRUE
)[[1]]
pheno <- Biobase::pData(geo)
require_columns(pheno, c("indiidual:ch1", "stimulation_status:ch1", "Sex:ch1",
                         "characteristics_ch1.5"), "GSE189148 metadata")

# 'indiidual:ch1' is the field spelling used in the original GEO-based report.
# Preserve its original NA/FA grouping rule; "NA" is a literal group label.
pheno$group <- ifelse(grepl("NA", pheno[["indiidual:ch1"]]), "NA", "FA")
pheno$status <- pheno[["stimulation_status:ch1"]]
pheno$individual <- pheno[["indiidual:ch1"]]
pheno$sex <- pheno[["Sex:ch1"]]
pheno$allergy_status <- sub(".*:\\s*", "", pheno[["characteristics_ch1.5"]])
pheno <- pheno[, c("group", "status", "individual", "sex", "allergy_status")]

# 3. Raw methylation intensities ---------------------------------------------
# Keep extended/force settings from the source report.
rgSet <- minfi::read.metharray.exp(base = idat_path, extended = TRUE, force = TRUE)
print(dim(rgSet))
print(length(Biobase::sampleNames(rgSet)))

# The source names start with ten-character GSM identifiers. Validate and
# explicitly reorder annotation rather than relying on coincidental row order.
sample_ids <- substr(Biobase::sampleNames(rgSet), 1, 10)
if (anyNA(sample_ids) || anyDuplicated(sample_ids) ||
    any(!grepl("^GSM[0-9]{7}$", sample_ids))) {
  stop("Expected unique ten-character GSM IDs at the start of IDAT sample names.")
}
if (!setequal(sample_ids, rownames(pheno))) {
  stop("IDAT sample IDs and GEO methylation annotation IDs do not match.")
}
pheno <- pheno[sample_ids, , drop = FALSE]
Biobase::sampleNames(rgSet) <- sample_ids
stopifnot(identical(Biobase::sampleNames(rgSet), rownames(pheno)))
Biobase::pData(rgSet) <- S4Vectors::DataFrame(pheno)
mSet <- minfi::preprocessRaw(rgSet)
print(class(mSet))
print(utils::head(as.data.frame(Biobase::pData(mSet))))
write.csv(pheno, file.path(output_dir, "methylation_sample_annotation.csv"))

# Original formulas, named to avoid masking minfi::getM.
# These helpers are defined but not evaluated automatically, as in the report.
calculate_m_values <- function(methyl_set, offset = 1) {
  log2((minfi::getMeth(methyl_set) + offset) /
         (minfi::getUnmeth(methyl_set) + offset))
}
calculate_beta_values <- function(methyl_set, offset = 100) {
  methylated <- minfi::getMeth(methyl_set)
  unmethylated <- minfi::getUnmeth(methyl_set)
  methylated / (methylated + unmethylated + offset)
}

# 4. RNA-seq count table ------------------------------------------------------
counts <- utils::read.table(counts_path, sep = "\t", header = TRUE,
                            check.names = FALSE, stringsAsFactors = FALSE)
if (ncol(counts) < 2L) stop("Expected a gene column and sample count columns.")
colnames(counts)[1] <- "GeneName"
if (anyNA(counts$GeneName) || any(!nzchar(as.character(counts$GeneName)))) {
  stop("Gene identifiers must be nonmissing and nonempty.")
}
# Preserve the report's rule: retain the first row for each gene name.
counts <- counts[!duplicated(counts$GeneName), , drop = FALSE]
if (!all(vapply(counts[, -1, drop = FALSE], is.numeric, logical(1)))) {
  stop("All sample count columns must be numeric.")
}
count_matrix <- as.matrix(counts[, -1, drop = FALSE])
if (any(!is.finite(count_matrix)) || any(count_matrix < 0)) {
  stop("Counts must be finite, nonnegative values.")
}
counts[, -1] <- round(counts[, -1, drop = FALSE], digits = 0)
rownames(counts) <- counts$GeneName
print(utils::head(counts))
print(dim(counts))
utils::write.table(counts, file.path(output_dir, "rna_counts_rounded.tsv"),
                   sep = "\t", quote = FALSE, row.names = FALSE)

# 5. RNA-seq sample annotation -----------------------------------------------
geo1 <- GEOquery::getGEO(
  GEO = "GSE189149", destdir = geo_cache, GSEMatrix = TRUE
)[[1]]
pheno1 <- Biobase::pData(geo1)
require_columns(pheno1, c("characteristics_ch1.4", "organism_ch1"),
                "GSE189149 metadata")
pheno1$allergy_status <- sub(".*:\\s*", "", pheno1[["characteristics_ch1.4"]])
pheno1$organism <- pheno1[["organism_ch1"]]
pheno1 <- pheno1[, c("organism", "allergy_status")]
print(utils::head(pheno1))
write.csv(pheno1, file.path(output_dir, "rna_sample_annotation.csv"))

# The source does not establish a mapping between count-column names and GEO
# IDs. Do not assume row order is a valid mapping for downstream analysis.
if (!setequal(colnames(counts)[-1], rownames(pheno1))) {
  warning("RNA count column names differ from GEO sample IDs. Establish an ",
          "explicit sample mapping before downstream analysis.", call. = FALSE)
}
capture.output(utils::sessionInfo(), file = file.path(output_dir, "sessionInfo.txt"))
message("Data overview complete. Outputs: ", normalizePath(output_dir))
