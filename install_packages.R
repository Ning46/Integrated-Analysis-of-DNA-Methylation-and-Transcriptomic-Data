# Optional dependency installation; run explicitly, never during an analysis.
if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")
BiocManager::install(c(
  "GEOquery", "Biobase", "S4Vectors", "minfi", "ENmix", "limma", "edgeR",
  "tximport", "AnnotationDbi", "org.Hs.eg.db", "missMethyl", "statmod",
  "IlluminaHumanMethylationEPICmanifest", "IlluminaHumanMethylationEPICanno.ilm10b4.hg19"
), ask = FALSE, update = FALSE)
