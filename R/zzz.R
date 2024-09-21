.onLoad <- function(libname, pkgname) {
  bioc_pkgs <- c("BiocParallel", "DESeq2", "S4Vectors", "SummarizedExperiment")
  for (pkg in bioc_pkgs) {
    if (!requireNamespace(pkg, quietly = TRUE)) {
      BiocManager::install(pkg)
    }
  }
}
