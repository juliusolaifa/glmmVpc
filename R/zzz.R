.onAttach <- function(libname, pkgname) {
  required_pkgs <- c("BiocManager", "BiocParallel", "DESeq2", "S4Vectors", "SummarizedExperiment")
  missing_pkgs <- required_pkgs[!sapply(required_pkgs, requireNamespace, quietly = TRUE)]

  if (length(missing_pkgs) > 0) {
    msg <- paste0(
      "Some required packages are not installed: ",
      paste(missing_pkgs, collapse = ", "),
      ". \nPlease install them using: \n",
      "if (!requireNamespace('BiocManager', quietly = TRUE)) install.packages('BiocManager')\n",
      "BiocManager::install(c('", paste(missing_pkgs, collapse = "', '"), "'))"
    )
    packageStartupMessage(msg)
  }
}
