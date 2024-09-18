#' Variance Stabilizing Transformation
#'
#' @param counts A matrix of integer count data.
#' @param blind Logical, whether to ignore the experimental design in the transformation.
#' @param num_cores Number of workers
#'
#' @return A matrix of transformed values.
#' @export
#'
#' @examples
#' ns <- rep(10, 10)
#' X <- rgen01(ns)
#' beta <- c(3, 5)
#' Sigma <- matrix(c(2, 1, 1, 2), nrow=2)
#' family <- "nb"
#' cov_values <- c(0, 1)
#' formula <- y ~ x + (1 | group)
#' datafrmMat <- batchGLMMData(num=5, X, beta=beta, Sigma=Sigma, ns=ns,
#' family = "negative_binomial", theta = 1.5)
#' ys <- datafrmMat[-1,]; X <- datafrmMat[1,]; group <- colnames(datafrmMat)
#' vstransform(ys, num_cores=4)

vstransform <- function(counts, blind = FALSE, num_cores=NULL) {
  counts <- counts[!apply(counts, 1, function(row) any(row >= .Machine$integer.max)), ]
  if (!is.matrix(counts) || any(is.na(counts))) {
    stop("The 'counts' parameter must be a numeric matrix without missing values.")
  }
  if (!all(counts == floor(counts))) {
    message("Rounding non-integer counts to nearest integers.")
    counts <- round(counts)
  }

  if (is.null(colnames(counts))) {
    stop("Error: 'counts' must have column names representing the sample groups.")
  }
  chk <- Sys.getenv("_R_CHECK_LIMIT_CORES_", "")
  if (nzchar(chk) && chk == "TRUE") {
    num_cores <- 2L
  } else {
    num_cores <- ifelse(!is.null(num_cores), num_cores,parallel::detectCores() - 1)
  }

  os <- Sys.info()["sysname"]
  if (os == "Windows") {
    param <- BiocParallel::SnowParam(workers = num_cores)
  } else {
    param <- BiocParallel::MulticoreParam(workers = num_cores)
  }
  if (is.null(colnames(counts))) {
    stop("The 'counts' matrix must have column names indicating the sample groups.")
  }

  condition <- factor(colnames(counts))
  dds <- DESeq2::DESeqDataSetFromMatrix(countData = counts,
                                        colData = S4Vectors::DataFrame(condition),
                                        design = ~ condition)
  dds <- tryCatch({
    DESeq2::estimateSizeFactors(dds)
  }, error = function(e) {
    if (grepl("cannot compute log geometric means", e$message)) {
      message("Error in size factor estimation: Switching to 'poscounts' method")
      DESeq2::estimateSizeFactors(dds, type = "poscounts")
    } else {
      stop(e)
    }
  })
  dds <- DESeq2::DESeq(dds, parallel = TRUE, BPPARAM = param)
  vst_data <- DESeq2::varianceStabilizingTransformation(dds, blind = blind)
  vst_data_s3 <- SummarizedExperiment::assay(vst_data)
  BiocParallel::bpstop(param)
  return(vst_data_s3)
}
