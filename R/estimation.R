compareGLMMFit <- function(params, covariates, X, ns, true_family,
                           false_family, formula, iter = 1, seed = NULL) {
  if (!is.null(seed)) {
    set.seed(seed)
  }

  datafrmMat <- parallelbatchGLMMData(params, X, ns, family = true_family,
                                      iter = iter)

  ys <- datafrmMat[-1, ]
  X <- datafrmMat[1, ]
  group <- colnames(datafrmMat)

  print("Tranforming Count")
  vst_ys <- tryCatch({
    vstransform(ys, num_cores=4)
  }, error = function(e) {
    message("Error in VST transformation: ", e)
    return(NULL)
  })

  print(paste("Fitting", true_family))
  fit_true <- batchGLMMFit(formula, ys, X, group, family = true_family,
                           cov_values = covariates)

  print(paste("Fitting", false_family))
  fit_false <- batchGLMMFit(formula, ys, X, group, family = false_family,
                            cov_values = covariates)

  if (!is.null(vst_ys)) {
    print("Fitting VST")
    fit_vst <- batchGLMMFit(formula, vst_ys, X, group, family = "ga", cov_values = covariates)
    print("coef VST")
    coefs_vst <- coef(fit_vst)
  } else {
    coefs_vst <- NULL
  }

  print(paste("Coef", true_family))
  coefs_true <- coef(fit_true)
  print(paste("Coef", false_family))
  coefs_false <- coef(fit_false)

  return(list(coefs_true = coefs_true, coefs_false = coefs_false, coefs_vst = coefs_vst))
}


b0 <- c(3,7)
b1 <- c(-5,5)

sigmas <- list(
  #Equal variances
  sig1 = c(2, 1, 2),
  sig2 = c(2, -1, 2),
  sig4 = c(2,0,2),
  # # #Unequal Variances
  sig3 = c(2, 1, 4),
  sig2 = c(2, -1, 4),
  sig4 = c(2,0,4)
)

phi <- c(seq(0.01,0.4,0.02),seq(0.41,1,0.03),seq(1.1,2,0.05), seq(2.1,3,0.1),
         seq(3,5,0.2), seq(5.1,13,0.5), seq(13.5,40,5))


params <- expand.grid(b0=b0, b1=b1, sigmas = sigmas, phi=phi)
params$sig11 <- sapply(params$sigmas, `[`, 1)
params$sig12 <- sapply(params$sigmas, `[`, 2)
params$sig22 <- sapply(params$sigmas, `[`, 3)

paramsnb <- params[c("b0","b1","phi","sig11","sig12","sig22")]
covariates <- c(0,1)
ns <-  rep(10,10)
X <- rgen01(ns)
formula <- y ~ x + (1 + x | group)

execution_timenb <- system.time({
  resultnb <- compareGLMMFit(paramsnb, covariates, X, ns, true_family="nb",
                             false_family="tw", formula=formula, iter=1, seed=123)
})

paramstw <- resultnb$coefs_false[c("b0","b1","phi","sig11","sig12","sig22", "power")]

execution_timetw <- system.time({
  resulttw <- compareGLMMFit(paramstw, covariates, X, ns, true_family="tw",
                             false_family="nb", formula=formula, iter=1, seed=123)
})

