test_that("batchGLMMData returns NULL and issues a warning for num_feat = 0", {
  X <- matrix(1:3, nrow = 3, byrow = TRUE)
  beta <- c(0.5, 1)
  Sigma <- matrix(c(1, 0.5, 0.5, 1), nrow = 2)
  ns <- c(2, 1)

  expect_warning(
    result <- batchGLMMData(num = 0, X = X, beta = beta, Sigma = Sigma, ns = ns, family = "gaussian"),
    "num_feat must be positive. Returning NULL."
  )

  expect_null(result)
})


# test_that("batchGLMMData returns NULL if all iterations fail", {
#   X <- matrix(1:3, nrow = 3, byrow = TRUE)
#   beta <- c(0.5, 1)
#   Sigma <- matrix(c(1, 0.5, 0.5, 1), nrow = 2)
#   ns <- c(2,1)
#   result <- batchGLMMData(num_feat = 5, X = X, beta = beta, Sigma = Sigma, ns = ns, family = "unsupported_family")
#   expect_null(result)
# })
