test_that("computeFixedEffects with design matrix X and beta", {
  X <- matrix(c(1, 2, 3, 4, 5), nrow = 5, ncol = 1)
  beta <- c(1, 0.5)  # Intercept and slope
  result <- computeFixedEffects(beta,X)
  expected <- c(1.5, 2, 2.5, 3, 3.5)
  expect_equal(as.vector(result), expected)
})

test_that("computeFixedEffects Intercept-only model", {
  result <- computeFixedEffects(5)
  expected <- 5
  expect_equal(result, expected)
})

test_that("computeFixedEffects adding intercept to X", {
  X <- matrix(c(2, 4, 6), nrow = 3, ncol = 1)  # No intercept column
  beta <- c(1, 0.5)  # Intercept and slope
  result <- computeFixedEffects(beta,X)
  expected <- c(1 + 2*0.5, 1 + 4*0.5, 1 + 6*0.5)
  expect_equal(as.vector(result), expected)
})

test_that("computeFixedEffects Intercept already included in X", {
  X <- matrix(c(1, 1, 1, 2, 4, 6), nrow = 3, ncol = 2)  # Already has intercept
  beta <- c(1, 0.5)  # Intercept and slope
  result <- computeFixedEffects(X=X, beta=beta, add_intercept = FALSE)
  expected <- c(1 + 2*0.5, 1 + 4*0.5, 1 + 6*0.5)
  expect_equal(as.vector(result), expected)
})

# test_that("computeFixedEffects Mismatch between X and beta dimensions", {
#   X <- matrix(c(1, 1, 1, 2, 4, 6), nrow = 3, ncol = 2)
#   beta <- c(1, 0.5, 0.2)  # Extra element
#   expect_error(computeFixedEffects(beta,X), "Number of columns in X")
# })

# test_that("computeFixedEffects Error when both X and n are NULL", {
#   expect_error(computeFixedEffects(NULL, 1), "Both X and n cannot be NULL.")
# })

test_that("computeRandomEffects Basic case with Z and U", {
  Z <- list(matrix(c(1, 1, 0, 0, 0, 1), nrow = 3, ncol = 2),
            matrix(c(0, 1, 0, 1, 1, 0), nrow = 3, ncol = 2))
  U <- matrix(c(1, 2, 0.5, 1), ncol = 2)  # Coefficients for random intercept and slope
  result <- computeRandomEffects(Z, U)
  expected <- c(2, 2.5, 2)
  expect_equal(as.vector(result), expected)
})

test_that("computeRandomEffects Error when U is not a matrix", {
  Z <- list(matrix(c(1, 0), nrow = 2))
  U <- c(1, 2)  # Not a matrix
  expect_error(computeRandomEffects(Z, U), "U must be a matrix")
})

test_that("computeRandomEffects throws error when number of matrices in Z is less than number of columns in U", {
  Z <- list(matrix(c(1, 2, 3, 4), nrow = 2))  # Only one matrix in Z
  U <- matrix(c(1, 2, 3, 4), ncol = 2)  # Two columns in U
  expect_error(computeRandomEffects(Z, U),
               "Number of random design matrices in Z can not be less than number of columns in U.")
})


test_that("computeRandomEffects Valid case with random intercept only", {
  Z <- list(matrix(c(1, 0, 0, 1), nrow = 2))  # Random intercept only
  U <- matrix(c(1, 0.5), ncol = 1)  # Only intercept
  result <- computeRandomEffects(Z, U)
  expected <- c(1, 0.5)
  expect_equal(as.vector(result), expected)
})

test_that("generateAndComputeRandomEffects works with basic input", {
  X <- matrix(1:8, nrow = 4, ncol = 2)
  ns <- c(2,2)
  Sigma <- matrix(c(3,2,1,2,2,2,1,2,4),nrow=3)

  result <- generateAndComputeRandomEffects(ns=ns, Sigma=Sigma, X=X)

  expect_length(result, 4)  # Should match the number of observations in ns
  expect_true(is.numeric(result))  # Ensure the output is numeric
})


test_that("generateAndComputeRandomEffects works with zero covariance matrix", {
  X <- matrix(c(1, 2, 3, 4), nrow = 2, ncol = 2)
  ns <- c(2)
  Sigma <- diag(3) * 0  # Zero covariance matrix

  result <- generateAndComputeRandomEffects(ns=ns, Sigma=Sigma, X=X)

  expect_equal(result, matrix(rep(0, length(ns)),nrow=ns))  # Expect zero random effects
})

