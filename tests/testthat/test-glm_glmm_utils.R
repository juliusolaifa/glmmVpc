# Test: Identity link without random effects (behaves like a GLM)
test_that("Identity link without random effects", {
  X <- matrix(1:10, nrow=5)
  beta <- c(1, 0.5, 2)  # Intercept and slope
  suppressWarnings(
    mu <- computeMuGLMM(beta, Sigma = NULL, ns = NULL, X, link = "identity")
  )
  X_with_intercept <- cbind(1,X)
  expected <- X_with_intercept %*% beta
  expect_equal(as.vector(mu), as.vector(expected))
})

# Test: Log link with random effects
test_that("Log link with random effects", {
  X <- matrix(1:3, nrow=3)
  beta <- c(0.5, 0.3)  # Fixed effects
  Sigma <- matrix(c(1, 0.5, 0.5, 1), nrow=2)  # Random effects covariance matrix
  ns <- c(3)
  mu <- computeMuGLMM(beta, Sigma, ns, X, link = "log")
  # Manually compute expected value for fixed + random effects
  fixed_effect <- X %*% beta
  expect_length(mu, 3)  # Length of response vector should match ns
  expect_true(all(mu > 0))  # Log link should ensure positive responses
})
#
# Test: Logit link with random effects
test_that("Logit link with random effects", {
  X <- matrix(1:3, nrow=3)
  beta <- c(0.2, 0.8)  # Fixed effects
  Sigma <- matrix(c(1, 0.2, 0.2, 1), nrow=2)  # Covariance matrix for random effects
  ns <- c(3)
  mu <- computeMuGLMM(beta, Sigma, ns,X, link = "logit")
  expect_length(mu, 3)
  expect_true(all(mu >= 0 & mu <= 1))  # Logit link should produce values in [0,1]
})
#
#Test: Inverse link with zero random effects
test_that("Inverse link with zero random effects", {
  X <- matrix(1:4, nrow=4)
  beta <- c(1, 0.5)
  Sigma <- diag(2)
  ns <- c(4)
  mu <- computeMuGLMM(beta, Sigma, ns, X, link = "inverse")
  expect_length(mu, 4)
  # expect_warning(computeMuGLMM(X, beta, Sigma, ns, link = "inverse"),
  #                "Zero values in 'eta' will result in Inf")
})
#
# Test: Unsupported link function
test_that("Unsupported link function", {
  X <- matrix(1:2, nrow=2)
  beta <- c(1, 0.5)
  Sigma <- diag(2)
  ns <- c(2)
  expect_error(computeMuGLMM(beta, Sigma, ns, X, link = "unsupported"),
               "Unsupported link function")
})

test_that("compute_mu_sig works with intercept-only (x = NULL)", {
  beta <- c(1.5)  # Intercept and slope
  Sigma <- 1  # Scalar variance

  result <- compute_mu_sig(beta, Sigma, x = NULL)

  expect_equal(result$mu, 1.5)  # Only intercept is used
  expect_equal(result$sig, 1)  # Sigma is scalar
})

test_that("compute_mu_sig works with x as a vector and Sigma as scalar", {
  beta <- c(1.5, 0.5)  # Intercept and slope
  Sigma <- 2  # Scalar variance
  x <- 2  # One covariate value
  result <- compute_mu_sig(beta, Sigma, x = x)
  expect_equal(result$mu, 1.5 + 0.5 * 2)  # mu = intercept + slope * x
  expect_equal(result$sig, 2)  # Sigma is scalar
})

test_that("compute_mu_sig works with x as vector and Sigma as matrix", {
  beta <- c(1.5, 0.5)  # Intercept and slope
  Sigma <- matrix(c(1, 0.5, 0.5, 2), nrow = 2)  # 2x2 covariance matrix
  x <- 2  # One covariate value
  result <- compute_mu_sig(beta, Sigma, x = x)
  expect_equal(result$mu, 1.5 + 0.5 * 2)  # mu = intercept + slope * x
  Z <- c(1, 2)  # Z = [1, x]
  expected_sigma <- as.numeric(Z %*% Sigma %*% Z)
  expect_equal(result$sig, expected_sigma)  # Computed sigma
})

test_that("compute_mu_sig throws error when length of X does not match beta", {
  beta <- c(1.5, 0.5, 0.2)  # Intercept and two slopes
  Sigma <- 1  # Scalar variance
  x <- c(2)  # Only one covariate value (should be two)
  expect_error(compute_mu_sig(beta, Sigma, x = x), "Length of X must match the length of beta.")
})

# test_that("compute_mu_sig throws error for non-numeric input", {
#   beta <- c(1.5, 0.5)
#   Sigma <- 1  # Scalar variance
#   x <- "non-numeric"  # Invalid input for x
#   expect_error(compute_mu_sig(beta, Sigma, x = x), "non-numeric argument")
# })


