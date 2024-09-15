library(testthat)

# Test: Identity link
test_that("Identity link returns eta unchanged", {
  eta <- c(1, 2, 3)
  result <- computeMuFromEta(eta, "identity")
  expect_equal(result, eta)
})

# Test: Log link
test_that("Log link returns exp(eta)", {
  eta <- c(0, 1, 2)
  result <- computeMuFromEta(eta, "log")
  expected <- exp(eta)
  expect_equal(result, expected)
})

# Test: Logit link
test_that("Logit link returns 1 / (1 + exp(-eta))", {
  eta <- c(0, 1, -1)
  result <- computeMuFromEta(eta, "logit")
  expected <- 1 / (1 + exp(-eta))
  expect_equal(result, expected)
})

# Test: Inverse link
test_that("Inverse link returns 1 / eta", {
  eta <- c(1, 2, 0.5)
  result <- computeMuFromEta(eta, "inverse")
  expected <- 1 / eta
  expect_equal(result, expected)
})

# Test: Inverse link warns about zero values
test_that("Inverse link warns about zero values", {
  eta <- c(0, 1, 2)
  expect_error(result <- computeMuFromEta(eta, "inverse"),
                 "Zero values in 'eta' are not allowed with the 'inverse' link function.")
  # expect_equal(result, c(Inf, 1, 0.5))
})

# Test: Unsupported link function throws error
test_that("Unsupported link throws an error", {
  eta <- c(1, 2, 3)
  expect_error(computeMuFromEta(eta, "unsupported_link"), "Unsupported link function: unsupported_link")
})

# Test: Non-numeric eta throws error
test_that("Non-numeric eta throws an error", {
  eta <- c("a", "b", "c")
  expect_error(computeMuFromEta(eta, "identity"), "Input 'eta' must be a numeric vector.")
})

# Test: NA values in eta issue warning
test_that("Warning for NA values in eta", {
  eta <- c(1, NA, 3)
  expect_warning(result <- computeMuFromEta(eta, "identity"),
                 "NA values detected in 'eta'. Results will contain NAs.")
  expect_equal(result, eta)
})

