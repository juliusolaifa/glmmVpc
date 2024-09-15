test_that("generateRandomInterceptMatrix with two groups", {
  ns <- c(3, 3)
  Z <- generateRandomInterceptMatrix(ns)
  expected_Z <- matrix(c(1,1,1,0,0,0,0,0,0,1,1,1), nrow = 6)
  expect_equal(Z, expected_Z)
})

test_that("generateRandomInterceptMatrix Single group", {
  ns <- c(4)
  Z <- generateRandomInterceptMatrix(ns)
  expected_Z <- matrix(1, nrow = 4, ncol = 1)
  expect_equal(Z, expected_Z)
})

test_that("generateRandomInterceptMatrix Multiple groups with different sizes", {
  ns <- c(2, 4, 3)
  Z <- generateRandomInterceptMatrix(ns)
  expected_Z <- matrix(c(1,1,0,0,0,0,0,0,0,
                         0,0,1,1,1,1,0,0,0,
                         0,0,0,0,0,0,1,1,1), nrow = 9)
  expect_equal(Z, expected_Z)
})

test_that("generateRandomDesignMatrices Intercept only", {
  ns <- c(3, 3)
  Z_list <- generateRandomDesignMatrices(ns)
  expected_Z0 <- matrix(c(1,1,1,0,0,0,0,0,0,1,1,1), nrow = 6)
  expect_equal(Z_list[[1]], expected_Z0)
  expect_equal(length(Z_list), 1)
})

test_that("generateRandomDesignMatrices Intercept and slope", {
  ns <- c(2, 3, 2)
  X <- matrix(c(1,1,1,1,1,1,1,
                1,2,3,4,5,6,7), nrow = 7)
  Z_list <- generateRandomDesignMatrices(ns=ns,X=X)

  expected_Z0 <- generateRandomInterceptMatrix(ns)
  expected_Z1 <- expected_Z0 * c(1,2,3,4,5,6,7)

  expect_equal(Z_list[[1]], expected_Z0)
  expect_equal(Z_list[[2]], expected_Z1)
})

test_that("No intercept in X", {
  ns <- c(2, 2)
  X <- matrix(c(5,6,7,8), nrow = 4)
  Z_list <- generateRandomDesignMatrices(ns=ns,X=X)

  expected_Z0 <- generateRandomInterceptMatrix(ns)
  expected_Z1 <- expected_Z0 * c(5,6,7,8)

  expect_equal(Z_list[[1]], expected_Z0)
  expect_equal(Z_list[[2]], expected_Z1)
})

test_that("generateRandomDesignMatrices mismatch between X and ns", {
  ns <- c(3, 2)
  X <- matrix(c(1,1,1,2,2,2,3), nrow = 7)
  expect_error(generateRandomDesignMatrices(ns=ns,X=X),
               paste("The number of rows in X must match the total number",
                    "of observations in ns.")
  )
})

test_that("generateRandomDesignMatrices Intercept-only matrix X", {
  ns <- c(2, 2)
  X <- matrix(rep(1, 4), nrow = 4)
  Z_list <- generateRandomDesignMatrices(ns=ns,X=X)

  expected_Z0 <- generateRandomInterceptMatrix(ns)

  expect_equal(Z_list[[1]], expected_Z0)
  expect_equal(length(Z_list), 1)  # No slope matrices generated
})
