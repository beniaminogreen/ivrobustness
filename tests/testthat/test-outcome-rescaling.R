testthat::test_that("SOO RV and bias accuracy are invariant to outcome units", {
  # Orthogonal, centered columns give SD(Z | W) = 1 and SD(Y | W,Z) = 1.
  E <- sqrt(3) / 2 * rbind(
    c(1, 1, 1), c(1, -1, -1), c(-1, 1, -1), c(-1, -1, 1)
  )
  w <- E[, 1]
  z <- 0.5 * w + E[, 2]
  y <- z + E[, 3]
  expected_rv <- (sqrt(5) - 1) / 2

  for (scale in c(1, 1e-8, 1e8)) {
    cm <- CovarianceMatrix$new(w, z, scale * y)
    result <- cm$soo_rv(scale)

    testthat::expect_true(is.finite(result$rv))
    testthat::expect_lt(abs(result$rv - expected_rv), 1e-4)
    testthat::expect_lte(abs(result$bias / scale - 1), 1.1e-5)
  }
})
