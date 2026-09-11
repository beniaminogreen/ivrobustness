test_that("Covariance Matrix Calculation Matches Results From R", {
  x <- CovarianceMatrix$new(
    iris$Sepal.Width,
    iris$Sepal.Length,
    iris$Petal.Width
  )

  diff <- x$to_r_mat() - unname(cov(iris[,c("Sepal.Width", "Sepal.Length", "Petal.Width")]))

  expect_true(max(abs(diff)) < 10^-7)
})


test_that("Extend Covariance Matrix is Self-Consistent", {
  x <- CovarianceMatrix$new(
    iris$Sepal.Width,
    iris$Sepal.Length,
    iris$Petal.Width
  )

  out_mat <- x$to_r_mat()
  extended_out_mat <- x$extend(c(0,0,0))$to_r_mat()

  expect_true(max(abs(out_mat - extended_out_mat[2:4,2:4])) < 10^6)

})


test_that("Extended Cov Mat Matches Reference Implimentation", {
  x <- CovarianceMatrix$new(
    iris$Sepal.Width,
    iris$Sepal.Length,
    iris$Petal.Width
  )

  S0 <- cov(iris[,c("Sepal.Width", "Sepal.Length", "Petal.Width")])
  for (i in 1:500) { 
    params <- runif(3)

    extended_out_mat <- x$extend(params)$to_r_mat()

    reference_extended_out_mat <- sigma_ext(S0, params)

    expect_true(
      max(abs(extended_out_mat - reference_extended_out_mat)) < 10^-6
    )
  }
})

test_that("Sensitivity Parameters are right on toy example", {
  x <- CovarianceMatrix$new(
    iris$Sepal.Width,
    iris$Sepal.Length,
    iris$Petal.Width
  )

  extended <- x$extend(c(0,0,0))
  expect_equal(extended$soo_y(), 0)
  expect_equal(extended$soo_z(), 0)
  expect_equal(extended$iv_z(), 0)
}) 

test_that("Robustness values match implimentations", {
  x <- CovarianceMatrix$new(
    iris$Sepal.Width,
    iris$Sepal.Length,
    iris$Petal.Width
  )

  S0 <- cov(iris[,c("Sepal.Width", "Sepal.Length", "Petal.Width")])

  for (i in 1:500) { 
      params <- runif(3)

      extended_out_mat <- x$extend(params)
      reference_extended_out_mat <- sigma_ext(S0, params)

      expect_equal(
        extended_out_mat$tau(),
        unname(tau(reference_extended_out_mat)), 
        tolerance = 10^-3
      )
      

      arrows_out <- arrows(reference_extended_out_mat)
      
      expect_equal(
        extended_out_mat$soo_z(),
        unname(arrows_out['soo_z']), 
        tolerance = 10^-3
      )

      expect_equal(
        extended_out_mat$soo_y(),
        unname(arrows_out['soo_y']),
        tolerance = 10^-3
      )

      expect_equal(
        extended_out_mat$iv_z(),
        unname(arrows_out['iv_z']),
        tolerance = 10^-3
      )

      expect_equal(
        extended_out_mat$iv_y(),
        unname(arrows_out['iv_y']),
        tolerance = 10^-3
      )
    }
})


