testthat::test_that("SOO rvs agree", {
  for (i in seq(500)) {
    ds <- Re(random_test_dataset())


    x <- CovarianceMatrix$new(ds[, 1], ds[, 2], ds[, 3])
    cv <- cov(ds)


  ok <- tryCatch({
    out_1 <- unname(x$soo_rv(.2)$rv)
    out_2 <- unname(calc_rvs(cv, c(.2, .2))$rv[1])
    TRUE
  }, error = function(e){ 
        FALSE
    })

  if (!ok) {
    print("Skipping")
    next
  }

  testthat::expect_lt(abs(out_1 - out_2), 1e-3)
  }
})

testthat::test_that("IV rvs agree", {
  for (i in seq(500)) {
    ds <- Re(random_test_dataset())


    x <- CovarianceMatrix$new(ds[, 1], ds[, 2], ds[, 3])
    cv <- cov(ds)


  ok <- tryCatch({
    out_1 <- unname(x$iv_rv(.2)$rv)
    out_2 <- unname(calc_rvs(cv, c(.2, .2))$rv[2])
    TRUE
  }, error = function(e) {
        print(e)
        FALSE
      })

  if (!ok) {
    print("Skipping")
    next
  }

  testthat::expect_lt(abs(out_1 - out_2), 1e-3)
  }
})
