random_cov_mat <- function() {
  p <- 3   # Number of variables
  df <- 500 # Integer >= p
  S <- rWishart(n = 1, df = df, Sigma = diag(p))[, , 1] / df

  S + diag(3)*10^-6
}

random_test_dataset <- function() {
  S <- random_cov_mat()
  n <- 100000
  p <- nrow(S)
  Z <- base::matrix(stats::rnorm(n * p), nrow = n, ncol = p)
  X <- Z %*% base::chol(S)
  colnames(X) <- letters[1:3]

  return(X)
}

