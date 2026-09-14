benchmark_covariates <- function(w, z, y, X) {
  X <- cbind(1,X)
  p <- ncol(X)

  rhos <- map(2:ncol(X), function(i) { 
    subset_x <- X[,-i]
    u <- X[,i]

    projection_matrix <- subset_x %*% solve(t(subset_x) %*% subset_x) %*% t(subset_x)
    
    residual_maker_matrix <- diag(ncol(projection_matrix)) - projection_matrix

    vars <- cbind(u,w,z,y)

    S <- cov(residual_maker_matrix %*% vars)

    pcor <- function(indices) {
        P <- solve(S[indices, indices, drop = FALSE])
        -P[1, 2] / sqrt(P[1, 1] * P[2, 2])
    }

    c(
      rho_1 = pcor(c(1, 2)),        # U ~ W
      rho_2 = pcor(c(1, 3, 2)),     # U ~ Z | W
      rho_3 = pcor(c(1, 4, 2, 3))   # U ~ Y | W, Z
    )
  })

  model <- CovarianceMatrix$new(w, z, y)

  original_tau <- model$extend(c(0,0,0))$tau()

  out <- rhos %>% map(
    function(rho) {
      bench <- model$extend(rho)
      benchmark_strength <- tibble(
        SOO_RV = max(bench$soo_z()^2, bench$soo_y()^2),
        IV_RV  = max(bench$iv_z()^2,  bench$iv_y()^2), 
        tau = bench$tau(), 
        change_in_tau = bench$tau() - original_tau, 
        rho_1 = rho[1],
        rho_2 = rho[2],
        rho_3 = rho[3]
      )

      return(benchmark_strength)
    }
  ) %>% 
    bind_rows()

  if (!is.null(colnames(X))) { 
    out$variable <- colnames(X)[-1] 
    out <- out %>% 
      relocate(variable)
  }

  return(out)
}
