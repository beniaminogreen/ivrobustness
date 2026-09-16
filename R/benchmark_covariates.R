#' Benchmark Robustness Values Based on Observed Covariates
#'
#' @param w the instrument vector, a numeric (double) vector with no NA values
#' @param z the treatment vector, a numeric (double) vector with no NA values
#' @param y the outcome vecctor, a numeric (double) vector with no NA values
#' @param X a matrix of control covaraites used for benchmarking
#'
#' @return a dataframe containing the robustness values associated with held out covariates
#'
#' @export
#' @importFrom dplyr %>% 
benchmark_covariates <- function(w, z, y, X) {
  p <- ncol(X)

  rhos <- purrr::map(seq(p), function(i) { 
    u <- X[,i]
    subset_x <- X[, -i, drop = FALSE]
    vars <- cbind(u = X[, i], w, z, y)

    residuals <- stats::lm.fit(
      x = subset_x,
      y = vars,
      singular.ok = TRUE
    )$residuals

    S <- stats::cov(residuals)

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

  out <- rhos %>% purrr::map(
    function(rho) {
      bench <- model$extend(rho)

      SOO_Z <- bench$soo_z()
      SOO_y <- bench$soo_y()
      IV_z <- bench$iv_z()
      IV_y <- bench$iv_y()

      SOO_RV = max(SOO_Z^2, SOO_Y^2)
      IV_RV  = max(IV_Z^2,  IV_Y^2)

      benchmark_strength <- tibble::tibble(
        SOO_RV = SOO_RV,
        IV_RV  = IV_RV,
        tau = bench$tau(), 
        change_in_tau = bench$tau() - original_tau, 
        SOO_Z = SOO_Z,
        SOO_Y = SOO_Y,
        IV_Z = IV_Z,
        IV_Y = IV_Y,
        rho_1 = rho[1],
        rho_2 = rho[2],
        rho_3 = rho[3]
      )

      return(benchmark_strength)
    }
  ) %>% 
    dplyr::bind_rows()

  if (!is.null(colnames(X))) { 
    out$variable <- colnames(X)
    out <- out %>% 
      dplyr::relocate("variable")
  }

  return(out)
}
