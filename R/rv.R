#' Calculate Robustness Values
#'
#' @param w the instrument vector, a numeric (double) vector with no NA values
#' @param z the treatment vector, a numeric (double) vector with no NA values
#' @param y the outcome vecctor, a numeric (double) vector with no NA values
#'
#' @return a dataframe containing robustness values with the following columns: 
#'  - *Method* The method for the associated Robustness Value: 
#'     - SOO - Selection on observables
#'     - IV - Instrumental variables
#'     - ExogenousIV - Instrumental variables, where exogeneity is known to be satisfied (as in a randomized experiment)
#'     - ExcludableIV - Instrumental variables, where the exclusion restriction is known to be satisfied
#'  - *RV* the robustness value 
#'  - *rho_1, rho_2, rho_3* The parameters used to achieve this robustness value 
#'  - *bias* The bias achieved at this combination of parameters
#'  - *assumption gap* the largest gap between partial correlations that are assumed to be 
#'  zero (only relevant for ExogenousIV and ExcludableIV) and the actual value found by the optimizer. 
#'  This should be negligible and close to the machine tolerance. 
#'
#' @export
robustness_values <- function(w, z, y) {
  stopifnot(
    "Instrument and treatment vectors must have same length" = length(w) == length(z)
  )

  stopifnot(
    "Treatment and outcome vectors must have same length" = length(z) == length(y)
  )

  stopifnot("Instrument must be a double vector" = is.double(w))
  stopifnot("Treatment must be a double vector" = is.double(z))
  stopifnot("Outcome must be a double vector" = is.double(y))

  cmat <- CovarianceMatrix$new(
    w, 
    z, 
    y
  )

  iv_estimate <- cmat$iv_estimate()
  soo_estimate <- cmat$soo_estimate()

  soo_rv <- cmat$soo_rv(soo_estimate)
  iv_rv <- cmat$iv_rv(iv_estimate)
  exog_iv_rv <- cmat$exogenous_iv_rv(iv_estimate)
  exclud_iv_rv <- cmat$excludable_iv_rv(iv_estimate)

  iv_df <- bind_rows(iv_rv, exog_iv_rv, exclud_iv_rv)
  iv_df$estimate <- iv_estimate

  soo_df <- as_tibble(soo_rv)
  soo_df$estimate <- soo_estimate 

  return(bind_rows(soo_df, iv_df))
}
