#' Calculate Robustness Values
#'
#' @param w the instrument vector, a numeric (double) vector with no NA values
#' @param z the treatment vector, a numeric (double) vector with no NA values
#' @param y the outcome vecctor, a numeric (double) vector with no NA values
#' @param X an optional matrix of covariates to residualize on
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
#' @examples 
#'create_simulation <- function(n=10^7, rho = .5) {
#'  w <- rnorm(n)
#'  z <- sqrt(rho)*w + sqrt(1-rho)*rnorm(n)
#'  y <- 2*z + rnorm(n)
#'
#'  return(list(
#'    w = w, 
#'    z = z, 
#'    y = y
#'  ))
#'}
#'
#'
#'strong_instrument_data <- create_simulation(rho = .90)
#'weak_instrument_data <- create_simulation(rho = .01)
#'
#'robustness_values(
#'  strong_instrument_data$w, 
#'  strong_instrument_data$z, 
#'  strong_instrument_data$y
#')
#'
#'robustness_values(
#'  weak_instrument_data$w, 
#'  weak_instrument_data$z, 
#'  weak_instrument_data$y
#')
robustness_values <- function(w, z, y, X = NULL) {
  stopifnot(
    "Instrument and treatment vectors must have same length" = length(w) == length(z)
  )

  stopifnot(
    "Treatment and outcome vectors must have same length" = length(z) == length(y)
  )

  stopifnot("Instrument must be a double vector" = is.double(w))
  stopifnot("Treatment must be a double vector" = is.double(z))
  stopifnot("Outcome must be a double vector" = is.double(y))


  if (!is.null(X)) {
    vars <- cbind(w, z, y)
    residuals <- stats::lm.fit(
        x = X,
        y = vars,
        singular.ok = TRUE
      )$residuals

    w <- residuals[,1]
    z <- residuals[,2]
    y <- residuals[,3]
  }

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

  iv_df <- dplyr::bind_rows(iv_rv, exog_iv_rv, exclud_iv_rv)
  iv_df$estimate <- iv_estimate

  soo_df <- tibble::as_tibble(soo_rv)
  soo_df$estimate <- soo_estimate 

  return(dplyr::bind_rows(soo_df, iv_df))
}


#' Interpret a Table Of Robustness Values
#' 
#'
#' Struggling to interpret the results of `robustness_values`? Look no further! 
#' `interpret_table` helps by providing a plain-text interpretation of the robustness 
#' values of a design. It interprets each of the four robustness values for you, 
#' and flags low robustness values. 
#'
#' @param df a datasdrame, as returned by `robustness_values`
#'
#' @return A string 
#'
#' @export
#' @examples
#'create_simulation <- function(n=10^7, rho = .5) {
#'  w <- rnorm(n)
#'  z <- sqrt(rho)*w + sqrt(1-rho)*rnorm(n)
#'  y <- 2*z + rnorm(n)
#'
#'  return(list(
#'    w = w, 
#'    z = z, 
#'    y = y
#'  ))
#'}
#'
#'
#'weak_instrument_data <- create_simulation(rho = .01)
#'
#'df <- robustness_values(
#'  weak_instrument_data$w, 
#'  weak_instrument_data$z, 
#'  weak_instrument_data$y
#')
#'cat(interpret_table(df))
interpret_table <- function(df) {
  paragraphs <- purrr::map(seq_len(nrow(df)), function(i) {
    rv <- signif(df$rv[i], 3)
    bias <- signif(df$bias[i], 3)
    pct <- signif(100 * df$rv[i], 3)

    warn <- (df$rv[i] < .05)

    out <- switch(tolower(as.character(df$method[i])),
      soo = stringr::str_glue(
        "The robustness value for selection on observables is {rv}. ",
        "A bias of {bias} requires an unobserved confounder to explain ",
        "at least {pct}% of residual variation in at least one of ",
        "the treatment or the outcome."
      ),
      iv = stringr::str_glue(
        "The robustness value for instrumental variables is {rv}. ",
        "A bias of {bias} requires either the confounder to explain ",
        "at least {pct}% of residual instrument variation, or the ",
        "instrument to explain at least {pct}% of outcome variation ",
        "remaining after adjusting for the treatment and confounder."
      ),
      exogenousiv = stringr::str_glue(
        "With instrument exogeneity imposed, the robustness value is {rv}. ",
        "A bias of {bias} requires the instrument to explain at least ",
        "{pct}% of outcome variation remaining after adjusting for ",
        "the treatment and confounder."
      ),
      excludableiv = stringr::str_glue(
        "With exclusion imposed, the robustness value is {rv}. ",
        "A bias of {bias} requires the confounder to explain at least ",
        "{pct}% of residual instrument variation."
      ),
      NULL
    )

    if (warn) {
      return(paste(out, 
        "This robustness value is below 5%. While there is no 'hard limit' on these values, a value this small means that relatively small violations of the corresponding identifying assumptions could drive the point estimate to zero."))
    } else {
      return(out)
    }
  })

  context <- ("Consider reading the associated paper for more information on interpreting robustness values (https://arxiv.org/pdf/2507.23743)")

  paste(c(unlist(paragraphs, use.names = FALSE), context), collapse = "\n\n")
}

