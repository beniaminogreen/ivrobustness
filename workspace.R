library(microbenchmark)
library(tidyverse)

# rextendr::vendor_crates()
devtools::document()
devtools::load_all()

create_simulation <- function(n=10^7, rho = .5) {
  w <- rnorm(n)
  z <- sqrt(rho)*w + sqrt(1-rho)*rnorm(n)
  y <- 2*z + rnorm(n)
  
  return(list(
    w = w, 
    z = z, 
    y = y
  ))
}


strong_instrument_data <- create_simulation(rho = .90)
weak_instrument_data <- create_simulation(rho = .01)

robustness_values(
  strong_instrument_data$w, 
  strong_instrument_data$z, 
  strong_instrument_data$y
)
df <- robustness_values(
  weak_instrument_data$w, 
  weak_instrument_data$z, 
  weak_instrument_data$y
)


interpret_table <- function(df) {
  paragraphs <- lapply(seq_len(nrow(df)), function(i) {
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

  paste(c(unlist(paragraphs, use.names = FALSE), context), collapse = "\n\n")
}


interpret_table(df)
