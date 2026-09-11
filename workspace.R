library(microbenchmark)
library(tidyverse)

rextendr::vendor_crates()

devtools::document()
# devtools::test()

calc_robustness_values(letters, letters, letters)
calc_robustness_values(1, letters, letters)
is.double(1L)

pkgdown::build_site(preview = TRUE)
