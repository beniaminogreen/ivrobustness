# Interpret a Table Of Robustness Values

Struggling to interpret the results of `robustness_values`? Look no
further! `interpret_table` helps by providing a plain-text
interpretation of the robustness values of a design. It interprets each
of the four robustness values for you, and flags low robustness values.

## Usage

``` r
interpret_table(df)
```

## Arguments

- df:

  a datasdrame, as returned by `robustness_values`

## Value

A string

## Examples

``` r
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


weak_instrument_data <- create_simulation(rho = .01)

df <- robustness_values(
 weak_instrument_data$w, 
 weak_instrument_data$z, 
 weak_instrument_data$y
)
cat(interpret_table(df))
#> The robustness value for selection on observables is 0.827. A bias of 2 requires an unobserved confounder to explain at least 82.7% of residual variation in at least one of the treatment or the outcome.
#> 
#> The robustness value for instrumental variables is 0.00458. A bias of 2 requires either the confounder to explain at least 0.458% of residual instrument variation, or the instrument to explain at least 0.458% of outcome variation remaining after adjusting for the treatment and confounder. This robustness value is below 5%. While there is no 'hard limit' on these values, a value this small means that relatively small violations of the corresponding identifying assumptions could drive the point estimate to zero.
#> 
#> With instrument exogeneity imposed, the robustness value is 0.0402. A bias of 2 requires the instrument to explain at least 4.02% of outcome variation remaining after adjusting for the treatment and confounder. This robustness value is below 5%. While there is no 'hard limit' on these values, a value this small means that relatively small violations of the corresponding identifying assumptions could drive the point estimate to zero.
#> 
#> With exclusion imposed, the robustness value is 0.00946. A bias of 2 requires the confounder to explain at least 0.946% of residual instrument variation. This robustness value is below 5%. While there is no 'hard limit' on these values, a value this small means that relatively small violations of the corresponding identifying assumptions could drive the point estimate to zero.
#> 
#> Consider reading the associated paper for more information on interpreting robustness values (https://arxiv.org/pdf/2507.23743)
```
