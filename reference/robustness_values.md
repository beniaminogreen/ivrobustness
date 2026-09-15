# Calculate Robustness Values

Calculate Robustness Values

## Usage

``` r
robustness_values(w, z, y)
```

## Arguments

- w:

  the instrument vector, a numeric (double) vector with no NA values

- z:

  the treatment vector, a numeric (double) vector with no NA values

- y:

  the outcome vecctor, a numeric (double) vector with no NA values

## Value

a dataframe containing robustness values with the following columns:

- *Method* The method for the associated Robustness Value:

  - SOO - Selection on observables

  - IV - Instrumental variables

  - ExogenousIV - Instrumental variables, where exogeneity is known to
    be satisfied (as in a randomized experiment)

  - ExcludableIV - Instrumental variables, where the exclusion
    restriction is known to be satisfied

- *RV* the robustness value

- *rho_1, rho_2, rho_3* The parameters used to achieve this robustness
  value

- *bias* The bias achieved at this combination of parameters

- *assumption gap* the largest gap between partial correlations that are
  assumed to be zero (only relevant for ExogenousIV and ExcludableIV)
  and the actual value found by the optimizer. This should be negligible
  and close to the machine tolerance.

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


strong_instrument_data <- create_simulation(rho = .90)
weak_instrument_data <- create_simulation(rho = .01)

robustness_values(
 strong_instrument_data$w, 
 strong_instrument_data$z, 
 strong_instrument_data$y
)
#> # A tibble: 4 × 8
#>   method           rv    rho_1 rho_2 rho_3  bias assumption_gap estimate
#>   <chr>         <dbl>    <dbl> <dbl> <dbl> <dbl>          <dbl>    <dbl>
#> 1 Soo          0.463  4.56e- 5 0.681 0.681  2.00       0            2.00
#> 2 Iv           0.0199 3.85e- 2 0.975 0.145  2.00       0            2.00
#> 3 ExogenousIV  0.0199 3.30e-12 0.975 0.145  2.00       3.30e-12     2.00
#> 4 ExcludableIV 0.720  8.48e- 1 0.535 1.000  2.00       1.38e-13     2.00

robustness_values(
 weak_instrument_data$w, 
 weak_instrument_data$z, 
 weak_instrument_data$y
)
#> # A tibble: 4 × 8
#>   method            rv     rho_1 rho_2 rho_3  bias assumption_gap estimate
#>   <chr>          <dbl>     <dbl> <dbl> <dbl> <dbl>          <dbl>    <dbl>
#> 1 Soo          0.827   -1.67e- 4 0.910 0.910  2.00       0            2.00
#> 2 Iv           0.00463  6.80e- 2 0.974 0.464  2.00       0            2.00
#> 3 ExogenousIV  0.0406   3.33e-12 0.975 0.457  2.00       3.33e-12     2.00
#> 4 ExcludableIV 0.00956  9.78e- 2 0.975 0.457  2.00       4.26e-12     2.00
```
