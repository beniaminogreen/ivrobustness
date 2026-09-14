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
