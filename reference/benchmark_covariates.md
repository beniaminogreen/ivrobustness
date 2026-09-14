# Benchmark Robustness Values Based on Observed Covariates

Benchmark Robustness Values Based on Observed Covariates

## Usage

``` r
benchmark_covariates(w, z, y, X)
```

## Arguments

- w:

  the instrument vector, a numeric (double) vector with no NA values

- z:

  the treatment vector, a numeric (double) vector with no NA values

- y:

  the outcome vecctor, a numeric (double) vector with no NA values

- X:

  a matrix of control covaraites used for benchmarking

## Value

a dataframe containing the robustness values associated with held out
covariates
