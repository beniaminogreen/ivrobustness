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
