
# IVrobustness

`IVrobustness` allows you to construct the robustness values for
comparing various identification strategies described by [Huang and
McCartan (2026)](https://arxiv.org/pdf/2507.23743) for selection on
observables and instrumental variables designs.

At the moment the following features are supported:

- Calculating Robustness Values (RVs) for selection-on-observables and
  IV designs.
- Benchmarking the robustness values using held-out covariates

Future versions will include:

- Standard errors and confidence intervals for Robustness Values based
  on

## Installation

To use the package, you must first install it. The core of the package
is written in [Rust](https://rust-lang.org/) for speed, which means you
must have a rust compiler installed if you are building the package from
scratch. If you don’t want to do this, you can use the versions hosted
on R-universe, which are available for all major platforms.

### Installing Compiled Version from R-Universe:

``` r
install.packages('ivrobustness', repos = c('https://beniaminogreen.r-universe.dev'))
```

### Installing Rust

If your operating system or version of R is not installed, you must have
the [Rust compiler](https://rust-lang.org/tools/install/) installed to
compile this package from sources. After the package is compiled, Rust
is no longer required, and can be safely uninstalled.

#### Installing Rust on Linux or Mac:

To install Rust on Linux or Mac, you can simply run the following
snippet in your terminal.

``` sh
curl --proto '=https' --tlsv1.2 -sSf https://sh.rustup.rs | sh
```

#### Installing Rust on Windows:

To install Rust on windows, you can use the Rust installation wizard,
`rustup-init.exe`, found [at this
site](https://forge.rust-lang.org/infra/other-installation-methods.html).
Depending on your version of Windows, you may see an error that looks
something like this:

    error: toolchain 'stable-x86_64-pc-windows-gnu' is not installed

In this case, you should run
`rustup install stable-x86_64-pc_windows-gnu` to install the missing
toolchain. If you’re missing another toolchain, simply type this in the
place of `stable-x86_64-pc_windows-gnu` in the command above.

### Installing Package from Github:

Once you have rust installed Rust, you should be able to install the
package with either the install.packages function as above, or using the
`install_github` function from the `devtools` package or with the
`pkg_install` function from the `pak` package.

``` r
## Install with devtools
# install.packages("devtools")
devtools::install_github("beniaminogreen/ivrobustness")

## Install with pak
# install.packages("pak")
pak::pkg_install("beniaminogreen/ivrobustness")
```

### Loading The Package

Once the package is installed, you can load it into memory as usual by
typing:

``` r
library(ivrobustness)
```
