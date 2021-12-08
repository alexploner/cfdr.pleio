# cfdr.pleio

The `cfdr.pleio` package implements conditional and conjunctional false discovery rates to leverage and study pleiotropic gene variants between pairs of phenotypic traits, as originally described by [Andreassen et al. (2013)](https://www.ncbi.nlm.nih.gov/pmc/articles/pmid/23375658/).

## Installation

You can install the development version of `cfdr.pleio` like so:

``` r
if ( !require("remotes", quietly = TRUE) ) install.packages("remotes")
remotes::install_github("alexploner/cfdr.pleio", build_vignettes = TRUE)
```

`cfdr.pleio` also needs a set of genetic reference data in a local directory; you can download a local copy like so:

``` r
download.file("https://zenodo.org/record/5750318/files/genref.zip",
               destfile = "/path/to/dir/genref.zip")
unzip("/path/to/dir/genref.zip", exdir = "/path/to/dir")               
```
You can delete the zip file after successful unpacking. 

## Example

For a basic example, install the package and read the introductory vignette:

``` r
vignette("Introduction", package = "cfdr.pleio")
```

## Credit

`cfdr.pleio` is a based on a fork (in February 2021) of the MATLAB package `pleiofdr` available from https://github.com/precimed/pleiofdr 





