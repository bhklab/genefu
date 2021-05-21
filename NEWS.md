# genefu Package Updates

## v1.25.1
- Loss of lazy loading broke functionality in downstream packages
- As a result, we added a check that a gene signature `exists` in the dependent
functions, and load the requisite data if it does not exist
- Also fixed a bug in `claudinLow` where a non-imported function `standardize`
was called if `std=TRUE`; `scale` is now used instead, which z-score normalizes
the data
- Removed LazyData: TRUE from DESCRIPTION to correct Bioconductor build WARNING
- Added data calls to examples and vignettes to compensate for loss of lazy loading

## v1.25.0
- Bioc 3.13 release!
- All documentation has been refactored to use Roxygen2 in markdown format