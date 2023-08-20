# tidyestimate 1.1.1

## Minor improvements and fixes

* Update R-CMD-Check GitHub actions
* Remove unused dependencies
* Add required "_PACKAGE" sentinel to package documentation

# tidyestimate 1.1.0

## Bugfixes

* When using `filter_common_genes`, supplying a dataframe that has a gene column that is NOT named `hgnc_symbol` or `entrezgene_id` now works as expected
* When using `filter_common_genes` and supplying entrez ids, the returned dataframe is now HGNC symbols to better flow into `estimate_score`
* Failure to supply `is_affymetrix` results in an earlier error, before calculating scores.

# tidyestimate 1.0.4

* CRAN release

# tidyestimate 1.0.3

# tidyestimate 1.0.2

* Added a `NEWS.md` file to track changes to the package.
