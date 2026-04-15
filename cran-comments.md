## Resubmission

This is a resubmission. Changes in this version (0.15.2):

* Fixed `\donttest{` spacing in `transform_pq()` examples that caused an Rd parse warning.
* Removed commented-out code from `transform_pq()` examples.
* Added `@return` documentation to `select_taxa` generic.
* Added new normalisation functions: `css_pq()`, `gmpr_pq()`, `mcknight_residuals_pq()`, `rarefy_pq()`, `srs_pq()`, `tmm_pq()`, `transform_pq()`, `vst_pq()`.
* Improved `hill_bar_pq()` with new parameters for error bars, point transparency, and letter placement.
* Fixed `tax_bar_pq()` OTU counting bug when using `nb_seq = FALSE` with a grouping factor.

## Test environments

* Local Linux (Pop!_OS 24.04), R 4.5.2

## R CMD check results

0 ERRORs | 0 WARNINGs | 0 NOTEs

## Known downstream dependency issue

* `ancombc_pq()` and `signif_ancombc()` fail when the 'CVXR' package API changes
  (symbol 'solve' removed from its namespace). This is an upstream issue in 'ANCOMBC2';
  a CVXR version constraint has been noted for the next release.

## Method References

Important features are described in Taudière A. (2023)
<doi:10.21105/joss.06038>. No additional published references describe
the methods in this package beyond the cited publication.
