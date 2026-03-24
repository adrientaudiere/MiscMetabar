## Resubmission

This is a resubmission. Changes in this version (0.15.1):

* Added `@examples` to all exported functions that were missing examples.
* Fixed commented-out code in `assign_sintax()` examples.
* Quoted 'R' in the DESCRIPTION field.
* Fixed `hill_tuckey_pq()` double-transpose bug causing "variable lengths differ" errors.
* Fixed `plot_refseq_extremity_pq()` example (replaced large dataset with smaller one).
* Added `divent` to NAMESPACE imports.

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
