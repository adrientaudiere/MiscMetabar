## Resubmission after archival

MiscMetabar was archived on CRAN on 2026-05-19 because the check time
exceeded CRAN's 10-minute budget. This resubmission addresses that issue
along with other improvements since 0.16.5.

Key changes (also in `NEWS.md`):

* Extensive example speed-ups: 26 functions switched from `data_fungi`
  (185 × 1420) to `data_fungi_mini` (137 × 45), cutting individual
  example times from seconds-to-minutes down to < 1 s.
* CPU-intensive examples (`hill_acc_pq()`, `iNEXT_pq()`,
  `format2dada2()`, `hill_test_rarperm_pq()`) moved from `\donttest{}`
  to `\dontrun{}` — they are inherently long-running (permutation ×
  rarefaction × q-loop) and are documented in the package vignettes.
* Hot loops vectorised: `divent_hill_matrix_pq()` and `circle_pq()`
  run ~10× faster, further reducing check overhead.
* `verify_tax_table()` ~10× faster on full-size taxonomy tables.
* `XVector` added to `Imports` (explicitly) and `importFrom(XVector,subseq)`
  declared in `NAMESPACE`. The package data files (`data_fungi*.rda`)
  store `DNAStringSet` (an `XVector`-based class); the explicit import
  ensures `XVector` is loaded before the data objects are accessed,
  regardless of load order.

## Test environments

* Local Linux (Pop!_OS 24.04), R 4.5.2 — 0 errors, 0 warnings, 3 NOTEs.

## NOTEs

1. **CRAN incoming feasibility** — "Package was archived on CRAN" and
   "Version contains large components (0.16.5.9000)". Both are addressed
   by this resubmission and the version bump to a release version.

2. **Future file timestamps** — "unable to verify current time". This
   is a transient note caused by an unreachable NTP server at check time;
   no action needed.

3. **Namespace in Imports field not imported from: 'XVector'** —
   Addressed in this version by adding `importFrom(XVector,subseq)` to
   NAMESPACE. `XVector` is declared in `Imports` because the package
   data files contain `DNAStringSet` objects (XVector-based classes)
   and `Biostrings` alone does not guarantee that `XVector` is loaded
   before those objects are accessed.

## URL notes

The DOI `https://doi.org/10.1111/j.1365-294X.2012.05542.x` returns
HTTP 403 for automated checkers (doi.org blocks bot requests), but
resolves correctly in a browser. No action needed.

## Package coverage (`covr`): 46.88 %

## Method References

Important features are described in Taudière A. (2023)
<doi:10.21105/joss.06038>. No additional published references describe
the methods in this package beyond the cited publication.
