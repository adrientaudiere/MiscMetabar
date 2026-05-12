## Submission

This is a resubmission of MiscMetabar following the failure of the
incoming-checks for 0.16.2:

> Status: 2 ERRORs, 1 WARNING (Windows); 2 ERRORs (Debian)
>
>   Error in .local(x, row.names, optional, ...) :
>     unused argument (validRN = FALSE)
>   Calls: write_pq ... as.data.frame.Vector -> as.data.frame -> as.data.frame

Changes since 0.16.2:

* `write_pq()` no longer passes the raw `DNAStringSet` `refseq` slot to
  `utils::write.table()`. It is now coerced to plain `character` via
  `as.character()` first. R-devel's `data.frame()` forwards an internal
  `validRN = FALSE` argument when converting its inputs; the
  `as.data.frame,XStringSet-method` defined in `Biostrings` re-dispatches
  through a `.local()` closure whose fixed signature does not accept
  `validRN`, which produced the example/test failure on win-builder.

* `Biostrings` moved from `Suggests` to `Imports`, so the `XVector` S4
  classes stored inside `data/data_fungi*.rda` are now in the package's
  recursive strong-dependency graph. This addresses the WARNING:

  > Data files with namespace references not in the recursive strong
  > package dependencies: data/data_fungi.rda: XVector …

* `NEWS.md` link to `https://www.indrapatil.com/ggstatsplot/` (which
  fails SSL handshake) replaced with the CRAN URL for `ggstatsplot`.

* Plus the minor `clean_pq()` enhancement listed in `NEWS.md`.

## Test environments

* Local Linux (Pop!_OS 24.04), R 4.5.2 — 0 errors, 0 warnings, 0 notes
* Local reproduction of R-devel's `data.frame(validRN = FALSE)` path,
  exercising `write_pq()` / `save_pq()` / `read_pq()` — passes.

## Method References

Important features are described in Taudière A. (2023)
<doi:10.21105/joss.06038>. No additional published references describe
the methods in this package beyond the cited publication.
