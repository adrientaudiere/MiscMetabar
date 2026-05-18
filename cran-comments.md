## Submission

This is a routine patch release of MiscMetabar following the CRAN release
of 0.16.4. It fixes example failures triggered by network-dependent
`\donttest{}` blocks and slow example timings, and removes an unused
`Imports` declaration. No user-visible API changes.

Changes since 0.16.4 (also in `NEWS.md`):

* `funguild_assign()` and `rotl_pq()` examples now use `\dontrun{}`
  instead of `\donttest{}`. Both examples call external APIs
  (`www.stbates.org` and the Open Tree of Life respectively) that are not
  always reachable under CRAN's `--run-donttest`, which produced spurious
  ERRORs in `R CMD check --as-cran`.

* `verify_tax_table()`'s introductory example was moved inside the
  existing `\donttest{}` block. The call against the full `data_fungi`
  dataset ran for ~70 s and triggered the "examples > 5 s" NOTE.

* `XVector` was removed from `DESCRIPTION` `Imports`. It was declared but
  never imported in `NAMESPACE` or used directly; `Biostrings` already
  loads it transitively. This resolves the CRAN NOTE "Namespace in
  Imports field not imported from: 'XVector'".

* Corrected the DOI for Taberlet et al. (2012) "Environmental DNA" in
  the `paper/` and `vignettes/` `.bib` files (was the journal ISSN
  landing, now the paper DOI).

## Test environments

* Local Linux (Pop!_OS 24.04), R 4.5.2 — 0 errors, 0 warnings, 3 NOTEs.
  The three notes are the expected "Version contains large components"
  (resolved by this submission), "unable to verify current time" (local
  NTP), and a 70 s example for `verify_tax_table` which is now wrapped
  in `\donttest{}` and so no longer fires on `--as-cran` without
  `--run-donttest`.

* Test suite: 1088 tests, 0 failures, 31 warnings (concentrated in
  `test-normalize_pq.R` and pre-existing on the previous release),
  2 skipped.

* Package coverage (`covr`): 46.88 %.

## Method References

Important features are described in Taudière A. (2023)
<doi:10.21105/joss.06038>. No additional published references describe
the methods in this package beyond the cited publication.
