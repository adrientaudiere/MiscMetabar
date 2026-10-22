## Update to version 0.18.0

This update gathers the changes of 0.17.0 (never submitted) and 0.18.0; see NEWS.md.

## Check time

MiscMetabar was archived on 2026-05-19 for exceeding the 10-minute check
budget, and 0.16.8 was accepted afterwards. This version keeps the check time
at the level of 0.16.8:

* Local `R CMD check --as-cran --run-donttest` (R 4.6.1, Linux, 2026-10-22):
  524 s in total for 0.18.0 versus 533 s for 0.16.8 on the same machine,
  of which examples 236 s and tests 75 s. Examples were shortened to one
  short call per function; longer demonstrations live in the vignettes.
* Examples that need external software (e.g. vsearch, mmseqs2) are wrapped
  in `@examplesIf` and do not run when the software is not installed.
* Heavy test files are skipped on CRAN with a file-level `skip_on_cran()`.

* win-builder R-devel (Windows Server 2022, 2026-10-22): examples 203 s,
  tests 116 s; 0 errors, 0 warnings, 2 notes (the dev version number, and
  an Rd link to the not-yet-published `bootpq`, since removed).

## R CMD check results

0 errors | 0 warnings | 1 note

* The NOTE is "Version contains large components (0.18.0.9000)" on the
  development version; it disappears with the release version 0.18.0.

## Known issue

* r-oldrel-macos-arm64: the ERROR "Package required but not available:
  'dada2'" comes from the Bioconductor binaries missing for that platform,
  not from MiscMetabar.
