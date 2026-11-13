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

## Examples in `\dontrun{}`

The 22 help pages that keep a `\dontrun{}` block, by reason:

* They download and install external software: `install_itsx()`,
  `install_itsxpress()`, `install_mmseqs2()`, `install_vsearch()`.
* They query a web service (FUNGuild, Open Tree of Life):
  `add_funguild_info()`, `funguild_assign()`, `get_funguild_db()`,
  `plot_guild_pq()`, `rotl_pq()`.
* They need files the user provides or that an external program produces:
  `lulu()`, `sam_data_matching_names()`.
* They take from 6 s to several minutes, which the 10-minute check budget
  cannot absorb: `add_new_taxonomy_pq()`, `ancombc_pq()`, `assign_dada2()`,
  `assign_idtaxa()`, `hill_acc_pq()`, `hill_test_rarperm_pq()` (deprecated,
  moved to bootpq), `iNEXT_pq()`, `learn_idtaxa()`, `plot_ancombc_pq()`,
  `plot_volcano_pq()`, `signif_ancombc()`.

Examples that only need an external program (BLAST, vsearch, MMseqs2, MAFFT,
Krona, fastp, cutadapt, mumu, ITSx, ITSxpress) run when it is installed, via
`@examplesIf` or `if ()`.

* win-builder R-devel (r90881 ucrt, 2026-11-13): examples 195 s, tests
  123 s; 0 errors, 0 warnings, 1 note (the dev version number).
* R-hub (2026-11-13), `R CMD check --as-cran` with `--run-donttest`:
  ubuntu-clang (R-devel r90876) examples 241 s + 397 s, tests 195 s;
  ubuntu-gcc12 (R-devel r89874) examples 235 s + 400 s, tests 176 s.
  Status OK on both.

## R CMD check results

0 errors | 0 warnings | 1 note

* The NOTE is "Version contains large components (0.18.0.9000)" on the
  development version; it disappears with the release version 0.18.0.

## Known issue

* r-oldrel-macos-arm64: the ERROR "Package required but not available:
  'dada2'" comes from the Bioconductor binaries missing for that platform,
  not from MiscMetabar.
