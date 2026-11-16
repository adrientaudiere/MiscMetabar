## Update to version 0.18.0

This update gathers the changes of 0.17.0 (never submitted) and 0.18.0; see NEWS.md.

## Resubmission

In reply to the incoming check of 2026-11-14 (r-devel-windows-x86_64:
examples 307 s, tests 198 s, overall 18 min), this resubmission shortens
both stages:

* Tests: the slowest tests now run on `data_fungi_mini` (or a 20-sample
  subset of it) instead of `data_fungi`, and the bootstrap test of
  `strassoc_pq()` uses 10 replicates instead of 99. Locally (Linux), the six
  test files concerned went from 82 s to 22 s. The 13 slowest remaining tests
  (about 29 s locally) are skipped on CRAN with `skip_on_cran()`; they still
  run in our own checks; locally, the suite run as on CRAN went from 86 s to
  53 s.
* Examples: the 14 slowest R-only examples no longer run on CRAN (see
  "Examples skipped on CRAN" below). Redundant calls were removed from the
  examples of `circle_pq()` and `graph_test_pq()`, and `tbl_sum_samdata()`
  uses `data_fungi_mini`.
* win-builder R-devel (2026-11-16), before the 13 `skip_on_cran()` above:
  360 s for examples and tests together instead of 505 s.

## Check time

MiscMetabar was archived on 2026-05-19 for exceeding the 10-minute check
budget, and 0.16.8 was accepted afterwards. Before this resubmission, the
check time of 0.18.0 was at the level of 0.16.8 (figures below predate the
changes listed under "Resubmission"):

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

## Examples skipped on CRAN

To fit the check budget, 14 examples are wrapped in
`@examplesIf identical(Sys.getenv("NOT_CRAN"), "true") || identical(Sys.getenv("IN_PKGDOWN"), "true")`.
They run in our own checks (`devtools::check()` sets `NOT_CRAN`) and on the
pkgdown website, but not on CRAN:

* They load heavy Bioconductor or Suggests packages (ALDEx2, ggstatsplot,
  lefser and mia, DESeq2) that are slow to load on Windows: `aldex_pq()`,
  `ggbetween_pq()`, `lefser_pq()`, `plot_deseq2_pq()`.
* They are the slowest of the remaining examples that need only R packages:
  `accu_plot()`, `accu_samp_threshold()`, `circle_pq()`, `glmutli_pq()`,
  `graph_test_pq()`, `LCBD_pq()`, `plot_LCBD_pq()`, `plot_SCBD_pq()`,
  `plot_var_part_pq()`, `var_par_rarperm_pq()`.

Examples that only need an external program (BLAST, vsearch, MMseqs2, MAFFT,
Krona, fastp, cutadapt, mumu, ITSx, ITSxpress) run when it is installed, via
`@examplesIf` or `if ()`.

## R CMD check results

0 errors | 0 warnings | 0 notes

## Known issue

* r-oldrel-macos-arm64: the ERROR "Package required but not available:
  'dada2'" comes from the Bioconductor binaries missing for that platform,
  not from MiscMetabar.
