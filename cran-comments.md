## Submission

This is a resubmission. Changes since last CRAN version (0.15.2):

### New functions
* `css_pq()`, `gmpr_pq()`, `mcknight_residuals_pq()`, `rarefy_pq()`, `srs_pq()`,
  `tmm_pq()`, `transform_pq()`, `vst_pq()` — normalisation and transformation helpers.
* `hill_bar_pq()` — Hill diversity bar charts with error bars and compact letter display.
* `profile_hill_pq()` — Hill diversity profiles via `divent::profile_hill()`.
* `ridges_sam_pq()` — sample-centric ridge plots.
* `divent_hill_matrix_pq()` — Hill numbers for all samples via `divent::div_hill()`.
* `find_vsearch()`, `install_vsearch()` — cross-platform vsearch discovery and install.
* `unwanted_tax_patterns` — exported named character vector of regex patterns for
  common problematic taxonomy values; used as default in `verify_tax_table()`.

### Enhancements
* `hill_pq()`, `hill_tuckey_pq()`, `ggbetween_pq()`, `compare_pairs_pq()`,
  `psmelt_samples_pq()`, `hill_acc_pq()` — switched to `divent::div_hill()` for
  Hill numbers; default estimator is now `"UnveilJ"` (bias-corrected).
* `tax_bar_pq()` — OTU counting bug fixed for `nb_seq = FALSE` with a grouping
  factor; `(n=X)` sample-count label now displayed below bars in a separate layer;
  gains `n_sample_text_size` parameter.
* `cutadapt_remove_primers()` — gains `cutadapt_args` and `verbose` parameters.
* `biplot_pq()` — gains `color_rank` and `taxa_names_rank` parameters.
* `search_exact_seq_pq()` — now emits a clear error when a plain character vector
  is passed instead of a `DNAStringSet`.
* Many functions accepting `fact` now handle single-level factors gracefully.

### Bug fixes
* `chimera_removal_vs()` — fixed matrix dimension drop with single-sample input.
* `format2sintax()` — fixed wrong internal name for `pattern_tax` parameter.
* `umap_pq()` — no longer emits a tibble `.name_repair` deprecation warning.

## Test environments

* Local Linux (Pop!_OS 24.04), R 4.5.2

## Method References

Important features are described in Taudière A. (2023)
<doi:10.21105/joss.06038>. No additional published references describe
the methods in this package beyond the cited publication.
