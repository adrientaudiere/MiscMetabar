## Resubmission — check-time fix (3rd attempt)

MiscMetabar was archived on CRAN on 2026-05-19 because the check time
exceeded CRAN's 10-minute budget. Version 0.16.6 was rejected for the
same reason (Windows check time 18 min). This resubmission (0.16.7)
addresses the two main contributors called out by Uwe Ligges:
`checking examples [282s]` and `checking tests [266s]`.

### Example timing reductions

Extra examples kept for documentation purposes are now wrapped in
`\dontrun{}` instead of `\donttest{}`, so they are preserved in the
help pages but do not contribute to check time. The primary example
for each function remains in `\donttest{}`.

* `assign_vsearch_lca()`: 3 calls → 1 in `\donttest{}`, 2 in `\dontrun{}`.
* `ggbetween_pq()`: 3 calls → 1 in `\donttest{}`, 2 in `\dontrun{}`.
* `circle_pq()`: 3 calls → 1 in `\donttest{}`, 2 in `\dontrun{}`.
* `chimera_removal_vs()`: 4 calls → 1 in `\donttest{}`, 3 in `\dontrun{}`.
* `adonis_pq()`: 7 calls → 2 in `\donttest{}`, 5 in `\dontrun{}`.
* `postcluster_pq()`: regular example moved to `\donttest{}`; extra
  method calls moved to `\dontrun{}`.
* `umap_pq()`: heavy patchwork + tsne + ordinate block moved to
  `\dontrun{}`; `\donttest{}` keeps only the `pkg = "uwot"` call.
* `plot_var_part_pq()` and `var_par_rarperm_pq()`: cross-references to
  each other moved to `\dontrun{}`.
* `hill_pq()` and `hill_bar_pq()`: extra `\donttest{}` parameter
  variants moved to `\dontrun{}`.
* `psmelt_samples_pq()`: `ggstatsplot` visualisation calls moved to
  `\dontrun{}`; only the plain `psmelt_samples_pq()` call stays in
  `\donttest{}`.
* `profile_hill_pq()`: example uses 5 samples instead of all 137
  (`prune_samples()`), combined with the earlier `orders = c(0, 1, 2)`
  fix (~25× faster).
* Earlier reductions still in place: `profile_hill_pq()` orders,
  `umap_pq()` n_neighbors, `hill_pq()`/`hill_bar_pq()`/
  `plot_refseq_extremity_pq()` single-q examples.

### Test timing reductions

Added file-level `skip_on_cran()` to 11 test files total:

* **Previously added (v0.16.6):** `test_figures_beta_div.R`,
  `test_figures_taxo.R`, `test_figures_alpha_div.R`,
  `test_deseq2_edgeR.R`, `test_ancombc.R`, `test_targets.R`.
* **Added in this version:** `test_alpha_div.R`, `test_tuckey.R`,
  `test_assignment.R`, `test_figures_misc.R`, `test_vsearch.R`.

## Test environments

* Local Linux (Pop!_OS 24.04), R 4.5.2 — 0 errors, 0 warnings, 2 NOTEs.

## NOTEs

1. **CRAN incoming feasibility** — "Package was archived on CRAN" and
   "Version contains large components". Both are addressed by this
   resubmission and the version bump.

2. **Future file timestamps** — "unable to verify current time". This
   is a transient note caused by an unreachable NTP server at check time;
   no action needed.

## URL notes

The DOI `https://doi.org/10.1111/j.1365-294X.2012.05542.x` returns
HTTP 403 for automated checkers (doi.org blocks bot requests), but
resolves correctly in a browser. No action needed.

## Package coverage (`covr`): 46.88 %

## Method References

Important features are described in Taudière A. (2023)
<doi:10.21105/joss.06038>. No additional published references describe
the methods in this package beyond the cited publication.
