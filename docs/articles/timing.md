# Function timing

``` r

library(MiscMetabar)
```

## Why this article exists

Metabarcoding datasets vary by two to three orders of magnitude in size
— a small pilot may have a few dozen samples and a few hundred OTUs,
while a regional-scale dataset can reach thousands of samples and tens
of thousands of OTUs. Functions in `MiscMetabar` were written with the
small / medium case in mind, and several of them have a cost that grows
non-linearly with the input.

This article reports wall-clock timings of the main exported functions
on two reference datasets shipped with the package, so you can:

- pick reasonable defaults for permutation- and bootstrap-based
  analyses,
- decide which calls are worth wrapping in
  [`targets::tar_target()`](https://docs.ropensci.org/targets/reference/tar_target.html)
  to avoid recomputing them,
- and know which functions to subset (samples or taxa) before running
  them on a larger dataset.

The two reference datasets are:

| Dataset | Samples | Taxa | Use |
|----|----|----|----|
| `data_fungi_mini` | 45 | 137 | Quick sanity checks; examples; unit tests |
| `data_fungi` | 185 | 1420 | Realistic medium-scale dataset; what most users will plug in |

The timings below were collected on the package’s CI machine
(`x86_64-pc-linux-gnu`, 4 cores, 16 GB RAM) with R 4.5.2 and one
replicate per call. Numbers should be read as orders of magnitude, not
as benchmarks.

## Timing table

|     | Function                   | data_fungi_mini (s) | data_fungi (s) |
|:----|:---------------------------|--------------------:|---------------:|
| 1   | hill_pq                    |                3.19 |           9.91 |
| 3   | hill_tuckey_pq             |                1.85 |           4.74 |
| 5   | profile_hill_pq            |               12.96 |          19.79 |
| 7   | hill_acc_pq\[sample,n=10\] |                3.44 |           5.70 |
| 9   | adonis_pq                  |                  NA |           0.26 |
| 11  | graph_test_pq              |                  NA |           0.78 |
| 13  | plot_tsne_pq               |                  NA |           0.26 |
| 15  | verify_pq                  |                0.00 |           0.01 |
| 17  | verify_tax_table           |                0.20 |           5.14 |
| 19  | summary_plot_pq            |                0.08 |           0.12 |
| 21  | ggvenn_pq                  |                0.30 |           0.37 |

Elapsed time (seconds) per function and dataset. One replicate. {.table}

## How to speed up your own analysis

### Subset before computing

For every function that ends in `_pq`, the cost is dominated either by
the number of samples (permutation tests, ordinations) or by the number
of taxa (taxonomy checks, plotting). Subset first:

``` r

# Keep only the most abundant taxa
data_fungi_subset <- subset_taxa_pq(data_fungi, taxa_sums(data_fungi) > 1000)

# Or merge samples by a grouping variable before plotting
data_fungi_merged <- merge_samples2(data_fungi, "Height")
```

### Cache permutation-based results

[`adonis_rarperm_pq()`](https://adrientaudiere.github.io/MiscMetabar/reference/adonis_rarperm_pq.md),
[`hill_test_rarperm_pq()`](https://adrientaudiere.github.io/MiscMetabar/reference/hill_test_rarperm_pq.md),
[`var_par_rarperm_pq()`](https://adrientaudiere.github.io/MiscMetabar/reference/var_par_rarperm_pq.md),
and `hill_acc_pq(type = "sample")` are the most expensive functions
because they run a full permutation loop. They are also the most worth
caching with `targets` or `memoise`:

``` r

library(memoise)
adonis_pq_m <- memoise(adonis_pq)
adonis_pq_m(data_fungi, "Height") # slow the first time
adonis_pq_m(data_fungi, "Height") # instant the second time
```

### Lower iteration counts for exploratory work

Most permutation-based functions take an iteration argument (`nperm`,
`n_permutations`, `nboot`, `n_simulations`). The defaults are tuned for
final, publishable analyses. For exploratory work, lowering them by an
order of magnitude is almost always fine:

``` r

# Default: 999 permutations
adonis_pq(data_fungi_mini, "Height")

# Exploration: 99 permutations
adonis_pq(data_fungi_mini, "Height", nperm = 99)
```

The function defaults stay at the published values for reproducibility;
the choice is yours to lower them in your own code.

### Use `divent::` and `iNEXT.4steps::` for very large datasets

For datasets larger than `data_fungi`, the underlying packages `divent`
(Hill numbers) and `iNEXT.4steps` (coverage-based rarefaction) expose
lower-level functions that skip phyloseq slot lookups and can be called
directly on the OTU matrix.

## Reproducing this table

The script that generates `inst/benchmark/function_timings.csv` is
shipped with the package:

``` r

system.file("benchmark", "function_timings.R", package = "MiscMetabar")
```

Run it from the package source root:

``` bash
Rscript inst/benchmark/function_timings.R
```

Re-rendering this article will then pick up the updated CSV.

## Internal speedups since 0.16.5

Recent versions of `MiscMetabar` include internal rewrites of the most
expensive helpers. The user-facing behaviour is unchanged; only the wall
time differs. Where the speedup is documented (`NEWS.md`), it is the
result of vectorising one or more inner loops or factoring out repeated
rarefaction work, **not** of changing the scientific defaults. See
`NEWS.md` for the version-by-version list.
