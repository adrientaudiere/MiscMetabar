# Calculate ecological distance among positive controls vs distance for all samples

[![lifecycle-experimental](https://img.shields.io/badge/lifecycle-experimental-orange)](https://adrientaudiere.github.io/MiscMetabar/articles/Rules.html#lifecycle)

Compute distance among positive controls, i.e. samples which are
duplicated to test for variation, for example in (i) a step in the
sampling, (ii) a step in the extraction, (iii) a step in the sequencing.

## Usage

``` r
dist_pos_control(physeq, samples_names, method = "bray")
```

## Arguments

- physeq:

  (required) a
  [`phyloseq-class`](https://rdrr.io/pkg/phyloseq/man/phyloseq-class.html)
  object obtained using the `phyloseq` package.

- samples_names:

  (required) a vector of names for samples with positives controls of
  the same samples having the same name

- method:

  (default: "bray") a method to calculate the distance, parsed to
  [`vegan::vegdist()`](https://vegandevs.github.io/vegan/reference/vegdist.html).
  See ?vegdist for a list of possible values.

## Value

A list of two data-frames with (i) the distance among positive controls
and (ii) the distance among all samples

## Author

Adrien Taudière

## Examples

``` r
data("enterotype")
sam_name_factice <- gsub("TS1_V2", "TS10_V2", sample_names(enterotype))
res_dist_cont <- dist_pos_control(enterotype, sam_name_factice)
hist(unlist(res_dist_cont$distAllSamples))
abline(
  v = mean(unlist(res_dist_cont$dist_controlontrolSamples), na.rm = TRUE),
  col = "red", lwd = 3
)
```
