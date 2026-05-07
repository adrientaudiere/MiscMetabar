# Reorder fill and color scales to maximize perceptual contrast between adjacent segments

[![lifecycle-experimental](https://img.shields.io/badge/lifecycle-experimental-orange)](https://adrientaudiere.github.io/MiscMetabar/articles/Rules.html#lifecycle)

In stacked bar plots, ggplot2's default discrete palette assigns colors
using level ordered (sometimes alphabetically), which often places
perceptually similar colors next to each other. This function reassigns
the **same set of colors** to factor levels so that visually adjacent
segments receive maximally different colors. Both the fill and color
scales are updated so that direct labels (e.g. from `label_taxa = TRUE`)
stay in sync with the bars.

## Usage

``` r
reorder_distinct_colors(
  p = NULL,
  alternate_lightness = FALSE,
  lightness_amount = 0.15,
  colorblind = FALSE
)
```

## Arguments

- p:

  A ggplot object that uses a discrete fill aesthetic. Can be omitted
  when using the `+` operator (e.g. `p + reorder_distinct_colors()`).

- alternate_lightness:

  (logical, default FALSE) If TRUE, darken every other level to add a
  luminance alternation cue on top of hue differences.

- lightness_amount:

  (numeric, default 0.15) Intensity of the lightness alternation
  (proportion to darken). Only used when `alternate_lightness = TRUE`.

- colorblind:

  (logical, default FALSE) If TRUE, compute perceptual distances under
  simulated deuteranopia so that the reordering optimizes contrast for
  colorblind viewers.

## Value

A new ggplot object with
[`ggplot2::scale_fill_manual()`](https://ggplot2.tidyverse.org/reference/scale_manual.html)
and (if a color scale is present)
[`ggplot2::scale_color_manual()`](https://ggplot2.tidyverse.org/reference/scale_manual.html)
replacing the original scales. When `p` is omitted, returns an object
that can be added to a ggplot with `+`.

## Author

Adrien Taudière

## Examples

``` r
p <- tax_bar_pq(data_fungi_mini, taxa = "Class", fact = "Time")
reorder_distinct_colors(p)

reorder_distinct_colors(p, colorblind = TRUE)

p + reorder_distinct_colors(alternate_lightness = TRUE)


tax_bar_pq(data_fungi_mini, fact = "Height", taxa = "Order",
 nb_seq = FALSE, percent_bar = TRUE, label_taxa = TRUE,
 add_ribbon = TRUE, value_size=7, ribbon_alpha = .6,
 show_values=TRUE, label_size = 4, top_label_size = 8,
 minimum_value_to_show=0.05) |>
 reorder_distinct_colors(alternate_lightness=TRUE)
```
