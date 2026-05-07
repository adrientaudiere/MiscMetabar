# Create a visualization table to describe taxa distribution across a modality

[![lifecycle-maturing](https://img.shields.io/badge/lifecycle-maturing-blue)](https://adrientaudiere.github.io/MiscMetabar/articles/Rules.html#lifecycle)

Allow to visualize a table with graphical input.

## Usage

``` r
formattable_pq(
  physeq,
  modality,
  taxonomic_levels = c("Phylum", "Order", "Family", "Genus"),
  min_nb_seq_taxa = 1000,
  log10trans = FALSE,
  void_style = FALSE,
  lev_col_taxa = "Phylum",
  arrange_by = "nb_seq",
  descending_order = TRUE,
  na_remove = TRUE,
  formattable_args = NULL
)
```

## Arguments

- physeq:

  (required) a
  [`phyloseq-class`](https://rdrr.io/pkg/phyloseq/man/phyloseq-class.html)
  object obtained using the `phyloseq` package.

- modality:

  (required) The name of a column present in the `@sam_data` slot of the
  physeq object. Must be a character vector or a factor.

- taxonomic_levels:

  (default = c("Phylum", "Order", "Family", "Genus")) The taxonomic
  levels (must be present in the `@sam_data` slot) you want to see
  and/or used (for example to compute a color) in the table.

- min_nb_seq_taxa:

  (default = 1000) filter out taxa with less than `min_nb_seq_taxa`
  sequences

- log10trans:

  (logical, default TRUE) Do sequences count is log10 transformed (using
  log10(x + 1) to allow 0)

- void_style:

  (logical, default FALSE) Do the default style is discard ?

- lev_col_taxa:

  Taxonomic level used to plot the background color of taxa names

- arrange_by:

  The column used to sort the table. Can take the values NULL,
  "proportion_samp", "nb_seq" (default), , "nb_sam" "OTU", or a column
  names from the levels of modality or from taxonomic levels

- descending_order:

  (logical, default TRUE) Do we use descending order when sort the table
  (if arrange_by is not NULL) ?

- na_remove:

  (logical, default TRUE) if TRUE remove all the samples with NA in the
  `split_by` variable of the `physeq@sam_data` slot

- formattable_args:

  Other args to the formattable function. See examples and
  [`formattable::formattable()`](https://renkun-ken.github.io/formattable/reference/formattable.html)

## Value

A datatable

## Details

This function is mainly a wrapper of the work of others. Please make a
reference to
[`formattable::formattable()`](https://renkun-ken.github.io/formattable/reference/formattable.html)
if you use this function.

## See also

[`formattable::formattable()`](https://renkun-ken.github.io/formattable/reference/formattable.html)

## Author

Adrien Taudière

## Examples

``` r
if (requireNamespace("formattable")) {
  ## Distribution of the nb of sequences per OTU across Height
  ##   modality (nb of sequences are log-transformed).
  ## Only OTU with more than 10000 sequences are taking into account
  ## The Phylum column is discarded
  formattable_pq(
    data_fungi,
    "Height",
    min_nb_seq_taxa = 10000,
    formattable_args = list("Phylum" = FALSE),
    log10trans = TRUE
  )

  ## Distribution of the nb of samples per OTU across Height modality
  ## Only OTU  present in more than 50 samples are taking into account
  formattable_pq(
    as_binary_otu_table(data_fungi),
    "Height",
    min_nb_seq_taxa = 50,
    formattable_args = list("nb_seq" = FALSE),
  )

  ## Distribution of the nb of sequences per OTU across Time modality
  ##  arranged by Family Name in ascending order.
  ## Only OTU with more than 10000 sequences are taking into account
  ## The Phylum column is discarded
  formattable_pq(
    data_fungi,
    "Time",
    min_nb_seq_taxa = 10000,
    taxonomic_levels = c("Order", "Family", "Genus", "Species"),
    formattable_args = list(
      Order = FALSE,
      Species = formattable::formatter(
        "span",
        style = x ~ formattable::style(
          "font-style" = "italic",
          `color` = ifelse(is.na(x), "white", "grey")
        )
      )
    ),
    arrange_by = "Family",
    descending_order = FALSE
  )
}
#> Loading required namespace: formattable
#> 54 samples were discarded due to NA in variable modality
#> Cleaning suppress 0 taxa (  ) and 22 sample(s) ( BE9-006-B_S27_MERGED.fastq.gz / BG7-010-H_S31_MERGED.fastq.gz / C21-NV1-M_S64_MERGED.fastq.gz / D9-027-B_S83_MERGED.fastq.gz / DJ2-008-B_S87_MERGED.fastq.gz / DJ2-008-H_S88_MERGED.fastq.gz / DY5-004-H_S97_MERGED.fastq.gz / DY5-004-M_S98_MERGED.fastq.gz / E9-009-B_S100_MERGED.fastq.gz / E9-009-H_S101_MERGED.fastq.gz / J18-004-B_S114_MERGED.fastq.gz / J18-004-H_S115_MERGED.fastq.gz / J18-004-M_S116_MERGED.fastq.gz / N22-001-B_S129_MERGED.fastq.gz / O20-X-B_S139_MERGED.fastq.gz / O21-007-M_S144_MERGED.fastq.gz / R28-008-H_S159_MERGED.fastq.gz / R28-008-M_S160_MERGED.fastq.gz / W26-001-H_S166_MERGED.fastq.gz / W26-001-M_S167_MERGED.fastq.gz / Y29-007-H_S182_MERGED.fastq.gz / Y29-007-M_S183_MERGED.fastq.gz ).
#> Number of non-matching ASV 0
#> Number of matching ASV 1420
#> Number of filtered-out ASV 1400
#> Number of kept ASV 20
#> Number of kept samples 109
#> Joining with `by = join_by(OTU)`
#> 54 samples were discarded due to NA in variable modality
#> Cleaning suppress 0 taxa (  ) and 35 sample(s) ( B18-006-B_S19_MERGED.fastq.gz / BE9-006-B_S27_MERGED.fastq.gz / BE9-006-H_S28_MERGED.fastq.gz / BE9-006-M_S29_MERGED.fastq.gz / BG7-010-B_S30_MERGED.fastq.gz / BG7-010-H_S31_MERGED.fastq.gz / BL7-006-B_S36_MERGED.fastq.gz / C21-NV1-M_S64_MERGED.fastq.gz / CB8-019-B_S69_MERGED.fastq.gz / CB8-019-H_S70_MERGED.fastq.gz / CB8-019-M_S71_MERGED.fastq.gz / D9-027-B_S83_MERGED.fastq.gz / DJ2-008-B_S87_MERGED.fastq.gz / DJ2-008-H_S88_MERGED.fastq.gz / DY5-004-H_S97_MERGED.fastq.gz / DY5-004-M_S98_MERGED.fastq.gz / E9-009-B_S100_MERGED.fastq.gz / E9-009-H_S101_MERGED.fastq.gz / J18-004-B_S114_MERGED.fastq.gz / J18-004-H_S115_MERGED.fastq.gz / J18-004-M_S116_MERGED.fastq.gz / N22-001-B_S129_MERGED.fastq.gz / N23-002-M_S132_MERGED.fastq.gz / O20-X-B_S139_MERGED.fastq.gz / O21-007-M_S144_MERGED.fastq.gz / R28-008-B_S158_MERGED.fastq.gz / R28-008-H_S159_MERGED.fastq.gz / R28-008-M_S160_MERGED.fastq.gz / T28-ABM602-B_S162_MERGED.fastq.gz / W26-001-H_S166_MERGED.fastq.gz / W26-001-M_S167_MERGED.fastq.gz / W9-025-M_S169_MERGED.fastq.gz / Y29-007-B_S181_MERGED.fastq.gz / Y29-007-H_S182_MERGED.fastq.gz / Y29-007-M_S183_MERGED.fastq.gz ).
#> Number of non-matching ASV 0
#> Number of matching ASV 1420
#> Number of filtered-out ASV 1417
#> Number of kept ASV 3
#> Number of kept samples 96
#> Joining with `by = join_by(OTU)`
#> 23 samples were discarded due to NA in variable modality
#> Cleaning suppress 0 taxa (  ) and 14 sample(s) ( BE9-006-B_S27_MERGED.fastq.gz / C21-NV1-M_S64_MERGED.fastq.gz / DJ2-008-B_S87_MERGED.fastq.gz / DY5-004-H_S97_MERGED.fastq.gz / DY5-004-M_S98_MERGED.fastq.gz / E9-009-B_S100_MERGED.fastq.gz / E9-009-H_S101_MERGED.fastq.gz / N22-001-B_S129_MERGED.fastq.gz / O21-007-M_S144_MERGED.fastq.gz / R28-008-H_S159_MERGED.fastq.gz / R28-008-M_S160_MERGED.fastq.gz / W26-001-M_S167_MERGED.fastq.gz / Y29-007-H_S182_MERGED.fastq.gz / Y29-007-M_S183_MERGED.fastq.gz ).
#> Number of non-matching ASV 0
#> Number of matching ASV 1420
#> Number of filtered-out ASV 1390
#> Number of kept ASV 30
#> Number of kept samples 148
#> Joining with `by = join_by(OTU)`
# \donttest{
if (requireNamespace("formattable")) {
  ## Distribution of the nb of sequences per OTU across Height modality
  ##  (nb of sequences are log-transformed).
  ## OTU name background is light gray for Basidiomycota
  ##  and dark grey otherwise (Ascomycota)
  ## A different color is defined for each modality level
  formattable_pq(
    data_fungi,
    "Height",
    taxonomic_levels = c("Phylum", "Family", "Genus"),
    void_style = TRUE,
    formattable_args = list(
      OTU = formattable::formatter(
        "span",
        style = ~ formattable::style(
          "display" = "block",
          `border-radius` = "5px",
          `background-color` = ifelse(Phylum == "Basidiomycota", transp("gray"), "gray")
        ),
        `padding-right` = "2px"
      ),
      High = formattable::formatter(
        "span",
        style = x ~ formattable::style(
          "font-size" = "80%",
          "display" = "inline-block",
          direction = "rtl",
          `border-radius` = "0px",
          `padding-right` = "2px",
          `background-color` = formattable::csscolor(formattable::gradient(
            as.numeric(x), transp("#1a91ff"), "#1a91ff"
          )),
          width = formattable::percent(formattable::proportion(as.numeric(x), na.rm = TRUE))
        )
      ),
      Low = formattable::formatter(
        "span",
        style = x ~ formattable::style(
          "font-size" = "80%",
          "display" = "inline-block",
          direction = "rtl",
          `border-radius` = "0px",
          `padding-right` = "2px",
          `background-color` = formattable::csscolor(formattable::gradient(
            as.numeric(x),
            transp("green"), "green"
          )),
          width = formattable::percent(formattable::proportion(as.numeric(x), na.rm = TRUE))
        )
      ),
      Middle = formattable::formatter(
        "span",
        style = x ~ formattable::style(
          "font-size" = "80%",
          "display" = "inline-block",
          direction = "rtl",
          `border-radius` = "0px",
          `padding-right` = "2px",
          `background-color` = formattable::csscolor(formattable::gradient(
            as.numeric(x), transp("orange"), "orange"
          )),
          width = formattable::percent(formattable::proportion(as.numeric(x), na.rm = TRUE))
        )
      )
    )
  )
}
#> 54 samples were discarded due to NA in variable modality
#> Cleaning suppress 0 taxa (  ) and 5 sample(s) ( DY5-004-M_S98_MERGED.fastq.gz / E9-009-B_S100_MERGED.fastq.gz / O21-007-M_S144_MERGED.fastq.gz / Y29-007-H_S182_MERGED.fastq.gz / Y29-007-M_S183_MERGED.fastq.gz ).
#> Number of non-matching ASV 0
#> Number of matching ASV 1420
#> Number of filtered-out ASV 1265
#> Number of kept ASV 155
#> Number of kept samples 126
#> Joining with `by = join_by(OTU)`
# }
```
