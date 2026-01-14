# Subset taxa using a conditional named boolean vector.

[![lifecycle-experimental](https://img.shields.io/badge/lifecycle-experimental-orange)](https://adrientaudiere.github.io/MiscMetabar/articles/Rules.html#lifecycle)

The main objective of this function is to complete the
[`phyloseq::subset_taxa()`](https://rdrr.io/pkg/phyloseq/man/subset_taxa-methods.html)
function by propose a more easy way of subset_taxa using a named boolean
vector. Names must match taxa_names.

## Usage

``` r
subset_taxa_pq(
  physeq,
  condition,
  verbose = TRUE,
  clean_pq = TRUE,
  taxa_names_from_physeq = FALSE
)
```

## Arguments

- physeq:

  (required) a
  [`phyloseq-class`](https://rdrr.io/pkg/phyloseq/man/phyloseq-class.html)
  object obtained using the `phyloseq` package.

- condition:

  A named boolean vector to subset taxa. Length must fit the number of
  taxa and names must match taxa_names. Can also be a condition using a
  column of the tax_table slot (see examples). If the order of condition
  is the same as taxa_names(physeq), you can use the parameter
  `taxa_names_from_physeq = TRUE`.

- verbose:

  (logical) Informations are printed

- clean_pq:

  (logical) If set to TRUE, empty samples are discarded after subsetting
  ASV

- taxa_names_from_physeq:

  (logical) If set to TRUE, rename the condition vector using
  taxa_names(physeq). Carefully check the result of this function if you
  use this parameter. No effect if the condition is of class
  `tax_table`.

## Value

a new phyloseq object

## Author

Adrien Taudière

## Examples

``` r
subset_taxa_pq(data_fungi, data_fungi@tax_table[, "Phylum"] == "Ascomycota")
#> Cleaning suppress 0 taxa (  ) and 0 sample(s) (  ).
#> Number of non-matching ASV 0
#> Number of matching ASV 1420
#> Number of filtered-out ASV 354
#> Number of kept ASV 1066
#> Number of kept samples 185
#> phyloseq-class experiment-level object
#> otu_table()   OTU Table:         [ 1066 taxa and 185 samples ]
#> sample_data() Sample Data:       [ 185 samples by 7 sample variables ]
#> tax_table()   Taxonomy Table:    [ 1066 taxa by 12 taxonomic ranks ]
#> refseq()      DNAStringSet:      [ 1066 reference sequences ]
subset_taxa_pq(data_fungi, taxa_sums(data_fungi) > 100)
#> Cleaning suppress 0 taxa (  ) and 0 sample(s) (  ).
#> Number of non-matching ASV 0
#> Number of matching ASV 1420
#> Number of filtered-out ASV 342
#> Number of kept ASV 1078
#> Number of kept samples 185
#> phyloseq-class experiment-level object
#> otu_table()   OTU Table:         [ 1078 taxa and 185 samples ]
#> sample_data() Sample Data:       [ 185 samples by 7 sample variables ]
#> tax_table()   Taxonomy Table:    [ 1078 taxa by 12 taxonomic ranks ]
#> refseq()      DNAStringSet:      [ 1078 reference sequences ]

cond_taxa <- grepl("Endophyte", data_fungi@tax_table[, "Guild"])
names(cond_taxa) <- taxa_names(data_fungi)
subset_taxa_pq(data_fungi, cond_taxa)
#> Cleaning suppress 0 taxa (  ) and 9 sample(s) ( A10-005-B_S188_MERGED.fastq.gz / A10-005-M_S190_MERGED.fastq.gz / BE9-006-H_S28_MERGED.fastq.gz / C21-NV1-M_S64_MERGED.fastq.gz / CA12-024_S66_MERGED.fastq.gz / H10-018-M_S110_MERGED.fastq.gz / L23-002-H_S123_MERGED.fastq.gz / O24-003-H_S146_MERGED.fastq.gz / O26-004-M_S150_MERGED.fastq.gz ).
#> Number of non-matching ASV 0
#> Number of matching ASV 1420
#> Number of filtered-out ASV 1292
#> Number of kept ASV 128
#> Number of kept samples 176
#> phyloseq-class experiment-level object
#> otu_table()   OTU Table:         [ 128 taxa and 176 samples ]
#> sample_data() Sample Data:       [ 176 samples by 7 sample variables ]
#> tax_table()   Taxonomy Table:    [ 128 taxa by 12 taxonomic ranks ]
#> refseq()      DNAStringSet:      [ 128 reference sequences ]

subset_taxa_pq(data_fungi, grepl("mycor", data_fungi@tax_table[, "Guild"]),
  taxa_names_from_physeq = TRUE
)
#> Cleaning suppress 0 taxa (  ) and 47 sample(s) ( A10-005-H_S189_MERGED.fastq.gz / A10-005-M_S190_MERGED.fastq.gz / A15-004_S3_MERGED.fastq.gz / A8-005_S4_MERGED.fastq.gz / BA16-036bis_S20_MERGED.fastq.gz / BB6-019-M_S25_MERGED.fastq.gz / BE9-006-B_S27_MERGED.fastq.gz / BE9-006-H_S28_MERGED.fastq.gz / BP11-001-H_S44_MERGED.fastq.gz / BP11-001-M_S45_MERGED.fastq.gz / BQ3-019_S48_MERGED.fastq.gz / BQ4-018-H_S50_MERGED.fastq.gz / C21-NV1-B_S62_MERGED.fastq.gz / C21-NV1-M_S64_MERGED.fastq.gz / C9-005_S65_MERGED.fastq.gz / CA9-027_S67_MERGED.fastq.gz / CB8-019-B_S69_MERGED.fastq.gz / CB8-019-M_S71_MERGED.fastq.gz / CC8-003_S74_MERGED.fastq.gz / D17-011_S77_MERGED.fastq.gz / D18-003-B_S78_MERGED.fastq.gz / D18-003-M_S80_MERGED.fastq.gz / DS1-ABM002-H_S92_MERGED.fastq.gz / DZ6-ABM-001_S99_MERGED.fastq.gz / F6-ABM-001_S105_MERGED.fastq.gz / F7-015-M_S106_MERGED.fastq.gz / H10-018-M_S110_MERGED.fastq.gz / N22-001-B_S129_MERGED.fastq.gz / N23-002-M_S132_MERGED.fastq.gz / N25-ABMX_S133_MERGED.fastq.gz / NVABM-0058_S134_MERGED.fastq.gz / NVABM0244-M_S137_MERGED.fastq.gz / O21-007-H_S143_MERGED.fastq.gz / O21-007-M_S144_MERGED.fastq.gz / O24-003-B_S145_MERGED.fastq.gz / O24-003-H_S146_MERGED.fastq.gz / O24-003-M_S147_MERGED.fastq.gz / O26-004-M_S150_MERGED.fastq.gz / P27-ABM001_S155_MERGED.fastq.gz / W26-001-M_S167_MERGED.fastq.gz / X24-009-B_S170_MERGED.fastq.gz / X24-009-H_S171_MERGED.fastq.gz / X24-009-M_S172_MERGED.fastq.gz / X24-010_S173_MERGED.fastq.gz / Y28-002-H_S179_MERGED.fastq.gz / Y28-002-M_S180_MERGED.fastq.gz / Z30-002_S186_MERGED.fastq.gz ).
#> Number of non-matching ASV 0
#> Number of matching ASV 1420
#> Number of filtered-out ASV 1376
#> Number of kept ASV 44
#> Number of kept samples 138
#> phyloseq-class experiment-level object
#> otu_table()   OTU Table:         [ 44 taxa and 138 samples ]
#> sample_data() Sample Data:       [ 138 samples by 7 sample variables ]
#> tax_table()   Taxonomy Table:    [ 44 taxa by 12 taxonomic ranks ]
#> refseq()      DNAStringSet:      [ 44 reference sequences ]
```
