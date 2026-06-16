# Rarefy a phyloseq object to even sequencing depth

[![lifecycle-stable](https://img.shields.io/badge/lifecycle-stable-brightgreen)](https://adrientaudiere.github.io/MiscMetabar/articles/Rules.html#lifecycle)

An R-version-robust drop-in replacement for
[`phyloseq::rarefy_even_depth()`](https://rdrr.io/pkg/phyloseq/man/rarefy_even_depth.html),
heavily inspired by that function.

## Usage

``` r
rarefy_even_depth_pq(
  physeq,
  sample_size = NULL,
  rngseed = FALSE,
  replace = TRUE,
  trimOTUs = TRUE
)
```

## Arguments

- physeq:

  (required) a
  [`phyloseq-class`](https://rdrr.io/pkg/phyloseq/man/phyloseq-class.html)
  object obtained using the `phyloseq` package.

- sample_size:

  (integer) the sequencing depth to rarefy to. If `NULL` (default),
  `min(sample_sums(physeq))` is used. Samples with fewer reads than
  `sample_size` are dropped.

- rngseed:

  (logical or integer, default `FALSE`) random seed. Set to an integer
  to seed the RNG before subsampling and restore the caller's global RNG
  state on exit (matching
  [`phyloseq::rarefy_even_depth()`](https://rdrr.io/pkg/phyloseq/man/rarefy_even_depth.html)
  behaviour). `FALSE` leaves the global RNG untouched.

- replace:

  (logical, default `TRUE`) sample with replacement? `TRUE` matches the
  [`phyloseq::rarefy_even_depth()`](https://rdrr.io/pkg/phyloseq/man/rarefy_even_depth.html)
  default.

- trimOTUs:

  (logical, default `TRUE`) if `TRUE`, taxa that are entirely emptied by
  subsampling are removed from the returned object.

## Value

A new
[`phyloseq-class`](https://rdrr.io/pkg/phyloseq/man/phyloseq-class.html)
object with a rarefied `otu_table`.

## See also

[`phyloseq::rarefy_even_depth()`](https://rdrr.io/pkg/phyloseq/man/rarefy_even_depth.html),
[`rarefy_pq()`](https://adrientaudiere.github.io/MiscMetabar/reference/rarefy_pq.md)

## Author

Adrien Taudière

## Examples

``` r
data_f_rar <- rarefy_even_depth_pq(data_fungi_mini, sample_size = 500)
sample_sums(data_f_rar)
#>    A10-005-B_S188_MERGED.fastq.gz    A10-005-H_S189_MERGED.fastq.gz 
#>                               500                               500 
#>    A10-005-M_S190_MERGED.fastq.gz      A12-007_S191_MERGED.fastq.gz 
#>                               500                               500 
#>      A12-007-B_S2_MERGED.fastq.gz         A8-005_S4_MERGED.fastq.gz 
#>                               500                               500 
#>    AB29-ABMX-H_S6_MERGED.fastq.gz       AC27-013_S7_MERGED.fastq.gz 
#>                               500                               500 
#>    AD26-005-H_S10_MERGED.fastq.gz    AD26-005-M_S11_MERGED.fastq.gz 
#>                               500                               500 
#>   AD30-ABMX-M_S12_MERGED.fastq.gz    AD32-007-M_S13_MERGED.fastq.gz 
#>                               500                               500 
#>    ADABM30X-B_S14_MERGED.fastq.gz    ADABM30X-H_S15_MERGED.fastq.gz 
#>                               500                               500 
#>    ADABM30X-M_S16_MERGED.fastq.gz   AE30-ABM507_S17_MERGED.fastq.gz 
#>                               500                               500 
#>       B17-014_S18_MERGED.fastq.gz   BA16-036bis_S20_MERGED.fastq.gz 
#>                               500                               500 
#>    BA17-050-B_S21_MERGED.fastq.gz    BB19-006-H_S22_MERGED.fastq.gz 
#>                               500                               500 
#>     BB6-019-M_S25_MERGED.fastq.gz      BD14-021_S26_MERGED.fastq.gz 
#>                               500                               500 
#>    BJ17-007-M_S34_MERGED.fastq.gz   BJ8-ABM-003_S35_MERGED.fastq.gz 
#>                               500                               500 
#>      BN11-041_S39_MERGED.fastq.gz       BO8-002_S41_MERGED.fastq.gz 
#>                               500                               500 
#>    BP11-001-B_S43_MERGED.fastq.gz    BP11-001-H_S44_MERGED.fastq.gz 
#>                               500                               500 
#>    BP11-001-M_S45_MERGED.fastq.gz      BP14-006_S47_MERGED.fastq.gz 
#>                               500                               500 
#>     BQ4-018-B_S49_MERGED.fastq.gz     BQ4-018-H_S50_MERGED.fastq.gz 
#>                               500                               500 
#>    BQ9ABM-002_S52_MERGED.fastq.gz      BS14-006_S54_MERGED.fastq.gz 
#>                               500                               500 
#>    BV11-002-B_S57_MERGED.fastq.gz    BV11-002-H_S58_MERGED.fastq.gz 
#>                               500                               500 
#>        C1-001_S61_MERGED.fastq.gz        C9-005_S65_MERGED.fastq.gz 
#>                               500                               500 
#>         CA9-X_S68_MERGED.fastq.gz       CB9-013_S72_MERGED.fastq.gz 
#>                               500                               500 
#>       CC3-044_S73_MERGED.fastq.gz       D17-011_S77_MERGED.fastq.gz 
#>                               500                               500 
#>     D18-003-B_S78_MERGED.fastq.gz    D22-NVABM1_S81_MERGED.fastq.gz 
#>                               500                               500 
#>     D61-010-B_S82_MERGED.fastq.gz   DBM-ABM-001_S86_MERGED.fastq.gz 
#>                               500                               500 
#>    DP4-ABM001_S90_MERGED.fastq.gz  DS1-ABM002-H_S92_MERGED.fastq.gz 
#>                               500                               500 
#>     DU3-045-B_S94_MERGED.fastq.gz       DW4-007_S95_MERGED.fastq.gz 
#>                               500                               500 
#>   DZ6-ABM-001_S99_MERGED.fastq.gz  EA5-ABM-001_S103_MERGED.fastq.gz 
#>                               500                               500 
#>    EC2-013-B_S104_MERGED.fastq.gz   F6-ABM-001_S105_MERGED.fastq.gz 
#>                               500                               500 
#>     F7-015-M_S106_MERGED.fastq.gz    FOMES19-H_S108_MERGED.fastq.gz 
#>                               500                               500 
#>    FOMES19-M_S109_MERGED.fastq.gz   K26-NVABM1_S118_MERGED.fastq.gz 
#>                               500                               500 
#>       L19X-H_S120_MERGED.fastq.gz       L19X-M_S121_MERGED.fastq.gz 
#>                               500                               500 
#>    L23-002-H_S123_MERGED.fastq.gz       N19X-H_S127_MERGED.fastq.gz 
#>                               500                               500 
#>       N19X-M_S128_MERGED.fastq.gz     N25-ABMX_S133_MERGED.fastq.gz 
#>                               500                               500 
#>   NVABM-0058_S134_MERGED.fastq.gz      O20-X-H_S140_MERGED.fastq.gz 
#>                               500                               500 
#>      O20-X-M_S141_MERGED.fastq.gz     O9-005-B_S152_MERGED.fastq.gz 
#>                               500                               500 
#>     R25-ABMX_S157_MERGED.fastq.gz      T28-011_S161_MERGED.fastq.gz 
#>                               500                               500 
#>   U27-ABM002_S163_MERGED.fastq.gz      W30-006_S168_MERGED.fastq.gz 
#>                               500                               500 
#>    X24-009-H_S171_MERGED.fastq.gz    X24-009-M_S172_MERGED.fastq.gz 
#>                               500                               500 
#>      X24-010_S173_MERGED.fastq.gz    X29-004-M_S176_MERGED.fastq.gz 
#>                               500                               500 
#> Y21-ABM484-H_S177_MERGED.fastq.gz Y31-ABM484-B_S184_MERGED.fastq.gz 
#>                               500                               500 
#> Z30-ABM560-M_S187_MERGED.fastq.gz 
#>                               500 

# \donttest{
data_f_rar_notrim <- rarefy_even_depth_pq(
  data_fungi_mini,
  sample_size = 500,
  trimOTUs = FALSE
)
# }
```
