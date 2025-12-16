# Show taxa which are present in only one given level of a modality

[![lifecycle-experimental](https://img.shields.io/badge/lifecycle-experimental-orange)](https://adrientaudiere.github.io/MiscMetabar/articles/Rules.html#lifecycle)

[![lifecycle-experimental](https://img.shields.io/badge/lifecycle-experimental-orange)](https://adrientaudiere.github.io/MiscMetabar/articles/Rules.html#lifecycle)

Given one modality name in sam_data and one level of the modality,
return the taxa strictly specific of this level.

## Usage

``` r
taxa_only_in_one_level(
  physeq,
  modality,
  level,
  min_nb_seq_taxa = 0,
  min_nb_samples_taxa = 0
)

taxa_only_in_one_level(
  physeq,
  modality,
  level,
  min_nb_seq_taxa = 0,
  min_nb_samples_taxa = 0
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

- level:

  (required) The level (must be present in modality) of interest

- min_nb_seq_taxa:

  (default 0 = no filter) The minimum number of sequences per taxa

- min_nb_samples_taxa:

  (default 0 = no filter) The minimum number of samples per taxa

## Value

A vector of taxa names

A vector of taxa names

## Author

Adrien Taudière

## Examples

``` r
data_fungi_mini_woNA4height <- subset_samples(
  data_fungi_mini,
  !is.na(data_fungi_mini@sam_data$Height)
)
taxa_only_in_one_level(data_fungi_mini_woNA4height, "Height", "High")
#> Cleaning suppress 3 taxa and 0 samples.
#> Cleaning suppress 0 taxa (  ) and 0 sample(s) (  ).
#> Number of non-matching ASV 0
#> Number of matching ASV 42
#> Number of filtered-out ASV 35
#> Number of kept ASV 7
#> Number of kept samples 3
#> Cleaning suppress 3 taxa and 0 samples.
#> [1] "ASV48" "ASV50" "ASV77" "ASV93"
#' # Taxa present only in low height samples
suppressMessages(suppressWarnings(
  taxa_only_in_one_level(data_fungi, "Height", "Low")
))
#>   [1] "ASV54"   "ASV84"   "ASV143"  "ASV190"  "ASV202"  "ASV229"  "ASV247" 
#>   [8] "ASV253"  "ASV254"  "ASV274"  "ASV307"  "ASV311"  "ASV318"  "ASV330" 
#>  [15] "ASV388"  "ASV394"  "ASV408"  "ASV412"  "ASV422"  "ASV433"  "ASV442" 
#>  [22] "ASV457"  "ASV476"  "ASV487"  "ASV491"  "ASV495"  "ASV549"  "ASV554" 
#>  [29] "ASV560"  "ASV571"  "ASV604"  "ASV659"  "ASV664"  "ASV672"  "ASV679" 
#>  [36] "ASV682"  "ASV683"  "ASV703"  "ASV710"  "ASV713"  "ASV716"  "ASV730" 
#>  [43] "ASV739"  "ASV751"  "ASV753"  "ASV777"  "ASV804"  "ASV808"  "ASV814" 
#>  [50] "ASV827"  "ASV829"  "ASV862"  "ASV867"  "ASV900"  "ASV943"  "ASV946" 
#>  [57] "ASV957"  "ASV993"  "ASV1033" "ASV1041" "ASV1045" "ASV1105" "ASV1115"
#>  [64] "ASV1116" "ASV1117" "ASV1124" "ASV1158" "ASV1166" "ASV1184" "ASV1197"
#>  [71] "ASV1199" "ASV1214" "ASV1218" "ASV1228" "ASV1242" "ASV1262" "ASV1265"
#>  [78] "ASV1298" "ASV1307" "ASV1364" "ASV1374" "ASV1378" "ASV1381" "ASV1384"
#>  [85] "ASV1405" "ASV1416" "ASV1430" "ASV1450" "ASV1454" "ASV1470" "ASV1471"
#>  [92] "ASV1480" "ASV1485" "ASV1498" "ASV1524" "ASV1527" "ASV1537" "ASV1542"
#>  [99] "ASV1559" "ASV1563" "ASV1569" "ASV1575" "ASV1576" "ASV1584" "ASV1594"
#> [106] "ASV1600" "ASV1611" "ASV1624" "ASV1627" "ASV1632" "ASV1645" "ASV1653"
#> [113] "ASV1654" "ASV1664" "ASV1665" "ASV1667" "ASV1670" "ASV1672" "ASV1680"
#> [120] "ASV1681" "ASV1704" "ASV1705" "ASV1716" "ASV1726"
# Number of taxa present only in sample of time equal to 15
suppressMessages(suppressWarnings(
  length(taxa_only_in_one_level(data_fungi, "Time", "15"))
))
#> [1] 126
data_fungi_mini_woNA4height <- subset_samples(
  data_fungi_mini,
  !is.na(data_fungi_mini@sam_data$Height)
)
taxa_only_in_one_level(data_fungi_mini_woNA4height, "Height", "High")
#> Cleaning suppress 3 taxa and 0 samples.
#> Cleaning suppress 0 taxa (  ) and 0 sample(s) (  ).
#> Number of non-matching ASV 0
#> Number of matching ASV 42
#> Number of filtered-out ASV 35
#> Number of kept ASV 7
#> Number of kept samples 3
#> Cleaning suppress 3 taxa and 0 samples.
#> [1] "ASV48" "ASV50" "ASV77" "ASV93"
#' # Taxa present only in low height samples
suppressMessages(suppressWarnings(
  taxa_only_in_one_level(data_fungi, "Height", "Low")
))
#>   [1] "ASV54"   "ASV84"   "ASV143"  "ASV190"  "ASV202"  "ASV229"  "ASV247" 
#>   [8] "ASV253"  "ASV254"  "ASV274"  "ASV307"  "ASV311"  "ASV318"  "ASV330" 
#>  [15] "ASV388"  "ASV394"  "ASV408"  "ASV412"  "ASV422"  "ASV433"  "ASV442" 
#>  [22] "ASV457"  "ASV476"  "ASV487"  "ASV491"  "ASV495"  "ASV549"  "ASV554" 
#>  [29] "ASV560"  "ASV571"  "ASV604"  "ASV659"  "ASV664"  "ASV672"  "ASV679" 
#>  [36] "ASV682"  "ASV683"  "ASV703"  "ASV710"  "ASV713"  "ASV716"  "ASV730" 
#>  [43] "ASV739"  "ASV751"  "ASV753"  "ASV777"  "ASV804"  "ASV808"  "ASV814" 
#>  [50] "ASV827"  "ASV829"  "ASV862"  "ASV867"  "ASV900"  "ASV943"  "ASV946" 
#>  [57] "ASV957"  "ASV993"  "ASV1033" "ASV1041" "ASV1045" "ASV1105" "ASV1115"
#>  [64] "ASV1116" "ASV1117" "ASV1124" "ASV1158" "ASV1166" "ASV1184" "ASV1197"
#>  [71] "ASV1199" "ASV1214" "ASV1218" "ASV1228" "ASV1242" "ASV1262" "ASV1265"
#>  [78] "ASV1298" "ASV1307" "ASV1364" "ASV1374" "ASV1378" "ASV1381" "ASV1384"
#>  [85] "ASV1405" "ASV1416" "ASV1430" "ASV1450" "ASV1454" "ASV1470" "ASV1471"
#>  [92] "ASV1480" "ASV1485" "ASV1498" "ASV1524" "ASV1527" "ASV1537" "ASV1542"
#>  [99] "ASV1559" "ASV1563" "ASV1569" "ASV1575" "ASV1576" "ASV1584" "ASV1594"
#> [106] "ASV1600" "ASV1611" "ASV1624" "ASV1627" "ASV1632" "ASV1645" "ASV1653"
#> [113] "ASV1654" "ASV1664" "ASV1665" "ASV1667" "ASV1670" "ASV1672" "ASV1680"
#> [120] "ASV1681" "ASV1704" "ASV1705" "ASV1716" "ASV1726"
# Number of taxa present only in sample of time equal to 15
suppressMessages(suppressWarnings(
  length(taxa_only_in_one_level(data_fungi, "Time", "15"))
))
#> [1] 126
```
