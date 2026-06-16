# Search for a list of sequence in an object to remove chimera taxa using [vsearch](https://github.com/torognes/vsearch)

[![lifecycle-experimental](https://img.shields.io/badge/lifecycle-experimental-orange)](https://adrientaudiere.github.io/MiscMetabar/articles/Rules.html#lifecycle)

Use the VSEARCH software.

## Usage

``` r
chimera_removal_vs(object, type = "Discard_only_chim", clean_pq = FALSE, ...)
```

## Arguments

- object:

  (required) A phyloseq-class object or one of dada, derep, data.frame
  or list coercible to sequences table using the function
  [`dada2::makeSequenceTable()`](https://rdrr.io/pkg/dada2/man/makeSequenceTable.html)

- type:

  (default "Discard_only_chim"). The type define the type of filtering.

  - "Discard_only_chim" will only discard taxa classify as chimera by
    vsearch

  - "Select_only_non_chim_seqlen_filtered" will only select taxa
    classify as non-chimera by vsearch(after filtering taxa based on
    their sequence length by the parameter `min_seq_length` from the
    [`chimera_detection_vs()`](https://adrientaudiere.github.io/MiscMetabar/reference/chimera_detection_vs.md)
    function)

  - "Select_only_chim" will only select taxa classify as chimera by
    vsearch (after filtering taxa based on their sequence length by the
    parameter `min_seq_length` from the
    [`chimera_detection_vs()`](https://adrientaudiere.github.io/MiscMetabar/reference/chimera_detection_vs.md)
    function)

- clean_pq:

  (logical; default FALSE) If TRUE, return the phyloseq object after
  cleaning using the default parameter of
  [`clean_pq()`](https://adrientaudiere.github.io/MiscMetabar/reference/clean_pq.md)
  function.

- ...:

  Additional arguments passed on to
  [`chimera_detection_vs()`](https://adrientaudiere.github.io/MiscMetabar/reference/chimera_detection_vs.md)
  function

## Value

- I/ a sequences tables if object is of class dada, derep, data.frame or
  list.

- II/ a phyloseq object without (or with if type = 'Select_only_chim')
  chimeric taxa

## Details

This function is mainly a wrapper of the work of others. Please cite
[vsearch](https://github.com/torognes/vsearch).

## See also

[`chimera_detection_vs()`](https://adrientaudiere.github.io/MiscMetabar/reference/chimera_detection_vs.md),
[`dada2::removeBimeraDenovo()`](https://rdrr.io/pkg/dada2/man/removeBimeraDenovo.html)

## Author

Adrien Taudière

## Examples

``` r
# \donttest{
data_fungi_nochim <- chimera_removal_vs(data_fungi)
#> Filtering for sequences under 100 bp remove a total of 0 ( 0 %) unique sequences for a total of 0 sequences removed ( 0 %)
#> Cleaning suppress 0 taxa (  ) and 0 sample(s) (  ).
#> Number of non-matching ASV 0
#> Number of matching ASV 1420
#> Number of filtered-out ASV 240
#> Number of kept ASV 1180
#> Number of kept samples 185
# }
if (FALSE) { # \dontrun{
# Adding a chimeric sequence for the example
data_fungi_with_chim <- data_fungi
data_fungi_with_chim@refseq["ASV1710"] <- Biostrings::xscat(
  Biostrings::subseq(data_fungi_with_chim@refseq[1], start = 1, end = 150),
  Biostrings::subseq(data_fungi_with_chim@refseq[4], start = 151, end = 300)
)
data_fungi_nochim <- chimera_removal_vs(data_fungi)

# Higher value of abskew parameter is less stringent
data_fungi_nochim_16 <- chimera_removal_vs(data_fungi,
  abskew = 16, min_seq_length = 10
)

# Potential Chimeric ASVs detected by vsearch
chim_asv <- taxa_names(data_fungi_with_chim)[!taxa_names(data_fungi_with_chim)
%in% taxa_names(data_fungi_nochim)]
"ASV1710" %in% chim_asv
track_wkflow(list(data_fungi_with_chim, data_fungi_nochim))

data_fungi_nochim2 <-
  chimera_removal_vs(data_fungi, type = "Select_only_non_chim_seqlen_filtered")
data_fungi_chimera <-
  chimera_removal_vs(data_fungi, type = "Select_only_chim")
} # }
```
