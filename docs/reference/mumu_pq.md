# MUMU reclustering of class `physeq`

[![lifecycle-experimental](https://img.shields.io/badge/lifecycle-experimental-orange)](https://adrientaudiere.github.io/MiscMetabar/articles/Rules.html#lifecycle)

See https://www.nature.com/articles/s41467-017-01312-x for more
information on the original method LULU. This is a wrapper of
[mumu](https://github.com/frederic-mahe/mumu) a C++ re-implementation of
LULU by Frédéric Mahé

## Usage

``` r
mumu_pq(
  physeq,
  nproc = 1,
  id = 0.84,
  vsearchpath = "vsearch",
  mumupath = "mumu",
  verbose = FALSE,
  clean_pq = TRUE,
  keep_temporary_files = FALSE
)
```

## Arguments

- physeq:

  (required) a
  [`phyloseq-class`](https://rdrr.io/pkg/phyloseq/man/phyloseq-class.html)
  object obtained using the `phyloseq` package.

- nproc:

  (default 1) Set to number of cpus/processors to use for the clustering

- id:

  (default: 0.84) id for –usearch_global.

- vsearchpath:

  (default: vsearch) path to vsearch.

- mumupath:

  path to mumu. See [mumu](https://github.com/frederic-mahe/mumu) for
  installation instruction

- verbose:

  (logical) If true, print some additional messages.

- clean_pq:

  (logical) If true, empty samples and empty ASV are discarded before
  clustering.

- keep_temporary_files:

  (logical, default: FALSE) Do we keep temporary files

## Value

a list of for object

- "new_physeq": The new phyloseq object (class physeq)

- "mumu_results": The log file of the mumu software. Run `man mumu` into
  bash to obtain details about columns' signification.

## Details

This function is mainly a wrapper of the work of others. Please cite
[mumu](https://github.com/frederic-mahe/mumu/blob/main/CITATION.cff) and
[lulu](https://www.nature.com/articles/s41467-017-01312-x) if you use
this function for your work.

## References

- MUMU: <https://github.com/frederic-mahe/mumu>

- VSEARCH can be downloaded from <https://github.com/torognes/vsearch>.

## See also

[`lulu_pq()`](https://adrientaudiere.github.io/MiscMetabar/reference/lulu_pq.md)

## Author

Frédéric Mahé & Adrien Taudière <adrien.taudiere@zaclys.net>

## Examples

``` r
if (FALSE) { # \dontrun{
mumu_pq(data_fungi_sp_known)
} # }
```
