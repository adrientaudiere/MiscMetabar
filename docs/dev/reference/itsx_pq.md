# Extract the ITS region of the sequences of a phyloseq object with ITSx

[![lifecycle-experimental](https://img.shields.io/badge/lifecycle-experimental-orange)](https://adrientaudiere.github.io/MiscMetabar/articles/Rules.html#lifecycle)

Run ITSx (Bengtsson-Palme et al. 2013) on the `refseq` slot of a
phyloseq object and return the object with its `refseq` replaced by the
extracted region. ITSx locates the conserved flanks (end of the SSU,
5.8S, start of the LSU) with HMM profiles, so the region is cut where
the ribosomal genes actually are, whatever the primers.

Removing the conserved flanks matters beyond tidiness: they are nearly
identical in every reference record, so an alignment-based assignment
such as
[`assign_blastn()`](https://adrientaudiere.github.io/MiscMetabar/dev/reference/assign_blastn.md)
aligns a query that still carries them against almost the whole
database. On fungal ITS2 amplicons, extracting the ITS2 made blastn more
than ten times faster.

ITSx also reports, for each sequence it detects, the organism group
whose profiles fit it best (its *putative origin*). `add_origin` keeps
it as a column of the `tax_table`, and `remove_other_origin` drops the
taxa of another origin, a cheap filter against non-target sequences.

## Usage

``` r
itsx_pq(
  physeq,
  region = c("ITS1", "ITS2", "full"),
  organism_groups = "all",
  keep_undetected = TRUE,
  add_origin = TRUE,
  origin_col = "ITSx_origin",
  add_detected = FALSE,
  detected_col = "ITSx_detected",
  remove_other_origin = FALSE,
  keep_origin = "F",
  duplicated_seqs = c("merge", "remove"),
  nproc = 1,
  itsxpath = find_itsx(),
  args_before_itsx = "",
  verbose = TRUE
)
```

## Arguments

- physeq:

  (required) a
  [`phyloseq-class`](https://rdrr.io/pkg/phyloseq/man/phyloseq-class.html)
  object obtained using the `phyloseq` package.

- region:

  (default: "ITS1") The region to keep: `"ITS1"`, `"ITS2"` or `"full"`
  (ITS1 + 5.8S + ITS2, only for the sequences in which ITSx finds all
  three).

- organism_groups:

  (default: "all") The ITSx profile sets to search, passed to `-t`:
  `"all"`, or comma-separated one-letter codes such as `"F"` (fungi) or
  `"F,T"`. See Details for the codes.

- keep_undetected:

  (logical, default TRUE) If TRUE, the taxa in which ITSx does not find
  `region` keep their full sequence, so the set of taxa is unchanged. If
  FALSE, they are removed.

- add_origin:

  (logical, default TRUE) If TRUE, a column `origin_col` is added to the
  `tax_table` with the putative origin ITSx gives to each taxon (e.g.
  `"Fungi"`), NA when ITSx detects nothing in the sequence.

- origin_col:

  (default: "ITSx_origin") Name of that column.

- add_detected:

  (logical, default FALSE) If TRUE, a column `detected_col` is added to
  the `tax_table`: `"TRUE"` when ITSx found `region` in the taxon (its
  sequence is the extracted region), `"FALSE"` when it did not (with
  `keep_undetected = TRUE`, its sequence is the original one). After
  merging, the value of the taxon kept.

- detected_col:

  (default: "ITSx_detected") Name of that column.

- remove_other_origin:

  (logical, default FALSE) If TRUE, the taxa whose putative origin is
  not in `keep_origin` are removed, including the taxa in which ITSx
  detects nothing.

- keep_origin:

  (default: "F") One-letter codes of the origins kept when
  `remove_other_origin` is TRUE.

- duplicated_seqs:

  (default: "merge") What to do with taxa whose extracted sequences are
  identical (they differed only in the removed flanks): `"merge"` merges
  them with
  [`merge_taxa_vec()`](https://adrientaudiere.github.io/MiscMetabar/dev/reference/merge_taxa_vec.md)
  into the most abundant one (its taxonomy is kept), `"remove"` keeps
  the most abundant and removes the others (ties: the first one).

- nproc:

  (default: 1) Number of CPUs given to ITSx (`--cpu`).

- itsxpath:

  (default:
  [`find_itsx()`](https://adrientaudiere.github.io/MiscMetabar/dev/reference/find_itsx.md))
  Path to ITSx.

- args_before_itsx:

  (String, default "") A one line bash command run before ITSx, e.g. the
  conda activation returned by
  [`install_itsx()`](https://adrientaudiere.github.io/MiscMetabar/dev/reference/install_itsx.md).

- verbose:

  (logical, default TRUE) If TRUE, report how many taxa were detected,
  removed or merged.

## Value

A new
[`phyloseq-class`](https://rdrr.io/pkg/phyloseq/man/phyloseq-class.html)
object whose `refseq` holds the extracted region.

## Details

Origin codes of ITSx 1.1.3: A Alveolata, B Bryophyta, C Bacillariophyta,
D Amoebozoa, E Euglenozoa, F Fungi, G Chlorophyta, H Rhodophyta, I
Phaeophyceae, L Marchantiophyta, M Metazoa, N Microsporidia, O Oomycota,
P Haptophyceae, Q Raphidophyceae, R Rhizaria, S Synurophyceae, T
Tracheophyta, U Eustigmatophyceae, X Apusozoa, Y Parabasalia.

The putative origin is the profile set with the best score; on short or
partial sequences (e.g. an ITS1 amplicon ending inside the 5.8S) it can
be wrong, so check it before using `remove_other_origin`.

## References

Bengtsson-Palme J, Ryberg M, Hartmann M, et al. (2013). Improved
software detection and extraction of ITS1 and ITS2 from ribosomal ITS
sequences of fungi and other eukaryotes for analysis of environmental
sequencing data. *Methods in Ecology and Evolution* 4: 914–919.
[doi:10.1111/2041-210X.12073](https://doi.org/10.1111/2041-210X.12073)

## See also

[`install_itsx()`](https://adrientaudiere.github.io/MiscMetabar/dev/reference/install_itsx.md),
[`is_itsx_installed()`](https://adrientaudiere.github.io/MiscMetabar/dev/reference/is_itsx_installed.md),
[`cutadapt_remove_primers()`](https://adrientaudiere.github.io/MiscMetabar/dev/reference/cutadapt_remove_primers.md)

## Author

Adrien Taudière

## Examples

``` r
if (FALSE) { # \dontrun{
prelude <- install_itsx()
d_its1 <- itsx_pq(data_fungi_mini, region = "ITS1", args_before_itsx = prelude)
table(d_its1@tax_table[, "ITSx_origin"], useNA = "ifany")
} # }
```
