# Align the reference sequences of a phyloseq object

[![lifecycle-experimental](https://img.shields.io/badge/lifecycle-experimental-orange)](https://adrientaudiere.github.io/MiscMetabar/articles/Rules.html#lifecycle)

Build a multiple sequence alignment from the `refseq` slot of a
[`phyloseq-class`](https://rdrr.io/pkg/phyloseq/man/phyloseq-class.html)
object, or from a set of sequences given directly. Two backends are
available:

- `"decipher"` (default), the pure-R
  [`DECIPHER::AlignSeqs()`](https://rdrr.io/pkg/DECIPHER/man/AlignSeqs.html),
  which needs no external program;

- `"mafft"`, the MAFFT command-line aligner driven through
  [`ips::mafft()`](https://rdrr.io/pkg/ips/man/mafft.html), which is
  considerably faster on large sets of sequences and is the usual choice
  for metabarcoding-scale data.

The alignment is returned as a `DNAStringSet`, **not** as a phyloseq
object. Gap characters are valid DNA letters, so an alignment can
technically be written back into a `refseq` slot, but the sequences of
that slot are consumed as unaligned elsewhere in the pqverse
(clustering, BLAST, primer trimming), so MiscMetabar never does it for
you. Feed the returned alignment to
[`build_phytree_pq()`](https://adrientaudiere.github.io/MiscMetabar/dev/reference/build_phytree_pq.md)
or
[`Biostrings::writeXStringSet()`](https://rdrr.io/pkg/Biostrings/man/XStringSet-io.html)
instead.

## Usage

``` r
align_pq(
  x,
  method = c("decipher", "mafft"),
  exec = NULL,
  mafft_method = "auto",
  thread = -1,
  force = FALSE,
  verbose = FALSE,
  ...
)
```

## Arguments

- x:

  (required) A
  [`phyloseq-class`](https://rdrr.io/pkg/phyloseq/man/phyloseq-class.html)
  object with a non-empty `refseq` slot, a `DNAStringSet` (or any
  `XStringSet`), or a `DNAbin` object.

- method:

  One of `"decipher"` (default) or `"mafft"`.

- exec:

  Path to the MAFFT executable. Only used when `method = "mafft"`.
  Default to NULL, in which case MAFFT is looked up with
  [`is_mafft_installed()`](https://adrientaudiere.github.io/MiscMetabar/dev/reference/is_mafft_installed.md):
  first the `MiscMetabar.mafftpath` option, then the system `PATH`.

- mafft_method:

  The MAFFT strategy passed to
  [`ips::mafft()`](https://rdrr.io/pkg/ips/man/mafft.html), e.g.
  `"auto"` (default, lets MAFFT pick from the size of the input),
  `"localpair"` (L-INS-i), `"globalpair"` (G-INS-i) or `"retree 2"`.

- thread:

  Number of threads given to MAFFT. Default to -1, i.e. MAFFT detects
  the number of available cores itself.

- force:

  Logical, if TRUE sequences that are already all of the same length are
  realigned anyway. Default to FALSE, in which case such sequences are
  assumed to be aligned already and returned untouched.

- verbose:

  Logical, if TRUE report the aligner used and the width of the
  resulting alignment. Default to FALSE.

- ...:

  Additional arguments passed on to
  [`DECIPHER::AlignSeqs()`](https://rdrr.io/pkg/DECIPHER/man/AlignSeqs.html)
  or [`ips::mafft()`](https://rdrr.io/pkg/ips/man/mafft.html), depending
  on `method`.

## Value

A `DNAStringSet` of aligned sequences, all of the same width, in the
order of the input.

## Installing MAFFT

MAFFT is an external program and is not installed by MiscMetabar. It is
available from <https://mafft.cbrc.jp/alignment/software/>, and from the
usual package managers:

- Debian / Ubuntu: `sudo apt install mafft`

- macOS (Homebrew): `brew install mafft`

- conda: `conda install -c bioconda mafft`

- Windows: use the installer from the MAFFT website, or run it under
  WSL.

If the executable is not on the `PATH`, either pass its path with `exec`
or set it once for the session with
`options(MiscMetabar.mafftpath = "/path/to/mafft")`.

## See also

[`is_mafft_installed()`](https://adrientaudiere.github.io/MiscMetabar/dev/reference/is_mafft_installed.md),
[`build_phytree_pq()`](https://adrientaudiere.github.io/MiscMetabar/dev/reference/build_phytree_pq.md)

## Author

Adrien Taudière

## Examples

``` r
# \donttest{
data(data_fungi_mini)
df <- subset_taxa_pq(data_fungi_mini, taxa_sums(data_fungi_mini) > 12000)
#> Cleaning suppress 0 taxa (  ) and 7 sample(s) ( AD26-005-H_S10_MERGED.fastq.gz / CB8-019-H_S70_MERGED.fastq.gz / DY5-004-H_S97_MERGED.fastq.gz / F7-015-M_S106_MERGED.fastq.gz / N23-002-B_S130_MERGED.fastq.gz / NVABM0244-M_S137_MERGED.fastq.gz / T28-ABM602-B_S162_MERGED.fastq.gz ).
#> Number of non-matching ASV 0
#> Number of matching ASV 45
#> Number of filtered-out ASV 32
#> Number of kept ASV 13
#> Number of kept samples 130

# Pure-R alignment, no external program required
ali <- align_pq(df, method = "decipher", verbose = TRUE)
#> Determining distance matrix based on shared 8-mers:
#> ================================================================================
#> 
#> Time difference of 0 secs
#> 
#> Clustering into groups by similarity:
#> ================================================================================
#> 
#> Time difference of 0 secs
#> 
#> Aligning Sequences:
#> ================================================================================
#> 
#> Time difference of 0.07 secs
#> 
#> Iteration 1 of 2:
#> 
#> Determining distance matrix based on alignment:
#> ================================================================================
#> 
#> Time difference of 0 secs
#> 
#> Reclustering into groups by similarity:
#> ================================================================================
#> 
#> Time difference of 0 secs
#> 
#> Realigning Sequences:
#> ================================================================================
#> 
#> Time difference of 0.05 secs
#> 
#> Iteration 2 of 2:
#> 
#> Determining distance matrix based on alignment:
#> ================================================================================
#> 
#> Time difference of 0 secs
#> 
#> Reclustering into groups by similarity:
#> ================================================================================
#> 
#> Time difference of 0 secs
#> 
#> Realigning Sequences:
#> ================================================================================
#> 
#> Time difference of 0.03 secs
#> 
#> Refining the alignment:
#> ================================================================================
#> 
#> Time difference of 0.02 secs
#> 
#> ✔ 13 sequences aligned with decipher over 452 positions.
ali
#> DNAStringSet object of length 13:
#>      width seq                                              names               
#>  [1]   452 AAATGCGATAAGTAATGTGAATT...TGGGACTACCCGCTGAACTTA- ASV7
#>  [2]   452 AAATGCGATAAGTAATGTGAATT...TGGGACTACCCGCTGAACTTA- ASV8
#>  [3]   452 AAATGCGATAAGTAATGTGAATT...TAGGACTACCCGCTGAACTTA- ASV12
#>  [4]   452 AAATGCGATAAGTAATGTGAATT...TGGGACTACCCGCTGAACTTA- ASV18
#>  [5]   452 AAATGCGATAAGTAATGTGAATT...TAGGACTACCCGCTGAACTTA- ASV25
#>  ...   ... ...
#>  [9]   452 AATTGCGATAAGTAATGTGAATT...TGGGACTACCCGCTGAACTTA- ASV32
#> [10]   452 AAATGCGATAAGTAATGTGAATT...CAGGACTACCCGCTGAACTTA- ASV34
#> [11]   452 AAATGCGATAAGTAATGTGAATT...TAGGACTACCCGCTGAACTTA- ASV35
#> [12]   452 AAATGCGATAAGTAATGTGAATT...TAGGACTACCCGCTGAACTTA- ASV41
#> [13]   452 AAATGCGATAAGTAATGTGAATT...TAGGACTACCCGCTGAACTTA- ASV42
# }

if (FALSE) { # \dontrun{
# The same alignment with MAFFT, much faster on large refseq slots
ali_mafft <- align_pq(df, method = "mafft", verbose = TRUE)

# A more accurate (and slower) MAFFT strategy
ali_linsi <- align_pq(df, method = "mafft", mafft_method = "localpair")

# An explicit path to the executable
ali <- align_pq(df, method = "mafft", exec = "/usr/local/bin/mafft")

# Sequences can also be given directly
align_pq(phyloseq::refseq(df), method = "mafft")

# Feed the alignment to a tree builder or write it to a FASTA file
Biostrings::writeXStringSet(ali_mafft, "refseq_aligned.fasta")
} # }
```
