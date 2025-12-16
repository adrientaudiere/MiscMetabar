# Format a fasta database in sintax format

[![lifecycle-experimental](https://img.shields.io/badge/lifecycle-experimental-orange)](https://adrientaudiere.github.io/MiscMetabar/articles/Rules.html#lifecycle)

Only tested with Unite and Eukaryome fasta file for the moment. Rely on
the presence of the pattern pattern_tax default "k\_\_" to format the
header.

A reference database in sintax format contain taxonomic information in
the header of each sequence in the form of a string starting with
";tax=" and followed by a comma-separated list of up to nine taxonomic
identifiers. Each taxonomic identifier must start with an indication of
the rank by one of the letters d (for domain) k (kingdom), p (phylum), c
(class), o (order), f (family), g (genus), s (species), or t (strain).
The letter is followed by a colon (:) and the name of that rank. Commas
and semicolons are not allowed in the name of the rank. Non-ascii
characters should be avoided in the names.

Example:

\\X80725_S000004313;tax=d:Bacteria,p:Proteobacteria,c:Gammaproteobacteria,o:Enterobacteriales,f:Enterobacteriaceae,g:Escherichia/Shigella,s:Escherichia_coli,t:str.\_K-12_substr.\_MG1655

## Usage

``` r
format2sintax(
  fasta_db = NULL,
  taxnames = NULL,
  pattern_tax = "k__",
  pattern_sintax = "tax=k:",
  output_path = NULL
)
```

## Arguments

- fasta_db:

  A link to a fasta files

- taxnames:

  A list of names to format. You must specify either fasta_db OR
  taxnames, not both.

- pattern_tax:

  (default "k\_\_") The pattern to replace by pattern_sintax.

- pattern_sintax:

  (default "tax=k:") Useless for most users. Sometimes you may want to
  replacte by "tax=d:" (d for domain instead of kingdom).

- output_path:

  (optional) A path to an output fasta files. Only used if fasta_db is
  set.

## Value

Either an object of class DNAStringSet or a vector of reformated names

## See also

[`format2dada2_species()`](https://adrientaudiere.github.io/MiscMetabar/reference/format2dada2_species.md),
[`format2dada2()`](https://adrientaudiere.github.io/MiscMetabar/reference/format2dada2.md)

## Author

Adrien Taudière
