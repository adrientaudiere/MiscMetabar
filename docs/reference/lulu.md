# Post Clustering Curation of Amplicon Data.

[![lifecycle-stable](https://img.shields.io/badge/lifecycle-stable-green)](https://adrientaudiere.github.io/MiscMetabar/articles/Rules.html#lifecycle)

The original function and documentation was written by Tobias Guldberg
Frøslev in the [lulu](https://github.com/tobiasgf/lulu) package.

This algorithm `lulu` consumes an OTU table and a matchlist, and
evaluates cooccurence of 'daughters' (potential analytical artefacts)
and their 'parents' (~= real biological species/OTUs). The algorithm
requires an OTU table (species/site matrix), and a match list. The OTU
table can be made with various r-packages (e.g. `DADA2`) or external
pipelines (`VSEARCH, USEARCH, QIIME`, etc.), and the match-list can be
made with external bioinformatic tools like `VSEARCH, USEARCH, BLASTN`
or another algorithm for pair-wise sequence matching.

## Usage

``` r
lulu(
  otu_table,
  matchlist,
  minimum_ratio_type = "min",
  minimum_ratio = 1,
  minimum_match = 84,
  minimum_relative_cooccurence = 0.95,
  progress_bar = TRUE,
  log_conserved = FALSE
)
```

## Arguments

- otu_table:

  a data.frame with with an OTU table that has sites/samples as columns
  and OTUs (unique OTU id's) as rows, and observations as read counts.

- matchlist:

  a data.frame containing three columns: (1) OTU id of potential
  child, (2) OTU id of potential parent, (3) match - % identiti between
  the sequences of the potential parent and potential child OTUs. **NB:
  The matchlist is the product of a mapping of OTU sequences against
  each other. This is currently carried out by an external script in
  e.g. Blastn or VSEARCH, prior to running lulu!**

- minimum_ratio_type:

  sets whether a potential error must have lower abundance than the
  parent in all samples `min` (default), or if an error just needs to
  have lower abundance on average `avg`. Choosing lower abundance on
  average over globally lower abundance will greatly increase the number
  of designated errors. This option was introduced to make it possible
  to account for non-sufficiently clustered intraspecific variation, but
  is not generally recommended, as it will also increase the potential
  of cluster well-separated, but co-occuring, sequence similar species.

- minimum_ratio:

  sets the minimim abundance ratio between a potential error and a
  potential parent to be identified as an error. If the
  `minimum_ratio_type` is set to `min` (default), the `minimum_ratio`
  applies to the lowest observed ration across the samples. If the
  `minimum_ratio_type` is set to `avg` (default), the `minimum_ratio`
  applies to the mean of observed ration across the samples.`avg`.
  (default is 1).

- minimum_match:

  minimum threshold of sequence similarity for considering any OTU as an
  error of another can be set (default 84%).

- minimum_relative_cooccurence:

  minimum co-occurrence rate, i.e. the lower rate of occurrence of the
  potential error explained by co-occurrence with the potential parent
  for considering error state.

- progress_bar:

  (Logical, default TRUE) print progress during the calculation or not.

- log_conserved:

  (Logical, default FALSE) conserved log files writed in the disk

## Value

Function `lulu` returns a list of results based on the input OTU table
and match list.

- `curated_table` - a curated OTU table with daughters merged with their
  matching parents.

- `curated_count` - number of curated (parent) OTUs.

- `curated_otus` - ids of the OTUs that were accepted as valid OTUs.

- `discarded_count` - number of discarded (merged with parent) OTUs.

- `discarded_otus` - ids of the OTUs that were identified as errors
  (daughters) and merged with respective parents.

- `runtime` - time used by the script.

- `minimum_match` - the id threshold (minimum match \\ by user).

- `minimum_relative_cooccurence` - minimum ratio of daughter-occurences
  explained by co-occurence with parent (set by user).

- `otu_map` - information of which daughters were mapped to which
  parents.

- `original_table` - original OTU table.

The matchlist is the product of a mapping of OTU sequences against each
other. This is currently carried out by an external script in e.g.
BLASTN or VSEARCH, prior to running `lulu`! Producing the match list
requires a file with all the OTU sequences (centroids) - e.g.
`OTUcentroids.fasta`. The matchlist can be produced by mapping all OTUs
against each other with an external algorithm like VSEARCH or BLASTN. In
`VSEARCH` a matchlist can be produced e.g. with the following command:
`vsearch --usearch_global OTUcentroids.fasta --db OTUcentroids.fasta --strand plus --self --id .80 --iddef 1 --userout matchlist.txt --userfields query+target+id --maxaccepts 0 --query_cov .9 --maxhits 10`.
In `BLASTN` a matchlist can be produces e.g. with the following
commands. First we produce a blast-database from the fasta file:
`makeblastdb -in OTUcentroids.fasta -parse_seqids -dbtype nucl`, then we
match the centroids against that database:
`blastn -db OTUcentoids.fasta -num_threads 10 -outfmt'6 qseqid sseqid pident' -out matchlist.txt -qcov_hsp_perc .90 -perc_identity .84 -query OTUcentroids.fasta`

## Details

Please cite the lulu original paper:
https://www.nature.com/articles/s41467-017-01312-x

## Author

Tobias Guldberg Frøslev (orcid:
[0000-0002-3530-013X](https://orcid.org/0000-0002-3530-013X)), modified
by Adrien Taudière
