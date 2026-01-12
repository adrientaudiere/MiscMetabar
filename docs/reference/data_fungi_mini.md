# Fungal OTU in phyloseq format

It is a subset of the data_fungi dataset including only Basidiomycota
with more than 5000 sequences.

## Usage

``` r
data(data_fungi_mini)

data(data_fungi_mini)
```

## Format

A physeq object containing 45 taxa with references sequences described
by 14 taxonomic ranks and 137 samples described by 7 sample variables:

- *X*: the name of the fastq-file

- *Sample_names*: the names of ... the samples

- *Treename*: the name of an tree

- *Sample_id*: identifier for each sample

- *Height*: height of the sample in the tree

- *Diameter*: diameter of the trunk

- *Time*: time since the dead of the tree

A physeq object containing 45 taxa with references sequences described
by 14 taxonomic ranks and 137 samples described by 7 sample variables:

- *X*: the name of the fastq-file

- *Sample_names*: the names of ... the samples

- *Treename*: the name of an tree

- *Sample_id*: identifier for each sample

- *Height*: height of the sample in the tree

- *Diameter*: diameter of the trunk

- *Time*: time since the dead of the tree

## Details

Obtain using
`data_fungi_mini <- subset_taxa(data_fungi, Phylum == "Basidiomycota")`
and then
`data_fungi_mini <- subset_taxa_pq(data_fungi_mini, colSums(data_fungi_mini@otu_table) > 5000)`
