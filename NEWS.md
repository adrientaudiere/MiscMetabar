# MiscMetabar 0.41  (in development)

- Add function `iNEXT_pq()` to calculate hill diversity using the [iNEXT](https://github.com/AnneChao/iNEXT) package.
- Add argument `paires` to `multi_biplot_pq()` in order to indicate all paires of samples we want to print.
- Improve `compare_pairs_pq()` with information about the number of shared sequences among paires
- Add function `upset_pq()` to plot upset of phyloseq object using the [ComplexUpset](https://krassowski.github.io/complex-upset/) package
- Add info (param `add_info`) in subtitle of the `hill_pq()` function


# MiscMetabar 0.40  

* Add function `multi_biplot_pq()` to visualize a collection of couples of samples for comparison through a list of `biplot_pq()`
* Add options `add_info`, `na_remove`, and `clean_pq` to `plot_tax_pq()` function
* Add options `vsearch_cluster_method` and `vsearch_args` to `otu2asv()` for more detailed control of the vsearch software
* Suppression of buggy function `MM_idtaxa()`
* Add a wrapper of `write_pq()` called `save_pq()` to save a *phyloseq* object in the three possible formats () at the same time
  * 4 separate tables
  * 1 table version 
  * 1 RData file
* Add a function `add_blast_info()` to add information from `blast_pq()` to the `tax_table` slot of a *phyloseq* object
* Add option `keep_temporary_files` in `asv2otu()` function
* Improve the documentation of `asv2otu()` and fix a little bug in the name of the conserved ASV after `asv2otu()`
* Test coverage largely improved leading to numerous minor bug fixes.
* Add function `search_exact_seq_pq()` to search for exact matching of sequences using complement, reverse and reverse-complement against a phyloseq object.
* Add function `add_new_taxonomy_pq()` to add new taxonomic rank to a phyloseq object. For exemple to add taxonomic assignment from a new database.
* Add a battery of test using `test_that` package and improve code compatibility with cran recommendations.

## BREAKING CHANGES
* `asv2otu()` with `method="vsearch"` change two default values (to repeat the precedent behavior, use `asv2otu(..., vsearch_cluster_method = "--cluster_fast", tax_adjust = 1)`): 
  * vsearch_cluster_method = "--cluster_size"
  * tax_adjust = 0


# MiscMetabar 0.34

* Add option `add_nb_samples` to `ggvenn_pq()` which add the number of samples to level name in the plot. Useful to see disequilibrium in the number of samples among the factor's levels.
* Add option `args_makedb` and `args_blastn` to funtions `blast_pq()`, `blast_to_phyloseq()`, `blast_to_derep()` and `filter_asv_blast()`. 
* Add option `rarefy_nb_seqs` to `ggven_pq()` in order to rarefy samples before plotting. 
* Add function `SRS_curve_pq()` to plot scaling with ranked subsampling (SRS) curves using the `SRS::SRS_curve()` function (see citation("SRS") for reference).
* Add option `nb_samples_info` to `biplot_pq()` in order to add the number of samples merged by level of factors.
* Add a message when two modalities differ greatly (more than x2) in their number of sequences in `biplot_pq()` and `ggvenn_pq()`
* Add options `na_remove`, `dist_method` (including Aitchinson and robust-Aitchinson distance), `correction_for_sample_size` and `rarefy_nb_seqs` options to `adonis_pq()` function.
* Add option `na_remove` to `graph_test_pq()` function.

# MiscMetabar 0.33

* Add function `plot_tax_pq()` to plot taxonomic distribution (nb of sequences or nb of ASV) across factor.
* Add option `add_points` and make better axis of `hill_pq()` function
* Add function `blast_to_derep()` in order to facilitate searching some fasta sequences in dereplicated sequences (obtained by `dada2::derepFastq`) 

|                       | Database (makeblastdb)                         | Sequences to blast (blastn)       |
|-----------------------|------------------------------------------------|-----------------------------------|
| `blast_to_phyloseq()` | Built from `ref_seq` slot(physeq-class)        | Custom fasta file                 |
| `blast_to_derep()`    | Built from dereplicate sequences (derep-class) | Custom fasta file                 |
| `blast_pq()`          | Custom database or custom fasta file           | `ref_seq` slot of a physeq object |

* Add functions `tsne_pq()` and `plot_tsne_pq()` to quickly visualize results of the t-SNE multidimensional analysis based on the `Rtsne::Rtsne()` function.


# MiscMetabar 0.32

* Add the possibility to select a folder in the function `count_seq()`
* Add functions `track_wkflow_samples()` and `select_one_sample()`
* Add option `sam_data_first` in function `write_pq()`
* Add option `reorder_asv` and `rename_asv` to in function `write_pq()` and `clean_pq`
* Add a function `rotl_pq()` to build a phylogenetic tree using the ASV binomial names of a physeq object and the Open Tree of Life tree.

# MiscMetabar 0.31

* Argument `split_by` to make multiple plot given a variable in `sam_data` slot (function `ggvenn_pq()`)  
* Argument `seq_names` in `asv2otu()` function allow to clusterize sequences from a character vector of DNA.
* Add a `blast_pq()` function to blast the sequences of the `@ref_seq` slot against a custom database
* Add a `filter_asv_blast()` function to filter ASV in *phyloseq* dataset using blast against a custom database
* Add a `subset_taxa_pq()` function to filter ASV based on a named conditional vector. Used in `filter_asv_blast()`.
* Add parameter `force_taxa_as_columns` (default FALSE) and `force_taxa_as_rows` (default FALSE) to `clean_pq()`.
* Add a first version of the function `count_fastq_seq()` to count sequences from fastq.gz files directly from R.
* Add taxonomic info to `track_wkflow()` function (parameter `taxonomy_rank`)

# MiscMetabar 0.3 

* Change some function names, mainly replacing `physeq` by `pk`.
\tabular{rl}{
  `graph_test_pq()` \tab now a synonym for `physeq_graph_test`\cr
  `adonis_pq()` \tab now a synonym for `adonis_phyloseq`\cr
  `clean_pq()` \tab now a synonym for `clean_physeq`\cr
  `lulu_pq()` \tab now a synonym for `lulu_phyloseq`\cr
  `circle_pq()` \tab now a synonym for `otu_circle`\cr
  `biplot_pq()` \tab now a synonym for `biplot_physeq`\cr
  `read_pq()` \tab now a synonym for `read_phyloseq`\cr
  `write_pq()` \tab now a synonym for `write_phyloseq`\cr
  `sankey_pq()` \tab now a synonym for `sankey_phyloseq`\cr
  `summary_plot_pq()` \tab now a synonym for `summary_plot_phyloseq`\cr
  `plot_edgeR_pq()` \tab now a synonym for `plot_edgeR_phyloseq`\cr
  `plot_deseq2_pq()` \tab now a synonym for `plot_deseq2_phyloseq`\cr
  `venn_pq()` \tab now a synonym for `venn_phyloseq`\cr
  `ggvenn_pq()` \tab now a synonym for `ggVenn_phyloseq`\cr
  `hill_tuckey_pq()` \tab now a synonym for `hill_tuckey_phyloseq`\cr
  `hill_pq()` \tab now a synonym for `hill_phyloseq`\cr
  `heat_tree_pq()` \tab now a synonym for `physeq_heat_tree`\cr
  `compare_pairs_pq()` \tab now a synonym for `multiple_share_bisamples`\cr
}
* Improve documentation using some rules documented in the Rules vignettes.
* Add a option `sam_names()` to `read_pq()`
* Correction of `data_fungi` and `data_fungi_sp_known` metadata


# MiscMetabar 0.24

* Add supplementary info in summary_plot_physeq()`
* Better arguments in biplot_physeq()`)
* Add merge_sample_by argument in biplot_physeq()`
* Better documentation with more example.
* For other minors bugs fixes and addition, see the list of commits

# MiscMetabar 0.23

* Adapt the function `asv2otu()` to *IdClusters* change in the DECIPHER package (commit 254100922f2093cc789d018c18a26752a3cda1e3). Then change the *IdClusters* function that was removed from DECIPHER to *Clusterize* function.

* Better functioning of `blast_to_phyloseq()` when none query sequences are founded.

* Add *tax_adjust* argument to `asv2otu()`function

* Add some functions useful for the targets package

* Add a `biplot_physeq()` function to visualize of two samples for comparison of physeq object

* Add an argument *modality* in the `tax_datatable()` function to split OTU abundancy by level of the sample modality

* Add a function `multiple_share_bisamples()` to help compare samples by pairs

* Add a new function (`ggVenn_phyloseq()`) for better venn diagramm but without area calculation (use `venn_phyloseq()` in this case).

* Add two functions helpful for beta-diversity analysis (`adonis_phyloseq()` and `physeq_graph_test()`)

# MiscMetabar 0.22

* Add badge to set the development lifecycle of each function
* Add the lulu_phyloseq function to make easy the reclustering of phyloseq object using the lulu algorithm (https://www.nature.com/articles/s41467-017-01312-x) from the [lulu package](https://github.com/adrientaudiere/lulu).


# MiscMetabar 0.21

* This is the first release of pkgdown.
