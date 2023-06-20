# MiscMetabar 0.32 (in development)

* Add the possibility to select a folder in the function `count_seq()`
* Add functions `track_wkflow_samples()` and `select_one_sample()`
* Add option `sam_data_first` in function `write_pq()`
* Add option `reorder_asv` and `rename_asv` to in function `write_pq()` and `clean_pq`
* Add a function `rotl_pq()` to build a phylogenetic tree using the ASV binomial names of a physeq object and the Open Tree of Life tree.

# MiscMetabar 0.31 (in development)

* Argument `split_by` to make multiple plot given a variable in `sam_data` slot (function `ggvenn_pq()`)  
* Argument `seq_names` in `asv2otu()` function allow to clusterize sequences from a character vector of DNA.
* Add a `blast_pq()` function to blast the sequences of the `@ref_seq` slot against a custom database
* Add a `filter_asv_blast()` function to filter ASV in phyloseq dataset using blast against a custom database
* Add a `subset_taxa_pq()` function to filter ASV based on a named conditional vector. Used in `filter_asv_blast()`.
* Add parameter `force_taxa_as_columns` (default FALSE) and `force_taxa_as_rows` (default FALSE) to `clean_pq()`.
* Add a first version of the function `count_fastq_seq()`  to count sequences from fastq.gz files directly from R.
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
