# MiscMetabar 0.7 (in development)

- Add functions `chimera_detection_vs()` and `chimera_removal_vs()` to process chimera detection and removal using [vsearch](https://github.com/torognes/vsearch) software 
- Add functions `filter_trim()`, `sample_data_with_new_names()` and `rename_samples()` to facilitate the use of [targets](https://books.ropensci.org/targets/) for bioinformatic pipeline.
- Add function `add_info_to_sam_data()` to expand sam_data slot using a data.frame and using nb_asv and nb_seq 
- Add functions `swarm_clustering()` and `vsearch_clustering()` and add `swarm` method in the function `asv2otu()`
- Add function `physeq_or_string_to_dna()` mostly for internal use
- Add function `cutadapt_remove_primers()` to remove primers using [cutadapt](https://github.com/marcelm/cutadapt/)
- Add internal functions `is_swarm_installed()`, `is_cutadapt_installed()`, `is_vsearch_installed()` and `is_falco_installed()` to test for the availability of external software in order to run examples and test from testthat.

- Submit to CRAN and change code to comply with their rules (patch 0.7.1 to 0.7.8)
- Numerous examples and tests are skipped on CRAN because it spends to much time to run. Rules vignettes is updated to details the strategy for this.


## BREAKING CHANGES

- Harmonization of parameters names:
  - `add_nb_sequences` -> `add_nb_seq` in `ggvenn_pq()`
  - `db` -> `db_url` in `get_funguild_db()`
  - `db` -> `db_funguild` in `get_funguild_db()`
  - `file` -> `file_path` in `get_file_extension()`
  - `n_seq` -> `nb_seq` in `subsample_fastq()`
  - `otutable` -> `otu_table` in `lulu()`
  - `alpha` -> `pval` in `plot_edgeR_pq()` and `plot_deseq2_pq()` and change default value from 0.01 to more classical 0.05
  - `sequences` -> `seq2search` in function `search_exact_seq_pq()`
  - `seq_names` -> `dna_seq` in function `asv2otu`


# MiscMetabar 0.6.0 

- Add function `ancombc_pq()` to simplify the call to `ANCOMBC::ancombc2()` : ANalysis of COmpositions of Microbiomes with Bias Correction 2 
- Add param `taxa_names_from_physeq` (default FALSE) to `subset_taxa_pq()` 
- Add param `rarefy_by_sample` (default FALSE) to function `ggbetween_pq()`
- Add function `are_modality_even_depth()` to test if samples depth significantly vary among the modalities of a factor
- Add functions `merge_taxa_vec()` and `merge_samples2()` from the [speedyseq](https://github.com/mikemc/speedyseq/) package into MiscMetabar to decrease package dependencies (Thanks to Mike R. Mclaren)
- Add function `reorder_taxa_pq()` in order to replace the unique call to package MicroViz to decrease package dependencies.
- Add functions `get_funguild_db()` and `funguild_assign()` from the [FUNGuildR](https://github.com/brendanf/FUNGuildR/) package into MiscMetabar to decrease package dependencies 
- Remove all dependencies from packages not available on CRAN or Bioconductor. Improve code using `goodpractice::gp`() and `devtools::check()` function
- Add messages in various cases (NA in samples data, low number of sequences in samples, low number of sequences by taxa) when using `verify_pq()` with args `verbose=TRUE`
- Fix a bug in `multitax_bar_pq()` when using `nb_seq = FALSE`

# MiscMetabar 0.52 

- Add function `ggbetween_pq()` to facilitate comparison of hill number using the power of `ggstatsplot::ggbetweenstats()`
- Add function `plot_SCBD_pq()` to plot species contributions to beta diversity (SCBD) of samples

# MiscMetabar 0.51 

- Add function `LCBD_pq()` and `plot_LCBD_pq()` to compute, test and plot local contributions to beta diversity (LCBD) of samples
- Add function `tbl_sum_samdata()` to summarize information from sample data in a table
- Add function `mumu_pq()` to use [mumu](https://github.com/frederic-mahe/mumu), a fast and robust C++ implementation of lulu.
- Add (a mostly internal) function `install_pkg_needed()` to install pkg (mostly for package list in *Suggest* in DESCRIPTION) if needed by a function.
- Add function `add_funguild_info()` and `plot_guild_pq()` to add and plot fungal guild information from taxonomy using `FUNGuild` package
- Add function `build_phytree_pq()` to build 3 phylogenetic trees (NJ, UPGMA and ML using `phangorn` R package) from the `refseq` slot of a `phyloseq` object, possibly with bootstrap values. See the vignettes [Tree visualization](https://adrientaudiere.github.io/MiscMetabar/articles/tree_visualization.html) for an introduction to tree visualization using `ggtree` R package.

# MiscMetabar 0.5

- Phyloseq object are converted in taxa_are_columns in the `ggvenn_pq()` thanks to issue #31

## BREAKING CHANGES

- Rename param `log_10` in function `biplot_pq()` into `log10trans`
- Rename param `log10transform` in function `circle_pq()` into `log10trans`
 
# MiscMetabar 0.42 

- Add argument `one_plot` (default FALSE, same behavior than before) to `hill_pq` function in order to return an unique ggplot2 object with the four plots inside.
- Add argument `correction_for_sample_size` (default TRUE, same behavior than before) to `hill_pq` and `hill_tuckey_pq` function to allow removing any correction for uneven sampling depth.
- Add function `multitax_bar_pq()` to plot 3 levels of taxonomy in function of samples attributes
- Add function `ridges_pq()` to plot ridges of one taxonomic level in function of samples attributes
- Add function `treemap_pq` to plot treemap of two taxonomic levels

# MiscMetabar 0.41

- Add function `iNEXT_pq()` to calculate hill diversity using the [iNEXT](https://github.com/AnneChao/iNEXT) package.
- Add argument `pairs` to `multi_biplot_pq()` in order to indicate all pairs of samples we want to print.
- Improve `compare_pairs_pq()` with information about the number of shared sequences among pairs.
- Add function `upset_pq()` to plot upset of phyloseq object using the [ComplexUpset](https://krassowski.github.io/complex-upset/) package.
- Add function `upset_test_pq` to test for differences between intersections (wrapper of `ComplexUpset::upset_test()` for `phyloseq-object`).
- Add info (param `add_info`) in subtitle of the `hill_pq()` function.
- Add argument `remove_space` to `simplify_taxo()` function.
- Add argument `simplify_taxo` to `clean_pq()` function.
- Change argument `rarefy_nb_seq` by `rarefy_before_merging` and add arguments `rarefy_after_merging` and `add_nb_seq` to `ggvenn_pq()` function.
- Add arguments `rarefy_after_merging` to `biplot_pq()` and `upset_pq()` functions.
- Add argument `taxa_fill` to `upset_pq()` function in order to fill the bar with taxonomic rank.
- Add a function `subsample_fastq()` to make subset of fastq files in order to test your pipeline with all samples but with a low number of reads.
- Add a function `accu_samp_threshold()` to compute the number of sequence to obtain a given proportion of ASV in accumulation curves (`accu_plot).
- Add a function `tax_bar_pq()` in order to plot taxonomic distribution across samples.

# MiscMetabar 0.40  

* Add function `multi_biplot_pq()` to visualize a collection of couples of samples for comparison through a list of `biplot_pq()`.
* Add options `add_info`, `na_remove`, and `clean_pq` to `plot_tax_pq()` function.
* Add options `vsearch_cluster_method` and `vsearch_args` to `otu2asv()` for more detailed control of the vsearch software.
* Suppression of buggy function `MM_idtaxa()`.
* Add a wrapper of `write_pq()` called `save_pq()` to save a *phyloseq* object in the three possible formats () at the same time
  * 4 separate tables
  * 1 table version 
  * 1 RData file
* Add a function `add_blast_info()` to add information from `blast_pq()` to the `tax_table` slot of a *phyloseq* object.
* Add option `keep_temporary_files` in `asv2otu()` function.
* Improve the documentation of `asv2otu()` and fix a little bug in the name of the conserved ASV after `asv2otu()`.
* Test coverage largely improved leading to numerous minor bug fixes.
* Add function `search_exact_seq_pq()` to search for exact matching of sequences using complement, reverse and reverse-complement against a phyloseq object.
* Add function `add_new_taxonomy_pq()` to add new taxonomic rank to a phyloseq object. For example to add taxonomic assignment from a new database.
* Add a battery of test using `test_that` package and improve code compatibility with cran recommendations.

## BREAKING CHANGES
* `asv2otu()` with `method="vsearch"` change two default values (to repeat the precedent behavior, use `asv2otu(..., vsearch_cluster_method = "--cluster_fast", tax_adjust = 1)`): 
  * vsearch_cluster_method = "--cluster_size"
  * tax_adjust = 0


# MiscMetabar 0.34

* Add option `add_nb_samples` to `ggvenn_pq()` which add the number of samples to level name in the plot. Useful to see disequilibrium in the number of samples among the factor's levels.
* Add option `args_makedb` and `args_blastn` to functions `blast_pq()`, `blast_to_phyloseq()`, `blast_to_derep()` and `filter_asv_blast()`.
* Add option `rarefy_nb_seqs` to `ggven_pq()` in order to rarefy samples before plotting. 
* Add function `SRS_curve_pq()` to plot scaling with ranked subsampling (SRS) curves using the `SRS::SRS_curve()` function (see citation("SRS") for reference).
* Add option `nb_samples_info` to `biplot_pq()` in order to add the number of samples merged by level of factors.
* Add a message when two modalities differ greatly (more than x2) in their number of sequences in `biplot_pq()` and `ggvenn_pq()`.
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

* Add a new function (`ggVenn_phyloseq()`) for better venn diagram but without area calculation (use `venn_phyloseq()` in this case).

* Add two functions helpful for beta-diversity analysis (`adonis_phyloseq()` and `physeq_graph_test()`)

# MiscMetabar 0.22

* Add badge to set the development lifecycle of each function
* Add the lulu_phyloseq function to make easy the reclustering of phyloseq object using the lulu algorithm (https://www.nature.com/articles/s41467-017-01312-x) from the [lulu package](https://github.com/adrientaudiere/lulu).


# MiscMetabar 0.21

* This is the first release of pkgdown.
