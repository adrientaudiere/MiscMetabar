# Changelog

## MiscMetabar 0.14.4 (in development)

CRAN release: 2025-09-30

### New features and improvements

- Add function
  [`plot_seq_ratio_pq()`](https://adrientaudiere.github.io/MiscMetabar/reference/plot_seq_ratio_pq.md)
  to explore the number of sequences per samples using difference ratio
  of the number of sequences per samples ordered by the number of
  sequences.

- Add params `discard_genus_alone`, `pattern_to_remove_tip` and
  `pattern_to_remove_node` to
  [`rotl_pq()`](https://adrientaudiere.github.io/MiscMetabar/reference/rotl_pq.md)
  to enhance the default naming of nodes and tips

- Improve documentation consistency following the style guide

- Allow `DNAStringSet` object as input of
  [`swarm_clustering()`](https://adrientaudiere.github.io/MiscMetabar/reference/swarm_clustering.md)
  and
  [`physeq_or_string_to_dna()`](https://adrientaudiere.github.io/MiscMetabar/reference/physeq_or_string_to_dna.md)

- Add param `rank_propagation` in
  [`merge_taxa_vec()`](https://adrientaudiere.github.io/MiscMetabar/reference/merge_taxa_vec.md)
  to dissable the rank propagation of NA when merging taxa. It is useful
  when merging taxa with informations in the tax_table slot that do not
  follow a strict taxonomic hierarchical structure (e.g. functional
  guilds).

### Bug fixes

- Fix a bug in
  [`subset_taxa_pq()`](https://adrientaudiere.github.io/MiscMetabar/reference/subset_taxa_pq.md)
  when the condition was TRUE only for one taxon

- Fix warnings in
  [`graph_test_pq()`](https://adrientaudiere.github.io/MiscMetabar/reference/graph_test_pq.md)
  with ggplot2 v.4.0.0

- Fix a bug in `upseq_pq()` when using `min_nb_seq` parameter.

- Fix a bug in blast function by allowing value to be equal (not
  strictly greater) to the threshold values `id_cut`, `bit_score_cut`,
  `min_cover_cut` and `e_value_cut`.

- Fix a bug in swarm associated functions
  (`swarm_clustering(), add_swarms_to_pq()`) to take into account the
  `d` parameter. Also add a parameter fastidious that is automatically
  set to FALSE is d is different from 1.

### BREAKING CHANGE

- Replace `species_colnames` by `taxonomic_ranks` in
  [`rotl_pq()`](https://adrientaudiere.github.io/MiscMetabar/reference/rotl_pq.md)
- Parameter name changes in
  [`plot_mt()`](https://adrientaudiere.github.io/MiscMetabar/reference/plot_mt.md)
  and
  [`krona()`](https://adrientaudiere.github.io/MiscMetabar/reference/krona.md)
  - [`plot_mt()`](https://adrientaudiere.github.io/MiscMetabar/reference/plot_mt.md):
    `alpha` → `pval` (aligns with existing pval pattern in other
    functions)
  - [`krona()`](https://adrientaudiere.github.io/MiscMetabar/reference/krona.md):
    `file` → `file_path` (aligns with existing file_path pattern)

## MiscMetabar 0.14.3

CRAN release: 2025-06-21

- Better message in
  [`subset_taxa_tax_control()`](https://adrientaudiere.github.io/MiscMetabar/reference/subset_taxa_tax_control.md)
- Add parameters `text_size` and `text_size_info` to expand or minimize
  text annotation in
  [`summary_plot_pq()`](https://adrientaudiere.github.io/MiscMetabar/reference/summary_plot_pq.md).
- Add function
  [`filt_taxa_wo_NA()`](https://adrientaudiere.github.io/MiscMetabar/reference/filt_taxa_wo_NA.md)
  to filter out taxa with NA values at given taxonomic rank(s)
- Fix a bug in
  [`format2dada2()`](https://adrientaudiere.github.io/MiscMetabar/reference/format2dada2.md)
  by adding semicolons to fill all the taxonomic levels if `from_sintax`
  is TRUE
- Fix a bug in
  [`adonis_pq()`](https://adrientaudiere.github.io/MiscMetabar/reference/adonis_pq.md)
  for method `aitchison` and `robust.aitchison`.

## MiscMetabar 0.14.2

CRAN release: 2025-03-20

- Minor bug fix for CRAN resubmission

## MiscMetabar 0.14.1

CRAN release: 2025-02-19

- Add the possibility to use to resolve conflict using
  [`resolve_vector_ranks()`](https://adrientaudiere.github.io/MiscMetabar/reference/resolve_vector_ranks.md)
  in the
  [`assign_sintax()`](https://adrientaudiere.github.io/MiscMetabar/reference/assign_sintax.md)
  function.
- Add numerous parameters to
  [`assign_sintax()`](https://adrientaudiere.github.io/MiscMetabar/reference/assign_sintax.md),
  in particular `vote_algorithm` to choose the algo resolving conflict.
- Add param `pattern_to_remove` in
  [`format2dada2()`](https://adrientaudiere.github.io/MiscMetabar/reference/format2dada2.md)

## MiscMetabar 0.14.0

- Better filter of parameters in
  [`add_new_taxonomy_pq()`](https://adrientaudiere.github.io/MiscMetabar/reference/add_new_taxonomy_pq.md).
  Only parameters used by the assign\_\* function corresponding to
  `method` are used.
- Add functions
  [`format2sintax()`](https://adrientaudiere.github.io/MiscMetabar/reference/format2sintax.md),
  [`format2dada2()`](https://adrientaudiere.github.io/MiscMetabar/reference/format2dada2.md)
  and `format2dada2_species` to format fasta database in sintax, dada2
  ([`dada2::assignTaxonomy()`](https://rdrr.io/pkg/dada2/man/assignTaxonomy.html))
  and dada2 Species
  ([`dada2::assignSpecies()`](https://rdrr.io/pkg/dada2/man/assignSpecies.html))
  format
- Add function
  [`assign_dada2()`](https://adrientaudiere.github.io/MiscMetabar/reference/assign_dada2.md)
  to assign Taxonomy (with missing ranks if needed) and to assign
  species using
  [`dada2::assignSpecies()`](https://rdrr.io/pkg/dada2/man/assignSpecies.html)
  with only one database input. Add method `dada2_2steps` in function
  [`add_new_taxonomy_pq()`](https://adrientaudiere.github.io/MiscMetabar/reference/add_new_taxonomy_pq.md)
  which use
  [`assign_dada2()`](https://adrientaudiere.github.io/MiscMetabar/reference/assign_dada2.md)
  function.

## MiscMetabar 0.13.0

- Add function
  [`assign_blastn()`](https://adrientaudiere.github.io/MiscMetabar/reference/assign_blastn.md)
  and add a method `blast` in the function
  [`add_new_taxonomy_pq()`](https://adrientaudiere.github.io/MiscMetabar/reference/add_new_taxonomy_pq.md).
- Add function
  [`resolve_vector_ranks()`](https://adrientaudiere.github.io/MiscMetabar/reference/resolve_vector_ranks.md)
  to resolve conflict in a vector of taxonomy values

## MiscMetabar 0.12.1

CRAN release: 2025-01-29

- Add parameter name `min_bootstrap` in
  [`add_new_taxonomy_pq()`](https://adrientaudiere.github.io/MiscMetabar/reference/add_new_taxonomy_pq.md)
- Bug fix in
  [`assign_idtaxa()`](https://adrientaudiere.github.io/MiscMetabar/reference/assign_idtaxa.md)
- Add parameters `pattern_to_remove` and `remove_NA` to
  [`simplify_taxo()`](https://adrientaudiere.github.io/MiscMetabar/reference/simplify_taxo.md)

## MiscMetabar 0.12.0

- Add function
  [`assign_idtaxa()`](https://adrientaudiere.github.io/MiscMetabar/reference/assign_idtaxa.md)
  and
  [`learn_idtaxa()`](https://adrientaudiere.github.io/MiscMetabar/reference/learn_idtaxa.md)
  to facilitate the taxonomic assignation using the idtaxa algorithm
  from the DECIPHER R package.
- Add option `idtaxa` to method in
  [`add_new_taxonomy_pq()`](https://adrientaudiere.github.io/MiscMetabar/reference/add_new_taxonomy_pq.md)
- Add function
  [`tbl_sum_taxtable()`](https://adrientaudiere.github.io/MiscMetabar/reference/tbl_sum_taxtable.md)
  to summarize tax_table from a phyloseq object
- In function
  [`assign_sintax()`](https://adrientaudiere.github.io/MiscMetabar/reference/assign_sintax.md),
  add params `too_few` (default value “align_start”) and `too_many`
  (default “merge”) to authorize db with variable numbers of rank and
  parenthesis in taxonomic name,

## MiscMetabar 0.11.1

- Add param `suffix` to
  [`add_blast_info()`](https://adrientaudiere.github.io/MiscMetabar/reference/add_blast_info.md)
  allowing multiple use of the function on the same phyloseq object
  (e.g. in order to used different database)
- Add param `return_DNAStringSet` to `write_temp_fasta()` function to
  return a DNAStringSet object in place of a temporary file.
- Add a vignette pkgnet-report.
- Add the possibility to send fasta.gz file to
  [`count_seq()`](https://adrientaudiere.github.io/MiscMetabar/reference/count_seq.md)

## MiscMetabar 0.11

- Add function
  [`filt_taxa_pq()`](https://adrientaudiere.github.io/MiscMetabar/reference/filt_taxa_pq.md)
  to filter taxa based on the number of sequences/occurences
- Add functions
  [`no_legend()`](https://adrientaudiere.github.io/MiscMetabar/reference/no_legend.md)
  and
  [`hill_curves_pq()`](https://adrientaudiere.github.io/MiscMetabar/reference/hill_curves_pq.md)
  to plot hill diversity accumulation curves for phyloseq
- Add function
  [`umap_pq()`](https://adrientaudiere.github.io/MiscMetabar/reference/umap_pq.md)
  to compute Dimensionality Reduction with UMAP
- Add function
  [`plot_complexity_pq()`](https://adrientaudiere.github.io/MiscMetabar/reference/plot_complexity_pq.md)
  to plot kmer complexity of references sequences of a phyloseq object
- Add param `type` to `ridge_pq()` to plot a cumulative version
  (type=“ecdf”) version of ridge
- Introduce the idea of a pq-verse: some other packages will complete
  the MiscMetabar packages to make package maintenance easier. The
  \`comparpq\](<https://github.com/adrientaudiere/comparpq>) package
  will facilitate the comparison of phyloseq object with different
  taxonomy, different clustering methods, different samples with same
  modality or different primers.
- Add functions
  [`assign_vsearch_lca()`](https://adrientaudiere.github.io/MiscMetabar/reference/assign_vsearch_lca.md),
  [`assign_sintax()`](https://adrientaudiere.github.io/MiscMetabar/reference/assign_sintax.md)
  and internal function `write_temp_fasta()`
- Add param `method` to
  [`add_new_taxonomy_pq()`](https://adrientaudiere.github.io/MiscMetabar/reference/add_new_taxonomy_pq.md)
  to allow the use of `dada2::assign_taxonomy()` (default, precedent
  only method available),
  [`assign_sintax()`](https://adrientaudiere.github.io/MiscMetabar/reference/assign_sintax.md)
  or
  [`assign_vsearch_lca()`](https://adrientaudiere.github.io/MiscMetabar/reference/assign_vsearch_lca.md)

## MiscMetabar 0.10.4

- Add functions
  [`plot_refseq_pq()`](https://adrientaudiere.github.io/MiscMetabar/reference/plot_refseq_pq.md)
  and
  [`plot_refseq_extremity_pq()`](https://adrientaudiere.github.io/MiscMetabar/reference/plot_refseq_extremity_pq.md)
  to plot the proportion of each nucleotide and the diversity of
  nucleotides from `@refseq` of a phyloseq object.

## MiscMetabar 0.10.3

- Add params `type`, `na_remove` and `verbose` to
  [`ggvenn_pq()`](https://adrientaudiere.github.io/MiscMetabar/reference/ggvenn_pq.md).
  The type = “nb_seq” allow to plot Venn diagram with the number of
  shared sequences instead of shared ASV.
- Add automatic report in json for the function
  [`cutadapt_remove_primers()`](https://adrientaudiere.github.io/MiscMetabar/reference/cutadapt_remove_primers.md).
- Add param `verbose` to
  [`track_wkflow()`](https://adrientaudiere.github.io/MiscMetabar/reference/track_wkflow.md)
  and improve examples for
  [`track_wkflow()`](https://adrientaudiere.github.io/MiscMetabar/reference/track_wkflow.md)
  and `list_fastq_files`

## MiscMetabar 0.10.2

- Improve code thanks to {lintr} package
- Add option `return_file_path` to
  [`cutadapt_remove_primers()`](https://adrientaudiere.github.io/MiscMetabar/reference/cutadapt_remove_primers.md)
  in order to facilitate targets pipeline
- Add function
  [`sam_data_matching_names()`](https://adrientaudiere.github.io/MiscMetabar/reference/sam_data_matching_names.md)
  to match and verify congruence between fastq files names and sample
  metadata (sam_data)

## MiscMetabar 0.10.1

CRAN release: 2024-10-07

> CRAN 2024-09-10

- Delete function `heat_tree_pq()` because {metacoder} package is
  archived from CRAN.

## MiscMetabar 0.9.4

- Set a seed in the example of `build_tree_pq` to resubmit to CRAN Add a
  param `return_a_vector` in function
  [`filter_trim()`](https://adrientaudiere.github.io/MiscMetabar/reference/filter_trim.md)
  to make possible to return a vector of path as it is useful when used
  with `targets::tar_targets(..., format="file")`)
- Make some storage amelioration by replacing
  [`list()`](https://rdrr.io/r/base/list.html) by `vector(list, ...)`

## MiscMetabar 0.9.3

CRAN release: 2024-09-09

> CRAN 2024-09-09

- Homogenize terminology replacing ASV by taxa/taxon in documentation
  and code
- Build an alias function
  [`filter_taxa_blast()`](https://adrientaudiere.github.io/MiscMetabar/reference/filter_asv_blast.md)
  for
  [`filter_asv_blast()`](https://adrientaudiere.github.io/MiscMetabar/reference/filter_asv_blast.md)
- Build an alias function
  [`postcluster_pq()`](https://adrientaudiere.github.io/MiscMetabar/reference/postcluster_pq.md)
  for
  [`asv2otu()`](https://adrientaudiere.github.io/MiscMetabar/reference/postcluster_pq.md)
- Add param `return_data_for_venn` in function `ggvenn_pq` in order to
  make more customizable plot following [ggVennDiagram
  tutorial](https://gaospecial.github.io/ggVennDiagram/articles/fully-customed.html)

### BREAKING CHANGES

- Replacing misnamed param `rename_asv` by `rename_taxons` in
  [`clean_pq()`](https://adrientaudiere.github.io/MiscMetabar/reference/clean_pq.md)
- Replacing misnamed param `reorder_asv` by `reorder_taxons` in
  [`clean_pq()`](https://adrientaudiere.github.io/MiscMetabar/reference/clean_pq.md)

## MiscMetabar 0.9.2

- Add param `default_fun` in function
  [`merge_samples2()`](https://adrientaudiere.github.io/MiscMetabar/reference/merge_samples2.md)
  in order to replace the default function that change the sample data
  in case of merging. A useful parameter is
  `default_fun=diff_fct_diff_class`.
- Add param `kruskal_test` to
  [`hill_pq()`](https://adrientaudiere.github.io/MiscMetabar/reference/hill_pq.md)
  function to prevent user to mis-interpret Tuckey HSD result (and
  letters) if the global effect of the tested factor on Hill diversity
  is non significant.
- Add param `vioplot` to hill_pq() function to allow violin plot instead
  of boxplot.
- Modify `rarefy_sample_count_by_modality` to debug the case of modality
  with level of length one.

## MiscMetabar 0.9.1

CRAN release: 2024-04-28

> CRAN 2024-04-28

### New functions

- Add functions
  [`taxa_as_rows()`](https://adrientaudiere.github.io/MiscMetabar/reference/taxa_as_rows.md)
  and
  [`taxa_as_columns()`](https://adrientaudiere.github.io/MiscMetabar/reference/taxa_as_columns.md)
  to replace verbose called to
  [`clean_pq()`](https://adrientaudiere.github.io/MiscMetabar/reference/clean_pq.md)
- Add function
  [`ggscatt_pq()`](https://adrientaudiere.github.io/MiscMetabar/reference/ggscatt_pq.md)
  to plot and test for effect of a numerical columns in sam_data on Hill
  number. Its the equivalent for numerical variables of
  [`ggbetween_pq()`](https://adrientaudiere.github.io/MiscMetabar/reference/ggbetween_pq.md)
  which focus on the effect of a factor.
- Add functions
  [`var_par_pq()`](https://adrientaudiere.github.io/MiscMetabar/reference/var_par_pq.md)
  ,
  [`var_par_rarperm_pq()`](https://adrientaudiere.github.io/MiscMetabar/reference/var_par_rarperm_pq.md)
  and
  [`plot_var_part_pq()`](https://adrientaudiere.github.io/MiscMetabar/reference/plot_var_part_pq.md)
  to compute the partition of the variation of community and plot it. It
  introduce the notion of `rarperm` part in the function name. It refers
  to the fact that this function compute permutation of samples depth
  rarefaction to measure the variation due to the random process in
  rarefaction.
- Add function
  [`hill_test_rarperm_pq()`](https://adrientaudiere.github.io/MiscMetabar/reference/hill_test_rarperm_pq.md)
  to test the effect of a factor on hill diversity accounting for the
  variation due to random nature of the rarefaction by sample depth.
- Add function
  [`rarefy_sample_count_by_modality()`](https://adrientaudiere.github.io/MiscMetabar/reference/rarefy_sample_count_by_modality.md)
  to equalize the number of samples for each levels of a modality
  (factor)
- Add function
  [`accu_plot_balanced_modality()`](https://adrientaudiere.github.io/MiscMetabar/reference/accu_plot_balanced_modality.md)
  to plot accumulation curves with balanced modality (same number of
  samples per level) and depth rarefaction (same number of sequences per
  sample)
- Add function
  [`adonis_rarperm_pq()`](https://adrientaudiere.github.io/MiscMetabar/reference/adonis_rarperm_pq.md)
  to compute multiple Permanova analyses on different sample depth
  rarefaction.
- Add function
  [`ggaluv_pq()`](https://adrientaudiere.github.io/MiscMetabar/reference/ggaluv_pq.md)
  to plot taxonomic distribution in alluvial fashion with ggplot2 (using
  the `ggalluvial` package)
- Add function
  [`glmutli_pq()`](https://adrientaudiere.github.io/MiscMetabar/reference/glmutli_pq.md)
  to use automated model selection and multimodel inference with (G)LMs
  for phyloseq object

### New parameters

- Add param `taxa_ranks` in function
  [`psmelt_samples_pq()`](https://adrientaudiere.github.io/MiscMetabar/reference/psmelt_samples_pq.md)
  to group results by samples AND taxonomic ranks.
- Add param `hill_scales` in functions
  [`hill_tuckey_pq()`](https://adrientaudiere.github.io/MiscMetabar/reference/hill_tuckey_pq.md)
  and `hill_p()` to choose the level of the hill number.
- Add param `na_remove` in function
  [`hill_pq()`](https://adrientaudiere.github.io/MiscMetabar/reference/hill_pq.md)
  to remove samples with NA in the factor fact.

## MiscMetabar 0.8.1

- Add param `plot_with_tuckey` to
  [`hill_pq()`](https://adrientaudiere.github.io/MiscMetabar/reference/hill_pq.md).,
- Add function
  [`formattable_pq()`](https://adrientaudiere.github.io/MiscMetabar/reference/formattable_pq.md)
  to make beautiful table of the distribution of taxa across a modality
  using visualization inside in the table.
- Add functions
  [`fac2col()`](https://adrientaudiere.github.io/MiscMetabar/reference/fac2col.md)
  and
  [`transp()`](https://adrientaudiere.github.io/MiscMetabar/reference/transp.md)
  to facilitate manipulation of colors, especially in function
  [`formattable_pq()`](https://adrientaudiere.github.io/MiscMetabar/reference/formattable_pq.md)
- Add functions
  [`signif_ancombc()`](https://adrientaudiere.github.io/MiscMetabar/reference/signif_ancombc.md)
  and
  [`plot_ancombc_pq()`](https://adrientaudiere.github.io/MiscMetabar/reference/plot_ancombc_pq.md)
  to plot significant results from
  [`ancombc_pq()`](https://adrientaudiere.github.io/MiscMetabar/reference/ancombc_pq.md)
  function
- Add function
  [`distri_1_taxa()`](https://adrientaudiere.github.io/MiscMetabar/reference/distri_1_taxa.md)
  to summarize the distribution of one given taxa across level of a
  modality
- Add function
  [`normalize_prop_pq()`](https://adrientaudiere.github.io/MiscMetabar/reference/normalize_prop_pq.md)
  to implement the method proposed by [McKnight et
  al. 2018](https://doi.org/10.5061/dryad.tn8qs35)
- Add function
  [`psmelt_samples_pq()`](https://adrientaudiere.github.io/MiscMetabar/reference/psmelt_samples_pq.md)
  to build data frame of samples information including the number of
  sequences (Abundance) and Hill diversity metrics. Useful to use with
  the [ggstatsplot](https://indrajeetpatil.github.io/ggstatsplot/)
  packages (see examples).
- Replace param `variable` by `fact` in function
  [`ggbetween_pq()`](https://adrientaudiere.github.io/MiscMetabar/reference/ggbetween_pq.md)
  and
  [`hill_pq()`](https://adrientaudiere.github.io/MiscMetabar/reference/hill_pq.md)
  (keeping the variable option in
  [`hill_pq()`](https://adrientaudiere.github.io/MiscMetabar/reference/hill_pq.md)
  for backward compatibility)
- Fix a bug in the class of the return object of function
  [`chimera_removal_vs()`](https://adrientaudiere.github.io/MiscMetabar/reference/chimera_removal_vs.md).
  Now it return a matrix to be able to be parsed on to
  [`dada2::getUniques()`](https://rdrr.io/pkg/dada2/man/getUniques.html)

## MiscMetabar 0.7

> CRAN 2024-03-08

- Add functions
  [`chimera_detection_vs()`](https://adrientaudiere.github.io/MiscMetabar/reference/chimera_detection_vs.md)
  and
  [`chimera_removal_vs()`](https://adrientaudiere.github.io/MiscMetabar/reference/chimera_removal_vs.md)
  to process chimera detection and removal using
  [vsearch](https://github.com/torognes/vsearch) software

- Add functions
  [`filter_trim()`](https://adrientaudiere.github.io/MiscMetabar/reference/filter_trim.md),
  [`sample_data_with_new_names()`](https://adrientaudiere.github.io/MiscMetabar/reference/sample_data_with_new_names.md)
  and
  [`rename_samples()`](https://adrientaudiere.github.io/MiscMetabar/reference/rename_samples.md)
  to facilitate the use of
  [targets](https://books.ropensci.org/targets/) for bioinformatic
  pipeline.

- Add function
  [`add_info_to_sam_data()`](https://adrientaudiere.github.io/MiscMetabar/reference/add_info_to_sam_data.md)
  to expand sam_data slot using a data.frame and using nb_asv and nb_seq

- Add functions
  [`swarm_clustering()`](https://adrientaudiere.github.io/MiscMetabar/reference/swarm_clustering.md)
  and
  [`vsearch_clustering()`](https://adrientaudiere.github.io/MiscMetabar/reference/vsearch_clustering.md)
  and add `swarm` method in the function
  [`asv2otu()`](https://adrientaudiere.github.io/MiscMetabar/reference/postcluster_pq.md)

- Add function
  [`physeq_or_string_to_dna()`](https://adrientaudiere.github.io/MiscMetabar/reference/physeq_or_string_to_dna.md)
  mostly for internal use

- Add function
  [`cutadapt_remove_primers()`](https://adrientaudiere.github.io/MiscMetabar/reference/cutadapt_remove_primers.md)
  to remove primers using
  [cutadapt](https://github.com/marcelm/cutadapt/)

- Add internal functions
  [`is_swarm_installed()`](https://adrientaudiere.github.io/MiscMetabar/reference/is_swarm_installed.md),
  [`is_cutadapt_installed()`](https://adrientaudiere.github.io/MiscMetabar/reference/is_cutadapt_installed.md),
  [`is_vsearch_installed()`](https://adrientaudiere.github.io/MiscMetabar/reference/is_vsearch_installed.md)
  and
  [`is_falco_installed()`](https://adrientaudiere.github.io/MiscMetabar/reference/is_falco_installed.md)
  to test for the availability of external software in order to run
  examples and test from testthat.

- Submit to CRAN and change code to comply with their rules (patch 0.7.1
  to 0.7.9)

- Numerous examples and tests are skipped on CRAN because it spends to
  much time to run. Rules vignettes is updated to details the strategy
  for this.

### BREAKING CHANGES

- Harmonization of parameters names:
  - `add_nb_sequences` -\> `add_nb_seq` in
    [`ggvenn_pq()`](https://adrientaudiere.github.io/MiscMetabar/reference/ggvenn_pq.md)
  - `db` -\> `db_url` in
    [`get_funguild_db()`](https://adrientaudiere.github.io/MiscMetabar/reference/get_funguild_db.md)
  - `db` -\> `db_funguild` in
    [`get_funguild_db()`](https://adrientaudiere.github.io/MiscMetabar/reference/get_funguild_db.md)
  - `file` -\> `file_path` in
    [`get_file_extension()`](https://adrientaudiere.github.io/MiscMetabar/reference/get_file_extension.md)
  - `n_seq` -\> `nb_seq` in
    [`subsample_fastq()`](https://adrientaudiere.github.io/MiscMetabar/reference/subsample_fastq.md)
  - `otutable` -\> `otu_table` in
    [`lulu()`](https://adrientaudiere.github.io/MiscMetabar/reference/lulu.md)
  - `alpha` -\> `pval` in
    [`plot_edgeR_pq()`](https://adrientaudiere.github.io/MiscMetabar/reference/plot_edgeR_pq.md)
    and
    [`plot_deseq2_pq()`](https://adrientaudiere.github.io/MiscMetabar/reference/plot_deseq2_pq.md)
    and change default value from 0.01 to more classical 0.05
  - `sequences` -\> `seq2search` in function
    [`search_exact_seq_pq()`](https://adrientaudiere.github.io/MiscMetabar/reference/search_exact_seq_pq.md)
  - `seq_names` -\> `dna_seq` in function `asv2otu`
- Removing the function `install_pkg_needed()` which do not comply with
  CRAN policies

## MiscMetabar 0.6.0

- Add function
  [`ancombc_pq()`](https://adrientaudiere.github.io/MiscMetabar/reference/ancombc_pq.md)
  to simplify the call to
  [`ANCOMBC::ancombc2()`](https://rdrr.io/pkg/ANCOMBC/man/ancombc2.html)
  : ANalysis of COmpositions of Microbiomes with Bias Correction 2
- Add param `taxa_names_from_physeq` (default FALSE) to
  [`subset_taxa_pq()`](https://adrientaudiere.github.io/MiscMetabar/reference/subset_taxa_pq.md)
- Add param `rarefy_by_sample` (default FALSE) to function
  [`ggbetween_pq()`](https://adrientaudiere.github.io/MiscMetabar/reference/ggbetween_pq.md)
- Add function
  [`are_modality_even_depth()`](https://adrientaudiere.github.io/MiscMetabar/reference/are_modality_even_depth.md)
  to test if samples depth significantly vary among the modalities of a
  factor
- Add functions
  [`merge_taxa_vec()`](https://adrientaudiere.github.io/MiscMetabar/reference/merge_taxa_vec.md)
  and
  [`merge_samples2()`](https://adrientaudiere.github.io/MiscMetabar/reference/merge_samples2.md)
  from the [speedyseq](https://github.com/mikemc/speedyseq/) package
  into MiscMetabar to decrease package dependencies (Thanks to Mike R.
  Mclaren)
- Add function
  [`reorder_taxa_pq()`](https://adrientaudiere.github.io/MiscMetabar/reference/reorder_taxa_pq.md)
  in order to replace the unique call to package MicroViz to decrease
  package dependencies.
- Add functions
  [`get_funguild_db()`](https://adrientaudiere.github.io/MiscMetabar/reference/get_funguild_db.md)
  and
  [`funguild_assign()`](https://adrientaudiere.github.io/MiscMetabar/reference/funguild_assign.md)
  from the [FUNGuildR](https://github.com/brendanf/FUNGuildR/) package
  into MiscMetabar to decrease package dependencies
- Remove all dependencies from packages not available on CRAN or
  Bioconductor. Improve code using `goodpractice::gp`() and
  [`devtools::check()`](https://devtools.r-lib.org/reference/check.html)
  function
- Add messages in various cases (NA in samples data, low number of
  sequences in samples, low number of sequences by taxa) when using
  [`verify_pq()`](https://adrientaudiere.github.io/MiscMetabar/reference/verify_pq.md)
  with args `verbose=TRUE`
- Fix a bug in
  [`multitax_bar_pq()`](https://adrientaudiere.github.io/MiscMetabar/reference/multitax_bar_pq.md)
  when using `nb_seq = FALSE`

## MiscMetabar 0.52

- Add function
  [`ggbetween_pq()`](https://adrientaudiere.github.io/MiscMetabar/reference/ggbetween_pq.md)
  to facilitate comparison of hill number using the power of
  [`ggstatsplot::ggbetweenstats()`](https://indrajeetpatil.github.io/ggstatsplot/reference/ggbetweenstats.html)
- Add function
  [`plot_SCBD_pq()`](https://adrientaudiere.github.io/MiscMetabar/reference/plot_SCBD_pq.md)
  to plot species contributions to beta diversity (SCBD) of samples

## MiscMetabar 0.51

- Add function
  [`LCBD_pq()`](https://adrientaudiere.github.io/MiscMetabar/reference/LCBD_pq.md)
  and
  [`plot_LCBD_pq()`](https://adrientaudiere.github.io/MiscMetabar/reference/plot_LCBD_pq.md)
  to compute, test and plot local contributions to beta diversity (LCBD)
  of samples
- Add function
  [`tbl_sum_samdata()`](https://adrientaudiere.github.io/MiscMetabar/reference/tbl_sum_samdata.md)
  to summarize information from sample data in a table
- Add function
  [`mumu_pq()`](https://adrientaudiere.github.io/MiscMetabar/reference/mumu_pq.md)
  to use [mumu](https://github.com/frederic-mahe/mumu), a fast and
  robust C++ implementation of lulu.
- Add (a mostly internal) function `install_pkg_needed()` to install pkg
  (mostly for package list in *Suggest* in DESCRIPTION) if needed by a
  function.
- Add function
  [`add_funguild_info()`](https://adrientaudiere.github.io/MiscMetabar/reference/add_funguild_info.md)
  and
  [`plot_guild_pq()`](https://adrientaudiere.github.io/MiscMetabar/reference/plot_guild_pq.md)
  to add and plot fungal guild information from taxonomy using
  `FUNGuild` package
- Add function
  [`build_phytree_pq()`](https://adrientaudiere.github.io/MiscMetabar/reference/build_phytree_pq.md)
  to build 3 phylogenetic trees (NJ, UPGMA and ML using `phangorn` R
  package) from the `refseq` slot of a `phyloseq` object, possibly with
  bootstrap values. See the vignettes [Tree
  visualization](https://adrientaudiere.github.io/MiscMetabar/articles/tree_visualization.html)
  for an introduction to tree visualization using `ggtree` R package.

## MiscMetabar 0.5

- Phyloseq object are converted in taxa_are_columns in the
  [`ggvenn_pq()`](https://adrientaudiere.github.io/MiscMetabar/reference/ggvenn_pq.md)
  thanks to issue
  [\#31](https://github.com/adrientaudiere/MiscMetabar/issues/31)

### BREAKING CHANGES

- Rename param `log_10` in function
  [`biplot_pq()`](https://adrientaudiere.github.io/MiscMetabar/reference/biplot_pq.md)
  into `log10trans`
- Rename param `log10transform` in function
  [`circle_pq()`](https://adrientaudiere.github.io/MiscMetabar/reference/circle_pq.md)
  into `log10trans`

## MiscMetabar 0.42

- Add argument `one_plot` (default FALSE, same behavior than before) to
  `hill_pq` function in order to return an unique ggplot2 object with
  the four plots inside.
- Add argument `correction_for_sample_size` (default TRUE, same behavior
  than before) to `hill_pq` and `hill_tuckey_pq` function to allow
  removing any correction for uneven sampling depth.
- Add function
  [`multitax_bar_pq()`](https://adrientaudiere.github.io/MiscMetabar/reference/multitax_bar_pq.md)
  to plot 3 levels of taxonomy in function of samples attributes
- Add function
  [`ridges_pq()`](https://adrientaudiere.github.io/MiscMetabar/reference/ridges_pq.md)
  to plot ridges of one taxonomic level in function of samples
  attributes
- Add function `treemap_pq` to plot treemap of two taxonomic levels

## MiscMetabar 0.41

- Add function
  [`iNEXT_pq()`](https://adrientaudiere.github.io/MiscMetabar/reference/iNEXT_pq.md)
  to calculate hill diversity using the
  [iNEXT](https://github.com/AnneChao/iNEXT) package.
- Add argument `pairs` to
  [`multi_biplot_pq()`](https://adrientaudiere.github.io/MiscMetabar/reference/multi_biplot_pq.md)
  in order to indicate all pairs of samples we want to print.
- Improve
  [`compare_pairs_pq()`](https://adrientaudiere.github.io/MiscMetabar/reference/compare_pairs_pq.md)
  with information about the number of shared sequences among pairs.
- Add function
  [`upset_pq()`](https://adrientaudiere.github.io/MiscMetabar/reference/upset_pq.md)
  to plot upset of phyloseq object using the
  [ComplexUpset](https://krassowski.github.io/complex-upset/) package.
- Add function `upset_test_pq` to test for differences between
  intersections (wrapper of
  [`ComplexUpset::upset_test()`](https://krassowski.github.io/complex-upset/reference/upset_test.html)
  for `phyloseq-object`).
- Add info (param `add_info`) in subtitle of the
  [`hill_pq()`](https://adrientaudiere.github.io/MiscMetabar/reference/hill_pq.md)
  function.
- Add argument `remove_space` to
  [`simplify_taxo()`](https://adrientaudiere.github.io/MiscMetabar/reference/simplify_taxo.md)
  function.
- Add argument `simplify_taxo` to
  [`clean_pq()`](https://adrientaudiere.github.io/MiscMetabar/reference/clean_pq.md)
  function.
- Change argument `rarefy_nb_seq` by `rarefy_before_merging` and add
  arguments `rarefy_after_merging` and `add_nb_seq` to
  [`ggvenn_pq()`](https://adrientaudiere.github.io/MiscMetabar/reference/ggvenn_pq.md)
  function.
- Add arguments `rarefy_after_merging` to
  [`biplot_pq()`](https://adrientaudiere.github.io/MiscMetabar/reference/biplot_pq.md)
  and
  [`upset_pq()`](https://adrientaudiere.github.io/MiscMetabar/reference/upset_pq.md)
  functions.
- Add argument `taxa_fill` to
  [`upset_pq()`](https://adrientaudiere.github.io/MiscMetabar/reference/upset_pq.md)
  function in order to fill the bar with taxonomic rank.
- Add a function
  [`subsample_fastq()`](https://adrientaudiere.github.io/MiscMetabar/reference/subsample_fastq.md)
  to make subset of fastq files in order to test your pipeline with all
  samples but with a low number of reads.
- Add a function
  [`accu_samp_threshold()`](https://adrientaudiere.github.io/MiscMetabar/reference/accu_samp_threshold.md)
  to compute the number of sequence to obtain a given proportion of ASV
  in accumulation curves (\`accu_plot).
- Add a function
  [`tax_bar_pq()`](https://adrientaudiere.github.io/MiscMetabar/reference/tax_bar_pq.md)
  in order to plot taxonomic distribution across samples.

## MiscMetabar 0.40

- Add function
  [`multi_biplot_pq()`](https://adrientaudiere.github.io/MiscMetabar/reference/multi_biplot_pq.md)
  to visualize a collection of couples of samples for comparison through
  a list of
  [`biplot_pq()`](https://adrientaudiere.github.io/MiscMetabar/reference/biplot_pq.md).
- Add options `add_info`, `na_remove`, and `clean_pq` to
  [`plot_tax_pq()`](https://adrientaudiere.github.io/MiscMetabar/reference/plot_tax_pq.md)
  function.
- Add options `vsearch_cluster_method` and `vsearch_args` to `otu2asv()`
  for more detailed control of the vsearch software.
- Suppression of buggy function `MM_idtaxa()`.
- Add a wrapper of
  [`write_pq()`](https://adrientaudiere.github.io/MiscMetabar/reference/write_pq.md)
  called
  [`save_pq()`](https://adrientaudiere.github.io/MiscMetabar/reference/save_pq.md)
  to save a *phyloseq* object in the three possible formats () at the
  same time
  - 4 separate tables
  - 1 table version
  - 1 RData file
- Add a function
  [`add_blast_info()`](https://adrientaudiere.github.io/MiscMetabar/reference/add_blast_info.md)
  to add information from
  [`blast_pq()`](https://adrientaudiere.github.io/MiscMetabar/reference/blast_pq.md)
  to the `tax_table` slot of a *phyloseq* object.
- Add option `keep_temporary_files` in
  [`asv2otu()`](https://adrientaudiere.github.io/MiscMetabar/reference/postcluster_pq.md)
  function.
- Improve the documentation of
  [`asv2otu()`](https://adrientaudiere.github.io/MiscMetabar/reference/postcluster_pq.md)
  and fix a little bug in the name of the conserved ASV after
  [`asv2otu()`](https://adrientaudiere.github.io/MiscMetabar/reference/postcluster_pq.md).
- Test coverage largely improved leading to numerous minor bug fixes.
- Add function
  [`search_exact_seq_pq()`](https://adrientaudiere.github.io/MiscMetabar/reference/search_exact_seq_pq.md)
  to search for exact matching of sequences using complement, reverse
  and reverse-complement against a phyloseq object.
- Add function
  [`add_new_taxonomy_pq()`](https://adrientaudiere.github.io/MiscMetabar/reference/add_new_taxonomy_pq.md)
  to add new taxonomic rank to a phyloseq object. For example to add
  taxonomic assignment from a new database.
- Add a battery of test using `test_that` package and improve code
  compatibility with cran recommendations.

### BREAKING CHANGES

- [`asv2otu()`](https://adrientaudiere.github.io/MiscMetabar/reference/postcluster_pq.md)
  with `method="vsearch"` change two default values (to repeat the
  precedent behavior, use
  `asv2otu(..., vsearch_cluster_method = "--cluster_fast", tax_adjust = 1)`):
  - vsearch_cluster_method = “–cluster_size”
  - tax_adjust = 0

## MiscMetabar 0.34

- Add option `add_nb_samples` to
  [`ggvenn_pq()`](https://adrientaudiere.github.io/MiscMetabar/reference/ggvenn_pq.md)
  which add the number of samples to level name in the plot. Useful to
  see disequilibrium in the number of samples among the factor’s levels.
- Add option `args_makedb` and `args_blastn` to functions
  [`blast_pq()`](https://adrientaudiere.github.io/MiscMetabar/reference/blast_pq.md),
  [`blast_to_phyloseq()`](https://adrientaudiere.github.io/MiscMetabar/reference/blast_to_phyloseq.md),
  [`blast_to_derep()`](https://adrientaudiere.github.io/MiscMetabar/reference/blast_to_derep.md)
  and
  [`filter_asv_blast()`](https://adrientaudiere.github.io/MiscMetabar/reference/filter_asv_blast.md).
- Add option `rarefy_nb_seqs` to `ggven_pq()` in order to rarefy samples
  before plotting.
- Add function
  [`SRS_curve_pq()`](https://adrientaudiere.github.io/MiscMetabar/reference/SRS_curve_pq.md)
  to plot scaling with ranked subsampling (SRS) curves using the
  `SRS::SRS_curve()` function (see citation(“SRS”) for reference).
- Add option `nb_samples_info` to
  [`biplot_pq()`](https://adrientaudiere.github.io/MiscMetabar/reference/biplot_pq.md)
  in order to add the number of samples merged by level of factors.
- Add a message when two modalities differ greatly (more than x2) in
  their number of sequences in
  [`biplot_pq()`](https://adrientaudiere.github.io/MiscMetabar/reference/biplot_pq.md)
  and
  [`ggvenn_pq()`](https://adrientaudiere.github.io/MiscMetabar/reference/ggvenn_pq.md).
- Add options `na_remove`, `dist_method` (including Aitchinson and
  robust-Aitchinson distance), `correction_for_sample_size` and
  `rarefy_nb_seqs` options to
  [`adonis_pq()`](https://adrientaudiere.github.io/MiscMetabar/reference/adonis_pq.md)
  function.
- Add option `na_remove` to
  [`graph_test_pq()`](https://adrientaudiere.github.io/MiscMetabar/reference/graph_test_pq.md)
  function.

## MiscMetabar 0.33

- Add function
  [`plot_tax_pq()`](https://adrientaudiere.github.io/MiscMetabar/reference/plot_tax_pq.md)
  to plot taxonomic distribution (nb of sequences or nb of ASV) across
  factor.
- Add option `add_points` and make better axis of
  [`hill_pq()`](https://adrientaudiere.github.io/MiscMetabar/reference/hill_pq.md)
  function
- Add function
  [`blast_to_derep()`](https://adrientaudiere.github.io/MiscMetabar/reference/blast_to_derep.md)
  in order to facilitate searching some fasta sequences in dereplicated
  sequences (obtained by
  [`dada2::derepFastq`](https://rdrr.io/pkg/dada2/man/derepFastq.html))

|  | Database (makeblastdb) | Sequences to blast (blastn) |
|----|----|----|
| [`blast_to_phyloseq()`](https://adrientaudiere.github.io/MiscMetabar/reference/blast_to_phyloseq.md) | Built from `ref_seq` slot(physeq-class) | Custom fasta file |
| [`blast_to_derep()`](https://adrientaudiere.github.io/MiscMetabar/reference/blast_to_derep.md) | Built from dereplicate sequences (derep-class) | Custom fasta file |
| [`blast_pq()`](https://adrientaudiere.github.io/MiscMetabar/reference/blast_pq.md) | Custom database or custom fasta file | `ref_seq` slot of a physeq object |

- Add functions
  [`tsne_pq()`](https://adrientaudiere.github.io/MiscMetabar/reference/tsne_pq.md)
  and
  [`plot_tsne_pq()`](https://adrientaudiere.github.io/MiscMetabar/reference/plot_tsne_pq.md)
  to quickly visualize results of the t-SNE multidimensional analysis
  based on the
  [`Rtsne::Rtsne()`](https://rdrr.io/pkg/Rtsne/man/Rtsne.html) function.

## MiscMetabar 0.32

- Add the possibility to select a folder in the function
  [`count_seq()`](https://adrientaudiere.github.io/MiscMetabar/reference/count_seq.md)
- Add functions
  [`track_wkflow_samples()`](https://adrientaudiere.github.io/MiscMetabar/reference/track_wkflow_samples.md)
  and
  [`select_one_sample()`](https://adrientaudiere.github.io/MiscMetabar/reference/select_one_sample.md)
- Add option `sam_data_first` in function
  [`write_pq()`](https://adrientaudiere.github.io/MiscMetabar/reference/write_pq.md)
- Add option `reorder_asv` and `rename_asv` to in function
  [`write_pq()`](https://adrientaudiere.github.io/MiscMetabar/reference/write_pq.md)
  and `clean_pq`
- Add a function
  [`rotl_pq()`](https://adrientaudiere.github.io/MiscMetabar/reference/rotl_pq.md)
  to build a phylogenetic tree using the ASV binomial names of a physeq
  object and the Open Tree of Life tree.

## MiscMetabar 0.31

- Argument `split_by` to make multiple plot given a variable in
  `sam_data` slot (function
  [`ggvenn_pq()`](https://adrientaudiere.github.io/MiscMetabar/reference/ggvenn_pq.md))  
- Argument `seq_names` in
  [`asv2otu()`](https://adrientaudiere.github.io/MiscMetabar/reference/postcluster_pq.md)
  function allow to clusterize sequences from a character vector of DNA.
- Add a
  [`blast_pq()`](https://adrientaudiere.github.io/MiscMetabar/reference/blast_pq.md)
  function to blast the sequences of the `@ref_seq` slot against a
  custom database
- Add a
  [`filter_asv_blast()`](https://adrientaudiere.github.io/MiscMetabar/reference/filter_asv_blast.md)
  function to filter ASV in *phyloseq* dataset using blast against a
  custom database
- Add a
  [`subset_taxa_pq()`](https://adrientaudiere.github.io/MiscMetabar/reference/subset_taxa_pq.md)
  function to filter ASV based on a named conditional vector. Used in
  [`filter_asv_blast()`](https://adrientaudiere.github.io/MiscMetabar/reference/filter_asv_blast.md).
- Add parameter `force_taxa_as_columns` (default FALSE) and
  `force_taxa_as_rows` (default FALSE) to
  [`clean_pq()`](https://adrientaudiere.github.io/MiscMetabar/reference/clean_pq.md).
- Add a first version of the function `count_fastq_seq()` to count
  sequences from fastq.gz files directly from R.
- Add taxonomic info to
  [`track_wkflow()`](https://adrientaudiere.github.io/MiscMetabar/reference/track_wkflow.md)
  function (parameter `taxonomy_rank`)

## MiscMetabar 0.3

- Change some function names, mainly replacing `physeq` by `pk`.
- Improve documentation using some rules documented in the Rules
  vignettes.
- Add a option `sam_names()` to
  [`read_pq()`](https://adrientaudiere.github.io/MiscMetabar/reference/read_pq.md)
- Correction of `data_fungi` and `data_fungi_sp_known` metadata

## MiscMetabar 0.24

- Add supplementary info in summary_plot_physeq()\`
- Better arguments in biplot_physeq()\`)
- Add merge_sample_by argument in biplot_physeq()\`
- Better documentation with more example.
- For other minors bugs fixes and addition, see the list of commits

## MiscMetabar 0.23

- Adapt the function
  [`asv2otu()`](https://adrientaudiere.github.io/MiscMetabar/reference/postcluster_pq.md)
  to *IdClusters* change in the DECIPHER package (commit
  254100922f2093cc789d018c18a26752a3cda1e3). Then change the
  *IdClusters* function that was removed from DECIPHER to *Clusterize*
  function.

- Better functioning of
  [`blast_to_phyloseq()`](https://adrientaudiere.github.io/MiscMetabar/reference/blast_to_phyloseq.md)
  when none query sequences are founded.

- Add *tax_adjust* argument to
  [`asv2otu()`](https://adrientaudiere.github.io/MiscMetabar/reference/postcluster_pq.md)function

- Add some functions useful for the targets package

- Add a
  [`biplot_physeq()`](https://adrientaudiere.github.io/MiscMetabar/reference/MiscMetabar-deprecated.md)
  function to visualize of two samples for comparison of physeq object

- Add an argument *modality* in the
  [`tax_datatable()`](https://adrientaudiere.github.io/MiscMetabar/reference/tax_datatable.md)
  function to split OTU abundancy by level of the sample modality

- Add a function `multiple_share_bisamples()` to help compare samples by
  pairs

- Add a new function
  ([`ggVenn_phyloseq()`](https://adrientaudiere.github.io/MiscMetabar/reference/MiscMetabar-deprecated.md))
  for better venn diagram but without area calculation (use
  [`venn_phyloseq()`](https://adrientaudiere.github.io/MiscMetabar/reference/MiscMetabar-deprecated.md)
  in this case).

- Add two functions helpful for beta-diversity analysis
  ([`adonis_phyloseq()`](https://adrientaudiere.github.io/MiscMetabar/reference/MiscMetabar-deprecated.md)
  and
  [`physeq_graph_test()`](https://adrientaudiere.github.io/MiscMetabar/reference/MiscMetabar-deprecated.md))

## MiscMetabar 0.22

- Add badge to set the development lifecycle of each function
- Add the lulu_phyloseq function to make easy the reclustering of
  phyloseq object using the lulu algorithm
  (<https://www.nature.com/articles/s41467-017-01312-x>) from the [lulu
  package](https://github.com/adrientaudiere/lulu).

## MiscMetabar 0.21

- This is the first release of pkgdown.
