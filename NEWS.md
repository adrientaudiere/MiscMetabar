# MiscMetabar 0.24 (in development)

* Add supplementary info in `summary_plot_physeq`
* Better arguments in biplot_physeq
* Add merge_sample_by argument in biplot_physeq
* Better documentation with more example.
* For other minors bugs fixes and addition, see the list of commits



# MiscMetabar 0.23

* Adapt the function `asv2otu` to *IdClusters* change in the DECIPHER package (commit 254100922f2093cc789d018c18a26752a3cda1e3). Then change the *IdClusters* function that was removed from DECIPHER to *Clusterize* function.

* Better functioning of `blast_to_phyloseq` when none query sequences are founded.

* Add *tax_adjust* argument to `asv2otu `function

* Add some functions useful for the targets package

* Add a `biplot_physeq` function to visualize of two samples for comparison of physeq object

* Add an argument *modality* in the `tax_datatable` function to split OTU abundancy by level of the sample modality

* Add a function `multiple_share_bisamples` to help compare samples by pairs

* Add a new function (`ggVenn_phyloseq`) for better venn diagramm but without area calculation (use `venn_phyloseq` in this case).

* Add two functions helpful for beta-diversity analysis (`adonis_phyloseq` and `physeq_graph_test`)

# MiscMetabar 0.22

* Add badge to set the development lifecycle of each function
* Add the lulu_phyloseq function to make easy the reclustering of phyloseq object using the lulu algorithm (https://www.nature.com/articles/s41467-017-01312-x) from the [lulu package](https://github.com/adrientaudiere/lulu).


# MiscMetabar 0.21

* This is the first release of pkgdown.
