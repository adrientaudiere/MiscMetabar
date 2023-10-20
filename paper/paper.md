---
title: 'MiscMetabar : a R packages to facilitate visualization and Reproducibility in metabarcoding analysis'
tags:
  - R
  - Bioinformatic
  - Metagenomics
  - Barcoding
  - Reproducibility
authors:
  - name: Adrien Taudière
    orcid: 0000-0003-1088-1182
    affiliation: "1" # (Multiple affiliations must be quoted)

affiliations:
 - name: IdEst, Saint-Bonnet-de-Salendrinque, 30460 France
   index: 1

date: 23/10/2023
bibliography: paper.bib
---

# Summary
Describing communities of living organisms increasingly relies on massive DNA sequencing from environmental samples (e-DNA). The analysis of these large amounts of sequences is well established in the R ecosystem, especially for metabarcoding, i.e. the massive sequencing of a given DNA region, called markers. The `MiscMetabar` package aims to facilitate the *description*, *transformation*, *exploration*, and *reproducibility* of metabarcoding analysis. Several tutorials are available [online](https://adrientaudiere.github.io/MiscMetabar/articles/).

# Statement of Need

Biological studies, especially in ecology, health sciences and taxonomy, need to describe the biological composition of samples. During the last twenty years, (i) the development of DNA sequencing, (ii) reference databases, (iii) high-throughput sequencing (HTS), and (iv) bioinformatics resources have allowed the description of biological communities through metabarcoding. Metabarcoding involves the sequencing of millions (*meta*-) of short regions of specific DNA (*-barcoding*, [@valentini2009]) often from environmental samples (eDNA, [@taberlet2012]) such as human stomach contents, lake water, soil and air.

The analysis of metabarcoding studies is difficult due to the complexity of the datasets (millions of sequences of variable quality, including chimeras and sequencing errors), often coupled with tricky sampling schemes and a wide variety of questions. Several plateforms (see Table 1 in [@tedersoo2022] for a list) such as QIIME2 ([@bolyen2019]), mothur ([@schloss2020]), and Galaxy ([@jalili2020]) allow complete analysis from raw fastq sequences to statistical analysis and visualization. However, the R ecosystem ([@rcran]), especially through the [bioconductor](http://bioconductor.org/) project, is very rich (see the "State of the field in R" section) and more flexible than these platforms.

`MiscMetabar` aims to facilitate the **description**, **transformation**, **exploration** and **reproducibility** of metabarcoding analysis using R. The development of `MiscMetabar` relies heavily on the R packages `dada2`, `phyloseq` and `targets`.


# State of the Field in R

The metabarcoding ecosystem in the R language is mature, well-constructed, and relies on a very active community in both the [bioconductor](https://www.bioconductor.org/) and [cran](https://cran.r-project.org/) projects. The [bioconductor](https://www.bioconductor.org/) even creates specific task views in [Metagenomics](http://bioconductor.org/packages/release/BiocViews.html#___Metagenomics) and [Microbiome](http://bioconductor.org/packages/release/BiocViews.html#___Microbiome).

R package `dada2` [@Callahan:2016] (http://bioconductor.org/packages/release/bioc/html/dada2.html) provides a highly cited and recommended clustering method [@Pauvert:]. `dada2` also provides tools to complete the metabarcoding analysis pipeline, including chimera detection and taxonomic assignment. `phyloseq` [@McMurdie:2013] (http://bioconductor.org/packages/release/bioc/html/phyloseq.html) facilitate metagenomics analysis by providing a way to store data (the `phyloseq` class) and both graphical and statistical functions.

The phyloseq package introduces the S4 class object (class physeq), which contains (i) an OTU sample matrix, (ii) a taxonomic table, (iii) a sample metadata table, and two optional slots for (iv) a phylogenetic tree and (v) reference sequences.

Some packages already extend the phyloseq packages. For example, the [`microbiome`](https://microbiome.github.io/) package collection [@Lahti2017] provides some scripts and functions for manipulating microbiome datasets. The `speedyseq` package [@mclaren2020] provides faster versions of phyloseq's plotting and taxonomic merging functions, some of which are used in `MiscMetabar`. The [phylosmith](https://schuyler-smith.github.io/phylosmith/) [Smith2019](https://joss.theoj.org/papers/10.21105/joss.01442) package already provides some functions to extend and simplify the use of the phyloseq packages.

Other packages ([`mia`](https://github.com/microbiome/mia/) forming the [`microbiome`](https://microbiome.github.io/) package collection and [`MicrobiotaProcess`](https://github.com/YuLab-SMU/MicrobiotaProcess) [@xu2023]) extend a new data structure using the comprehensive Bioconductor ecosystem of the `SummarizedExperiment` family.

`MiscMetabar` enriches this R ecosystem by providing functions to (i) **describe** your dataset visually, (ii) **transform** your data, (iii) **explore** biological diversity (alpha, beta, and taxonomic diversity), and (iv) simplify **reproducibility**. `MiscMetabar` is designed to complement and not compete with other R packages mentioned above. For example. The `mia` package is recommended for studies focusing on phylogenetic trees, and `phylosmith` allows easy visualization of co-occurrence networks. Using the `MicrobiotaProcess::as.MPSE` function, most of the utilities in the `MicrobiotaProcess` package are available with functions from the `MiscMetabar`.

We do not try to reinvent the wheel and prefer to rely on existing packages and classes rather than building a new framework. `MiscMetabar` is based on the phyloseq class from phyloseq, the most cited package in metagenomics ([@wen2023]). For a description and comparison of these integrated packages competing with phyloseq (e.g. [microeco](https://github.com/ChiLiubio/microeco); [@liu2020], [EasyAmplicon](https://github.com/YongxinLiu/EasyAmplicon) ([@liu2023]) and [MicrobiomeAnalystR](https://www.microbiomeanalyst.ca) [@lu2023]) see [@wen2023]. Note that some limitations of the phyloseq packages are circumvented thanks to [phylosmith](https://schuyler-smith.github.io/phylosmith/) [Smith2019](https://joss.theoj.org/papers/10.21105/joss.01442), [`microViz`](https://david-barnett.github.io/microViz/) ([@Barnett2021]) and [`MiscMetabar`](https://adrientaudiere.github.io/MiscMetabar/).

Some packages provide an interactive interface (such as [animalcules](https://github.com/compbiomed/animalcules) ([@zhao2021]) and [`microViz`](https://david-barnett.github.io/microViz/) ([@Barnett2021]) or web-based platform ([MicrobiomeAnalystR](https://www.microbiomeanalyst.ca) [@lu2023]) useful for rapid exploration and for code-beginner biologists.

# Features

## Description

A quick graphical representation of the phyloseq object is available using the `summary_plot_pq()` function. This plot allows the novice to understand the structure of a phyloseq object and contains useful information  The functions `krona()` and `tax_datatable()` describe the taxonomy of organisms using krona interactive pie chart ([@ondov2011]) and [datatable](https://datatables.net/) libraries, respectively.

## Transformation

## Cleaning and filtering

The function `clean_pq()` validates a phyloseq object, mainly by removing empty taxa and samples, and checking for discrepancies between taxa and sample names in different slots.  

The filter functions `subset_samples_pq()` and `subset_taxa_pq()` complement `phyloseq::subset_samples()` and `phyloseq::subset_taxa()`, allowing the use of a boolean vector to filter samples or taxa from a phyloseq object. <!-- See Table 2 for a list of available functions to manipulate metabarcoding data. -->

We also implement a function to filter taxa based on their blast to a custom database (`filter_asv_blast()`). This function uses the blastn software ([@altschul1990]) to compare ASV sequences to a database and filter out species that are below a given threshold of e-value and/or bit-score.

<!-- ![Table 2: List of available functions to manipulate metabarcoding data](figures_svg/table_2.svg) -->.

### ASV (or ESV) to OTU

**ASV** (stands for *Amplicon Sequence Variant*; also called **ESV** for Exact Amplicon Variant) is a DNA sequence obtained from high-throughput analysis of marker genes. **OTU* are a group of closely related individuals created by clustering sequences based on a threshold of similarity. An ASV is a special case of an OTU with a similarity threshold of 100%.

The choice between ASV and OTU is important because they lead to different results ([@joos2020], Box 2 in [@tedersoo2022]). Most articles recommend making a choice depending on the question. For example, ASV may be better than OTU for describing a group of very closely related species. In addition, ASV are comparable across different datasets (obtained using identical marker genes). On the other hand, [@tedersoo2022] report that ASV approaches overestimate the richness of common fungal species (due to haplotype variation), but underestimate the richness of rare species.  They therefore recommend the use of OTUs in metabarcoding analyses of fungal communities. Finally, [@kauserud2023] argues that the ASV term falls within the original OTU term and recommends adopting  only the OTU terms, but with a concise and clear report on how the OTUs were generated.

Recent articles ([@forster2019],[@antich2021]) propose to use both approach together. They recommend (i) using ASV to denoise the dataset and (ii) for some questions, clustering the ASV sequences into OTUs. This is the goal of the function `asv2otu()`, using either the `DECIPHER::Clusterize` function from R or the [vsearch](https://github.com/torognes/vsearch) software. Another transformation method is implemented in `lulu_pq()`, which uses [@froslev2017]'s method for post-clustering curation of DNA amplicon data.


## Exploration 

`MiscMetabar` provides a large number of facilities to explore the biological diversity in a phyloseq object. In most functions, a parameter enables the effect of the number of reads (sampling depth) to be controlled by rarefaction or other statistical methods, depending on the function. For example, the alpha diversity analysis (function `hill_pq()`) uses the HSD-Tuckey test on a linear model that includes the square roots of the number of reads as the first explanatory variable.

The effect of an environmental variable (beta-diversity) on a biological organism can be explored by Venn diagram (`ggvenn_pq()`), upset plot (`pset_pq()`), and circle plot (`circle_pq()`). This effect can be tested with PERMANOVA (`adonis_pq()`) and the network test (`graph_test_pq()`). If only two modalities are compared, `biplot_pq()` is very useful. Differential abundance analysis can be performed directly using the `plot_deseq2_pq()` function. 

To illustrate the effect of sample variables on the taxonomy, `MiscMetabar` provides the functions `treemap_pq()`, `multitax_bar_pq()` and `heat_tree_pq()`. 

## Reproducibility

The targets R package ([@Landau2021]) improves the efficiency and reproducibility of the pipeline in R. It orchestrates the stage of the pipeline and stores the objects to skip tasks that are already up to date. Given the complexity, runtime, and parameter sensitivity of bioinformatic analysis, the use of targets is particularly relevant for metabarcoding. We develop functions to list fastq files in a directory (`list_fastq_files()`) and to track the number of sequences, clusters and samples through the pipeline (`track_wkflow()`) for a variety of objects (fastq files, OTU matrix, dada-class, derep-class and phyloseq-class). Function `write_pq()` save an object of class phyloseq and `read_pq()` read a phyloseq object from files. 

# Already cited package

`MiscMetabar` is already used by the scientific community in several teams ([@Vieira2021], [@Pleic2022], [@McCauley2022], [@McCauley2023], [@bouilloud2023], [@vieira2023]).

# Acknowledgements

I thank Will Landau, Paul J. McMurdie, and Benjamin Callahan for their excellent R packages on which `MiscMetabar` rests. This text were grammatically corrected using [Language Tools](https://languagetool.org/) and [Deepl write](https://www.deepl.com/write) softwares. 

# References