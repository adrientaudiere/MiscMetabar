---
title: 'MiscMetabar : a R packages to facilitate visualization and reproductibility in metabarcoding analysis'
tags:
  - R
  - Bioinformatic
  - Metagenomics
  - Barcoding
  - Reproductibility
authors:
  - name: Adrien Taudière
    orcid: 0000-0000-0000-0000
    affiliation: "1" # (Multiple affiliations must be quoted)

affiliations:
 - name: IdEst, Saint-Bonnet-de-Salendrinque, 30460 France
   index: 1

date: TODO
bibliography: paper.bib
---




TODO


# Summary

Metabarcoding is 


# Statement of need

The Metabarcoding ecosystem in R language is mature, well constructed and 
rely on a very active communities in both [bioconductor](https://www.bioconductor.org/) 
and [cran](https://cran.r-project.org/) projects. `dada2` [@Callahan:2016] (http://bioconductor.org/packages/release/bioc/html/dada2.html)
provides a clustering method highly cited and recommended [@Pauvert:]. `dada2` also
provides tools to complete pipeline of metabarcoding analysis. `phyloseq` [@McMurdie:2013] (http://bioconductor.org/packages/release/bioc/html/phyloseq.html) facilitate metagenomics analysis by providing a way to stock data (the `phyloseq` class) and graphical/statistical functions. 

Some packages already extend the phyloseq packages. For instance the `microbiome` package [@Lahti2017] (https://microbiome.github.io/) propose some scripts and functions to manipulate microbiome data sets. The `speedyseq` package [@mclaren2020] provides faster versions of phyloseq's plotting and taxonomic merging functions, some of which are used in MiscMetabar.

Others packages (e.g. the `mia` package) extend a recent data structure using the `TreeSummarizedExperiment` package [@]  to facilitate phylogenetic tree analysis and to benefits from the comprehensive Bioconductor ecosystem of the `SummarizedExperiment` family.

MiscMetabar enrich this R ecosystem by proposing (i) new bioinformatics functions (e.g. `asv2otu`, ),


Bioconductor task views [Metagenomics](http://bioconductor.org/packages/release/BiocViews.html#___Metagenomics) and 
[Microbiome](http://bioconductor.org/packages/release/BiocViews.html#___Microbiome)


# Finding fungus in a needle stack showcase




# Already used packages

MiscMetabar is already used by scientific communities in different teams ([@Pleic2022], [@Vieira2021], [@McCauley2022]). 











# Citations

Citations to entries in paper.bib should be in
[rMarkdown](http://rmarkdown.rstudio.com/authoring_bibliographies_and_citations.html)
format.

If you want to cite a software repository URL (e.g. something on GitHub without a preferred
citation) then you can do it with the example BibTeX entry below for @fidgit.

For a quick reference, the following citation commands can be used:
- `@author:2001`  ->  "Author et al. (2001)"
- `[@author:2001]` -> "(Author et al., 2001)"
- `[@author1:2001; @author2:2001]` -> "(Author1 et al., 2001; Author2 et al., 2002)"

# Figures

<--!

Figures can be included like this:
![Caption for example figure.\label{fig:example}](figure.png)
and referenced from text using \autoref{fig:example}.

Figure sizes can be customized by adding an optional second parameter:
![Caption for example figure.](figure.png){ width=20% }

-->

# Acknowledgements

We acknowledge contributions from Brigitta Sipocz, Syrtis Major, and Semyeong
Oh, and support from Kathryn Johnston during the genesis of this project.

# References 


Huang R, Soneson C, Ernst FGM et al. TreeSummarizedExperiment: a S4 class for data with hierarchical structure [version 2; peer review: 3 approved]. F1000Research 2021, 9:1246 (https://doi.org/10.12688/f1000research.26669.2) 