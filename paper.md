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

MiscMetabar aims at facilitate the *description*, the *transformation*, the *exploration* and the *reproductibility* of metabarcoding analysis using R.


# Statement of need

Metabarcoding is 



# States of the field in R

The Metabarcoding ecosystem in R language is mature, well constructed and 
rely on a very active communities in both [bioconductor](https://www.bioconductor.org/) 
and [cran](https://cran.r-project.org/) projects. [bioconductor](https://www.bioconductor.org/) even create specific task views in [Metagenomics](http://bioconductor.org/packages/release/BiocViews.html#___Metagenomics) and 
[Microbiome](http://bioconductor.org/packages/release/BiocViews.html#___Microbiome). 

`dada2` [@Callahan:2016] (http://bioconductor.org/packages/release/bioc/html/dada2.html) provides a clustering method highly cited and recommended [@Pauvert:]. `dada2` also supply tools to complete pipeline of metabarcoding analysis including chimera detection and taxonomic assignation. `phyloseq` [@McMurdie:2013] (http://bioconductor.org/packages/release/bioc/html/phyloseq.html) facilitate metagenomics analysis by providing a way to stock data (the `phyloseq` class) and both graphical and statistical functions. 

Some packages already extend the phyloseq packages. For instance the `microbiome` package [@Lahti2017] (https://microbiome.github.io/) propose some scripts and functions to manipulate microbiome data sets. The `speedyseq` package [@mclaren2020] provides faster versions of phyloseq's plotting and taxonomic merging functions, some of which are used in MiscMetabar.

Others packages (e.g. the `mia` package) extend a recent data structure using the `TreeSummarizedExperiment` package [@]  to facilitate phylogenetic tree analysis and to benefits from the comprehensive Bioconductor ecosystem of the `SummarizedExperiment` family.

MiscMetabar enrich this R ecosystem by proposing functions to (i) describe your dataset visually, (ii) transform your data, (iii) explore the biological diversity (alpha, beta and taxonomic diversity) and (iv) simplify reproductibility.

# Features

# Finding fungus in a needle stack showcase




# Already used packages

MiscMetabar is already used by scientific communities in different teams ([@Pleic2022], [@Vieira2021], [@McCauley2022]). 




# Acknowledgements




















For a quick reference, the following citation commands can be used:
- `@author:2001`  ->  "Author et al. (2001)"
- `[@author:2001]` -> "(Author et al., 2001)"
- `[@author1:2001; @author2:2001]` -> "(Author1 et al., 2001; Author2 et al., 2002)"




