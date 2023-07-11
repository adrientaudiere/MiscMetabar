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
    orcid: 0000-0003-1088-1182
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

The phyloseq package introduce the S4-class object (class physeq).  S4-class object that contains  (i) a OTU-samples matrix,  (ii) a taxonomic table,  (iii) a sample metadata table as well as a two optionnal slot for (iv) phylogenetic tree and (v) reference sequences. 

Some packages already extend the phyloseq packages. For instance the `microbiome` package [@Lahti2017] (https://microbiome.github.io/) propose some scripts and functions to manipulate microbiome data sets. The `speedyseq` package [@mclaren2020] provides faster versions of phyloseq's plotting and taxonomic merging functions, some of which are used in MiscMetabar. The package [phylosmith](https://schuyler-smith.github.io/phylosmith/) [Smith2019](https://joss.theoj.org/papers/10.21105/joss.01442) already offers some functions to extend and facilitate the use of the phyloseq packages. 

Others packages ([mia](https://github.com/microbiome/mia/) [@] and [MicrobiotaProcess](https://github.com/YuLab-SMU/MicrobiotaProcess) [@xu2023] packages) extend a recent data structure using the comprehensive Bioconductor ecosystem of the `SummarizedExperiment` family. 

MiscMetabar enrich this R ecosystem by proposing functions to (i) **describe** your dataset visually, (ii) **transform** your data, (iii) **explore** the biological diversity (alpha, beta and taxonomic diversity) and (iv) simplify **reproductibility**. MiscMetabar is designed to complete and not compete others R packages previuosly cited. For ex. `mia` package is recommended for studies focusing on phylogenetic tree and `phylosmith` allow to visualize easily co-occurrence networks. Using as.MPSE function, most of utilities in `MicrobiotaProcess` package are available in conjonction with functions form `MiscMetabar`.  


# Features

# Finding fungus in a needle stack showcase




# Already cited package

MiscMetabar is already used by scientific communities in different teams ([@Pleic2022], [@Vieira2021], [@McCauley2022]). 




# Acknowledgements





















For a quick reference, the following citation commands can be used:
- `@author:2001`  ->  "Author et al. (2001)"
- `[@author:2001]` -> "(Author et al., 2001)"
- `[@author1:2001; @author2:2001]` -> "(Author1 et al., 2001; Author2 et al., 2002)"




