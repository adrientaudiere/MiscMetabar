# R ecosystem for metabarcoding

This is a short introduction to other R packages in the field of
metabarcoding analysis.

## State of the Field in R

The metabarcoding ecosystem in the R language is mature,
well-constructed, and relies on a very active community in both the
[bioconductor](https://www.bioconductor.org/) and
[cran](https://cran.r-project.org/) projects. The
[bioconductor](https://www.bioconductor.org/) even creates specific task
views in
[Metagenomics](https://bioconductor.org/packages/release/BiocViews.html#___Metagenomics)
and
[Microbiome](https://bioconductor.org/packages/release/BiocViews.html#___Microbiome).

R package
[`dada2`](https://bioconductor.org/packages/release/bioc/html/dada2.html)
(Callahan et al. 2016) provides a highly cited and recommended
clustering method (Pauvert et al. 2019). `dada2` also provides tools to
complete the metabarcoding analysis pipeline, including chimera
detection and taxonomic assignment. `phyloseq` (McMurdie and Holmes
2013)
(<https://bioconductor.org/packages/release/bioc/html/phyloseq.html>)
facilitate metagenomics analysis by providing a way to store data (the
`phyloseq` class) and both graphical and statistical functions.

The phyloseq package introduces the S4 class object (class physeq),
which contains (i) an OTU sample matrix, (ii) a taxonomic table, (iii) a
sample metadata table, and two optional slots for (iv) a phylogenetic
tree and (v) reference sequences.

Some packages already extend the phyloseq packages. For example, the
[`microbiome`](https://microbiome.github.io/) package collection (Ernst
et al. 2023) provides some scripts and functions for manipulating
microbiome datasets.The `speedyseq` package (McLaren 2020) provides
faster versions of phyloseq’s plotting and taxonomic merging functions,
some of which (\[merge_samples2()\] and \[merge_taxa_vec()\]) are
integrated in `MiscMetabar` (thanks to Mike. R. McLaren). The
[phylosmith](https://schuyler-smith.github.io/phylosmith/) [Smith
(2023)](https://joss.theoj.org/papers/10.21105/joss.01442) package
already provides some functions to extend and simplify the use of the
phyloseq packages.

Other packages ([`mia`](https://github.com/microbiome/mia/) forming the
[`microbiome`](https://microbiome.github.io/) package collection and
[`MicrobiotaProcess`](https://github.com/YuLab-SMU/MicrobiotaProcess)
(Xu et al. 2023)) extend a new data structure using the comprehensive
Bioconductor ecosystem of the `SummarizedExperiment` family.

`MiscMetabar` enriches this R ecosystem by providing functions to (i)
**describe** your dataset visually, (ii) **transform** your data, (iii)
**explore** biological diversity (alpha, beta, and taxonomic diversity),
and (iv) simplify **reproducibility**. `MiscMetabar` is designed to
complement and not compete with other R packages mentioned above. For
example. The `mia` package is recommended for studies focusing on
phylogenetic trees, and `phylosmith` allows easy visualization of
co-occurrence networks. Using the
[`MicrobiotaProcess::as.MPSE`](https://rdrr.io/pkg/MicrobiotaProcess/man/as.MPSE.html)
function, most of the utilities in the `MicrobiotaProcess` package are
available with functions from the `MiscMetabar`.

I do not try to reinvent the wheel and prefer to rely on existing
packages and classes rather than building a new framework. `MiscMetabar`
is based on the phyloseq class from phyloseq, the most cited package in
metagenomics (Wen et al. 2023). For a description and comparison of
these integrated packages competing with phyloseq
(e.g. [microeco](https://github.com/ChiLiubio/microeco) by C. Liu et al.
(2020), [EasyAmplicon](https://github.com/YongxinLiu/EasyAmplicon) by
Y.-X. Liu et al. (2023) and
[MicrobiomeAnalystR](https://www.microbiomeanalyst.ca) by Lu et al.
(2023)) see Wen et al. (2023). Note that some limitations of the
phyloseq packages are circumvented thanks to
[phylosmith](https://schuyler-smith.github.io/phylosmith/) (Smith 2023),
[`microViz`](https://david-barnett.github.io/microViz/) ((Barnett, Arts,
and Penders 2021)) and
[`MiscMetabar`](https://adrientaudiere.github.io/MiscMetabar/).

Some packages provide an interactive interface useful for rapid
exploration and for code-beginner biologists.
[Animalcules](https://github.com/wejlab/animalcules) (Zhao et al. 2021)
and [`microViz`](https://david-barnett.github.io/microViz/) (Barnett,
Arts, and Penders 2021) provides shiny interactive interface whereas
[MicrobiomeAnalystR](https://www.microbiomeanalyst.ca) (Lu et al. 2023)
is a web-based platform.

## Session information

``` r
sessionInfo()
```

    ## R version 4.5.1 (2025-06-13)
    ## Platform: x86_64-pc-linux-gnu
    ## Running under: Kali GNU/Linux Rolling
    ## 
    ## Matrix products: default
    ## BLAS:   /usr/lib/x86_64-linux-gnu/openblas-pthread/libblas.so.3 
    ## LAPACK: /usr/lib/x86_64-linux-gnu/openblas-pthread/libopenblasp-r0.3.29.so;  LAPACK version 3.12.0
    ## 
    ## locale:
    ##  [1] LC_CTYPE=fr_FR.UTF-8       LC_NUMERIC=C              
    ##  [3] LC_TIME=fr_FR.UTF-8        LC_COLLATE=fr_FR.UTF-8    
    ##  [5] LC_MONETARY=fr_FR.UTF-8    LC_MESSAGES=fr_FR.UTF-8   
    ##  [7] LC_PAPER=fr_FR.UTF-8       LC_NAME=C                 
    ##  [9] LC_ADDRESS=C               LC_TELEPHONE=C            
    ## [11] LC_MEASUREMENT=fr_FR.UTF-8 LC_IDENTIFICATION=C       
    ## 
    ## time zone: Europe/Paris
    ## tzcode source: system (glibc)
    ## 
    ## attached base packages:
    ## [1] stats     graphics  grDevices utils     datasets  methods   base     
    ## 
    ## loaded via a namespace (and not attached):
    ##  [1] digest_0.6.37     desc_1.4.3        R6_2.6.1          fastmap_1.2.0    
    ##  [5] xfun_0.53         cachem_1.1.0      knitr_1.50        htmltools_0.5.8.1
    ##  [9] rmarkdown_2.29    lifecycle_1.0.4   cli_3.6.5         sass_0.4.10      
    ## [13] pkgdown_2.1.3     textshaping_1.0.3 jquerylib_0.1.4   systemfonts_1.2.3
    ## [17] compiler_4.5.1    tools_4.5.1       ragg_1.5.0        evaluate_1.0.5   
    ## [21] bslib_0.9.0       yaml_2.3.10       jsonlite_2.0.0    rlang_1.1.6      
    ## [25] fs_1.6.6          htmlwidgets_1.6.4

## References

Barnett, David J. M., Ilja C. W. Arts, and John Penders. 2021.
“microViz: An r Package for Microbiome Data Visualization and
Statistics.” *Journal of Open Source Software* 6 (63): 3201.
<https://doi.org/10.21105/joss.03201>.

Callahan, Benjamin J, Paul J McMurdie, Michael J Rosen, Andrew W Han,
Amy Jo A Johnson, and Susan P Holmes. 2016. “DADA2: High-Resolution
Sample Inference from Illumina Amplicon Data.” *Nature Methods* 13 (7):
581–83. <https://doi.org/10.1038/nmeth.3869>.

Ernst, Felix G. M., Sudarshan A. Shetty, Tuomas Borman, and Leo Lahti.
2023. *Mia: Microbiome Analysis*.
<https://doi.org/10.18129/B9.bioc.mia>.

Liu, Chi, Yaoming Cui, Xiangzhen Li, and Minjie Yao. 2020. “microeco: an
R package for data mining in microbial community ecology.” *FEMS
Microbiology Ecology* 97 (2): fiaa255.
<https://doi.org/10.1093/femsec/fiaa255>.

Liu, Yong-Xin, Lei Chen, Tengfei Ma, Xiaofang Li, Maosheng Zheng, Xin
Zhou, Liang Chen, et al. 2023. “EasyAmplicon: An Easy-to-Use,
Open-Source, Reproducible, and Community-Based Pipeline for Amplicon
Data Analysis in Microbiome Research.” *iMeta* 2 (1): e83.

Lu, Yao, Guangyan Zhou, Jessica Ewald, Zhiqiang Pang, Tanisha Shiri, and
Jianguo Xia. 2023. “MicrobiomeAnalyst 2.0: comprehensive statistical,
functional and integrative analysis of microbiome data.” *Nucleic Acids
Research* 51 (W1): W310–18. <https://doi.org/10.1093/nar/gkad407>.

McLaren, Michael. 2020. “Mikemc/Speedyseq: Speedyseq V0.2.0.” Zenodo.
<https://doi.org/10.5281/zenodo.3923184>.

McMurdie, Paul J., and Susan Holmes. 2013. “Phyloseq: An r Package for
Reproducible Interactive Analysis and Graphics of Microbiome Census
Data.” *PLoS ONE* 8 (4): e61217.
<https://doi.org/10.1371/journal.pone.0061217>.

Pauvert, Charlie, Marc Buée, Valérie Laval, Véronique Edel-Hermann,
Laure Fauchery, Angélique Gautier, Isabelle Lesur, Jessica Vallance, and
Corinne Vacher. 2019. “Bioinformatics Matters: The Accuracy of Plant and
Soil Fungal Community Data Is Highly Dependent on the Metabarcoding
Pipeline.” *Fungal Ecology* 41: 23–33.
<https://doi.org/10.1016/j.funeco.2019.03.005>.

Smith, Schuyler. 2023. *Phylosmith: Functions to Help Analyze Data as
Phyloseq Objects*. <https://schuyler-smith.github.io/phylosmith/>.

Wen, Tao, Guoqing Niu, Tong Chen, Qirong Shen, Jun Yuan, and Yong-Xin
Liu. 2023. “The Best Practice for Microbiome Analysis Using r.” *Protein
& Cell*, pwad024.

Xu, Shuangbin, Li Zhan, Wenli Tang, Qianwen Wang, Zehan Dai, Lang Zhou,
Tingze Feng, et al. 2023. “MicrobiotaProcess: A Comprehensive r Package
for Deep Mining Microbiome.” *The Innovation* 4 (2).
<https://doi.org/10.1016/j.xinn.2023.100388>.

Zhao, Yue, Anthony Federico, Tyler Faits, Solaiappan Manimaran, Daniel
Segrè, Stefano Monti, and W Evan Johnson. 2021. “Animalcules:
Interactive Microbiome Analytics and Visualization in r.” *Microbiome* 9
(1): 1–16. <https://doi.org/10.21203/rs.3.rs-29649/v2>.
