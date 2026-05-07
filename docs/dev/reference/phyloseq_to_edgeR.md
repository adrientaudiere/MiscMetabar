# Convert phyloseq OTU count data into DGEList for edgeR package

Convert phyloseq OTU count data into DGEList for edgeR package

## Usage

``` r
phyloseq_to_edgeR(physeq, group, method = "RLE", remove_na = TRUE, ...)
```

## Arguments

- physeq:

  (required) a
  [`phyloseq-class`](https://rdrr.io/pkg/phyloseq/man/phyloseq-class.html)
  object obtained using the `phyloseq` package.

- group:

  (required) A character vector or factor giving the experimental
  group/condition for each sample/library. Alternatively, you may
  provide the name of a sample variable. This name should be among the
  output of `sample_variables(physeq)`, in which case
  `get_variable(physeq, group)` would return either a character vector
  or factor. This is passed on to
  [`DGEList`](https://rdrr.io/pkg/edgeR/man/DGEList.html), and you may
  find further details or examples in its documentation.

- method:

  The label of the edgeR-implemented normalization to use. See
  [`calcNormFactors`](https://rdrr.io/pkg/edgeR/man/calcNormFactors.html)
  for supported options and details. The default option is `"RLE"`,
  which is a scaling factor method proposed by Anders and Huber (2010).
  At time of writing, the
  [edgeR](https://rdrr.io/pkg/edgeR/man/edgeR-package.html) package
  supported the following options to the `method` argument:

  `c("TMM", "RLE", "upperquartile", "none")`.

- remove_na:

  (logical) If TRUE, samples with NA values in the group variable will
  be removed.

- ...:

  Additional arguments passed on to
  [`DGEList`](https://rdrr.io/pkg/edgeR/man/DGEList.html)

## Value

A DGEList object. See
[`edgeR::estimateTagwiseDisp()`](https://rdrr.io/pkg/edgeR/man/estimateTagwiseDisp.html)
for more details.

## Examples

``` r
# \donttest{
if (requireNamespace("edgeR")) {
  phyloseq_to_edgeR(data_fungi_mini, group = "Height")
}
#> Loading required namespace: edgeR
#> An object of class "DGEList"
#> $counts
#>       A10-005-B_S188_MERGED.fastq.gz A10-005-H_S189_MERGED.fastq.gz
#> ASV7                             768                          10188
#> ASV8                             579                             10
#> ASV12                              1                              2
#> ASV18                             31                              1
#> ASV25                              1                              1
#>       A10-005-M_S190_MERGED.fastq.gz A12-007-B_S2_MERGED.fastq.gz
#> ASV7                             264                            1
#> ASV8                             400                          830
#> ASV12                              1                            1
#> ASV18                              4                         1041
#> ASV25                              1                            1
#>       AB29-ABMX-H_S6_MERGED.fastq.gz AD26-005-B_S9_MERGED.fastq.gz
#> ASV7                               1                             1
#> ASV8                               2                             1
#> ASV12                              2                             1
#> ASV18                              1                             1
#> ASV25                              1                             1
#>       AD26-005-H_S10_MERGED.fastq.gz AD26-005-M_S11_MERGED.fastq.gz
#> ASV7                               1                              2
#> ASV8                               1                              1
#> ASV12                              1                           1339
#> ASV18                              1                              1
#> ASV25                              1                              1
#>       AD30-ABMX-M_S12_MERGED.fastq.gz AD32-007-M_S13_MERGED.fastq.gz
#> ASV7                                1                              1
#> ASV8                                2                              8
#> ASV12                               1                              1
#> ASV18                               1                              1
#> ASV25                               1                              1
#>       ADABM30X-B_S14_MERGED.fastq.gz ADABM30X-H_S15_MERGED.fastq.gz
#> ASV7                               1                            540
#> ASV8                               3                              3
#> ASV12                              1                              1
#> ASV18                              1                              1
#> ASV25                              1                              3
#>       ADABM30X-M_S16_MERGED.fastq.gz B18-006-B_S19_MERGED.fastq.gz
#> ASV7                               1                             1
#> ASV8                               1                             1
#> ASV12                              1                             1
#> ASV18                              1                             1
#> ASV25                              1                             1
#>       BA17-050-B_S21_MERGED.fastq.gz BB19-006-H_S22_MERGED.fastq.gz
#> ASV7                            2885                              1
#> ASV8                               5                           9742
#> ASV12                              1                              1
#> ASV18                              1                          11964
#> ASV25                              1                              1
#>       BB6-019-M_S25_MERGED.fastq.gz BG7-010-H_S31_MERGED.fastq.gz
#> ASV7                            874                             1
#> ASV8                            337                             1
#> ASV12                             1                             1
#> ASV18                             8                             1
#> ASV25                             1                             1
#>       BJ17-007-M_S34_MERGED.fastq.gz BL7-006-H_S37_MERGED.fastq.gz
#> ASV7                               1                             1
#> ASV8                               7                             3
#> ASV12                              7                             1
#> ASV18                              1                             1
#> ASV25                              2                             1
#>       BP11-001-B_S43_MERGED.fastq.gz BP11-001-H_S44_MERGED.fastq.gz
#> ASV7                               1                            995
#> ASV8                               2                              9
#> ASV12                              2                              1
#> ASV18                              1                           1846
#> ASV25                              1                              1
#>       BP11-001-M_S45_MERGED.fastq.gz BP12-025-B_S46_MERGED.fastq.gz
#> ASV7                               7                              1
#> ASV8                               7                              2
#> ASV12                              1                              1
#> ASV18                              1                              1
#> ASV25                              1                              1
#>       BQ4-018-B_S49_MERGED.fastq.gz BQ4-018-H_S50_MERGED.fastq.gz
#> ASV7                              1                             1
#> ASV8                              2                             2
#> ASV12                             1                             1
#> ASV18                             1                             1
#> ASV25                             1                             2
#>       BQ4-018-M_S51_MERGED.fastq.gz BT-006-M_S55_MERGED.fastq.gz
#> ASV7                              1                            1
#> ASV8                              1                            1
#> ASV12                             1                            1
#> ASV18                             1                            1
#> ASV25                             1                            1
#>       BV11-002-B_S57_MERGED.fastq.gz BV11-002-H_S58_MERGED.fastq.gz
#> ASV7                               1                              1
#> ASV8                               4                              3
#> ASV12                           7551                              1
#> ASV18                              1                              1
#> ASV25                              1                             67
#>       BV11-002-M_S59_MERGED.fastq.gz C21-NV1-B_S62_MERGED.fastq.gz
#> ASV7                               1                             1
#> ASV8                               1                             1
#> ASV12                              1                             1
#> ASV18                              1                             1
#> ASV25                              1                             1
#>       CB8-019-B_S69_MERGED.fastq.gz CB8-019-H_S70_MERGED.fastq.gz
#> ASV7                              1                             1
#> ASV8                              1                             1
#> ASV12                             2                             1
#> ASV18                             1                             1
#> ASV25                             1                             1
#>       CB8-019-M_S71_MERGED.fastq.gz D18-003-B_S78_MERGED.fastq.gz
#> ASV7                              1                            40
#> ASV8                              1                            32
#> ASV12                             3                             1
#> ASV18                             1                             1
#> ASV25                             1                             1
#>       D18-003-M_S80_MERGED.fastq.gz D61-010-B_S82_MERGED.fastq.gz
#> ASV7                              1                             1
#> ASV8                              1                             6
#> ASV12                             2                             2
#> ASV18                             1                             1
#> ASV25                             1                             5
#>       D9-027-H_S84_MERGED.fastq.gz D9-027-M_S85_MERGED.fastq.gz
#> ASV7                             1                            1
#> ASV8                             4                            2
#> ASV12                            1                            1
#> ASV18                            1                            1
#> ASV25                            1                            1
#>       DJ2-008-H_S88_MERGED.fastq.gz DS1-ABM002-B_S91_MERGED.fastq.gz
#> ASV7                              1                                1
#> ASV8                              1                                1
#> ASV12                             1                                1
#> ASV18                             1                                1
#> ASV25                             1                                1
#>       DS1-ABM002-H_S92_MERGED.fastq.gz DS1-ABM002-M_S93_MERGED.fastq.gz
#> ASV7                                 1                                1
#> ASV8                                 1                                1
#> ASV12                                1                                1
#> ASV18                                1                                1
#> ASV25                                2                                4
#>       DU3-045-B_S94_MERGED.fastq.gz DY5-004-B_S96_MERGED.fastq.gz
#> ASV7                             34                             1
#> ASV8                           1109                             1
#> ASV12                             8                             1
#> ASV18                             8                             1
#> ASV25                             1                             1
#>       DY5-004-H_S97_MERGED.fastq.gz EC2-013-B_S104_MERGED.fastq.gz
#> ASV7                              1                              1
#> ASV8                              1                              2
#> ASV12                             1                              1
#> ASV18                             1                              1
#> ASV25                             1                              1
#>       F7-015-M_S106_MERGED.fastq.gz FOMES19-H_S108_MERGED.fastq.gz
#> ASV7                              1                              2
#> ASV8                              1                              1
#> ASV12                             1                              1
#> ASV18                             1                              1
#> ASV25                             1                              1
#>       FOMES19-M_S109_MERGED.fastq.gz H10-018-M_S110_MERGED.fastq.gz
#> ASV7                               1                              1
#> ASV8                               1                              2
#> ASV12                              1                              1
#> ASV18                              1                              1
#> ASV25                              1                              1
#>       H24-NVABM1-H_S111_MERGED.fastq.gz J18-004-B_S114_MERGED.fastq.gz
#> ASV7                                105                              1
#> ASV8                                  1                              1
#> ASV12                                 3                              1
#> ASV18                                 1                              1
#> ASV25                                10                              2
#>       J18-004-M_S116_MERGED.fastq.gz K18-002-H_S117_MERGED.fastq.gz
#> ASV7                               1                             36
#> ASV8                               1                              1
#> ASV12                              1                              1
#> ASV18                              1                             34
#> ASV25                              1                              1
#>       L19X-B_S119_MERGED.fastq.gz L19X-H_S120_MERGED.fastq.gz
#> ASV7                            1                         336
#> ASV8                          183                          37
#> ASV12                           1                           1
#> ASV18                           1                          10
#> ASV25                           1                           1
#>       L19X-M_S121_MERGED.fastq.gz L23-002-B_S122_MERGED.fastq.gz
#> ASV7                         2801                              7
#> ASV8                          445                             34
#> ASV12                           1                              1
#> ASV18                          22                             28
#> ASV25                           1                              1
#>       L23-002-H_S123_MERGED.fastq.gz L23-002-M_S124_MERGED.fastq.gz
#> ASV7                            2025                              1
#> ASV8                            2111                            152
#> ASV12                              1                              1
#> ASV18                            104                             38
#> ASV25                              1                              1
#>       N19X-H_S127_MERGED.fastq.gz N19X-M_S128_MERGED.fastq.gz
#> ASV7                         1484                         172
#> ASV8                         2318                         282
#> ASV12                           1                           1
#> ASV18                         167                          12
#> ASV25                           1                           1
#>       N23-002-B_S130_MERGED.fastq.gz N23-002-H_S131_MERGED.fastq.gz
#> ASV7                               1                              1
#> ASV8                               1                              2
#> ASV12                              1                              1
#> ASV18                              1                              1
#> ASV25                              1                            437
#>       NVABM-0163-H_S135_MERGED.fastq.gz NVABM0244-M_S137_MERGED.fastq.gz
#> ASV7                                  1                                1
#> ASV8                                  1                                1
#> ASV12                                 1                                1
#> ASV18                                 1                                1
#> ASV25                                 1                                1
#>       O20-X-H_S140_MERGED.fastq.gz O20-X-M_S141_MERGED.fastq.gz
#> ASV7                           102                         2111
#> ASV8                            23                          883
#> ASV12                         1579                            4
#> ASV18                            4                           70
#> ASV25                            1                            1
#>       O24-003-B_S145_MERGED.fastq.gz O26-004-B_S148_MERGED.fastq.gz
#> ASV7                               1                              1
#> ASV8                               3                              1
#> ASV12                              1                              1
#> ASV18                              1                              1
#> ASV25                              1                              1
#>       O9-005-B_S152_MERGED.fastq.gz P19-023-M_S153_MERGED.fastq.gz
#> ASV7                              1                              2
#> ASV8                              2                             27
#> ASV12                             2                              1
#> ASV18                             1                              1
#> ASV25                             1                              1
#>       P27-015-M_S154_MERGED.fastq.gz Q27-ABM003-B_S156_MERGED.fastq.gz
#> ASV7                               1                                 1
#> ASV8                             193                               472
#> ASV12                              8                                 1
#> ASV18                            277                                 1
#> ASV25                              1                                 1
#>       T28-ABM602-B_S162_MERGED.fastq.gz W26-001-B_S165_MERGED.fastq.gz
#> ASV7                                  1                              1
#> ASV8                                  1                              1
#> ASV12                                 1                              1
#> ASV18                                 1                              1
#> ASV25                                 1                              1
#>       W9-025-M_S169_MERGED.fastq.gz X24-009-B_S170_MERGED.fastq.gz
#> ASV7                              1                              1
#> ASV8                              1                              1
#> ASV12                             1                              1
#> ASV18                             1                              1
#> ASV25                             1                              1
#>       X24-009-H_S171_MERGED.fastq.gz X24-009-M_S172_MERGED.fastq.gz
#> ASV7                               2                              1
#> ASV8                               6                              8
#> ASV12                              1                              1
#> ASV18                              1                              1
#> ASV25                              1                              1
#>       X29-004-B_S174_MERGED.fastq.gz X29-004-H_S175_MERGED.fastq.gz
#> ASV7                               1                            108
#> ASV8                               1                            113
#> ASV12                              1                              1
#> ASV18                              1                              2
#> ASV25                              1                              1
#>       X29-004-M_S176_MERGED.fastq.gz Y21-ABM484-H_S177_MERGED.fastq.gz
#> ASV7                               1                                 1
#> ASV8                               3                                 2
#> ASV12                              1                               250
#> ASV18                              1                                 1
#> ASV25                              1                                 1
#>       Y28-002-B_S178_MERGED.fastq.gz Y31-ABM484-B_S184_MERGED.fastq.gz
#> ASV7                               1                                90
#> ASV8                               1                                52
#> ASV12                              1                               630
#> ASV18                              1                                 3
#> ASV25                              1                                 1
#>       Z29-001-H_S185_MERGED.fastq.gz Z30-ABM560-M_S187_MERGED.fastq.gz
#> ASV7                               1                               277
#> ASV8                               1                              9389
#> ASV12                              1                                 1
#> ASV18                              1                               130
#> ASV25                              1                                 1
#> 40 more rows ...
#> 
#> $samples
#>                                 group lib.size norm.factors
#> A10-005-B_S188_MERGED.fastq.gz    Low     1930   0.32203429
#> A10-005-H_S189_MERGED.fastq.gz   High    10246   0.06021446
#> A10-005-M_S190_MERGED.fastq.gz Middle     1179   0.52759695
#> A12-007-B_S2_MERGED.fastq.gz      Low     1914   0.30630862
#> AB29-ABMX-H_S6_MERGED.fastq.gz   High     1851   0.33312036
#> 85 more rows ...
#> 
#> $genes
#>       Domain        Phylum          Class           Order         Family
#> ASV7   Fungi Basidiomycota Agaricomycetes      Russulales     Stereaceae
#> ASV8   Fungi Basidiomycota Agaricomycetes      Russulales     Stereaceae
#> ASV12  Fungi Basidiomycota Agaricomycetes Hymenochaetales Schizoporaceae
#> ASV18  Fungi Basidiomycota Agaricomycetes      Russulales     Stereaceae
#> ASV25  Fungi Basidiomycota Agaricomycetes      Agaricales  Lyophyllaceae
#>            Genus    Species Trophic.Mode                                Guild
#> ASV7        <NA>       <NA>   Saprotroph Wood Saprotroph-Undefined Saprotroph
#> ASV8     Stereum     ostrea   Saprotroph                 Undefined Saprotroph
#> ASV12    Xylodon raduloides   Saprotroph                 Undefined Saprotroph
#> ASV18    Stereum     ostrea   Saprotroph                 Undefined Saprotroph
#> ASV25 Ossicaulis  lachnopus   Saprotroph                      Wood Saprotroph
#>           Trait Confidence.Ranking        Genus_species
#> ASV7       NULL           Probable                NA_NA
#> ASV8  White Rot           Probable       Stereum_ostrea
#> ASV12 White Rot           Probable   Xylodon_raduloides
#> ASV18 White Rot           Probable       Stereum_ostrea
#> ASV25 Brown Rot           Probable Ossicaulis_lachnopus
#> 40 more rows ...
#> 
#> $common.dispersion
#> [1] 3.269437
#> 
#> $pseudo.counts
#>       A10-005-B_S188_MERGED.fastq.gz A10-005-H_S189_MERGED.fastq.gz
#> ASV7                     739.9845140                   9888.7370815
#> ASV8                     557.8797556                      9.7036208
#> ASV12                      0.9605827                      1.9385836
#> ASV18                     29.8644385                      0.9683215
#> ASV25                      0.9603193                      0.9678363
#>       A10-005-M_S190_MERGED.fastq.gz A12-007-B_S2_MERGED.fastq.gz
#> ASV7                     254.1481708                     1.023169
#> ASV8                     385.0731375                   847.747007
#> ASV12                      0.9594440                     1.023116
#> ASV18                      3.8464905                  1063.167359
#> ASV25                      0.9595064                     1.023198
#>       AB29-ABMX-H_S6_MERGED.fastq.gz AD26-005-B_S9_MERGED.fastq.gz
#> ASV7                       0.9689209                      1.023169
#> ASV8                       1.9400514                      1.023180
#> ASV12                      1.9397400                      1.023116
#> ASV18                      0.9689180                      1.023324
#> ASV25                      0.9684418                      1.023198
#>       AD26-005-H_S10_MERGED.fastq.gz AD26-005-M_S11_MERGED.fastq.gz
#> ASV7                        1.023085                      1.9399666
#> ASV8                        1.023087                      0.9689068
#> ASV12                       1.023234                   1300.5607631
#> ASV18                       1.023087                      0.9684377
#> ASV25                       1.023452                      0.9687233
#>       AD30-ABMX-M_S12_MERGED.fastq.gz AD32-007-M_S13_MERGED.fastq.gz
#> ASV7                         1.025971                      0.9682712
#> ASV8                         2.050025                      7.7623786
#> ASV12                        1.026152                      0.9680578
#> ASV18                        1.026342                      0.9678321
#> ASV25                        1.026014                      0.9681221
#>       ADABM30X-B_S14_MERGED.fastq.gz ADABM30X-H_S15_MERGED.fastq.gz
#> ASV7                       0.9604967                    524.1303861
#> ASV8                       2.8871709                      2.9094842
#> ASV12                      0.9605827                      0.9681280
#> ASV18                      0.9602419                      0.9683215
#> ASV25                      0.9603193                      2.9086060
#>       ADABM30X-M_S16_MERGED.fastq.gz B18-006-B_S19_MERGED.fastq.gz
#> ASV7                        1.020353                      1.028159
#> ASV8                        1.020327                      1.028173
#> ASV12                       1.020494                      1.028094
#> ASV18                       1.020643                      1.028348
#> ASV25                       1.020393                      1.028186
#>       BA17-050-B_S21_MERGED.fastq.gz BB19-006-H_S22_MERGED.fastq.gz
#> ASV7                    2800.3646390                       1.020317
#> ASV8                       4.8502866                    9925.530722
#> ASV12                      0.9682836                       1.020448
#> ASV18                      0.9680080                   12189.387541
#> ASV25                      0.9680837                       1.020640
#>       BB6-019-M_S25_MERGED.fastq.gz BG7-010-H_S31_MERGED.fastq.gz
#> ASV7                    848.3304311                      1.028057
#> ASV8                    327.0949300                      1.028059
#> ASV12                     0.9680578                      1.028238
#> ASV18                     7.7614518                      1.028059
#> ASV25                     0.9681221                      1.028504
#>       BJ17-007-M_S34_MERGED.fastq.gz BL7-006-H_S37_MERGED.fastq.gz
#> ASV7                       0.9513405                      1.028057
#> ASV8                       6.6806803                      3.080174
#> ASV12                      6.6795323                      1.028238
#> ASV18                      0.9506736                      1.028059
#> ASV25                      1.9092417                      1.028504
#>       BP11-001-B_S43_MERGED.fastq.gz BP11-001-H_S44_MERGED.fastq.gz
#> ASV7                        1.023169                    958.6644070
#> ASV8                        2.044684                      8.6682003
#> ASV12                       2.044579                      0.9603903
#> ASV18                       1.023324                   1778.5974801
#> ASV25                       1.023198                      0.9600293
#>       BP11-001-M_S45_MERGED.fastq.gz BP12-025-B_S46_MERGED.fastq.gz
#> ASV7                       6.7955611                       1.028159
#> ASV8                       6.7957113                       2.054307
#> ASV12                      0.9686592                       1.028094
#> ASV18                      0.9684377                       1.028348
#> ASV25                      0.9687233                       1.028186
#>       BQ4-018-B_S49_MERGED.fastq.gz BQ4-018-H_S50_MERGED.fastq.gz
#> ASV7                       1.028159                      1.023085
#> ASV8                       2.054307                      2.044529
#> ASV12                      1.028094                      1.023234
#> ASV18                      1.028348                      1.023087
#> ASV25                      1.028186                      2.045074
#>       BQ4-018-M_S51_MERGED.fastq.gz BT-006-M_S55_MERGED.fastq.gz
#> ASV7                       1.023126                     1.025971
#> ASV8                       1.023096                     1.025938
#> ASV12                      1.023287                     1.026152
#> ASV18                      1.023455                     1.026342
#> ASV25                      1.023168                     1.026014
#>       BV11-002-B_S57_MERGED.fastq.gz BV11-002-H_S58_MERGED.fastq.gz
#> ASV7                       0.9688125                      0.9689209
#> ASV8                       3.8820095                      2.9111887
#> ASV12                   7333.4712076                      0.9687281
#> ASV18                      0.9686103                      0.9689180
#> ASV25                      0.9686856                     65.0759416
#>       BV11-002-M_S59_MERGED.fastq.gz C21-NV1-B_S62_MERGED.fastq.gz
#> ASV7                        1.025971                      1.023169
#> ASV8                        1.025938                      1.023180
#> ASV12                       1.026152                      1.023116
#> ASV18                       1.026342                      1.023324
#> ASV25                       1.026014                      1.023198
#>       CB8-019-B_S69_MERGED.fastq.gz CB8-019-H_S70_MERGED.fastq.gz
#> ASV7                       1.026019                      1.025925
#> ASV8                       1.026032                      1.025928
#> ASV12                      2.050063                      1.026092
#> ASV18                      1.026194                      1.025928
#> ASV25                      1.026048                      1.026338
#>       CB8-019-M_S71_MERGED.fastq.gz D18-003-B_S78_MERGED.fastq.gz
#> ASV7                       1.025971                    38.5348734
#> ASV8                       1.025938                    30.8270362
#> ASV12                      3.074539                     0.9605827
#> ASV18                      1.026342                     0.9602419
#> ASV25                      1.026014                     0.9603193
#>       D18-003-M_S80_MERGED.fastq.gz D61-010-B_S82_MERGED.fastq.gz
#> ASV7                       1.025971                     0.9688125
#> ASV8                       1.025938                     5.8242097
#> ASV12                      2.050366                     1.9399871
#> ASV18                      1.026342                     0.9686103
#> ASV25                      1.026014                     4.8647768
#>       D9-027-H_S84_MERGED.fastq.gz D9-027-M_S85_MERGED.fastq.gz
#> ASV7                      1.028057                     1.028106
#> ASV8                      4.106226                     2.054137
#> ASV12                     1.028238                     1.028302
#> ASV18                     1.028059                     1.028508
#> ASV25                     1.028504                     1.028149
#>       DJ2-008-H_S88_MERGED.fastq.gz DS1-ABM002-B_S91_MERGED.fastq.gz
#> ASV7                       1.020317                         1.026019
#> ASV8                       1.020319                         1.026032
#> ASV12                      1.020448                         1.025960
#> ASV18                      1.020319                         1.026194
#> ASV25                      1.020640                         1.026048
#>       DS1-ABM002-H_S92_MERGED.fastq.gz DS1-ABM002-M_S93_MERGED.fastq.gz
#> ASV7                         0.9683244                        0.9682712
#> ASV8                         0.9683218                        0.9683101
#> ASV12                        0.9681280                        0.9680578
#> ASV18                        0.9683215                        0.9678321
#> ASV25                        1.9381709                        3.8883517
#>       DU3-045-B_S94_MERGED.fastq.gz DY5-004-B_S96_MERGED.fastq.gz
#> ASV7                     32.7270259                      1.023169
#> ASV8                   1067.6859062                      1.023180
#> ASV12                     7.6981393                      1.023116
#> ASV18                     7.6971939                      1.023324
#> ASV25                     0.9594583                      1.023198
#>       DY5-004-H_S97_MERGED.fastq.gz EC2-013-B_S104_MERGED.fastq.gz
#> ASV7                       1.025925                      0.9688125
#> ASV8                       1.025928                      1.9398500
#> ASV12                      1.026092                      0.9688808
#> ASV18                      1.025928                      0.9686103
#> ASV25                      1.026338                      0.9686856
#>       F7-015-M_S106_MERGED.fastq.gz FOMES19-H_S108_MERGED.fastq.gz
#> ASV7                       1.025971                       2.039187
#> ASV8                       1.025938                       1.020319
#> ASV12                      1.026152                       1.020448
#> ASV18                      1.026342                       1.020319
#> ASV25                      1.026014                       1.020640
#>       FOMES19-M_S109_MERGED.fastq.gz H10-018-M_S110_MERGED.fastq.gz
#> ASV7                        1.023126                       1.028106
#> ASV8                        1.023096                       2.054137
#> ASV12                       1.023287                       1.028302
#> ASV18                       1.023455                       1.028508
#> ASV25                       1.023168                       1.028149
#>       H24-NVABM1-H_S111_MERGED.fastq.gz J18-004-B_S114_MERGED.fastq.gz
#> ASV7                         106.981018                       1.026019
#> ASV8                           1.020319                       1.026032
#> ASV12                          3.058330                       1.025960
#> ASV18                          1.020319                       1.026194
#> ASV25                         10.190628                       2.048230
#>       J18-004-M_S116_MERGED.fastq.gz K18-002-H_S117_MERGED.fastq.gz
#> ASV7                        1.028106                     34.9592516
#> ASV8                        1.028070                      0.9689183
#> ASV12                       1.028302                      0.9687281
#> ASV18                       1.028508                     33.0169042
#> ASV25                       1.028149                      0.9684418
#>       L19X-B_S119_MERGED.fastq.gz L19X-H_S120_MERGED.fastq.gz
#> ASV7                     1.026019                 326.1238049
#> ASV8                   187.400182                  35.9098544
#> ASV12                    1.025960                   0.9681280
#> ASV18                    1.026194                   9.7036189
#> ASV25                    1.026048                   0.9678363
#>       L19X-M_S121_MERGED.fastq.gz L23-002-B_S122_MERGED.fastq.gz
#> ASV7                 2698.7914099                      6.7407990
#> ASV8                  428.7466466                     32.7539900
#> ASV12                   0.9603035                      0.9605827
#> ASV18                  21.1940575                     26.9737587
#> ASV25                   0.9603664                      0.9603193
#>       L23-002-H_S123_MERGED.fastq.gz L23-002-M_S124_MERGED.fastq.gz
#> ASV7                    1949.4596380                       1.023126
#> ASV8                    2032.2538022                     155.257416
#> ASV12                      0.9595327                       1.023287
#> ASV18                    100.1154865                      38.812658
#> ASV25                      0.9591640                       1.023168
#>       N19X-H_S127_MERGED.fastq.gz N19X-M_S128_MERGED.fastq.gz
#> ASV7                 1440.3991066                 166.9431303
#> ASV8                 2249.9023873                 273.7107127
#> ASV12                   0.9681280                   0.9680578
#> ASV18                 162.0894130                  11.6441341
#> ASV25                   0.9678363                   0.9681221
#>       N23-002-B_S130_MERGED.fastq.gz N23-002-H_S131_MERGED.fastq.gz
#> ASV7                        1.026019                       1.023085
#> ASV8                        1.026032                       2.044529
#> ASV12                       1.025960                       1.023234
#> ASV18                       1.026194                       1.023087
#> ASV25                       1.026048                     446.265835
#>       NVABM-0163-H_S135_MERGED.fastq.gz NVABM0244-M_S137_MERGED.fastq.gz
#> ASV7                           1.025925                         1.025971
#> ASV8                           1.025928                         1.025938
#> ASV12                          1.026092                         1.026152
#> ASV18                          1.025928                         1.026342
#> ASV25                          1.026338                         1.026014
#>       O20-X-H_S140_MERGED.fastq.gz O20-X-M_S141_MERGED.fastq.gz
#> ASV7                     98.190154                 2049.0165024
#> ASV8                     22.138363                  857.0565881
#> ASV12                  1520.257340                    3.8794388
#> ASV18                     3.847720                   67.9518495
#> ASV25                     0.959164                    0.9681221
#>       O24-003-B_S145_MERGED.fastq.gz O26-004-B_S148_MERGED.fastq.gz
#> ASV7                        1.028159                       1.023169
#> ASV8                        3.080422                       1.023180
#> ASV12                       1.028094                       1.023116
#> ASV18                       1.028348                       1.023324
#> ASV25                       1.028186                       1.023198
#>       O9-005-B_S152_MERGED.fastq.gz P19-023-M_S153_MERGED.fastq.gz
#> ASV7                      0.9682140                       2.044594
#> ASV8                      1.9386956                      27.580317
#> ASV12                     1.9388353                       1.023287
#> ASV18                     0.9680080                       1.023455
#> ASV25                     0.9680837                       1.023168
#>       P27-015-M_S154_MERGED.fastq.gz Q27-ABM003-B_S156_MERGED.fastq.gz
#> ASV7                        1.020353                          1.026019
#> ASV8                      196.639498                        483.338087
#> ASV12                       8.152919                          1.025960
#> ASV18                     282.172711                          1.026194
#> ASV25                       1.020393                          1.026048
#>       T28-ABM602-B_S162_MERGED.fastq.gz W26-001-B_S165_MERGED.fastq.gz
#> ASV7                           1.023169                       1.028159
#> ASV8                           1.023180                       1.028173
#> ASV12                          1.023116                       1.028094
#> ASV18                          1.023324                       1.028348
#> ASV25                          1.023198                       1.028186
#>       W9-025-M_S169_MERGED.fastq.gz X24-009-B_S170_MERGED.fastq.gz
#> ASV7                       1.028106                      0.9688125
#> ASV8                       1.028070                      0.9687977
#> ASV12                      1.028302                      0.9688808
#> ASV18                      1.028508                      0.9686103
#> ASV25                      1.028149                      0.9686856
#>       X24-009-H_S171_MERGED.fastq.gz X24-009-M_S172_MERGED.fastq.gz
#> ASV7                       1.9389054                      0.9605674
#> ASV8                       5.8212488                      7.7046802
#> ASV12                      0.9681280                      0.9603035
#> ASV18                      0.9683215                      0.9600240
#> ASV25                      0.9678363                      0.9603664
#>       X29-004-B_S174_MERGED.fastq.gz X29-004-H_S175_MERGED.fastq.gz
#> ASV7                        1.028159                    104.0517756
#> ASV8                        1.028173                    108.8691295
#> ASV12                       1.028094                      0.9603903
#> ASV18                       1.028348                      1.9240643
#> ASV25                       1.028186                      0.9600293
#>       X29-004-M_S176_MERGED.fastq.gz Y21-ABM484-H_S177_MERGED.fastq.gz
#> ASV7                        1.023126                          1.020317
#> ASV8                        3.065988                          2.039190
#> ASV12                       1.023287                        254.704323
#> ASV18                       1.023455                          1.020319
#> ASV25                       1.023168                          1.020640
#>       Y28-002-B_S178_MERGED.fastq.gz Y31-ABM484-B_S184_MERGED.fastq.gz
#> ASV7                        1.028159                        86.6383224
#> ASV8                        1.028173                        50.0555269
#> ASV12                       1.028094                       606.4994165
#> ASV18                       1.028348                         2.8842769
#> ASV25                       1.028186                         0.9594583
#>       Z29-001-H_S185_MERGED.fastq.gz Z30-ABM560-M_S187_MERGED.fastq.gz
#> ASV7                        1.028057                       269.0126340
#> ASV8                        1.028059                      9118.4411470
#> ASV12                       1.028238                         0.9686592
#> ASV18                       1.028059                       126.2769656
#> ASV25                       1.028504                         0.9687233
#> 40 more rows ...
#> 
#> $pseudo.lib.size
#> [1] 598.8289
#> 
#> $AveLogCPM
#> [1] 18.85258 19.06397 17.66999 18.18946 13.84983
#> 40 more elements ...
#> 
#> $prior.df
#> [1] 10
#> 
#> $prior.n
#> [1] 0.1149425
#> 
#> $tagwise.dispersion
#> [1] 5.655012 5.359018 5.694498 4.969208 1.650035
#> 40 more elements ...
#> 
#> $span
#> [1] 0.3
#> 
# }
```
