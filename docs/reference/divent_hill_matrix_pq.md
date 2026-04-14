# Compute Hill diversity numbers for all samples in an OTU table

Iterates over all samples in an OTU table and computes Hill diversity
numbers using
[`divent::div_hill()`](https://ericmarcon.github.io/divent/reference/div_hill.html).

## Usage

``` r
divent_hill_matrix_pq(comm, q, ...)
```

## Arguments

- comm:

  (data.frame or matrix) OTU table with samples as rows and taxa as
  columns.

- q:

  (numeric vector) Hill diversity orders to compute. Hill numbers are
  more appropriate in DNA metabarcoding studies when `q > 0` (Alberdi &
  Gilbert, 2019; Calderón-Sanou et al., 2019).

- ...:

  Additional arguments passed to
  [`divent::div_hill()`](https://ericmarcon.github.io/divent/reference/div_hill.html)
  (e.g. `estimator = "naive"` to reproduce vegan-equivalent results).

## Value

A data.frame with one row per sample and one column per value in `q`.
Column names are the string representation of the `q` values. Row names
match the input row names.

## References

Alberdi, A., & Gilbert, M. T. P. (2019). A guide to the application of
Hill numbers to DNA-based diversity analyses. *Molecular Ecology
Resources*.
[doi:10.1111/1755-0998.13014](https://doi.org/10.1111/1755-0998.13014)

Calderón-Sanou, I., Münkemüller, T., Boyer, F., Zinger, L., & Thuiller,
W. (2019). From environmental DNA sequences to ecological conclusions:
How strong is the influence of methodological choices? *Journal of
Biogeography*, 47.
[doi:10.1111/jbi.13681](https://doi.org/10.1111/jbi.13681)

## See also

[`divent::div_hill()`](https://ericmarcon.github.io/divent/reference/div_hill.html),
[`hill_pq()`](https://adrientaudiere.github.io/MiscMetabar/reference/hill_pq.md),
[`hill_tuckey_pq()`](https://adrientaudiere.github.io/MiscMetabar/reference/hill_tuckey_pq.md)

## Examples

``` r
data("data_fungi_mini", package = "MiscMetabar")
otu <- as.data.frame(phyloseq::otu_table(
  taxa_as_columns(data_fungi_mini)
))
#> Taxa are now in columns.
divent_hill_matrix_pq(otu, q = c(0, 1, 2))
#> ! Sample coverage is 0, most estimators will return `NaN`.
#> ! Sample coverage is 0, most estimators will return `NaN`.
#> ! Sample coverage is 0, most estimators will return `NaN`.
#> ! Sample coverage is 0, most estimators will return `NaN`.
#> ! Sample coverage is 0, most estimators will return `NaN`.
#> ! Sample coverage is 0, most estimators will return `NaN`.
#>                                    0        1        2
#> A10-005-B_S188_MERGED.fastq.gz     6 3.778681 3.321022
#> A10-005-H_S189_MERGED.fastq.gz     7 1.012331 1.002750
#> A10-005-M_S190_MERGED.fastq.gz     7 4.151787 3.664966
#> A12-007_S191_MERGED.fastq.gz      10 1.658119 1.226407
#> A12-007-B_S2_MERGED.fastq.gz       2 1.987268 1.974830
#> A15-004_S3_MERGED.fastq.gz         3 2.233885 1.840764
#> A8-005_S4_MERGED.fastq.gz          6 2.458596 1.934101
#> AB29-ABMX-H_S6_MERGED.fastq.gz     8 2.178496 1.848196
#> AC27-013_S7_MERGED.fastq.gz        5 1.048577 1.014673
#> AC29033_S8_MERGED.fastq.gz         5 1.209420 1.078992
#> AD26-005-B_S9_MERGED.fastq.gz      5 1.617010 1.272737
#> AD26-005-H_S10_MERGED.fastq.gz     2 1.060478 1.021415
#> AD26-005-M_S11_MERGED.fastq.gz     7 2.824963 2.685002
#> AD30-ABMX-M_S12_MERGED.fastq.gz    4 1.037957 1.010218
#> AD32-007-M_S13_MERGED.fastq.gz     5 1.910754 1.602471
#> ADABM30X-B_S14_MERGED.fastq.gz     8 2.033137 1.850398
#> ADABM30X-H_S15_MERGED.fastq.gz    10 2.474791 1.864346
#> ADABM30X-M_S16_MERGED.fastq.gz     3 1.939813 1.816966
#> AE30-ABM507_S17_MERGED.fastq.gz    5 2.306849 1.753408
#> B17-014_S18_MERGED.fastq.gz        7 1.263748 1.094259
#> B18-006-B_S19_MERGED.fastq.gz      2 1.643473 1.463907
#> BA16-036bis_S20_MERGED.fastq.gz    7 3.950072 3.881949
#> BA17-050-B_S21_MERGED.fastq.gz     7 2.556756 2.269958
#> BB19-006-H_S22_MERGED.fastq.gz     3 2.065895 2.007785
#> BB6-019-M_S25_MERGED.fastq.gz      5 2.242332 1.879950
#> BD14-021_S26_MERGED.fastq.gz       9 2.843831 2.114546
#> BG7-010-H_S31_MERGED.fastq.gz      1 1.000000 1.000000
#> BH9-021_S33_MERGED.fastq.gz        1 1.000000 1.000000
#> BJ17-007-M_S34_MERGED.fastq.gz    24 2.432108 2.234763
#> BJ8-ABM-003_S35_MERGED.fastq.gz    4 2.090369 1.808210
#> BL7-006-H_S37_MERGED.fastq.gz      1 1.000000 1.000000
#> BN11-041_S39_MERGED.fastq.gz       7 2.214998 2.098953
#> BO8-002_S41_MERGED.fastq.gz        8 1.415531 1.155352
#> BO8-005_S42_MERGED.fastq.gz        1 1.000000 1.000000
#> BP11-001-B_S43_MERGED.fastq.gz     6 1.148648 1.062169
#> BP11-001-H_S44_MERGED.fastq.gz     7 3.834000 3.351706
#> BP11-001-M_S45_MERGED.fastq.gz     5 1.374991 1.206796
#> BP12-025-B_S46_MERGED.fastq.gz     1 1.000000 1.000000
#> BP14-006_S47_MERGED.fastq.gz       8 6.583211 6.173399
#> BQ3-019_S48_MERGED.fastq.gz        1 1.000000 1.000000
#> BQ4-018-B_S49_MERGED.fastq.gz      7 1.003077 1.000606
#> BQ4-018-H_S50_MERGED.fastq.gz      5 1.002235 1.000448
#> BQ4-018-M_S51_MERGED.fastq.gz      6 1.647788 1.410376
#> BQ9ABM-002_S52_MERGED.fastq.gz     3 1.130692 1.053861
#> BR8-005_S53_MERGED.fastq.gz        1 1.000000 1.000000
#> BS14-006_S54_MERGED.fastq.gz       7 2.323191 2.119591
#> BT-006-M_S55_MERGED.fastq.gz       2 1.134309 1.056558
#> BT7-006_S56_MERGED.fastq.gz        1 1.000000 1.000000
#> BV11-002-B_S57_MERGED.fastq.gz     7 2.469983 2.219008
#> BV11-002-H_S58_MERGED.fastq.gz     6 3.006670 2.554421
#> BV11-002-M_S59_MERGED.fastq.gz     3 1.790660 1.433604
#> BW8-003_S60_MERGED.fastq.gz        7 2.531708 1.865068
#> C1-001_S61_MERGED.fastq.gz         7 2.854908 2.421389
#> C21-NV1-B_S62_MERGED.fastq.gz      3 1.297132 1.130844
#> C9-005_S65_MERGED.fastq.gz         7 2.302059 1.647721
#> CA9-X_S68_MERGED.fastq.gz          4 1.298521 1.127462
#> CB8-019-B_S69_MERGED.fastq.gz      1 1.000000 1.000000
#> CB8-019-H_S70_MERGED.fastq.gz      1 1.000000 1.000000
#> CB8-019-M_S71_MERGED.fastq.gz      1 1.000000 1.000000
#> CB9-013_S72_MERGED.fastq.gz       10 3.734605 2.628309
#> CC3-044_S73_MERGED.fastq.gz        6 3.929486 3.445494
#> D17-011_S77_MERGED.fastq.gz        7 1.514946 1.213532
#> D18-003-B_S78_MERGED.fastq.gz      6 1.422399 1.159295
#> D18-003-M_S80_MERGED.fastq.gz      7 2.476772 1.915471
#> D22-NVABM1_S81_MERGED.fastq.gz     5 1.471269 1.249903
#> D61-010-B_S82_MERGED.fastq.gz      4 1.005035 1.001059
#> D9-027-H_S84_MERGED.fastq.gz       1 1.000000 1.000000
#> D9-027-M_S85_MERGED.fastq.gz       1 1.000000 1.000000
#> DBM-ABM-001_S86_MERGED.fastq.gz    7 2.833076 2.403005
#> DJ2-008-H_S88_MERGED.fastq.gz      3 1.472166 1.221700
#> DP4-ABM001_S90_MERGED.fastq.gz     2 1.005814 1.001404
#> DS1-ABM002-B_S91_MERGED.fastq.gz   1 1.000000 1.000000
#> DS1-ABM002-H_S92_MERGED.fastq.gz   6 1.886016 1.741623
#> DS1-ABM002-M_S93_MERGED.fastq.gz   5 3.583387 2.825210
#> DU3-045-B_S94_MERGED.fastq.gz      7 1.629168 1.352777
#> DW4-007_S95_MERGED.fastq.gz        3 1.005002 1.001084
#> DY5-004-B_S96_MERGED.fastq.gz      2 1.979626 1.960000
#> DY5-004-H_S97_MERGED.fastq.gz      1 1.000000 1.000000
#> DZ6-ABM-001_S99_MERGED.fastq.gz    7 1.011354 1.002605
#> EA5-ABM-001_S103_MERGED.fastq.gz   3 1.017959 1.004920
#> EC2-013-B_S104_MERGED.fastq.gz     6 1.076665 1.026479
#> F6-ABM-001_S105_MERGED.fastq.gz    6 1.129525 1.043131
#> F7-015-M_S106_MERGED.fastq.gz      1 1.000000 1.000000
#> FOMES19-H_S108_MERGED.fastq.gz     4 1.225468 1.102917
#> FOMES19-M_S109_MERGED.fastq.gz     3 1.005605 1.001251
#> H10-018-M_S110_MERGED.fastq.gz     1 1.000000 1.000000
#> H24-NVABM1-H_S111_MERGED.fastq.gz  3 1.434442 1.213191
#> J18-004-B_S114_MERGED.fastq.gz     1 1.000000 1.000000
#> J18-004-M_S116_MERGED.fastq.gz     1 1.000000 1.000000
#> K18-002-H_S117_MERGED.fastq.gz     5 3.414645 2.910891
#> K26-NVABM1_S118_MERGED.fastq.gz    9 2.200540 2.017670
#> L19X-B_S119_MERGED.fastq.gz        1 1.000000 1.000000
#> L19X-H_S120_MERGED.fastq.gz        5 3.038439 2.368687
#> L19X-M_S121_MERGED.fastq.gz        6 2.938505 2.493204
#> L23-002-B_S122_MERGED.fastq.gz     6 1.982746 1.463123
#> L23-002-H_S123_MERGED.fastq.gz     7 3.066403 2.672528
#> L23-002-M_S124_MERGED.fastq.gz     3 1.692714 1.477845
#> N19X-H_S127_MERGED.fastq.gz        5 2.996599 2.609700
#> N19X-M_S128_MERGED.fastq.gz        5 2.927734 2.567688
#> N23-002-B_S130_MERGED.fastq.gz     3 2.667090 2.417902
#> N23-002-H_S131_MERGED.fastq.gz     6 1.081911 1.023059
#> N25-ABMX_S133_MERGED.fastq.gz      4 1.004941 1.001027
#> NVABM-0058_S134_MERGED.fastq.gz    4 1.033902 1.009856
#> NVABM-0163-H_S135_MERGED.fastq.gz  2 1.035493 1.011363
#> NVABM-0397_S138_MERGED.fastq.gz    2 1.499248 1.317176
#> NVABM0216_S136_MERGED.fastq.gz     3 2.377973 2.008688
#> NVABM0244-M_S137_MERGED.fastq.gz   1 1.000000 1.000000
#> O20-X-H_S140_MERGED.fastq.gz       7 2.812120 2.400615
#> O20-X-M_S141_MERGED.fastq.gz       5 2.729337 2.284681
#> O24-003-B_S145_MERGED.fastq.gz     1 1.000000 1.000000
#> O26-004-B_S148_MERGED.fastq.gz     2 1.369564 1.208219
#> O27-012_S151_MERGED.fastq.gz       1 1.000000 1.000000
#> O9-005-B_S152_MERGED.fastq.gz      8 2.170883 1.785130
#> P19-023-M_S153_MERGED.fastq.gz     3 2.152862 2.060241
#> P27-015-M_S154_MERGED.fastq.gz     3 2.103820 1.995110
#> P27-ABM001_S155_MERGED.fastq.gz    2 2.000000 2.000000
#> Q27-ABM003-B_S156_MERGED.fastq.gz  1 1.000000 1.000000
#> R25-ABMX_S157_MERGED.fastq.gz      7 3.614251 2.618286
#> T28-011_S161_MERGED.fastq.gz       5 1.427775 1.251109
#> T28-ABM602-B_S162_MERGED.fastq.gz  2 1.832713 1.710059
#> U27-ABM002_S163_MERGED.fastq.gz    7 3.010872 2.503639
#> W25-ABMX_S164_MERGED.fastq.gz      4 1.484866 1.222420
#> W26-001-B_S165_MERGED.fastq.gz     1 1.000000 1.000000
#> W30-006_S168_MERGED.fastq.gz      10 3.584496 3.019785
#> W9-025-M_S169_MERGED.fastq.gz      2 2.000000 2.000000
#> X24-009-B_S170_MERGED.fastq.gz     4 3.989820 3.979741
#> X24-009-H_S171_MERGED.fastq.gz     8 4.508687 4.245029
#> X24-009-M_S172_MERGED.fastq.gz     6 3.805046 3.561964
#> X24-010_S173_MERGED.fastq.gz       2 1.000956 1.000186
#> X29-004-B_S174_MERGED.fastq.gz     1 1.000000 1.000000
#> X29-004-H_S175_MERGED.fastq.gz     9 4.781716 4.427453
#> X29-004-M_S176_MERGED.fastq.gz     5 2.037983 1.961554
#> Y21-ABM484-H_S177_MERGED.fastq.gz  5 2.542504 2.130989
#> Y28-002-B_S178_MERGED.fastq.gz     2 1.417411 1.246154
#> Y31-ABM484-B_S184_MERGED.fastq.gz  7 2.148479 1.566315
#> Z29-001-H_S185_MERGED.fastq.gz     3 2.667090 2.417902
#> Z30-ABM560-M_S187_MERGED.fastq.gz  4 1.339906 1.131557
```
