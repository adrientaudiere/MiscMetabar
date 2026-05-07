# Geometric Mean of Pairwise Ratios (GMPR) normalization of a phyloseq object

[![lifecycle-experimental](https://img.shields.io/badge/lifecycle-experimental-orange)](https://adrientaudiere.github.io/MiscMetabar/articles/Rules.html#lifecycle)

Pure-R implementation of the Geometric Mean of Pairwise Ratios
normalization (Chen et al. 2018,
[doi:10.7717/peerj.4600](https://doi.org/10.7717/peerj.4600) ) tailored
for zero-inflated count tables such as microbial OTU tables. Returns
counts divided by the per-sample GMPR size factors.

## Usage

``` r
gmpr_pq(physeq, intersect_no = 4, ct_min = 2)
```

## Arguments

- physeq:

  (required) a
  [`phyloseq-class`](https://rdrr.io/pkg/phyloseq/man/phyloseq-class.html)
  object obtained using the `phyloseq` package.

- intersect_no:

  (integer, default 4) minimum number of shared taxa between two samples
  for the pairwise ratio to be computed.

- ct_min:

  (integer, default 2) minimum count for a taxon to be considered
  "shared" between two samples.

## Value

A new
[`phyloseq-class`](https://rdrr.io/pkg/phyloseq/man/phyloseq-class.html)
object with a GMPR-normalised `otu_table`. Size factors are stored as an
attribute `"gmpr_size_factors"` on the otu_table.

## References

Chen L. et al. (2018) GMPR: a robust normalization method for
zero-inflated count data with application to microbiome sequencing data.
PeerJ 6:e4600.
[doi:10.7717/peerj.4600](https://doi.org/10.7717/peerj.4600)

## Author

Adrien Taudière

## Examples

``` r
data_f_gmpr <- gmpr_pq(data_fungi_mini)
#> Warning: GMPR size factors could not be computed for 90 sample(s); these samples are left unscaled.
sample_sums(data_f_gmpr)
#>    A10-005-B_S188_MERGED.fastq.gz    A10-005-H_S189_MERGED.fastq.gz 
#>                         1506.3028                        10201.0000 
#>    A10-005-M_S190_MERGED.fastq.gz      A12-007_S191_MERGED.fastq.gz 
#>                         2163.9241                        10267.8240 
#>      A12-007-B_S2_MERGED.fastq.gz        A15-004_S3_MERGED.fastq.gz 
#>                         1869.0000                           17.0000 
#>         A8-005_S4_MERGED.fastq.gz    AB29-ABMX-H_S6_MERGED.fastq.gz 
#>                         1736.0591                         1806.0000 
#>       AC27-013_S7_MERGED.fastq.gz        AC29033_S8_MERGED.fastq.gz 
#>                         8247.0000                           54.0000 
#>     AD26-005-B_S9_MERGED.fastq.gz    AD26-005-H_S10_MERGED.fastq.gz 
#>                           18.0000                         4719.0000 
#>    AD26-005-M_S11_MERGED.fastq.gz   AD30-ABMX-M_S12_MERGED.fastq.gz 
#>                         3615.0000                          592.0000 
#>    AD32-007-M_S13_MERGED.fastq.gz    ADABM30X-B_S14_MERGED.fastq.gz 
#>                          796.3391                         6982.0863 
#>    ADABM30X-H_S15_MERGED.fastq.gz    ADABM30X-M_S16_MERGED.fastq.gz 
#>                         1187.6686                          894.0000 
#>   AE30-ABM507_S17_MERGED.fastq.gz       B17-014_S18_MERGED.fastq.gz 
#>                         1655.0000                         6670.9711 
#>     B18-006-B_S19_MERGED.fastq.gz   BA16-036bis_S20_MERGED.fastq.gz 
#>                          233.0000                        38335.0000 
#>    BA17-050-B_S21_MERGED.fastq.gz    BB19-006-H_S22_MERGED.fastq.gz 
#>                         3165.6960                        21861.0000 
#>     BB6-019-M_S25_MERGED.fastq.gz      BD14-021_S26_MERGED.fastq.gz 
#>                         2206.2171                          906.6284 
#>     BG7-010-H_S31_MERGED.fastq.gz       BH9-021_S33_MERGED.fastq.gz 
#>                            1.0000                            1.0000 
#>    BJ17-007-M_S34_MERGED.fastq.gz   BJ8-ABM-003_S35_MERGED.fastq.gz 
#>                        61090.7241                         1360.0000 
#>     BL7-006-H_S37_MERGED.fastq.gz      BN11-041_S39_MERGED.fastq.gz 
#>                            2.0000                        22649.5899 
#>       BO8-002_S41_MERGED.fastq.gz       BO8-005_S42_MERGED.fastq.gz 
#>                          695.0000                            1.0000 
#>    BP11-001-B_S43_MERGED.fastq.gz    BP11-001-H_S44_MERGED.fastq.gz 
#>                         3316.0000                         1508.8342 
#>    BP11-001-M_S45_MERGED.fastq.gz    BP12-025-B_S46_MERGED.fastq.gz 
#>                        15602.0000                            1.0000 
#>      BP14-006_S47_MERGED.fastq.gz       BQ3-019_S48_MERGED.fastq.gz 
#>                         1666.5322                            1.0000 
#>     BQ4-018-B_S49_MERGED.fastq.gz     BQ4-018-H_S50_MERGED.fastq.gz 
#>                         9898.0000                         8929.0000 
#>     BQ4-018-M_S51_MERGED.fastq.gz    BQ9ABM-002_S52_MERGED.fastq.gz 
#>                          440.0000                        15818.0000 
#>       BR8-005_S53_MERGED.fastq.gz      BS14-006_S54_MERGED.fastq.gz 
#>                            1.0000                         5999.7248 
#>      BT-006-M_S55_MERGED.fastq.gz       BT7-006_S56_MERGED.fastq.gz 
#>                          109.0000                           43.0000 
#>    BV11-002-B_S57_MERGED.fastq.gz    BV11-002-H_S58_MERGED.fastq.gz 
#>                         3018.0074                          681.0978 
#>    BV11-002-M_S59_MERGED.fastq.gz       BW8-003_S60_MERGED.fastq.gz 
#>                           23.0000                         4189.6647 
#>        C1-001_S61_MERGED.fastq.gz     C21-NV1-B_S62_MERGED.fastq.gz 
#>                         3342.0391                           82.0000 
#>        C9-005_S65_MERGED.fastq.gz         CA9-X_S68_MERGED.fastq.gz 
#>                         1491.0000                         1207.2588 
#>     CB8-019-B_S69_MERGED.fastq.gz     CB8-019-H_S70_MERGED.fastq.gz 
#>                            1.0000                            1.0000 
#>     CB8-019-M_S71_MERGED.fastq.gz       CB9-013_S72_MERGED.fastq.gz 
#>                            2.0000                         5236.2514 
#>       CC3-044_S73_MERGED.fastq.gz       D17-011_S77_MERGED.fastq.gz 
#>                         1193.5960                        20487.0585 
#>     D18-003-B_S78_MERGED.fastq.gz     D18-003-M_S80_MERGED.fastq.gz 
#>                         3395.2792                           68.0000 
#>    D22-NVABM1_S81_MERGED.fastq.gz     D61-010-B_S82_MERGED.fastq.gz 
#>                        13780.8455                        18890.0000 
#>      D9-027-H_S84_MERGED.fastq.gz      D9-027-M_S85_MERGED.fastq.gz 
#>                            3.0000                            1.0000 
#>   DBM-ABM-001_S86_MERGED.fastq.gz     DJ2-008-H_S88_MERGED.fastq.gz 
#>                         1123.4650                           51.0000 
#>    DP4-ABM001_S90_MERGED.fastq.gz  DS1-ABM002-B_S91_MERGED.fastq.gz 
#>                         7125.0000                           13.0000 
#>  DS1-ABM002-H_S92_MERGED.fastq.gz  DS1-ABM002-M_S93_MERGED.fastq.gz 
#>                         2434.0000                           41.0000 
#>     DU3-045-B_S94_MERGED.fastq.gz       DW4-007_S95_MERGED.fastq.gz 
#>                        55241.2897                        11080.0000 
#>     DY5-004-B_S96_MERGED.fastq.gz     DY5-004-H_S97_MERGED.fastq.gz 
#>                            7.0000                            1.0000 
#>   DZ6-ABM-001_S99_MERGED.fastq.gz  EA5-ABM-001_S103_MERGED.fastq.gz 
#>                         2307.0000                         6113.0000 
#>    EC2-013-B_S104_MERGED.fastq.gz   F6-ABM-001_S105_MERGED.fastq.gz 
#>                         8118.0000                          571.0000 
#>     F7-015-M_S106_MERGED.fastq.gz    FOMES19-H_S108_MERGED.fastq.gz 
#>                          520.0000                         2552.0000 
#>    FOMES19-M_S109_MERGED.fastq.gz    H10-018-M_S110_MERGED.fastq.gz 
#>                         6398.0000                            1.0000 
#> H24-NVABM1-H_S111_MERGED.fastq.gz    J18-004-B_S114_MERGED.fastq.gz 
#>                          115.0000                            1.0000 
#>    J18-004-M_S116_MERGED.fastq.gz    K18-002-H_S117_MERGED.fastq.gz 
#>                            4.0000                          106.2525 
#>   K26-NVABM1_S118_MERGED.fastq.gz       L19X-B_S119_MERGED.fastq.gz 
#>                         1875.5830                          182.0000 
#>       L19X-H_S120_MERGED.fastq.gz       L19X-M_S121_MERGED.fastq.gz 
#>                         1641.2297                         3608.3115 
#>    L23-002-B_S122_MERGED.fastq.gz    L23-002-H_S123_MERGED.fastq.gz 
#>                         5491.7380                          994.0011 
#>    L23-002-M_S124_MERGED.fastq.gz       N19X-H_S127_MERGED.fastq.gz 
#>                          189.0000                          797.6245 
#>       N19X-M_S128_MERGED.fastq.gz    N23-002-B_S130_MERGED.fastq.gz 
#>                         1113.5429                            3.0000 
#>    N23-002-H_S131_MERGED.fastq.gz     N25-ABMX_S133_MERGED.fastq.gz 
#>                          441.0000                        19484.0000 
#>   NVABM-0058_S134_MERGED.fastq.gz NVABM-0163-H_S135_MERGED.fastq.gz 
#>                        10611.0000                          177.0000 
#>   NVABM-0397_S138_MERGED.fastq.gz    NVABM0216_S136_MERGED.fastq.gz 
#>                           50.0000                          136.0000 
#>  NVABM0244-M_S137_MERGED.fastq.gz      O20-X-H_S140_MERGED.fastq.gz 
#>                            2.0000                        21324.6058 
#>      O20-X-M_S141_MERGED.fastq.gz    O24-003-B_S145_MERGED.fastq.gz 
#>                         1234.8654                            2.0000 
#>    O26-004-B_S148_MERGED.fastq.gz      O27-012_S151_MERGED.fastq.gz 
#>                           21.0000                            1.0000 
#>     O9-005-B_S152_MERGED.fastq.gz    P19-023-M_S153_MERGED.fastq.gz 
#>                         3888.0000                           57.0000 
#>    P27-015-M_S154_MERGED.fastq.gz   P27-ABM001_S155_MERGED.fastq.gz 
#>                          475.0000                            2.0000 
#> Q27-ABM003-B_S156_MERGED.fastq.gz     R25-ABMX_S157_MERGED.fastq.gz 
#>                          471.0000                         1122.6977 
#>      T28-011_S161_MERGED.fastq.gz T28-ABM602-B_S162_MERGED.fastq.gz 
#>                         9390.0000                           17.0000 
#>   U27-ABM002_S163_MERGED.fastq.gz     W25-ABMX_S164_MERGED.fastq.gz 
#>                         2244.7856                          203.0000 
#>    W26-001-B_S165_MERGED.fastq.gz      W30-006_S168_MERGED.fastq.gz 
#>                            1.0000                        29329.3969 
#>     W9-025-M_S169_MERGED.fastq.gz    X24-009-B_S170_MERGED.fastq.gz 
#>                            2.0000                        12526.4824 
#>    X24-009-H_S171_MERGED.fastq.gz    X24-009-M_S172_MERGED.fastq.gz 
#>                         2270.0862                         2036.6323 
#>      X24-010_S173_MERGED.fastq.gz    X29-004-B_S174_MERGED.fastq.gz 
#>                        10756.0000                            2.0000 
#>    X29-004-H_S175_MERGED.fastq.gz    X29-004-M_S176_MERGED.fastq.gz 
#>                         1462.7702                        12564.0000 
#> Y21-ABM484-H_S177_MERGED.fastq.gz    Y28-002-B_S178_MERGED.fastq.gz 
#>                         1194.0000                            9.0000 
#> Y31-ABM484-B_S184_MERGED.fastq.gz    Z29-001-H_S185_MERGED.fastq.gz 
#>                         5674.1376                            3.0000 
#> Z30-ABM560-M_S187_MERGED.fastq.gz 
#>                         3760.1777 
```
