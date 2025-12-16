# Compute and test local contributions to beta diversity (LCBD) of samples

[![lifecycle-experimental](https://img.shields.io/badge/lifecycle-experimental-orange)](https://adrientaudiere.github.io/MiscMetabar/articles/Rules.html#lifecycle)

A wrapper for the
[`adespatial::beta.div()`](http://adeverse.github.io/adespatial/reference/beta.div.md)
function in the case of `physeq` object.

## Usage

``` r
LCBD_pq(physeq, p_adjust_method = "BH", ...)
```

## Arguments

- physeq:

  (required) a
  [`phyloseq-class`](https://rdrr.io/pkg/phyloseq/man/phyloseq-class.html)
  object obtained using the `phyloseq` package.

- p_adjust_method:

  (chr, default "BH"): the method used to adjust p-value

- ...:

  Additional arguments passed on to
  [`adespatial::beta.div()`](http://adeverse.github.io/adespatial/reference/beta.div.md)
  function

## Value

An object of class `beta.div` see
[`adespatial::beta.div()`](http://adeverse.github.io/adespatial/reference/beta.div.md)
function for more information

## See also

[plot_LCBD_pq](https://adrientaudiere.github.io/MiscMetabar/reference/plot_LCBD_pq.md),
[`adespatial::beta.div()`](http://adeverse.github.io/adespatial/reference/beta.div.md)

## Author

Adrien Taudière This function is mainly a wrapper of the work of others.
Please make a reference to
[`adespatial::beta.div()`](http://adeverse.github.io/adespatial/reference/beta.div.md)
if you use this function.

## Examples

``` r
if (requireNamespace("adespatial")) {
  res <- LCBD_pq(data_fungi_sp_known, nperm = 5)
  str(res)
  length(res$LCBD)
  length(res$SCBD)
}
#> List of 8
#>  $ beta  : Named num [1:2] 173.703 0.944
#>   ..- attr(*, "names")= chr [1:2] "SStotal" "BDtotal"
#>  $ SCBD  : Named num [1:651] 0.0533 0.0307 0.0206 0.0123 0.0122 ...
#>   ..- attr(*, "names")= chr [1:651] "ASV2" "ASV8" "ASV12" "ASV18" ...
#>  $ LCBD  : Named num [1:185] 0.00522 0.00559 0.00545 0.00516 0.00488 ...
#>   ..- attr(*, "names")= chr [1:185] "A10-005-B_S188_MERGED.fastq.gz" "A10-005-H_S189_MERGED.fastq.gz" "A10-005-M_S190_MERGED.fastq.gz" "A12-007_S191_MERGED.fastq.gz" ...
#>  $ p.LCBD: Named num [1:185] 0.833 0.333 0.167 1 1 ...
#>   ..- attr(*, "names")= chr [1:185] "A10-005-B_S188_MERGED.fastq.gz" "A10-005-H_S189_MERGED.fastq.gz" "A10-005-M_S190_MERGED.fastq.gz" "A12-007_S191_MERGED.fastq.gz" ...
#>  $ p.adj : Named num [1:185] 1 0.709 0.434 1 1 ...
#>   ..- attr(*, "names")= chr [1:185] "A10-005-B_S188_MERGED.fastq.gz" "A10-005-H_S189_MERGED.fastq.gz" "A10-005-M_S190_MERGED.fastq.gz" "A12-007_S191_MERGED.fastq.gz" ...
#>  $ method: chr [1:2] "hellinger" NA
#>  $ note  : chr "Info -- This coefficient is Euclidean"
#>  $ D     : logi NA
#>  - attr(*, "class")= chr "beta.div"
#> [1] 651
# \donttest{
if (requireNamespace("adespatial")) {
  LCBD_pq(data_fungi_sp_known, nperm = 5, method = "jaccard")
}
#> $beta
#>    SStotal    BDtotal 
#> 86.5283720  0.4702629 
#> 
#> $SCBD
#> [1] NA
#> 
#> $LCBD
#>    A10-005-B_S188_MERGED.fastq.gz    A10-005-H_S189_MERGED.fastq.gz 
#>                       0.005398284                       0.005550668 
#>    A10-005-M_S190_MERGED.fastq.gz      A12-007_S191_MERGED.fastq.gz 
#>                       0.005444845                       0.005441464 
#>      A12-007-B_S2_MERGED.fastq.gz        A15-004_S3_MERGED.fastq.gz 
#>                       0.005461737                       0.005078194 
#>         A8-005_S4_MERGED.fastq.gz         A9-012_S5_MERGED.fastq.gz 
#>                       0.005360522                       0.005128489 
#>    AB29-ABMX-H_S6_MERGED.fastq.gz       AC27-013_S7_MERGED.fastq.gz 
#>                       0.005563478                       0.005330179 
#>        AC29033_S8_MERGED.fastq.gz     AD26-005-B_S9_MERGED.fastq.gz 
#>                       0.005402056                       0.005555455 
#>    AD26-005-H_S10_MERGED.fastq.gz    AD26-005-M_S11_MERGED.fastq.gz 
#>                       0.005506894                       0.005566377 
#>   AD30-ABMX-M_S12_MERGED.fastq.gz    AD32-007-M_S13_MERGED.fastq.gz 
#>                       0.005466833                       0.005382835 
#>    ADABM30X-B_S14_MERGED.fastq.gz    ADABM30X-H_S15_MERGED.fastq.gz 
#>                       0.005487170                       0.005481407 
#>    ADABM30X-M_S16_MERGED.fastq.gz   AE30-ABM507_S17_MERGED.fastq.gz 
#>                       0.005628821                       0.005384814 
#>       B17-014_S18_MERGED.fastq.gz     B18-006-B_S19_MERGED.fastq.gz 
#>                       0.005192392                       0.005094756 
#>   BA16-036bis_S20_MERGED.fastq.gz    BA17-050-B_S21_MERGED.fastq.gz 
#>                       0.005561127                       0.005350902 
#>    BB19-006-H_S22_MERGED.fastq.gz     BB6-019-B_S23_MERGED.fastq.gz 
#>                       0.005446338                       0.005439750 
#>     BB6-019-H_S24_MERGED.fastq.gz     BB6-019-M_S25_MERGED.fastq.gz 
#>                       0.005482589                       0.005373188 
#>      BD14-021_S26_MERGED.fastq.gz     BE9-006-B_S27_MERGED.fastq.gz 
#>                       0.005358611                       0.005383989 
#>     BE9-006-H_S28_MERGED.fastq.gz     BE9-006-M_S29_MERGED.fastq.gz 
#>                       0.005355821                       0.005395364 
#>     BG7-010-B_S30_MERGED.fastq.gz     BG7-010-H_S31_MERGED.fastq.gz 
#>                       0.004943752                       0.005220289 
#>     BG7-010-M_S32_MERGED.fastq.gz       BH9-021_S33_MERGED.fastq.gz 
#>                       0.004970198                       0.005257841 
#>    BJ17-007-M_S34_MERGED.fastq.gz   BJ8-ABM-003_S35_MERGED.fastq.gz 
#>                       0.005377047                       0.005577087 
#>     BL7-006-B_S36_MERGED.fastq.gz     BL7-006-H_S37_MERGED.fastq.gz 
#>                       0.004999042                       0.004963713 
#>     BL7-006-M_S38_MERGED.fastq.gz      BN11-041_S39_MERGED.fastq.gz 
#>                       0.004781180                       0.005429140 
#>       BN9-002_S40_MERGED.fastq.gz       BO8-002_S41_MERGED.fastq.gz 
#>                       0.005665252                       0.005517642 
#>       BO8-005_S42_MERGED.fastq.gz    BP11-001-B_S43_MERGED.fastq.gz 
#>                       0.005573845                       0.005494164 
#>    BP11-001-H_S44_MERGED.fastq.gz    BP11-001-M_S45_MERGED.fastq.gz 
#>                       0.005338535                       0.005495800 
#>    BP12-025-B_S46_MERGED.fastq.gz      BP14-006_S47_MERGED.fastq.gz 
#>                       0.005556655                       0.005326438 
#>       BQ3-019_S48_MERGED.fastq.gz     BQ4-018-B_S49_MERGED.fastq.gz 
#>                       0.005665118                       0.005356941 
#>     BQ4-018-H_S50_MERGED.fastq.gz     BQ4-018-M_S51_MERGED.fastq.gz 
#>                       0.005491212                       0.005571424 
#>    BQ9ABM-002_S52_MERGED.fastq.gz       BR8-005_S53_MERGED.fastq.gz 
#>                       0.005497845                       0.005491499 
#>      BS14-006_S54_MERGED.fastq.gz      BT-006-M_S55_MERGED.fastq.gz 
#>                       0.005384537                       0.005168301 
#>       BT7-006_S56_MERGED.fastq.gz    BV11-002-B_S57_MERGED.fastq.gz 
#>                       0.005084315                       0.005458463 
#>    BV11-002-H_S58_MERGED.fastq.gz    BV11-002-M_S59_MERGED.fastq.gz 
#>                       0.005438616                       0.005596859 
#>       BW8-003_S60_MERGED.fastq.gz        C1-001_S61_MERGED.fastq.gz 
#>                       0.005346791                       0.005276780 
#>     C21-NV1-B_S62_MERGED.fastq.gz     C21-NV1-H_S63_MERGED.fastq.gz 
#>                       0.005511602                       0.005477861 
#>     C21-NV1-M_S64_MERGED.fastq.gz        C9-005_S65_MERGED.fastq.gz 
#>                       0.005707347                       0.005590082 
#>      CA12-024_S66_MERGED.fastq.gz       CA9-027_S67_MERGED.fastq.gz 
#>                       0.005552383                       0.005706414 
#>         CA9-X_S68_MERGED.fastq.gz     CB8-019-B_S69_MERGED.fastq.gz 
#>                       0.005426847                       0.005571891 
#>     CB8-019-H_S70_MERGED.fastq.gz     CB8-019-M_S71_MERGED.fastq.gz 
#>                       0.005502441                       0.005530143 
#>       CB9-013_S72_MERGED.fastq.gz       CC3-044_S73_MERGED.fastq.gz 
#>                       0.005490287                       0.005457918 
#>       CC8-003_S74_MERGED.fastq.gz       D17-011_S77_MERGED.fastq.gz 
#>                       0.005601528                       0.005358224 
#>     D18-003-B_S78_MERGED.fastq.gz     D18-003-H_S79_MERGED.fastq.gz 
#>                       0.005457737                       0.005522691 
#>     D18-003-M_S80_MERGED.fastq.gz    D22-NVABM1_S81_MERGED.fastq.gz 
#>                       0.005577862                       0.005366414 
#>     D61-010-B_S82_MERGED.fastq.gz      D9-027-B_S83_MERGED.fastq.gz 
#>                       0.005491082                       0.005212816 
#>      D9-027-H_S84_MERGED.fastq.gz      D9-027-M_S85_MERGED.fastq.gz 
#>                       0.005028374                       0.005239751 
#>   DBM-ABM-001_S86_MERGED.fastq.gz     DJ2-008-B_S87_MERGED.fastq.gz 
#>                       0.005312083                       0.005399521 
#>     DJ2-008-H_S88_MERGED.fastq.gz     DJ2-008-M_S89_MERGED.fastq.gz 
#>                       0.005631249                       0.005280451 
#>    DP4-ABM001_S90_MERGED.fastq.gz  DS1-ABM002-B_S91_MERGED.fastq.gz 
#>                       0.005565384                       0.005478764 
#>  DS1-ABM002-H_S92_MERGED.fastq.gz  DS1-ABM002-M_S93_MERGED.fastq.gz 
#>                       0.005653426                       0.005492775 
#>     DU3-045-B_S94_MERGED.fastq.gz       DW4-007_S95_MERGED.fastq.gz 
#>                       0.005348383                       0.005466873 
#>     DY5-004-B_S96_MERGED.fastq.gz     DY5-004-H_S97_MERGED.fastq.gz 
#>                       0.005152746                       0.005059661 
#>     DY5-004-M_S98_MERGED.fastq.gz   DZ6-ABM-001_S99_MERGED.fastq.gz 
#>                       0.005173344                       0.005558530 
#>     E9-009-B_S100_MERGED.fastq.gz     E9-009-H_S101_MERGED.fastq.gz 
#>                       0.005064167                       0.005023463 
#>     E9-009-M_S102_MERGED.fastq.gz  EA5-ABM-001_S103_MERGED.fastq.gz 
#>                       0.004953764                       0.005481721 
#>    EC2-013-B_S104_MERGED.fastq.gz   F6-ABM-001_S105_MERGED.fastq.gz 
#>                       0.005547660                       0.005730904 
#>     F7-015-M_S106_MERGED.fastq.gz    FOMES19-H_S108_MERGED.fastq.gz 
#>                       0.005720220                       0.005089936 
#>    FOMES19-M_S109_MERGED.fastq.gz    H10-018-M_S110_MERGED.fastq.gz 
#>                       0.005206987                       0.005538178 
#> H24-NVABM1-H_S111_MERGED.fastq.gz    J18-004-B_S114_MERGED.fastq.gz 
#>                       0.005239803                       0.005151929 
#>    J18-004-H_S115_MERGED.fastq.gz    J18-004-M_S116_MERGED.fastq.gz 
#>                       0.005762547                       0.005134888 
#>    K18-002-H_S117_MERGED.fastq.gz   K26-NVABM1_S118_MERGED.fastq.gz 
#>                       0.005491317                       0.005502090 
#>       L19X-B_S119_MERGED.fastq.gz       L19X-H_S120_MERGED.fastq.gz 
#>                       0.005739303                       0.005516046 
#>       L19X-M_S121_MERGED.fastq.gz    L23-002-B_S122_MERGED.fastq.gz 
#>                       0.005365199                       0.005143392 
#>    L23-002-H_S123_MERGED.fastq.gz    L23-002-M_S124_MERGED.fastq.gz 
#>                       0.005572630                       0.005257109 
#>      M22-001_S125_MERGED.fastq.gz       N19X-B_S126_MERGED.fastq.gz 
#>                       0.005584895                       0.005041025 
#>       N19X-H_S127_MERGED.fastq.gz       N19X-M_S128_MERGED.fastq.gz 
#>                       0.005264857                       0.005377916 
#>    N22-001-B_S129_MERGED.fastq.gz    N23-002-B_S130_MERGED.fastq.gz 
#>                       0.005760180                       0.005470941 
#>    N23-002-H_S131_MERGED.fastq.gz    N23-002-M_S132_MERGED.fastq.gz 
#>                       0.005471282                       0.006089028 
#>     N25-ABMX_S133_MERGED.fastq.gz   NVABM-0058_S134_MERGED.fastq.gz 
#>                       0.005535886                       0.005464185 
#> NVABM-0163-H_S135_MERGED.fastq.gz   NVABM-0397_S138_MERGED.fastq.gz 
#>                       0.005517807                       0.005507368 
#>    NVABM0216_S136_MERGED.fastq.gz  NVABM0244-M_S137_MERGED.fastq.gz 
#>                       0.005634269                       0.005716539 
#>      O20-X-B_S139_MERGED.fastq.gz      O20-X-H_S140_MERGED.fastq.gz 
#>                       0.005550458                       0.005216059 
#>      O20-X-M_S141_MERGED.fastq.gz    O21-007-B_S142_MERGED.fastq.gz 
#>                       0.005206341                       0.005466052 
#>    O21-007-H_S143_MERGED.fastq.gz    O21-007-M_S144_MERGED.fastq.gz 
#>                       0.005573792                       0.005933239 
#>    O24-003-B_S145_MERGED.fastq.gz    O24-003-H_S146_MERGED.fastq.gz 
#>                       0.005314380                       0.005468160 
#>    O24-003-M_S147_MERGED.fastq.gz    O26-004-B_S148_MERGED.fastq.gz 
#>                       0.005296319                       0.005563424 
#>    O26-004-H_S149_MERGED.fastq.gz    O26-004-M_S150_MERGED.fastq.gz 
#>                       0.005538888                       0.005737268 
#>      O27-012_S151_MERGED.fastq.gz     O9-005-B_S152_MERGED.fastq.gz 
#>                       0.005585947                       0.005091717 
#>    P19-023-M_S153_MERGED.fastq.gz    P27-015-M_S154_MERGED.fastq.gz 
#>                       0.005550866                       0.005619124 
#>   P27-ABM001_S155_MERGED.fastq.gz Q27-ABM003-B_S156_MERGED.fastq.gz 
#>                       0.005636609                       0.005457268 
#>     R25-ABMX_S157_MERGED.fastq.gz    R28-008-B_S158_MERGED.fastq.gz 
#>                       0.005415096                       0.005101819 
#>    R28-008-H_S159_MERGED.fastq.gz    R28-008-M_S160_MERGED.fastq.gz 
#>                       0.004987584                       0.005077265 
#>      T28-011_S161_MERGED.fastq.gz T28-ABM602-B_S162_MERGED.fastq.gz 
#>                       0.005478888                       0.005776477 
#>   U27-ABM002_S163_MERGED.fastq.gz     W25-ABMX_S164_MERGED.fastq.gz 
#>                       0.005370920                       0.005453257 
#>    W26-001-B_S165_MERGED.fastq.gz    W26-001-H_S166_MERGED.fastq.gz 
#>                       0.004758101                       0.004908588 
#>    W26-001-M_S167_MERGED.fastq.gz      W30-006_S168_MERGED.fastq.gz 
#>                       0.005065333                       0.005394839 
#>     W9-025-M_S169_MERGED.fastq.gz    X24-009-B_S170_MERGED.fastq.gz 
#>                       0.005728282                       0.005267850 
#>    X24-009-H_S171_MERGED.fastq.gz    X24-009-M_S172_MERGED.fastq.gz 
#>                       0.005586239                       0.005548480 
#>      X24-010_S173_MERGED.fastq.gz    X29-004-B_S174_MERGED.fastq.gz 
#>                       0.005564110                       0.005657518 
#>    X29-004-H_S175_MERGED.fastq.gz    X29-004-M_S176_MERGED.fastq.gz 
#>                       0.005273738                       0.005466336 
#> Y21-ABM484-H_S177_MERGED.fastq.gz    Y28-002-B_S178_MERGED.fastq.gz 
#>                       0.005373932                       0.005279263 
#>    Y28-002-H_S179_MERGED.fastq.gz    Y28-002-M_S180_MERGED.fastq.gz 
#>                       0.005125520                       0.005339353 
#>    Y29-007-B_S181_MERGED.fastq.gz    Y29-007-H_S182_MERGED.fastq.gz 
#>                       0.005127480                       0.005094909 
#>    Y29-007-M_S183_MERGED.fastq.gz Y31-ABM484-B_S184_MERGED.fastq.gz 
#>                       0.005505642                       0.005582052 
#>    Z29-001-H_S185_MERGED.fastq.gz      Z30-002_S186_MERGED.fastq.gz 
#>                       0.005431650                       0.005088953 
#> Z30-ABM560-M_S187_MERGED.fastq.gz 
#>                       0.005368628 
#> 
#> $p.LCBD
#>    A10-005-B_S188_MERGED.fastq.gz    A10-005-H_S189_MERGED.fastq.gz 
#>                         0.8333333                         0.1666667 
#>    A10-005-M_S190_MERGED.fastq.gz      A12-007_S191_MERGED.fastq.gz 
#>                         0.3333333                         0.5000000 
#>      A12-007-B_S2_MERGED.fastq.gz        A15-004_S3_MERGED.fastq.gz 
#>                         0.3333333                         1.0000000 
#>         A8-005_S4_MERGED.fastq.gz         A9-012_S5_MERGED.fastq.gz 
#>                         1.0000000                         1.0000000 
#>    AB29-ABMX-H_S6_MERGED.fastq.gz       AC27-013_S7_MERGED.fastq.gz 
#>                         0.1666667                         0.6666667 
#>        AC29033_S8_MERGED.fastq.gz     AD26-005-B_S9_MERGED.fastq.gz 
#>                         0.1666667                         0.1666667 
#>    AD26-005-H_S10_MERGED.fastq.gz    AD26-005-M_S11_MERGED.fastq.gz 
#>                         0.1666667                         0.1666667 
#>   AD30-ABMX-M_S12_MERGED.fastq.gz    AD32-007-M_S13_MERGED.fastq.gz 
#>                         0.1666667                         0.8333333 
#>    ADABM30X-B_S14_MERGED.fastq.gz    ADABM30X-H_S15_MERGED.fastq.gz 
#>                         0.1666667                         0.6666667 
#>    ADABM30X-M_S16_MERGED.fastq.gz   AE30-ABM507_S17_MERGED.fastq.gz 
#>                         0.1666667                         0.8333333 
#>       B17-014_S18_MERGED.fastq.gz     B18-006-B_S19_MERGED.fastq.gz 
#>                         1.0000000                         1.0000000 
#>   BA16-036bis_S20_MERGED.fastq.gz    BA17-050-B_S21_MERGED.fastq.gz 
#>                         0.1666667                         0.5000000 
#>    BB19-006-H_S22_MERGED.fastq.gz     BB6-019-B_S23_MERGED.fastq.gz 
#>                         0.3333333                         0.5000000 
#>     BB6-019-H_S24_MERGED.fastq.gz     BB6-019-M_S25_MERGED.fastq.gz 
#>                         0.5000000                         0.5000000 
#>      BD14-021_S26_MERGED.fastq.gz     BE9-006-B_S27_MERGED.fastq.gz 
#>                         0.8333333                         0.6666667 
#>     BE9-006-H_S28_MERGED.fastq.gz     BE9-006-M_S29_MERGED.fastq.gz 
#>                         0.5000000                         0.6666667 
#>     BG7-010-B_S30_MERGED.fastq.gz     BG7-010-H_S31_MERGED.fastq.gz 
#>                         1.0000000                         0.8333333 
#>     BG7-010-M_S32_MERGED.fastq.gz       BH9-021_S33_MERGED.fastq.gz 
#>                         1.0000000                         1.0000000 
#>    BJ17-007-M_S34_MERGED.fastq.gz   BJ8-ABM-003_S35_MERGED.fastq.gz 
#>                         0.3333333                         0.1666667 
#>     BL7-006-B_S36_MERGED.fastq.gz     BL7-006-H_S37_MERGED.fastq.gz 
#>                         1.0000000                         1.0000000 
#>     BL7-006-M_S38_MERGED.fastq.gz      BN11-041_S39_MERGED.fastq.gz 
#>                         1.0000000                         0.6666667 
#>       BN9-002_S40_MERGED.fastq.gz       BO8-002_S41_MERGED.fastq.gz 
#>                         0.1666667                         0.1666667 
#>       BO8-005_S42_MERGED.fastq.gz    BP11-001-B_S43_MERGED.fastq.gz 
#>                         0.1666667                         0.1666667 
#>    BP11-001-H_S44_MERGED.fastq.gz    BP11-001-M_S45_MERGED.fastq.gz 
#>                         1.0000000                         0.1666667 
#>    BP12-025-B_S46_MERGED.fastq.gz      BP14-006_S47_MERGED.fastq.gz 
#>                         0.1666667                         0.8333333 
#>       BQ3-019_S48_MERGED.fastq.gz     BQ4-018-B_S49_MERGED.fastq.gz 
#>                         0.1666667                         0.8333333 
#>     BQ4-018-H_S50_MERGED.fastq.gz     BQ4-018-M_S51_MERGED.fastq.gz 
#>                         0.3333333                         0.1666667 
#>    BQ9ABM-002_S52_MERGED.fastq.gz       BR8-005_S53_MERGED.fastq.gz 
#>                         0.3333333                         0.3333333 
#>      BS14-006_S54_MERGED.fastq.gz      BT-006-M_S55_MERGED.fastq.gz 
#>                         0.6666667                         1.0000000 
#>       BT7-006_S56_MERGED.fastq.gz    BV11-002-B_S57_MERGED.fastq.gz 
#>                         1.0000000                         0.3333333 
#>    BV11-002-H_S58_MERGED.fastq.gz    BV11-002-M_S59_MERGED.fastq.gz 
#>                         0.3333333                         0.1666667 
#>       BW8-003_S60_MERGED.fastq.gz        C1-001_S61_MERGED.fastq.gz 
#>                         0.8333333                         1.0000000 
#>     C21-NV1-B_S62_MERGED.fastq.gz     C21-NV1-H_S63_MERGED.fastq.gz 
#>                         0.3333333                         0.1666667 
#>     C21-NV1-M_S64_MERGED.fastq.gz        C9-005_S65_MERGED.fastq.gz 
#>                         0.1666667                         0.1666667 
#>      CA12-024_S66_MERGED.fastq.gz       CA9-027_S67_MERGED.fastq.gz 
#>                         0.1666667                         0.1666667 
#>         CA9-X_S68_MERGED.fastq.gz     CB8-019-B_S69_MERGED.fastq.gz 
#>                         0.6666667                         0.5000000 
#>     CB8-019-H_S70_MERGED.fastq.gz     CB8-019-M_S71_MERGED.fastq.gz 
#>                         0.3333333                         0.1666667 
#>       CB9-013_S72_MERGED.fastq.gz       CC3-044_S73_MERGED.fastq.gz 
#>                         0.3333333                         0.5000000 
#>       CC8-003_S74_MERGED.fastq.gz       D17-011_S77_MERGED.fastq.gz 
#>                         0.1666667                         0.8333333 
#>     D18-003-B_S78_MERGED.fastq.gz     D18-003-H_S79_MERGED.fastq.gz 
#>                         0.5000000                         0.3333333 
#>     D18-003-M_S80_MERGED.fastq.gz    D22-NVABM1_S81_MERGED.fastq.gz 
#>                         0.1666667                         0.8333333 
#>     D61-010-B_S82_MERGED.fastq.gz      D9-027-B_S83_MERGED.fastq.gz 
#>                         0.1666667                         1.0000000 
#>      D9-027-H_S84_MERGED.fastq.gz      D9-027-M_S85_MERGED.fastq.gz 
#>                         1.0000000                         1.0000000 
#>   DBM-ABM-001_S86_MERGED.fastq.gz     DJ2-008-B_S87_MERGED.fastq.gz 
#>                         0.6666667                         0.6666667 
#>     DJ2-008-H_S88_MERGED.fastq.gz     DJ2-008-M_S89_MERGED.fastq.gz 
#>                         0.1666667                         1.0000000 
#>    DP4-ABM001_S90_MERGED.fastq.gz  DS1-ABM002-B_S91_MERGED.fastq.gz 
#>                         0.3333333                         0.3333333 
#>  DS1-ABM002-H_S92_MERGED.fastq.gz  DS1-ABM002-M_S93_MERGED.fastq.gz 
#>                         0.1666667                         0.1666667 
#>     DU3-045-B_S94_MERGED.fastq.gz       DW4-007_S95_MERGED.fastq.gz 
#>                         0.8333333                         0.1666667 
#>     DY5-004-B_S96_MERGED.fastq.gz     DY5-004-H_S97_MERGED.fastq.gz 
#>                         1.0000000                         1.0000000 
#>     DY5-004-M_S98_MERGED.fastq.gz   DZ6-ABM-001_S99_MERGED.fastq.gz 
#>                         1.0000000                         0.3333333 
#>     E9-009-B_S100_MERGED.fastq.gz     E9-009-H_S101_MERGED.fastq.gz 
#>                         1.0000000                         1.0000000 
#>     E9-009-M_S102_MERGED.fastq.gz  EA5-ABM-001_S103_MERGED.fastq.gz 
#>                         1.0000000                         0.3333333 
#>    EC2-013-B_S104_MERGED.fastq.gz   F6-ABM-001_S105_MERGED.fastq.gz 
#>                         0.1666667                         0.1666667 
#>     F7-015-M_S106_MERGED.fastq.gz    FOMES19-H_S108_MERGED.fastq.gz 
#>                         0.1666667                         1.0000000 
#>    FOMES19-M_S109_MERGED.fastq.gz    H10-018-M_S110_MERGED.fastq.gz 
#>                         1.0000000                         0.1666667 
#> H24-NVABM1-H_S111_MERGED.fastq.gz    J18-004-B_S114_MERGED.fastq.gz 
#>                         0.8333333                         1.0000000 
#>    J18-004-H_S115_MERGED.fastq.gz    J18-004-M_S116_MERGED.fastq.gz 
#>                         0.1666667                         1.0000000 
#>    K18-002-H_S117_MERGED.fastq.gz   K26-NVABM1_S118_MERGED.fastq.gz 
#>                         0.3333333                         0.3333333 
#>       L19X-B_S119_MERGED.fastq.gz       L19X-H_S120_MERGED.fastq.gz 
#>                         0.1666667                         0.1666667 
#>       L19X-M_S121_MERGED.fastq.gz    L23-002-B_S122_MERGED.fastq.gz 
#>                         0.8333333                         1.0000000 
#>    L23-002-H_S123_MERGED.fastq.gz    L23-002-M_S124_MERGED.fastq.gz 
#>                         0.1666667                         0.8333333 
#>      M22-001_S125_MERGED.fastq.gz       N19X-B_S126_MERGED.fastq.gz 
#>                         0.3333333                         1.0000000 
#>       N19X-H_S127_MERGED.fastq.gz       N19X-M_S128_MERGED.fastq.gz 
#>                         1.0000000                         0.5000000 
#>    N22-001-B_S129_MERGED.fastq.gz    N23-002-B_S130_MERGED.fastq.gz 
#>                         0.1666667                         0.3333333 
#>    N23-002-H_S131_MERGED.fastq.gz    N23-002-M_S132_MERGED.fastq.gz 
#>                         0.5000000                         0.1666667 
#>     N25-ABMX_S133_MERGED.fastq.gz   NVABM-0058_S134_MERGED.fastq.gz 
#>                         0.3333333                         0.1666667 
#> NVABM-0163-H_S135_MERGED.fastq.gz   NVABM-0397_S138_MERGED.fastq.gz 
#>                         0.5000000                         0.6666667 
#>    NVABM0216_S136_MERGED.fastq.gz  NVABM0244-M_S137_MERGED.fastq.gz 
#>                         0.1666667                         0.1666667 
#>      O20-X-B_S139_MERGED.fastq.gz      O20-X-H_S140_MERGED.fastq.gz 
#>                         0.1666667                         1.0000000 
#>      O20-X-M_S141_MERGED.fastq.gz    O21-007-B_S142_MERGED.fastq.gz 
#>                         0.8333333                         0.6666667 
#>    O21-007-H_S143_MERGED.fastq.gz    O21-007-M_S144_MERGED.fastq.gz 
#>                         0.3333333                         0.1666667 
#>    O24-003-B_S145_MERGED.fastq.gz    O24-003-H_S146_MERGED.fastq.gz 
#>                         0.5000000                         0.5000000 
#>    O24-003-M_S147_MERGED.fastq.gz    O26-004-B_S148_MERGED.fastq.gz 
#>                         1.0000000                         0.1666667 
#>    O26-004-H_S149_MERGED.fastq.gz    O26-004-M_S150_MERGED.fastq.gz 
#>                         0.1666667                         0.1666667 
#>      O27-012_S151_MERGED.fastq.gz     O9-005-B_S152_MERGED.fastq.gz 
#>                         0.3333333                         1.0000000 
#>    P19-023-M_S153_MERGED.fastq.gz    P27-015-M_S154_MERGED.fastq.gz 
#>                         0.5000000                         0.3333333 
#>   P27-ABM001_S155_MERGED.fastq.gz Q27-ABM003-B_S156_MERGED.fastq.gz 
#>                         0.1666667                         0.1666667 
#>     R25-ABMX_S157_MERGED.fastq.gz    R28-008-B_S158_MERGED.fastq.gz 
#>                         0.3333333                         1.0000000 
#>    R28-008-H_S159_MERGED.fastq.gz    R28-008-M_S160_MERGED.fastq.gz 
#>                         1.0000000                         1.0000000 
#>      T28-011_S161_MERGED.fastq.gz T28-ABM602-B_S162_MERGED.fastq.gz 
#>                         0.3333333                         0.1666667 
#>   U27-ABM002_S163_MERGED.fastq.gz     W25-ABMX_S164_MERGED.fastq.gz 
#>                         1.0000000                         0.5000000 
#>    W26-001-B_S165_MERGED.fastq.gz    W26-001-H_S166_MERGED.fastq.gz 
#>                         1.0000000                         1.0000000 
#>    W26-001-M_S167_MERGED.fastq.gz      W30-006_S168_MERGED.fastq.gz 
#>                         1.0000000                         0.3333333 
#>     W9-025-M_S169_MERGED.fastq.gz    X24-009-B_S170_MERGED.fastq.gz 
#>                         0.1666667                         1.0000000 
#>    X24-009-H_S171_MERGED.fastq.gz    X24-009-M_S172_MERGED.fastq.gz 
#>                         0.1666667                         0.1666667 
#>      X24-010_S173_MERGED.fastq.gz    X29-004-B_S174_MERGED.fastq.gz 
#>                         0.1666667                         0.3333333 
#>    X29-004-H_S175_MERGED.fastq.gz    X29-004-M_S176_MERGED.fastq.gz 
#>                         1.0000000                         0.1666667 
#> Y21-ABM484-H_S177_MERGED.fastq.gz    Y28-002-B_S178_MERGED.fastq.gz 
#>                         0.6666667                         0.8333333 
#>    Y28-002-H_S179_MERGED.fastq.gz    Y28-002-M_S180_MERGED.fastq.gz 
#>                         1.0000000                         0.8333333 
#>    Y29-007-B_S181_MERGED.fastq.gz    Y29-007-H_S182_MERGED.fastq.gz 
#>                         1.0000000                         1.0000000 
#>    Y29-007-M_S183_MERGED.fastq.gz Y31-ABM484-B_S184_MERGED.fastq.gz 
#>                         0.1666667                         0.3333333 
#>    Z29-001-H_S185_MERGED.fastq.gz      Z30-002_S186_MERGED.fastq.gz 
#>                         0.6666667                         1.0000000 
#> Z30-ABM560-M_S187_MERGED.fastq.gz 
#>                         0.6666667 
#> 
#> $p.adj
#>    A10-005-B_S188_MERGED.fastq.gz    A10-005-H_S189_MERGED.fastq.gz 
#>                         1.0000000                         0.5138889 
#>    A10-005-M_S190_MERGED.fastq.gz      A12-007_S191_MERGED.fastq.gz 
#>                         0.6851852                         0.8726415 
#>      A12-007-B_S2_MERGED.fastq.gz        A15-004_S3_MERGED.fastq.gz 
#>                         0.6851852                         1.0000000 
#>         A8-005_S4_MERGED.fastq.gz         A9-012_S5_MERGED.fastq.gz 
#>                         1.0000000                         1.0000000 
#>    AB29-ABMX-H_S6_MERGED.fastq.gz       AC27-013_S7_MERGED.fastq.gz 
#>                         0.5138889                         1.0000000 
#>        AC29033_S8_MERGED.fastq.gz     AD26-005-B_S9_MERGED.fastq.gz 
#>                         0.5138889                         0.5138889 
#>    AD26-005-H_S10_MERGED.fastq.gz    AD26-005-M_S11_MERGED.fastq.gz 
#>                         0.5138889                         0.5138889 
#>   AD30-ABMX-M_S12_MERGED.fastq.gz    AD32-007-M_S13_MERGED.fastq.gz 
#>                         0.5138889                         1.0000000 
#>    ADABM30X-B_S14_MERGED.fastq.gz    ADABM30X-H_S15_MERGED.fastq.gz 
#>                         0.5138889                         1.0000000 
#>    ADABM30X-M_S16_MERGED.fastq.gz   AE30-ABM507_S17_MERGED.fastq.gz 
#>                         0.5138889                         1.0000000 
#>       B17-014_S18_MERGED.fastq.gz     B18-006-B_S19_MERGED.fastq.gz 
#>                         1.0000000                         1.0000000 
#>   BA16-036bis_S20_MERGED.fastq.gz    BA17-050-B_S21_MERGED.fastq.gz 
#>                         0.5138889                         0.8726415 
#>    BB19-006-H_S22_MERGED.fastq.gz     BB6-019-B_S23_MERGED.fastq.gz 
#>                         0.6851852                         0.8726415 
#>     BB6-019-H_S24_MERGED.fastq.gz     BB6-019-M_S25_MERGED.fastq.gz 
#>                         0.8726415                         0.8726415 
#>      BD14-021_S26_MERGED.fastq.gz     BE9-006-B_S27_MERGED.fastq.gz 
#>                         1.0000000                         1.0000000 
#>     BE9-006-H_S28_MERGED.fastq.gz     BE9-006-M_S29_MERGED.fastq.gz 
#>                         0.8726415                         1.0000000 
#>     BG7-010-B_S30_MERGED.fastq.gz     BG7-010-H_S31_MERGED.fastq.gz 
#>                         1.0000000                         1.0000000 
#>     BG7-010-M_S32_MERGED.fastq.gz       BH9-021_S33_MERGED.fastq.gz 
#>                         1.0000000                         1.0000000 
#>    BJ17-007-M_S34_MERGED.fastq.gz   BJ8-ABM-003_S35_MERGED.fastq.gz 
#>                         0.6851852                         0.5138889 
#>     BL7-006-B_S36_MERGED.fastq.gz     BL7-006-H_S37_MERGED.fastq.gz 
#>                         1.0000000                         1.0000000 
#>     BL7-006-M_S38_MERGED.fastq.gz      BN11-041_S39_MERGED.fastq.gz 
#>                         1.0000000                         1.0000000 
#>       BN9-002_S40_MERGED.fastq.gz       BO8-002_S41_MERGED.fastq.gz 
#>                         0.5138889                         0.5138889 
#>       BO8-005_S42_MERGED.fastq.gz    BP11-001-B_S43_MERGED.fastq.gz 
#>                         0.5138889                         0.5138889 
#>    BP11-001-H_S44_MERGED.fastq.gz    BP11-001-M_S45_MERGED.fastq.gz 
#>                         1.0000000                         0.5138889 
#>    BP12-025-B_S46_MERGED.fastq.gz      BP14-006_S47_MERGED.fastq.gz 
#>                         0.5138889                         1.0000000 
#>       BQ3-019_S48_MERGED.fastq.gz     BQ4-018-B_S49_MERGED.fastq.gz 
#>                         0.5138889                         1.0000000 
#>     BQ4-018-H_S50_MERGED.fastq.gz     BQ4-018-M_S51_MERGED.fastq.gz 
#>                         0.6851852                         0.5138889 
#>    BQ9ABM-002_S52_MERGED.fastq.gz       BR8-005_S53_MERGED.fastq.gz 
#>                         0.6851852                         0.6851852 
#>      BS14-006_S54_MERGED.fastq.gz      BT-006-M_S55_MERGED.fastq.gz 
#>                         1.0000000                         1.0000000 
#>       BT7-006_S56_MERGED.fastq.gz    BV11-002-B_S57_MERGED.fastq.gz 
#>                         1.0000000                         0.6851852 
#>    BV11-002-H_S58_MERGED.fastq.gz    BV11-002-M_S59_MERGED.fastq.gz 
#>                         0.6851852                         0.5138889 
#>       BW8-003_S60_MERGED.fastq.gz        C1-001_S61_MERGED.fastq.gz 
#>                         1.0000000                         1.0000000 
#>     C21-NV1-B_S62_MERGED.fastq.gz     C21-NV1-H_S63_MERGED.fastq.gz 
#>                         0.6851852                         0.5138889 
#>     C21-NV1-M_S64_MERGED.fastq.gz        C9-005_S65_MERGED.fastq.gz 
#>                         0.5138889                         0.5138889 
#>      CA12-024_S66_MERGED.fastq.gz       CA9-027_S67_MERGED.fastq.gz 
#>                         0.5138889                         0.5138889 
#>         CA9-X_S68_MERGED.fastq.gz     CB8-019-B_S69_MERGED.fastq.gz 
#>                         1.0000000                         0.8726415 
#>     CB8-019-H_S70_MERGED.fastq.gz     CB8-019-M_S71_MERGED.fastq.gz 
#>                         0.6851852                         0.5138889 
#>       CB9-013_S72_MERGED.fastq.gz       CC3-044_S73_MERGED.fastq.gz 
#>                         0.6851852                         0.8726415 
#>       CC8-003_S74_MERGED.fastq.gz       D17-011_S77_MERGED.fastq.gz 
#>                         0.5138889                         1.0000000 
#>     D18-003-B_S78_MERGED.fastq.gz     D18-003-H_S79_MERGED.fastq.gz 
#>                         0.8726415                         0.6851852 
#>     D18-003-M_S80_MERGED.fastq.gz    D22-NVABM1_S81_MERGED.fastq.gz 
#>                         0.5138889                         1.0000000 
#>     D61-010-B_S82_MERGED.fastq.gz      D9-027-B_S83_MERGED.fastq.gz 
#>                         0.5138889                         1.0000000 
#>      D9-027-H_S84_MERGED.fastq.gz      D9-027-M_S85_MERGED.fastq.gz 
#>                         1.0000000                         1.0000000 
#>   DBM-ABM-001_S86_MERGED.fastq.gz     DJ2-008-B_S87_MERGED.fastq.gz 
#>                         1.0000000                         1.0000000 
#>     DJ2-008-H_S88_MERGED.fastq.gz     DJ2-008-M_S89_MERGED.fastq.gz 
#>                         0.5138889                         1.0000000 
#>    DP4-ABM001_S90_MERGED.fastq.gz  DS1-ABM002-B_S91_MERGED.fastq.gz 
#>                         0.6851852                         0.6851852 
#>  DS1-ABM002-H_S92_MERGED.fastq.gz  DS1-ABM002-M_S93_MERGED.fastq.gz 
#>                         0.5138889                         0.5138889 
#>     DU3-045-B_S94_MERGED.fastq.gz       DW4-007_S95_MERGED.fastq.gz 
#>                         1.0000000                         0.5138889 
#>     DY5-004-B_S96_MERGED.fastq.gz     DY5-004-H_S97_MERGED.fastq.gz 
#>                         1.0000000                         1.0000000 
#>     DY5-004-M_S98_MERGED.fastq.gz   DZ6-ABM-001_S99_MERGED.fastq.gz 
#>                         1.0000000                         0.6851852 
#>     E9-009-B_S100_MERGED.fastq.gz     E9-009-H_S101_MERGED.fastq.gz 
#>                         1.0000000                         1.0000000 
#>     E9-009-M_S102_MERGED.fastq.gz  EA5-ABM-001_S103_MERGED.fastq.gz 
#>                         1.0000000                         0.6851852 
#>    EC2-013-B_S104_MERGED.fastq.gz   F6-ABM-001_S105_MERGED.fastq.gz 
#>                         0.5138889                         0.5138889 
#>     F7-015-M_S106_MERGED.fastq.gz    FOMES19-H_S108_MERGED.fastq.gz 
#>                         0.5138889                         1.0000000 
#>    FOMES19-M_S109_MERGED.fastq.gz    H10-018-M_S110_MERGED.fastq.gz 
#>                         1.0000000                         0.5138889 
#> H24-NVABM1-H_S111_MERGED.fastq.gz    J18-004-B_S114_MERGED.fastq.gz 
#>                         1.0000000                         1.0000000 
#>    J18-004-H_S115_MERGED.fastq.gz    J18-004-M_S116_MERGED.fastq.gz 
#>                         0.5138889                         1.0000000 
#>    K18-002-H_S117_MERGED.fastq.gz   K26-NVABM1_S118_MERGED.fastq.gz 
#>                         0.6851852                         0.6851852 
#>       L19X-B_S119_MERGED.fastq.gz       L19X-H_S120_MERGED.fastq.gz 
#>                         0.5138889                         0.5138889 
#>       L19X-M_S121_MERGED.fastq.gz    L23-002-B_S122_MERGED.fastq.gz 
#>                         1.0000000                         1.0000000 
#>    L23-002-H_S123_MERGED.fastq.gz    L23-002-M_S124_MERGED.fastq.gz 
#>                         0.5138889                         1.0000000 
#>      M22-001_S125_MERGED.fastq.gz       N19X-B_S126_MERGED.fastq.gz 
#>                         0.6851852                         1.0000000 
#>       N19X-H_S127_MERGED.fastq.gz       N19X-M_S128_MERGED.fastq.gz 
#>                         1.0000000                         0.8726415 
#>    N22-001-B_S129_MERGED.fastq.gz    N23-002-B_S130_MERGED.fastq.gz 
#>                         0.5138889                         0.6851852 
#>    N23-002-H_S131_MERGED.fastq.gz    N23-002-M_S132_MERGED.fastq.gz 
#>                         0.8726415                         0.5138889 
#>     N25-ABMX_S133_MERGED.fastq.gz   NVABM-0058_S134_MERGED.fastq.gz 
#>                         0.6851852                         0.5138889 
#> NVABM-0163-H_S135_MERGED.fastq.gz   NVABM-0397_S138_MERGED.fastq.gz 
#>                         0.8726415                         1.0000000 
#>    NVABM0216_S136_MERGED.fastq.gz  NVABM0244-M_S137_MERGED.fastq.gz 
#>                         0.5138889                         0.5138889 
#>      O20-X-B_S139_MERGED.fastq.gz      O20-X-H_S140_MERGED.fastq.gz 
#>                         0.5138889                         1.0000000 
#>      O20-X-M_S141_MERGED.fastq.gz    O21-007-B_S142_MERGED.fastq.gz 
#>                         1.0000000                         1.0000000 
#>    O21-007-H_S143_MERGED.fastq.gz    O21-007-M_S144_MERGED.fastq.gz 
#>                         0.6851852                         0.5138889 
#>    O24-003-B_S145_MERGED.fastq.gz    O24-003-H_S146_MERGED.fastq.gz 
#>                         0.8726415                         0.8726415 
#>    O24-003-M_S147_MERGED.fastq.gz    O26-004-B_S148_MERGED.fastq.gz 
#>                         1.0000000                         0.5138889 
#>    O26-004-H_S149_MERGED.fastq.gz    O26-004-M_S150_MERGED.fastq.gz 
#>                         0.5138889                         0.5138889 
#>      O27-012_S151_MERGED.fastq.gz     O9-005-B_S152_MERGED.fastq.gz 
#>                         0.6851852                         1.0000000 
#>    P19-023-M_S153_MERGED.fastq.gz    P27-015-M_S154_MERGED.fastq.gz 
#>                         0.8726415                         0.6851852 
#>   P27-ABM001_S155_MERGED.fastq.gz Q27-ABM003-B_S156_MERGED.fastq.gz 
#>                         0.5138889                         0.5138889 
#>     R25-ABMX_S157_MERGED.fastq.gz    R28-008-B_S158_MERGED.fastq.gz 
#>                         0.6851852                         1.0000000 
#>    R28-008-H_S159_MERGED.fastq.gz    R28-008-M_S160_MERGED.fastq.gz 
#>                         1.0000000                         1.0000000 
#>      T28-011_S161_MERGED.fastq.gz T28-ABM602-B_S162_MERGED.fastq.gz 
#>                         0.6851852                         0.5138889 
#>   U27-ABM002_S163_MERGED.fastq.gz     W25-ABMX_S164_MERGED.fastq.gz 
#>                         1.0000000                         0.8726415 
#>    W26-001-B_S165_MERGED.fastq.gz    W26-001-H_S166_MERGED.fastq.gz 
#>                         1.0000000                         1.0000000 
#>    W26-001-M_S167_MERGED.fastq.gz      W30-006_S168_MERGED.fastq.gz 
#>                         1.0000000                         0.6851852 
#>     W9-025-M_S169_MERGED.fastq.gz    X24-009-B_S170_MERGED.fastq.gz 
#>                         0.5138889                         1.0000000 
#>    X24-009-H_S171_MERGED.fastq.gz    X24-009-M_S172_MERGED.fastq.gz 
#>                         0.5138889                         0.5138889 
#>      X24-010_S173_MERGED.fastq.gz    X29-004-B_S174_MERGED.fastq.gz 
#>                         0.5138889                         0.6851852 
#>    X29-004-H_S175_MERGED.fastq.gz    X29-004-M_S176_MERGED.fastq.gz 
#>                         1.0000000                         0.5138889 
#> Y21-ABM484-H_S177_MERGED.fastq.gz    Y28-002-B_S178_MERGED.fastq.gz 
#>                         1.0000000                         1.0000000 
#>    Y28-002-H_S179_MERGED.fastq.gz    Y28-002-M_S180_MERGED.fastq.gz 
#>                         1.0000000                         1.0000000 
#>    Y29-007-B_S181_MERGED.fastq.gz    Y29-007-H_S182_MERGED.fastq.gz 
#>                         1.0000000                         1.0000000 
#>    Y29-007-M_S183_MERGED.fastq.gz Y31-ABM484-B_S184_MERGED.fastq.gz 
#>                         0.5138889                         0.6851852 
#>    Z29-001-H_S185_MERGED.fastq.gz      Z30-002_S186_MERGED.fastq.gz 
#>                         1.0000000                         1.0000000 
#> Z30-ABM560-M_S187_MERGED.fastq.gz 
#>                         1.0000000 
#> 
#> $method
#> [1] "jaccard"      "sqrt.D=FALSE"
#> 
#> $note
#> [1] "Info -- D is Euclidean because beta.div outputs D[jk] = sqrt(1-S[jk])"
#> [2] "For this D functions, use beta.div with option sqrt.D=FALSE"          
#> 
#> $D
#> [1] NA
#> 
#> attr(,"class")
#> [1] "beta.div"
# }
```
