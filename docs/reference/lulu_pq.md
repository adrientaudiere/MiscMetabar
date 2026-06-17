# Lulu reclustering of class `physeq`

[![lifecycle-experimental](https://img.shields.io/badge/lifecycle-experimental-orange)](https://adrientaudiere.github.io/MiscMetabar/articles/Rules.html#lifecycle)

See https://www.nature.com/articles/s41467-017-01312-x for more
information on the method.

## Usage

``` r
lulu_pq(
  physeq,
  nproc = 1,
  id = 0.84,
  vsearchpath = find_vsearch(),
  verbose = FALSE,
  clean_pq = FALSE,
  keep_temporary_files = FALSE,
  ...
)
```

## Arguments

- physeq:

  (required) a
  [`phyloseq-class`](https://rdrr.io/pkg/phyloseq/man/phyloseq-class.html)
  object obtained using the `phyloseq` package.

- nproc:

  (default 1) Set to number of cpus/processors to use for the clustering

- id:

  (default: 0.84) id for –usearch_global.

- vsearchpath:

  (default: vsearch) path to vsearch.

- verbose:

  (logical) If true, print some additional messages.

- clean_pq:

  (logical) If true, empty samples and empty ASV are discarded before
  clustering.

- keep_temporary_files:

  (logical, default: FALSE) Do we keep temporary files

- ...:

  Additional arguments passed on to function
  [`lulu()`](https://adrientaudiere.github.io/MiscMetabar/reference/lulu.md)

## Value

a list of for object

- "new_physeq": The new phyloseq object (class physeq)

- "discrepancy_vector": A vector of discrepancy showing for each
  taxonomic level the proportion of identic value before and after lulu
  reclustering. A value of 0.6 stands for 60% of ASV before
  re-clustering have identical value after re-clustering. In other word,
  40% of ASV are assigned to a different taxonomic value. NA value are
  not counted as discrepancy.

- "res_lulu": A list of the result from the lulu function

- "merged_ASV": the data.frame used to merged ASV

## Details

The version of LULU is a fork of Adrien Taudière
(<https://github.com/adrientaudiere/lulu>) from
<https://github.com/tobiasgf/lulu>

## References

- LULU : <https://github.com/adrientaudiere/lulu> forked from
  <https://github.com/tobiasgf/lulu>.

- VSEARCH can be downloaded from <https://github.com/torognes/vsearch>.

## See also

[`mumu_pq()`](https://adrientaudiere.github.io/MiscMetabar/reference/mumu_pq.md)

## Author

Tobias Guldberg Frøslev <tobiasgf@snm.ku.dk> & Adrien Taudière
<adrien.taudiere@zaclys.net>

## Examples

``` r
# \donttest{
data_f <- clean_pq(prune_samples(
  sample_names(data_fungi_sp_known)[1:20],
  data_fungi_sp_known
))
#> Cleaning suppress 240 taxa and 0 samples.
lulu_pq(data_f)
#> Start Vsearch usearch_global
#> Lulu algorithm
#>   |                                                          |                                                  |   0%
#> 
#> ####processing: ASV2 #####4
#> 
#> ---hits: ASV1027
#> ---hits: ASV1082
#> ---hits: ASV962
#> ---hits: ASV1247
#> ---hits: ASV11024
#> 
#> ---potential parent: 4
#> 
#> No parent found!
#> 4
#> 
#> ####processing: ASV8 #####4
#> 
#> ---hits: ASV18
#> ---hits: ASV666
#> ---hits: ASV170
#> ---hits: ASV493
#> ---hits: ASV831
#> ---hits: ASV261
#> ---hits: ASV474
#> ---hits: ASV209
#> ---hits: ASV94
#> ---hits: ASV3484
#> 
#> ---potential parent: 4
#> 
#> No parent found!
#> 4
#>   |                                                          |                                                  |   1%
#> 
#> ####processing: ASV38 #####4
#> 
#> ---hits: ASV466
#> ---hits: ASV989
#> ---hits: ASV853
#> ---hits: ASV816
#> ---hits: ASV975
#> ---hits: ASV1526
#> ---hits: ASV254
#> ---hits: ASV500
#> ---hits: ASV592
#> ---hits: ASV11074
#> 
#> ---potential parent: 4
#> 
#> No parent found!
#> 4
#> 
#> ####processing: ASV756 #####4
#> 
#> ---hits: ASV904
#> ---hits: ASV1035
#> ---hits: ASV1074
#> ---hits: ASV906
#> ---hits: ASV1278
#> ---hits: ASV1007
#> ---hits: ASV1192
#> ---hits: ASV1359
#> ---hits: ASV828
#> ---hits: ASV12244
#> 
#> ---potential parent: 4
#> 
#> No parent found!
#> 4
#>   |                                                          |=                                                 |   1%
#> 
#> ####processing: ASV12 #####4
#> 
#> ---hits: ASV67
#> ---hits: ASV107
#> ---hits: ASV1542
#> ---hits: ASV624
#> ---hits: ASV1262
#> ---hits: ASV847
#> ---hits: ASV1054
#> 
#> ---potential parent: 4
#> 
#> No parent found!
#> 4
#> 
#> ####processing: ASV19 #####4
#> 
#> ---hits: ASV333
#> ---hits: ASV178
#> ---hits: ASV594
#> ---hits: ASV1131
#> ---hits: ASV1010
#> ---hits: ASV858
#> ---hits: ASV1467
#> ---hits: ASV1078
#> ---hits: ASV986
#> ---hits: ASV9794
#> 
#> ---potential parent: 4
#> 
#> No parent found!
#> 4
#>   |                                                          |=                                                 |   2%
#> 
#> ####processing: ASV175 #####4
#> 
#> ---hits: ASV694
#> ---hits: ASV46
#> ---hits: ASV1624
#> 
#> ---potential parent: ASV6944
#> 
#> ------checking: ASV6944
#> 
#> ------relative cooccurence: 14
#>  which is sufficient!4
#> 
#> ------min avg abundance: 0.03030303030303034
#> 
#> No parent found!
#> 4
#> 
#> ####processing: ASV18 #####4
#> 
#> ---hits: ASV831
#> ---hits: ASV8
#> ---hits: ASV474
#> ---hits: ASV666
#> ---hits: ASV170
#> ---hits: ASV261
#> ---hits: ASV493
#> ---hits: ASV94
#> ---hits: ASV209
#> ---hits: ASV3484
#> 
#> ---potential parent: ASV84
#> 
#> ------checking: ASV84
#> 
#> ------relative cooccurence: 14
#>  which is sufficient!4
#> 
#> ------min avg abundance: 0.7971153846153854
#> 
#> No parent found!
#> 4
#> 
#> ####processing: ASV694 #####4
#> 
#> ---hits: ASV175
#> ---hits: ASV46
#> ---hits: ASV1624
#> 
#> ---potential parent: ASV1754
#> 
#> ------checking: ASV1754
#> 
#> ------relative cooccurence: 14
#>  which is sufficient!4
#> 
#> ------min avg abundance: 14
#> 
#> No parent found!
#> 4
#> 
#> ####processing: ASV113 #####4
#> 
#> ---hits: ASV415
#> ---hits: ASV60
#> ---hits: ASV1380
#> ---hits: ASV984
#> ---hits: ASV796
#> ---hits: ASV1066
#> ---hits: ASV1020
#> ---hits: ASV1011
#> ---hits: ASV1036
#> ---hits: ASV15324
#> 
#> ---potential parent: 4
#> 
#> No parent found!
#> 4
#>   |                                                          |=                                                 |   3%
#> 
#> ####processing: ASV378 #####4
#> 
#> ---hits: ASV1303
#> ---hits: ASV1332
#> ---hits: ASV1562
#> ---hits: ASV1230
#> ---hits: ASV746
#> ---hits: ASV650
#> ---hits: ASV1030
#> ---hits: ASV1200
#> ---hits: ASV8324
#> 
#> ---potential parent: 4
#> 
#> No parent found!
#> 4
#> 
#> ####processing: ASV717 #####4
#> 
#> ---hits: ASV1058
#> ---hits: ASV1059
#> ---hits: ASV926
#> ---hits: ASV1493
#> ---hits: ASV1022
#> ---hits: ASV344
#> ---hits: ASV722
#> ---hits: ASV1459
#> ---hits: ASV89
#> ---hits: ASV4624
#> 
#> ---potential parent: 4
#> 
#> No parent found!
#> 4
#>   |                                                          |==                                                |   3%
#> 
#> ####processing: ASV41 #####4
#> 
#> ---hits: ASV672
#> ---hits: ASV251
#> ---hits: ASV10924
#> 
#> ---potential parent: 4
#> 
#> No parent found!
#> 4
#> 
#> ####processing: ASV65 #####4
#> 
#> ---hits: ASV322
#> ---hits: ASV3654
#> 
#> ---potential parent: 4
#> 
#> No parent found!
#> 4
#>   |                                                          |==                                                |   4%
#> 
#> ####processing: ASV89 #####4
#> 
#> ---hits: ASV462
#> ---hits: ASV534
#> ---hits: ASV313
#> ---hits: ASV1101
#> ---hits: ASV1557
#> ---hits: ASV1058
#> ---hits: ASV1459
#> ---hits: ASV1548
#> ---hits: ASV926
#> ---hits: ASV10594
#> 
#> ---potential parent: 4
#> 
#> No parent found!
#> 4
#> 
#> ####processing: ASV618 #####4
#> 
#> ---hits: ASV953
#> ---hits: ASV346
#> ---hits: ASV592
#> ---hits: ASV1261
#> ---hits: ASV906
#> ---hits: ASV38
#> ---hits: ASV1278
#> ---hits: ASV466
#> ---hits: ASV828
#> ---hits: ASV10074
#> 
#> ---potential parent: ASV38
#> ---potential parent: ASV9534
#> 
#> ------checking: ASV384
#> 
#> ------relative cooccurence: 0.64
#> 
#> ------checking: ASV9534
#> 
#> ------relative cooccurence: 14
#>  which is sufficient!4
#> 
#> ------min avg abundance: 0.4045801526717564
#> 
#> No parent found!
#> 4
#> 
#> ####processing: ASV170 #####4
#> 
#> ---hits: ASV474
#> ---hits: ASV8
#> ---hits: ASV209
#> ---hits: ASV18
#> ---hits: ASV666
#> ---hits: ASV94
#> ---hits: ASV493
#> ---hits: ASV831
#> ---hits: ASV261
#> ---hits: ASV3484
#> 
#> ---potential parent: ASV8
#> ---potential parent: ASV18
#> ---potential parent: ASV6664
#> 
#> ------checking: ASV84
#> 
#> ------relative cooccurence: 0.84
#> 
#> ------checking: ASV184
#> 
#> ------relative cooccurence: 0.84
#> 
#> ------checking: ASV6664
#> 
#> ------relative cooccurence: 0.64
#> 
#> No parent found!
#> 4
#> 
#> ####processing: ASV746 #####4
#> 
#> ---hits: ASV1303
#> ---hits: ASV1332
#> ---hits: ASV1230
#> ---hits: ASV378
#> ---hits: ASV1200
#> ---hits: ASV1030
#> ---hits: ASV1562
#> ---hits: ASV1458
#> ---hits: ASV693
#> ---hits: ASV6504
#> 
#> ---potential parent: ASV3784
#> 
#> ------checking: ASV3784
#> 
#> ------relative cooccurence: 0.44
#> 
#> No parent found!
#> 4
#>   |                                                          |==                                                |   5%
#> 
#> ####processing: ASV953 #####4
#> 
#> ---hits: ASV618
#> ---hits: ASV592
#> ---hits: ASV346
#> ---hits: ASV1261
#> ---hits: ASV906
#> ---hits: ASV38
#> ---hits: ASV1278
#> ---hits: ASV466
#> ---hits: ASV828
#> ---hits: ASV7564
#> 
#> ---potential parent: ASV38
#> ---potential parent: ASV756
#> ---potential parent: ASV6184
#> 
#> ------checking: ASV384
#> 
#> ------relative cooccurence: 0.64
#> 
#> ------checking: ASV7564
#> 
#> ------relative cooccurence: 0.64
#> 
#> ------checking: ASV6184
#> 
#> ------relative cooccurence: 14
#>  which is sufficient!4
#> 
#> ------min avg abundance: 0.6254
#> 
#> No parent found!
#> 4
#> 
#> ####processing: ASV28 #####4
#> 
#> ---hits: ASV3674
#> 
#> ---potential parent: 4
#> 
#> No parent found!
#> 4
#>   |                                                          |===                                               |   5%
#> 
#> ####processing: ASV310 #####4
#> 
#> ---hits: ASV543
#> ---hits: ASV685
#> ---hits: ASV199
#> ---hits: ASV692
#> ---hits: ASV10374
#> 
#> ---potential parent: 4
#> 
#> No parent found!
#> 4
#> 
#> ####processing: ASV666 #####4
#> 
#> ---hits: ASV831
#> ---hits: ASV8
#> ---hits: ASV18
#> ---hits: ASV170
#> ---hits: ASV493
#> ---hits: ASV261
#> ---hits: ASV474
#> ---hits: ASV209
#> ---hits: ASV94
#> ---hits: ASV3484
#> 
#> ---potential parent: ASV8
#> ---potential parent: ASV18
#> ---potential parent: ASV1704
#> 
#> ------checking: ASV84
#> 
#> ------relative cooccurence: 14
#>  which is sufficient!4
#> 
#> ------min avg abundance: 18.06254
#>  which is OK!4
#> 
#> SETTING ASV666 to be an ERROR of ASV8
#> 4
#> 
#> ------checking: ASV184
#> 
#> ------checking: ASV1704
#>   |                                                          |===                                               |   6%
#> 
#> ####processing: ASV1022 #####4
#> 
#> ---hits: ASV1493
#> ---hits: ASV1058
#> ---hits: ASV1059
#> ---hits: ASV717
#> ---hits: ASV926
#> ---hits: ASV344
#> ---hits: ASV722
#> ---hits: ASV14594
#> 
#> ---potential parent: ASV7174
#> 
#> ------checking: ASV7174
#> 
#> ------relative cooccurence: 14
#>  which is sufficient!4
#> 
#> ------min avg abundance: 0.2142857142857144
#> 
#> No parent found!
#> 4
#> 
#> ####processing: ASV1192 #####4
#> 
#> ---hits: ASV1074
#> ---hits: ASV1278
#> ---hits: ASV906
#> ---hits: ASV1007
#> ---hits: ASV1359
#> ---hits: ASV828
#> ---hits: ASV1224
#> ---hits: ASV1319
#> ---hits: ASV756
#> ---hits: ASV9504
#> 
#> ---potential parent: ASV756
#> ---potential parent: ASV10744
#> 
#> ------checking: ASV7564
#> 
#> ------relative cooccurence: 0.44
#> 
#> ------checking: ASV10744
#> 
#> ------relative cooccurence: 0.64
#> 
#> No parent found!
#> 4
#> 
#> ####processing: ASV727 #####4
#> 
#> ---hits: ASV546
#> ---hits: ASV569
#> ---hits: ASV1396
#> ---hits: ASV1233
#> ---hits: ASV13114
#> 
#> ---potential parent: 4
#> 
#> No parent found!
#> 4
#> 
#> ####processing: ASV989 #####4
#> 
#> ---hits: ASV38
#> ---hits: ASV816
#> ---hits: ASV466
#> ---hits: ASV853
#> ---hits: ASV975
#> ---hits: ASV1526
#> ---hits: ASV254
#> ---hits: ASV500
#> ---hits: ASV1107
#> ---hits: ASV1164
#> 
#> ---potential parent: ASV384
#> 
#> ------checking: ASV384
#> 
#> ------relative cooccurence: 14
#>  which is sufficient!4
#> 
#> ------min avg abundance: 6.333333333333334
#>  which is OK!4
#> 
#> SETTING ASV989 to be an ERROR of ASV38
#> 4
#>   |                                                          |===                                               |   7%
#> 
#> ####processing: ASV635 #####4
#> 
#> ---hits: ASV1356
#> ---hits: ASV12794
#> 
#> ---potential parent: 4
#> 
#> No parent found!
#> 4
#> 
#> ####processing: ASV1074 #####4
#> 
#> ---hits: ASV906
#> ---hits: ASV1278
#> ---hits: ASV1359
#> ---hits: ASV1007
#> ---hits: ASV1192
#> ---hits: ASV1224
#> ---hits: ASV828
#> ---hits: ASV756
#> ---hits: ASV1319
#> ---hits: ASV9504
#> 
#> ---potential parent: ASV756
#> ---potential parent: ASV11924
#> 
#> ------checking: ASV7564
#> 
#> ------relative cooccurence: 0.44
#> 
#> ------checking: ASV11924
#> 
#> ------relative cooccurence: 0.64
#> 
#> No parent found!
#> 4
#>   |                                                          |====                                              |   7%
#> 
#> ####processing: ASV1108 #####4
#> 
#> ---hits: ASV1332
#> ---hits: ASV7464
#> 
#> ---potential parent: ASV7464
#> 
#> ------checking: ASV7464
#> 
#> ------relative cooccurence: 0.84
#> 
#> No parent found!
#> 4
#> 
#> ####processing: ASV1483 #####4
#> 
#> ---hits: 4
#> 
#> ---potential parent: 4
#> 
#> No parent found!
#> 4
#>   |                                                          |====                                              |   8%
#> 
#> ####processing: ASV43 #####4
#> 
#> ---hits: ASV201
#> ---hits: ASV507
#> ---hits: ASV542
#> ---hits: ASV999
#> ---hits: ASV1387
#> ---hits: ASV915
#> ---hits: ASV1117
#> ---hits: ASV139
#> ---hits: ASV244
#> ---hits: ASV14194
#> 
#> ---potential parent: 4
#> 
#> No parent found!
#> 4
#> 
#> ####processing: ASV46 #####4
#> 
#> ---hits: ASV175
#> ---hits: ASV162
#> ---hits: ASV6944
#> 
#> ---potential parent: ASV175
#> ---potential parent: ASV6944
#> 
#> ------checking: ASV1754
#> 
#> ------relative cooccurence: 14
#>  which is sufficient!4
#> 
#> ------min avg abundance: 0.2801358234295424
#> 
#> ------checking: ASV6944
#> 
#> ------relative cooccurence: 14
#>  which is sufficient!4
#> 
#> ------min avg abundance: 0.008488964346349754
#> 
#> No parent found!
#> 4
#> 
#> ####processing: ASV365 #####4
#> 
#> ---hits: ASV65
#> ---hits: ASV3224
#> 
#> ---potential parent: ASV654
#> 
#> ------checking: ASV654
#> 
#> ------relative cooccurence: 14
#>  which is sufficient!4
#> 
#> ------min avg abundance: 14
#> 
#> No parent found!
#> 4
#> 
#> ####processing: ASV94 #####4
#> 
#> ---hits: ASV209
#> ---hits: ASV261
#> ---hits: ASV493
#> ---hits: ASV170
#> ---hits: ASV474
#> ---hits: ASV348
#> ---hits: ASV1365
#> ---hits: ASV18
#> ---hits: ASV8
#> ---hits: ASV8314
#> 
#> ---potential parent: ASV8
#> ---potential parent: ASV18
#> ---potential parent: ASV170
#> ---potential parent: ASV209
#> ---potential parent: ASV2614
#> 
#> ------checking: ASV84
#> 
#> ------relative cooccurence: 14
#>  which is sufficient!4
#> 
#> ------min avg abundance: 2.829113924050634
#>  which is OK!4
#> 
#> SETTING ASV94 to be an ERROR of ASV8
#> 4
#> 
#> ------checking: ASV184
#> 
#> ------checking: ASV1704
#> 
#> ------checking: ASV2094
#> 
#> ------checking: ASV2614
#>   |                                                          |====                                              |   9%
#> 
#> ####processing: ASV602 #####4
#> 
#> ---hits: 4
#> 
#> ---potential parent: 4
#> 
#> No parent found!
#> 4
#> 
#> ####processing: ASV178 #####4
#> 
#> ---hits: ASV1131
#> ---hits: ASV19
#> ---hits: ASV333
#> ---hits: ASV594
#> ---hits: ASV1467
#> ---hits: ASV1078
#> ---hits: ASV1010
#> ---hits: ASV858
#> ---hits: ASV986
#> ---hits: ASV9794
#> 
#> ---potential parent: ASV194
#> 
#> ------checking: ASV194
#> 
#> ------relative cooccurence: 14
#>  which is sufficient!4
#> 
#> ------min avg abundance: 1.254
#>  which is OK!4
#> 
#> SETTING ASV178 to be an ERROR of ASV19
#> 4
#>   |                                                          |=====                                             |   9%
#> 
#> ####processing: ASV172 #####4
#> 
#> ---hits: ASV963
#> ---hits: ASV981
#> ---hits: ASV715
#> ---hits: ASV1369
#> ---hits: ASV7374
#> 
#> ---potential parent: 4
#> 
#> No parent found!
#> 4
#> 
#> ####processing: ASV209 #####4
#> 
#> ---hits: ASV94
#> ---hits: ASV493
#> ---hits: ASV170
#> ---hits: ASV261
#> ---hits: ASV348
#> ---hits: ASV8
#> ---hits: ASV474
#> ---hits: ASV1365
#> ---hits: ASV666
#> ---hits: ASV184
#> 
#> ---potential parent: ASV8
#> ---potential parent: ASV18
#> ---potential parent: ASV170
#> ---potential parent: ASV666
#> ---potential parent: ASV94
#> ---potential parent: ASV2614
#> 
#> ------checking: ASV84
#> 
#> ------relative cooccurence: 14
#>  which is sufficient!4
#> 
#> ------min avg abundance: 6.807106598984774
#>  which is OK!4
#> 
#> SETTING ASV209 to be an ERROR of ASV8
#> 4
#> 
#> ------checking: ASV184
#> 
#> ------checking: ASV1704
#> 
#> ------checking: ASV6664
#> 
#> ------checking: ASV944
#> 
#> ------checking: ASV2614
#> 
#> ####processing: ASV741 #####4
#> 
#> ---hits: 4
#> 
#> ---potential parent: 4
#> 
#> No parent found!
#> 4
#>   |                                                          |=====                                             |  10%
#> 
#> ####processing: ASV1058 #####4
#> 
#> ---hits: ASV1059
#> ---hits: ASV717
#> ---hits: ASV1493
#> ---hits: ASV926
#> ---hits: ASV1022
#> ---hits: ASV344
#> ---hits: ASV722
#> ---hits: ASV1459
#> ---hits: ASV89
#> ---hits: ASV4624
#> 
#> ---potential parent: ASV717
#> ---potential parent: ASV89
#> ---potential parent: ASV10224
#> 
#> ------checking: ASV7174
#> 
#> ------relative cooccurence: 14
#>  which is sufficient!4
#> 
#> ------min avg abundance: 0.03092783505154644
#> 
#> ------checking: ASV894
#> 
#> ------relative cooccurence: 14
#>  which is sufficient!4
#> 
#> ------min avg abundance: 1.484536082474234
#>  which is OK!4
#> 
#> SETTING ASV1058 to be an ERROR of ASV89
#> 4
#> 
#> ------checking: ASV10224
#> 
#> ####processing: ASV261 #####4
#> 
#> ---hits: ASV493
#> ---hits: ASV94
#> ---hits: ASV209
#> ---hits: ASV18
#> ---hits: ASV8
#> ---hits: ASV831
#> ---hits: ASV666
#> ---hits: ASV170
#> ---hits: ASV474
#> ---hits: ASV3484
#> 
#> ---potential parent: ASV8
#> ---potential parent: ASV18
#> ---potential parent: ASV170
#> ---potential parent: ASV666
#> ---potential parent: ASV94
#> ---potential parent: ASV2094
#> 
#> ------checking: ASV84
#> 
#> ------relative cooccurence: 14
#>  which is sufficient!4
#> 
#> ------min avg abundance: 18.8873239436624
#>  which is OK!4
#> 
#> SETTING ASV261 to be an ERROR of ASV8
#> 4
#> 
#> ------checking: ASV184
#> 
#> ------checking: ASV1704
#> 
#> ------checking: ASV6664
#> 
#> ------checking: ASV944
#> 
#> ------checking: ASV2094
#> 
#> ####processing: ASV329 #####4
#> 
#> ---hits: 4
#> 
#> ---potential parent: 4
#> 
#> No parent found!
#> 4
#> 
#> ####processing: ASV906 #####4
#> 
#> ---hits: ASV1074
#> ---hits: ASV1278
#> ---hits: ASV828
#> ---hits: ASV1007
#> ---hits: ASV1359
#> ---hits: ASV1192
#> ---hits: ASV1319
#> ---hits: ASV950
#> ---hits: ASV756
#> ---hits: ASV12244
#> 
#> ---potential parent: ASV756
#> ---potential parent: ASV1192
#> ---potential parent: ASV10744
#> 
#> ------checking: ASV7564
#> 
#> ------relative cooccurence: 0.54
#> 
#> ------checking: ASV11924
#> 
#> ------relative cooccurence: 0.54
#> 
#> ------checking: ASV10744
#> 
#> ------relative cooccurence: 14
#>  which is sufficient!4
#> 
#> ------min avg abundance: 0.04166666666666674
#> 
#> No parent found!
#> 4
#>   |                                                          |=====                                             |  11%
#> 
#> ####processing: ASV546 #####4
#> 
#> ---hits: ASV727
#> ---hits: ASV569
#> ---hits: ASV1396
#> ---hits: ASV1233
#> ---hits: ASV13114
#> 
#> ---potential parent: ASV7274
#> 
#> ------checking: ASV7274
#> 
#> ------relative cooccurence: 14
#>  which is sufficient!4
#> 
#> ------min avg abundance: 0.14
#> 
#> No parent found!
#> 4
#> 
#> ####processing: ASV1253 #####4
#> 
#> ---hits: 4
#> 
#> ---potential parent: 4
#> 
#> No parent found!
#> 4
#>   |                                                          |======                                            |  11%
#> 
#> ####processing: ASV24 #####4
#> 
#> ---hits: ASV736
#> ---hits: ASV266
#> ---hits: ASV979
#> ---hits: ASV731
#> ---hits: ASV986
#> ---hits: ASV858
#> ---hits: ASV1010
#> ---hits: ASV19
#> ---hits: ASV333
#> ---hits: ASV1784
#> 
#> ---potential parent: ASV19
#> ---potential parent: ASV178
#> ---potential parent: ASV2664
#> 
#> ------checking: ASV194
#> 
#> ------relative cooccurence: 0.54
#> 
#> ------checking: ASV1784
#> 
#> ------relative cooccurence: 0.54
#> 
#> ------checking: ASV2664
#> 
#> ------relative cooccurence: 0.754
#> 
#> No parent found!
#> 4
#> 
#> ####processing: ASV1421 #####4
#> 
#> ---hits: ASV784
#> ---hits: ASV1315
#> ---hits: ASV1303
#> ---hits: ASV693
#> ---hits: ASV12304
#> 
#> ---potential parent: 4
#> 
#> No parent found!
#> 4
#>   |                                                          |======                                            |  12%
#> 
#> ####processing: ASV816 #####4
#> 
#> ---hits: ASV989
#> ---hits: ASV38
#> ---hits: ASV466
#> ---hits: ASV975
#> ---hits: ASV853
#> ---hits: ASV1526
#> ---hits: ASV254
#> ---hits: ASV500
#> ---hits: ASV1107
#> ---hits: ASV1164
#> 
#> ---potential parent: ASV38
#> ---potential parent: ASV9894
#> 
#> ------checking: ASV384
#> 
#> ------relative cooccurence: 14
#>  which is sufficient!4
#> 
#> ------min avg abundance: 194
#>  which is OK!4
#> 
#> SETTING ASV816 to be an ERROR of ASV38
#> 4
#> 
#> ------checking: ASV9894
#> 
#> ####processing: ASV662 #####4
#> 
#> ---hits: ASV724
#> ---hits: ASV420
#> ---hits: ASV505
#> ---hits: ASV509
#> ---hits: ASV790
#> ---hits: ASV1198
#> ---hits: ASV1052
#> ---hits: ASV1410
#> ---hits: ASV1427
#> ---hits: ASV7494
#> 
#> ---potential parent: 4
#> 
#> No parent found!
#> 4
#> 
#> ####processing: ASV266 #####4
#> 
#> ---hits: ASV24
#> ---hits: ASV731
#> ---hits: ASV736
#> ---hits: ASV979
#> ---hits: ASV986
#> ---hits: ASV858
#> ---hits: ASV1010
#> ---hits: ASV178
#> ---hits: ASV1467
#> ---hits: ASV11314
#> 
#> ---potential parent: ASV178
#> ---potential parent: ASV244
#> 
#> ------checking: ASV1784
#> 
#> ------relative cooccurence: 0.54
#> 
#> ------checking: ASV244
#> 
#> ------relative cooccurence: 0.754
#> 
#> No parent found!
#> 4
#> 
#> ####processing: ASV162 #####4
#> 
#> ---hits: ASV175
#> ---hits: ASV46
#> ---hits: ASV6944
#> 
#> ---potential parent: ASV175
#> ---potential parent: ASV694
#> ---potential parent: ASV464
#> 
#> ------checking: ASV1754
#> 
#> ------relative cooccurence: 14
#>  which is sufficient!4
#> 
#> ------min avg abundance: 0.167512690355334
#> 
#> ------checking: ASV6944
#> 
#> ------relative cooccurence: 14
#>  which is sufficient!4
#> 
#> ------min avg abundance: 0.00507614213197974
#> 
#> ------checking: ASV464
#> 
#> ------relative cooccurence: 14
#>  which is sufficient!4
#> 
#> ------min avg abundance: 0.5979695431472084
#> 
#> No parent found!
#> 4
#>   |                                                          |======                                            |  13%
#> 
#> ####processing: ASV139 #####4
#> 
#> ---hits: ASV244
#> ---hits: ASV1268
#> ---hits: ASV766
#> ---hits: ASV744
#> ---hits: ASV461
#> ---hits: ASV43
#> ---hits: ASV201
#> ---hits: ASV999
#> ---hits: ASV5074
#> 
#> ---potential parent: ASV43
#> ---potential parent: ASV201
#> ---potential parent: ASV244
#> ---potential parent: ASV5074
#> 
#> ------checking: ASV434
#> 
#> ------relative cooccurence: 14
#>  which is sufficient!4
#> 
#> ------min avg abundance: 0.01641079936474324
#> 
#> ------checking: ASV2014
#> 
#> ------relative cooccurence: 0.6666666666666674
#> 
#> ------checking: ASV2444
#> 
#> ------relative cooccurence: 14
#>  which is sufficient!4
#> 
#> ------min avg abundance: 0.09219858156028374
#> 
#> ------checking: ASV5074
#> 
#> ------relative cooccurence: 14
#>  which is sufficient!4
#> 
#> ------min avg abundance: 0.03546099290780144
#> 
#> No parent found!
#> 4
#> 
#> ####processing: ASV219 #####4
#> 
#> ---hits: ASV8994
#> 
#> ---potential parent: ASV8994
#> 
#> ------checking: ASV8994
#> 
#> ------relative cooccurence: 0.6666666666666674
#> 
#> No parent found!
#> 4
#>   |                                                          |=======                                           |  13%
#> 
#> ####processing: ASV201 #####4
#> 
#> ---hits: ASV43
#> ---hits: ASV542
#> ---hits: ASV1387
#> ---hits: ASV507
#> ---hits: ASV915
#> ---hits: ASV999
#> ---hits: ASV1117
#> ---hits: ASV244
#> ---hits: ASV1394
#> 
#> ---potential parent: ASV43
#> ---potential parent: ASV139
#> ---potential parent: ASV244
#> ---potential parent: ASV507
#> ---potential parent: ASV5424
#> 
#> ------checking: ASV434
#> 
#> ------relative cooccurence: 14
#>  which is sufficient!4
#> 
#> ------min avg abundance: 14
#> 
#> ------checking: ASV1394
#> 
#> ------relative cooccurence: 0.6666666666666674
#> 
#> ------checking: ASV2444
#> 
#> ------relative cooccurence: 0.6666666666666674
#> 
#> ------checking: ASV5074
#> 
#> ------relative cooccurence: 0.6666666666666674
#> 
#> ------checking: ASV5424
#> 
#> ------relative cooccurence: 0.6666666666666674
#> 
#> No parent found!
#> 4
#> 
#> ####processing: ASV251 #####4
#> 
#> ---hits: ASV1092
#> ---hits: ASV41
#> ---hits: ASV6724
#> 
#> ---potential parent: ASV414
#> 
#> ------checking: ASV414
#> 
#> ------relative cooccurence: 14
#>  which is sufficient!4
#> 
#> ------min avg abundance: 0.008038585209003214
#> 
#> No parent found!
#> 4
#>   |                                                          |=======                                           |  14%
#> 
#> ####processing: ASV244 #####4
#> 
#> ---hits: ASV139
#> ---hits: ASV1268
#> ---hits: ASV766
#> ---hits: ASV744
#> ---hits: ASV461
#> ---hits: ASV201
#> ---hits: ASV43
#> ---hits: ASV5424
#> 
#> ---potential parent: ASV43
#> ---potential parent: ASV139
#> ---potential parent: ASV201
#> ---potential parent: ASV5424
#> 
#> ------checking: ASV434
#> 
#> ------relative cooccurence: 14
#>  which is sufficient!4
#> 
#> ------min avg abundance: 0.0322580645161294
#> 
#> ------checking: ASV1394
#> 
#> ------relative cooccurence: 14
#>  which is sufficient!4
#> 
#> ------min avg abundance: 0.8387096774193554
#> 
#> ------checking: ASV2014
#> 
#> ------relative cooccurence: 0.6666666666666674
#> 
#> ------checking: ASV5424
#> 
#> ------relative cooccurence: 14
#>  which is sufficient!4
#> 
#> ------min avg abundance: 0.100519930675914
#> 
#> No parent found!
#> 4
#> 
#> ####processing: ASV507 #####4
#> 
#> ---hits: ASV542
#> ---hits: ASV999
#> ---hits: ASV43
#> ---hits: ASV201
#> ---hits: ASV1387
#> ---hits: ASV915
#> ---hits: ASV1117
#> ---hits: ASV1394
#> 
#> ---potential parent: ASV43
#> ---potential parent: ASV139
#> ---potential parent: ASV201
#> ---potential parent: ASV5424
#> 
#> ------checking: ASV434
#> 
#> ------relative cooccurence: 14
#>  which is sufficient!4
#> 
#> ------min avg abundance: 0.002463054187192124
#> 
#> ------checking: ASV1394
#> 
#> ------relative cooccurence: 14
#>  which is sufficient!4
#> 
#> ------min avg abundance: 0.06403940886699514
#> 
#> ------checking: ASV2014
#> 
#> ------relative cooccurence: 0.6666666666666674
#> 
#> ------checking: ASV5424
#> 
#> ------relative cooccurence: 14
#>  which is sufficient!4
#> 
#> ------min avg abundance: 0.5321100917431194
#> 
#> No parent found!
#> 4
#> 
#> ####processing: ASV483 #####4
#> 
#> ---hits: ASV854
#> 
#> ---potential parent: ASV854
#> 
#> ------checking: ASV854
#> 
#> ------relative cooccurence: 14
#>  which is sufficient!4
#> 
#> ------min avg abundance: 0.696202531645574
#> 
#> No parent found!
#> 4
#> 
#> ####processing: ASV500 #####4
#> 
#> ---hits: ASV254
#> ---hits: ASV55
#> ---hits: ASV816
#> ---hits: ASV989
#> ---hits: ASV1526
#> ---hits: ASV1107
#> ---hits: ASV116
#> ---hits: ASV38
#> ---hits: ASV466
#> ---hits: ASV8534
#> 
#> ---potential parent: ASV38
#> ---potential parent: ASV989
#> ---potential parent: ASV8164
#> 
#> ------checking: ASV384
#> 
#> ------relative cooccurence: 0.3333333333333334
#> 
#> ------checking: ASV9894
#> 
#> ------relative cooccurence: 04
#> 
#> ------checking: ASV8164
#> 
#> ------relative cooccurence: 04
#> 
#> No parent found!
#> 4
#>   |                                                          |=======                                           |  15%
#> 
#> ####processing: ASV577 #####4
#> 
#> ---hits: ASV8704
#> 
#> ---potential parent: ASV8704
#> 
#> ------checking: ASV8704
#> 
#> ------relative cooccurence: 14
#>  which is sufficient!4
#> 
#> ------min avg abundance: 0.4012738853503184
#> 
#> No parent found!
#> 4
#> 
#> ####processing: ASV344 #####4
#> 
#> ---hits: ASV722
#> ---hits: ASV926
#> ---hits: ASV717
#> ---hits: ASV1058
#> ---hits: ASV1059
#> ---hits: ASV1493
#> ---hits: ASV1022
#> ---hits: ASV1459
#> ---hits: ASV89
#> ---hits: ASV5344
#> 
#> ---potential parent: ASV717
#> ---potential parent: ASV89
#> ---potential parent: ASV1022
#> ---potential parent: ASV1058
#> ---potential parent: ASV926
#> ---potential parent: ASV1459
#> ---potential parent: ASV10594
#> 
#> ------checking: ASV7174
#> 
#> ------relative cooccurence: 14
#>  which is sufficient!4
#> 
#> ------min avg abundance: 0.01179941002949854
#> 
#> ------checking: ASV894
#> 
#> ------relative cooccurence: 0.6666666666666674
#> 
#> ------checking: ASV10224
#> 
#> ------relative cooccurence: 14
#>  which is sufficient!4
#> 
#> ------min avg abundance: 0.005899705014749264
#> 
#> ------checking: ASV10584
#> 
#> ------relative cooccurence: 0.6666666666666674
#> 
#> ------checking: ASV9264
#> 
#> ------relative cooccurence: 14
#>  which is sufficient!4
#> 
#> ------min avg abundance: 0.1710914454277294
#> 
#> ------checking: ASV14594
#> 
#> ------relative cooccurence: 0.3333333333333334
#> 
#> ------checking: ASV10594
#> 
#> ------relative cooccurence: 0.6666666666666674
#> 
#> No parent found!
#> 4
#>   |                                                          |========                                          |  15%
#> 
#> ####processing: ASV542 #####4
#> 
#> ---hits: ASV507
#> ---hits: ASV201
#> ---hits: ASV999
#> ---hits: ASV43
#> ---hits: ASV1387
#> ---hits: ASV915
#> ---hits: ASV1117
#> ---hits: ASV2444
#> 
#> ---potential parent: ASV43
#> ---potential parent: ASV201
#> ---potential parent: ASV244
#> ---potential parent: ASV5074
#> 
#> ------checking: ASV434
#> 
#> ------relative cooccurence: 14
#>  which is sufficient!4
#> 
#> ------min avg abundance: 0.002724795640326984
#> 
#> ------checking: ASV2014
#> 
#> ------relative cooccurence: 0.6666666666666674
#> 
#> ------checking: ASV2444
#> 
#> ------relative cooccurence: 14
#>  which is sufficient!4
#> 
#> ------min avg abundance: 0.08446866485013624
#> 
#> ------checking: ASV5074
#> 
#> ------relative cooccurence: 14
#>  which is sufficient!4
#> 
#> ------min avg abundance: 1.106267029972754
#>  which is OK!4
#> 
#> SETTING ASV542 to be an ERROR of ASV507
#> 4
#> 
#> ####processing: ASV85 #####4
#> 
#> ---hits: ASV4834
#> 
#> ---potential parent: ASV4834
#> 
#> ------checking: ASV4834
#> 
#> ------relative cooccurence: 14
#>  which is sufficient!4
#> 
#> ------min avg abundance: 1.315789473684214
#>  which is OK!4
#> 
#> SETTING ASV85 to be an ERROR of ASV483
#> 4
#>   |                                                          |========                                          |  16%
#> 
#> ####processing: ASV870 #####4
#> 
#> ---hits: ASV5774
#> 
#> ---potential parent: ASV5774
#> 
#> ------checking: ASV5774
#> 
#> ------relative cooccurence: 14
#>  which is sufficient!4
#> 
#> ------min avg abundance: 1.111111111111114
#>  which is OK!4
#> 
#> SETTING ASV870 to be an ERROR of ASV577
#> 4
#> 
#> ####processing: ASV899 #####4
#> 
#> ---hits: ASV2194
#> 
#> ---potential parent: ASV2194
#> 
#> ------checking: ASV2194
#> 
#> ------relative cooccurence: 0.6666666666666674
#> 
#> No parent found!
#> 4
#> 
#> ####processing: ASV52 #####4
#> 
#> ---hits: ASV159
#> ---hits: ASV880
#> ---hits: ASV1276
#> ---hits: ASV10704
#> 
#> ---potential parent: ASV1594
#> 
#> ------checking: ASV1594
#> 
#> ------relative cooccurence: 14
#>  which is sufficient!4
#> 
#> ------min avg abundance: 0.3809523809523814
#> 
#> No parent found!
#> 4
#> 
#> ####processing: ASV926 #####4
#> 
#> ---hits: ASV1058
#> ---hits: ASV1059
#> ---hits: ASV717
#> ---hits: ASV1493
#> ---hits: ASV1022
#> ---hits: ASV344
#> ---hits: ASV722
#> ---hits: ASV1459
#> ---hits: ASV89
#> ---hits: ASV4624
#> 
#> ---potential parent: ASV717
#> ---potential parent: ASV89
#> ---potential parent: ASV1022
#> ---potential parent: ASV1058
#> ---potential parent: ASV344
#> ---potential parent: ASV1459
#> ---potential parent: ASV10594
#> 
#> ------checking: ASV7174
#> 
#> ------relative cooccurence: 14
#>  which is sufficient!4
#> 
#> ------min avg abundance: 0.06896551724137934
#> 
#> ------checking: ASV894
#> 
#> ------relative cooccurence: 0.6666666666666674
#> 
#> ------checking: ASV10224
#> 
#> ------relative cooccurence: 14
#>  which is sufficient!4
#> 
#> ------min avg abundance: 0.03448275862068974
#> 
#> ------checking: ASV10584
#> 
#> ------relative cooccurence: 0.6666666666666674
#> 
#> ------checking: ASV3444
#> 
#> ------relative cooccurence: 14
#>  which is sufficient!4
#> 
#> ------min avg abundance: 0.6666666666666674
#> 
#> ------checking: ASV14594
#> 
#> ------relative cooccurence: 0.3333333333333334
#> 
#> ------checking: ASV10594
#> 
#> ------relative cooccurence: 0.6666666666666674
#> 
#> No parent found!
#> 4
#>   |                                                          |========                                          |  17%
#> 
#> ####processing: ASV159 #####4
#> 
#> ---hits: ASV52
#> ---hits: ASV880
#> ---hits: ASV1276
#> ---hits: ASV10704
#> 
#> ---potential parent: ASV524
#> 
#> ------checking: ASV524
#> 
#> ------relative cooccurence: 14
#>  which is sufficient!4
#> 
#> ------min avg abundance: 0.4285714285714294
#> 
#> No parent found!
#> 4
#> 
#> ####processing: ASV333 #####4
#> 
#> ---hits: ASV19
#> ---hits: ASV178
#> ---hits: ASV1010
#> ---hits: ASV594
#> ---hits: ASV1131
#> ---hits: ASV858
#> ---hits: ASV1467
#> ---hits: ASV1078
#> ---hits: ASV986
#> ---hits: ASV9794
#> 
#> ---potential parent: ASV19
#> ---potential parent: ASV178
#> ---potential parent: ASV10784
#> 
#> ------checking: ASV194
#> 
#> ------relative cooccurence: 14
#>  which is sufficient!4
#> 
#> ------min avg abundance: 44
#>  which is OK!4
#> 
#> SETTING ASV333 to be an ERROR of ASV19
#> 4
#> 
#> ------checking: ASV1784
#> 
#> ------checking: ASV10784
#>   |                                                          |=========                                         |  17%
#> 
#> ####processing: ASV631 #####4
#> 
#> ---hits: ASV1934
#> 
#> ---potential parent: 4
#> 
#> No parent found!
#> 4
#> 
#> ####processing: ASV348 #####4
#> 
#> ---hits: ASV1365
#> ---hits: ASV209
#> ---hits: ASV94
#> ---hits: ASV493
#> ---hits: ASV170
#> ---hits: ASV261
#> ---hits: ASV1182
#> ---hits: ASV8
#> ---hits: ASV26
#> ---hits: ASV4744
#> 
#> ---potential parent: ASV8
#> ---potential parent: ASV170
#> ---potential parent: ASV94
#> ---potential parent: ASV209
#> ---potential parent: ASV261
#> ---potential parent: ASV4934
#> 
#> ------checking: ASV84
#> 
#> ------relative cooccurence: 14
#>  which is sufficient!4
#> 
#> ------min avg abundance: 16.15662650602414
#>  which is OK!4
#> 
#> SETTING ASV348 to be an ERROR of ASV8
#> 4
#> 
#> ------checking: ASV1704
#> 
#> ------checking: ASV944
#> 
#> ------checking: ASV2094
#> 
#> ------checking: ASV2614
#> 
#> ------checking: ASV4934
#>   |                                                          |=========                                         |  18%
#> 
#> ####processing: ASV817 #####4
#> 
#> ---hits: 4
#> 
#> ---potential parent: 4
#> 
#> No parent found!
#> 4
#> 
#> ####processing: ASV1007 #####4
#> 
#> ---hits: ASV1359
#> ---hits: ASV1278
#> ---hits: ASV1074
#> ---hits: ASV906
#> ---hits: ASV1224
#> ---hits: ASV828
#> ---hits: ASV1192
#> ---hits: ASV1035
#> ---hits: ASV950
#> ---hits: ASV7564
#> 
#> ---potential parent: ASV756
#> ---potential parent: ASV1192
#> ---potential parent: ASV1074
#> ---potential parent: ASV906
#> ---potential parent: ASV1359
#> ---potential parent: ASV12784
#> 
#> ------checking: ASV7564
#> 
#> ------relative cooccurence: 0.3333333333333334
#> 
#> ------checking: ASV11924
#> 
#> ------relative cooccurence: 0.6666666666666674
#> 
#> ------checking: ASV10744
#> 
#> ------relative cooccurence: 14
#>  which is sufficient!4
#> 
#> ------min avg abundance: 0.02941176470588244
#> 
#> ------checking: ASV9064
#> 
#> ------relative cooccurence: 0.6666666666666674
#> 
#> ------checking: ASV13594
#> 
#> ------relative cooccurence: 0.6666666666666674
#> 
#> ------checking: ASV12784
#> 
#> ------relative cooccurence: 0.3333333333333334
#> 
#> No parent found!
#> 4
#> 
#> ####processing: ASV543 #####4
#> 
#> ---hits: ASV199
#> ---hits: ASV692
#> ---hits: ASV310
#> ---hits: ASV685
#> ---hits: ASV10374
#> 
#> ---potential parent: ASV3104
#> 
#> ------checking: ASV3104
#> 
#> ------relative cooccurence: 0.6666666666666674
#> 
#> No parent found!
#> 4
#> 
#> ####processing: ASV60 #####4
#> 
#> ---hits: ASV984
#> ---hits: ASV113
#> ---hits: ASV796
#> ---hits: ASV1020
#> ---hits: ASV415
#> ---hits: ASV1380
#> ---hits: ASV1011
#> ---hits: ASV1036
#> ---hits: ASV1066
#> ---hits: ASV15324
#> 
#> ---potential parent: ASV1134
#> 
#> ------checking: ASV1134
#> 
#> ------relative cooccurence: 0.6666666666666674
#> 
#> No parent found!
#> 4
#> 
#> ####processing: ASV559 #####4
#> 
#> ---hits: 4
#> 
#> ---potential parent: 4
#> 
#> No parent found!
#> 4
#>   |                                                          |=========                                         |  19%
#> 
#> ####processing: ASV1370 #####4
#> 
#> ---hits: 4
#> 
#> ---potential parent: 4
#> 
#> No parent found!
#> 4
#> 
#> ####processing: ASV1459 #####4
#> 
#> ---hits: ASV1059
#> ---hits: ASV1058
#> ---hits: ASV926
#> ---hits: ASV717
#> ---hits: ASV1493
#> ---hits: ASV462
#> ---hits: ASV1022
#> ---hits: ASV89
#> ---hits: ASV313
#> ---hits: ASV5344
#> 
#> ---potential parent: ASV717
#> ---potential parent: ASV89
#> ---potential parent: ASV1022
#> ---potential parent: ASV1058
#> ---potential parent: ASV926
#> ---potential parent: ASV10594
#> 
#> ------checking: ASV7174
#> 
#> ------relative cooccurence: 14
#>  which is sufficient!4
#> 
#> ------min avg abundance: 0.02083333333333334
#> 
#> ------checking: ASV894
#> 
#> ------relative cooccurence: 0.6666666666666674
#> 
#> ------checking: ASV10224
#> 
#> ------relative cooccurence: 0.6666666666666674
#> 
#> ------checking: ASV10584
#> 
#> ------relative cooccurence: 0.6666666666666674
#> 
#> ------checking: ASV9264
#> 
#> ------relative cooccurence: 0.3333333333333334
#> 
#> ------checking: ASV10594
#> 
#> ------relative cooccurence: 0.6666666666666674
#> 
#> No parent found!
#> 4
#>   |                                                          |==========                                        |  19%
#> 
#> ####processing: ASV1359 #####4
#> 
#> ---hits: ASV1278
#> ---hits: ASV1007
#> ---hits: ASV1224
#> ---hits: ASV1074
#> ---hits: ASV906
#> ---hits: ASV828
#> ---hits: ASV1192
#> ---hits: ASV950
#> ---hits: ASV756
#> ---hits: ASV13194
#> 
#> ---potential parent: ASV756
#> ---potential parent: ASV1192
#> ---potential parent: ASV1074
#> ---potential parent: ASV906
#> ---potential parent: ASV1007
#> ---potential parent: ASV12784
#> 
#> ------checking: ASV7564
#> 
#> ------relative cooccurence: 0.6666666666666674
#> 
#> ------checking: ASV11924
#> 
#> ------relative cooccurence: 0.3333333333333334
#> 
#> ------checking: ASV10744
#> 
#> ------relative cooccurence: 0.6666666666666674
#> 
#> ------checking: ASV9064
#> 
#> ------relative cooccurence: 0.3333333333333334
#> 
#> ------checking: ASV10074
#> 
#> ------relative cooccurence: 0.6666666666666674
#> 
#> ------checking: ASV12784
#> 
#> ------relative cooccurence: 0.3333333333333334
#> 
#> No parent found!
#> 4
#> 
#> ####processing: ASV292 #####4
#> 
#> ---hits: ASV368
#> ---hits: ASV98
#> ---hits: ASV961
#> ---hits: ASV443
#> ---hits: ASV566
#> ---hits: ASV1111
#> ---hits: ASV11764
#> 
#> ---potential parent: ASV443
#> ---potential parent: ASV98
#> ---potential parent: ASV368
#> ---potential parent: ASV5664
#> 
#> ------checking: ASV4434
#> 
#> ------relative cooccurence: 14
#>  which is sufficient!4
#> 
#> ------min avg abundance: 0.7727272727272734
#> 
#> ------checking: ASV984
#> 
#> ------relative cooccurence: 14
#>  which is sufficient!4
#> 
#> ------min avg abundance: 0.5909090909090914
#> 
#> ------checking: ASV3684
#> 
#> ------relative cooccurence: 0.6666666666666674
#> 
#> ------checking: ASV5664
#> 
#> ------relative cooccurence: 14
#>  which is sufficient!4
#> 
#> ------min avg abundance: 0.24
#> 
#> No parent found!
#> 4
#>   |                                                          |==========                                        |  20%
#> 
#> ####processing: ASV1059 #####4
#> 
#> ---hits: ASV1058
#> ---hits: ASV717
#> ---hits: ASV1493
#> ---hits: ASV926
#> ---hits: ASV1022
#> ---hits: ASV344
#> ---hits: ASV722
#> ---hits: ASV1459
#> ---hits: ASV89
#> ---hits: ASV4624
#> 
#> ---potential parent: ASV717
#> ---potential parent: ASV89
#> ---potential parent: ASV1022
#> ---potential parent: ASV1058
#> ---potential parent: ASV344
#> ---potential parent: ASV926
#> ---potential parent: ASV14594
#> 
#> ------checking: ASV7174
#> 
#> ------relative cooccurence: 14
#>  which is sufficient!4
#> 
#> ------min avg abundance: 0.2142857142857144
#> 
#> ------checking: ASV894
#> 
#> ------relative cooccurence: 0.6666666666666674
#> 
#> ------checking: ASV10224
#> 
#> ------relative cooccurence: 14
#>  which is sufficient!4
#> 
#> ------min avg abundance: 0.1428571428571434
#> 
#> ------checking: ASV10584
#> 
#> ------relative cooccurence: 0.6666666666666674
#> 
#> ------checking: ASV3444
#> 
#> ------relative cooccurence: 0.6666666666666674
#> 
#> ------checking: ASV9264
#> 
#> ------relative cooccurence: 0.6666666666666674
#> 
#> ------checking: ASV14594
#> 
#> ------relative cooccurence: 0.6666666666666674
#> 
#> No parent found!
#> 4
#> 
#> ####processing: ASV443 #####4
#> 
#> ---hits: ASV368
#> ---hits: ASV1111
#> ---hits: ASV1176
#> ---hits: ASV292
#> ---hits: ASV566
#> ---hits: ASV98
#> ---hits: ASV9614
#> 
#> ---potential parent: ASV292
#> ---potential parent: ASV98
#> ---potential parent: ASV368
#> ---potential parent: ASV5664
#> 
#> ------checking: ASV2924
#> 
#> ------relative cooccurence: 14
#>  which is sufficient!4
#> 
#> ------min avg abundance: 0.54
#> 
#> ------checking: ASV984
#> 
#> ------relative cooccurence: 14
#>  which is sufficient!4
#> 
#> ------min avg abundance: 0.54
#> 
#> ------checking: ASV3684
#> 
#> ------relative cooccurence: 0.6666666666666674
#> 
#> ------checking: ASV5664
#> 
#> ------relative cooccurence: 14
#>  which is sufficient!4
#> 
#> ------min avg abundance: 0.1538461538461544
#> 
#> No parent found!
#> 4
#> 
#> ####processing: ASV505 #####4
#> 
#> ---hits: ASV724
#> ---hits: ASV420
#> ---hits: ASV509
#> ---hits: ASV790
#> ---hits: ASV662
#> ---hits: ASV1052
#> ---hits: ASV1410
#> ---hits: ASV1427
#> ---hits: ASV1198
#> ---hits: ASV6414
#> 
#> ---potential parent: ASV662
#> ---potential parent: ASV509
#> ---potential parent: ASV7244
#> 
#> ------checking: ASV6624
#> 
#> ------relative cooccurence: 0.6666666666666674
#> 
#> ------checking: ASV5094
#> 
#> ------relative cooccurence: 14
#>  which is sufficient!4
#> 
#> ------min avg abundance: 0.1754
#> 
#> ------checking: ASV7244
#> 
#> ------relative cooccurence: 14
#>  which is sufficient!4
#> 
#> ------min avg abundance: 0.054
#> 
#> No parent found!
#> 4
#> 
#> ####processing: ASV509 #####4
#> 
#> ---hits: ASV420
#> ---hits: ASV505
#> ---hits: ASV790
#> ---hits: ASV662
#> ---hits: ASV1052
#> ---hits: ASV1410
#> ---hits: ASV1427
#> ---hits: ASV724
#> ---hits: ASV1198
#> ---hits: ASV6414
#> 
#> ---potential parent: ASV662
#> ---potential parent: ASV505
#> ---potential parent: ASV7244
#> 
#> ------checking: ASV6624
#> 
#> ------relative cooccurence: 0.6666666666666674
#> 
#> ------checking: ASV5054
#> 
#> ------relative cooccurence: 14
#>  which is sufficient!4
#> 
#> ------min avg abundance: 0.084
#> 
#> ------checking: ASV7244
#> 
#> ------relative cooccurence: 14
#>  which is sufficient!4
#> 
#> ------min avg abundance: 0.2857142857142864
#> 
#> No parent found!
#> 4
#>   |                                                          |==========                                        |  21%
#> 
#> ####processing: ASV1624 #####4
#> 
#> ---hits: ASV1644
#> ---hits: ASV13004
#> 
#> ---potential parent: ASV16444
#> 
#> ------checking: ASV16444
#> 
#> ------relative cooccurence: 0.3333333333333334
#> 
#> No parent found!
#> 4
#> 
#> ####processing: ASV98 #####4
#> 
#> ---hits: ASV566
#> ---hits: ASV961
#> ---hits: ASV292
#> ---hits: ASV368
#> ---hits: ASV443
#> ---hits: ASV1111
#> ---hits: ASV11764
#> 
#> ---potential parent: ASV292
#> ---potential parent: ASV443
#> ---potential parent: ASV368
#> ---potential parent: ASV5664
#> 
#> ------checking: ASV2924
#> 
#> ------relative cooccurence: 14
#>  which is sufficient!4
#> 
#> ------min avg abundance: 0.8333333333333334
#> 
#> ------checking: ASV4434
#> 
#> ------relative cooccurence: 14
#>  which is sufficient!4
#> 
#> ------min avg abundance: 1.083333333333334
#>  which is OK!4
#> 
#> SETTING ASV98 to be an ERROR of ASV443
#> 4
#> 
#> ------checking: ASV3684
#> 
#> ------checking: ASV5664
#>   |                                                          |===========                                       |  21%
#> 
#> ####processing: ASV368 #####4
#> 
#> ---hits: ASV292
#> ---hits: ASV443
#> ---hits: ASV566
#> ---hits: ASV961
#> ---hits: ASV98
#> ---hits: ASV1111
#> ---hits: ASV11764
#> 
#> ---potential parent: ASV292
#> ---potential parent: ASV443
#> ---potential parent: ASV98
#> ---potential parent: ASV5664
#> 
#> ------checking: ASV2924
#> 
#> ------relative cooccurence: 0.6666666666666674
#> 
#> ------checking: ASV4434
#> 
#> ------relative cooccurence: 0.6666666666666674
#> 
#> ------checking: ASV984
#> 
#> ------relative cooccurence: 0.6666666666666674
#> 
#> ------checking: ASV5664
#> 
#> ------relative cooccurence: 0.6666666666666674
#> 
#> No parent found!
#> 4
#> 
#> ####processing: ASV493 #####4
#> 
#> ---hits: ASV261
#> ---hits: ASV209
#> ---hits: ASV8
#> ---hits: ASV94
#> ---hits: ASV666
#> ---hits: ASV18
#> ---hits: ASV170
#> ---hits: ASV348
#> ---hits: ASV831
#> ---hits: ASV4744
#> 
#> ---potential parent: ASV8
#> ---potential parent: ASV18
#> ---potential parent: ASV170
#> ---potential parent: ASV666
#> ---potential parent: ASV94
#> ---potential parent: ASV209
#> ---potential parent: ASV261
#> ---potential parent: ASV348
#> ---potential parent: ASV8314
#> 
#> ------checking: ASV84
#> 
#> ------relative cooccurence: 14
#>  which is sufficient!4
#> 
#> ------min avg abundance: 47.89285714285714
#>  which is OK!4
#> 
#> SETTING ASV493 to be an ERROR of ASV8
#> 4
#> 
#> ------checking: ASV184
#> 
#> ------checking: ASV1704
#> 
#> ------checking: ASV6664
#> 
#> ------checking: ASV944
#> 
#> ------checking: ASV2094
#> 
#> ------checking: ASV2614
#> 
#> ------checking: ASV3484
#> 
#> ------checking: ASV8314
#>   |                                                          |===========                                       |  22%
#> 
#> ####processing: ASV1279 #####4
#> 
#> ---hits: ASV635
#> ---hits: ASV13564
#> 
#> ---potential parent: ASV635
#> ---potential parent: ASV13564
#> 
#> ------checking: ASV6354
#> 
#> ------relative cooccurence: 14
#>  which is sufficient!4
#> 
#> ------min avg abundance: 0.07692307692307694
#> 
#> ------checking: ASV13564
#> 
#> ------relative cooccurence: 14
#>  which is sufficient!4
#> 
#> ------min avg abundance: 0.07692307692307694
#> 
#> No parent found!
#> 4
#> 
#> ####processing: ASV904 #####4
#> 
#> ---hits: ASV756
#> ---hits: ASV1035
#> ---hits: ASV1074
#> ---hits: ASV906
#> ---hits: ASV1278
#> ---hits: ASV828
#> ---hits: ASV950
#> ---hits: ASV1007
#> ---hits: ASV1192
#> ---hits: ASV13594
#> 
#> ---potential parent: ASV756
#> ---potential parent: ASV1192
#> ---potential parent: ASV1074
#> ---potential parent: ASV906
#> ---potential parent: ASV1007
#> ---potential parent: ASV1359
#> ---potential parent: ASV12784
#> 
#> ------checking: ASV7564
#> 
#> ------relative cooccurence: 14
#>  which is sufficient!4
#> 
#> ------min avg abundance: 0.1212121212121214
#> 
#> ------checking: ASV11924
#> 
#> ------relative cooccurence: 04
#> 
#> ------checking: ASV10744
#> 
#> ------relative cooccurence: 0.3333333333333334
#> 
#> ------checking: ASV9064
#> 
#> ------relative cooccurence: 0.3333333333333334
#> 
#> ------checking: ASV10074
#> 
#> ------relative cooccurence: 0.3333333333333334
#> 
#> ------checking: ASV13594
#> 
#> ------relative cooccurence: 0.3333333333333334
#> 
#> ------checking: ASV12784
#> 
#> ------relative cooccurence: 0.3333333333333334
#> 
#> No parent found!
#> 4
#> 
#> ####processing: ASV563 #####4
#> 
#> ---hits: ASV431
#> ---hits: ASV1420
#> ---hits: ASV673
#> ---hits: ASV1218
#> ---hits: ASV358
#> ---hits: ASV405
#> ---hits: ASV694
#> 
#> ---potential parent: ASV694
#> 
#> ------checking: ASV694
#> 
#> ------relative cooccurence: 04
#> 
#> No parent found!
#> 4
#> 
#> ####processing: ASV1644 #####4
#> 
#> ---hits: ASV1624
#> ---hits: ASV13004
#> 
#> ---potential parent: ASV16244
#> 
#> ------checking: ASV16244
#> 
#> ------relative cooccurence: 0.3333333333333334
#> 
#> No parent found!
#> 4
#>   |                                                          |===========                                       |  23%
#> 
#> ####processing: ASV1369 #####4
#> 
#> ---hits: ASV737
#> ---hits: ASV715
#> ---hits: ASV981
#> ---hits: ASV963
#> ---hits: ASV1724
#> 
#> ---potential parent: ASV1724
#> 
#> ------checking: ASV1724
#> 
#> ------relative cooccurence: 0.3333333333333334
#> 
#> No parent found!
#> 4
#> 
#> ####processing: ASV1409 #####4
#> 
#> ---hits: ASV772
#> ---hits: ASV504
#> ---hits: ASV860
#> ---hits: ASV832
#> ---hits: ASV693
#> ---hits: ASV1332
#> ---hits: ASV1200
#> ---hits: ASV1458
#> ---hits: ASV1303
#> ---hits: ASV10304
#> 
#> ---potential parent: 4
#> 
#> No parent found!
#> 4
#>   |                                                          |============                                      |  23%
#> 
#> ####processing: ASV724 #####4
#> 
#> ---hits: ASV505
#> ---hits: ASV662
#> ---hits: ASV1052
#> ---hits: ASV420
#> ---hits: ASV509
#> ---hits: ASV790
#> ---hits: ASV1198
#> ---hits: ASV859
#> ---hits: ASV1410
#> ---hits: ASV14274
#> 
#> ---potential parent: ASV662
#> ---potential parent: ASV505
#> ---potential parent: ASV509
#> ---potential parent: ASV8594
#> 
#> ------checking: ASV6624
#> 
#> ------relative cooccurence: 0.6666666666666674
#> 
#> ------checking: ASV5054
#> 
#> ------relative cooccurence: 14
#>  which is sufficient!4
#> 
#> ------min avg abundance: 0.1666666666666674
#> 
#> ------checking: ASV5094
#> 
#> ------relative cooccurence: 14
#>  which is sufficient!4
#> 
#> ------min avg abundance: 24
#>  which is OK!4
#> 
#> SETTING ASV724 to be an ERROR of ASV509
#> 4
#> 
#> ------checking: ASV8594
#> 
#> ####processing: ASV255 #####4
#> 
#> ---hits: ASV854
#> ---hits: ASV5704
#> 
#> ---potential parent: 4
#> 
#> No parent found!
#> 4
#>   |                                                          |============                                      |  24%
#> 
#> ####processing: ASV69 #####4
#> 
#> ---hits: ASV405
#> ---hits: ASV358
#> ---hits: ASV914
#> ---hits: ASV930
#> ---hits: ASV221
#> ---hits: ASV833
#> ---hits: ASV1420
#> ---hits: ASV673
#> ---hits: ASV1218
#> ---hits: ASV4314
#> 
#> ---potential parent: 4
#> 
#> No parent found!
#> 4
#> 
#> ####processing: ASV151 #####4
#> 
#> ---hits: ASV355
#> ---hits: ASV8344
#> 
#> ---potential parent: 4
#> 
#> No parent found!
#> 4
#> 
#> ####processing: ASV566 #####4
#> 
#> ---hits: ASV98
#> ---hits: ASV368
#> ---hits: ASV961
#> ---hits: ASV292
#> ---hits: ASV443
#> ---hits: ASV1111
#> ---hits: ASV11764
#> 
#> ---potential parent: ASV292
#> ---potential parent: ASV443
#> ---potential parent: ASV98
#> ---potential parent: ASV3684
#> 
#> ------checking: ASV2924
#> 
#> ------relative cooccurence: 14
#>  which is sufficient!4
#> 
#> ------min avg abundance: 14
#> 
#> ------checking: ASV4434
#> 
#> ------relative cooccurence: 14
#>  which is sufficient!4
#> 
#> ------min avg abundance: 24
#>  which is OK!4
#> 
#> SETTING ASV566 to be an ERROR of ASV443
#> 4
#> 
#> ------checking: ASV984
#> 
#> ------checking: ASV3684
#> 
#> ####processing: ASV831 #####4
#> 
#> ---hits: ASV18
#> ---hits: ASV666
#> ---hits: ASV8
#> ---hits: ASV474
#> ---hits: ASV170
#> ---hits: ASV261
#> ---hits: ASV493
#> ---hits: ASV94
#> ---hits: ASV209
#> ---hits: ASV3484
#> 
#> ---potential parent: ASV8
#> ---potential parent: ASV18
#> ---potential parent: ASV170
#> ---potential parent: ASV666
#> ---potential parent: ASV94
#> ---potential parent: ASV209
#> ---potential parent: ASV261
#> ---potential parent: ASV348
#> ---potential parent: ASV4934
#> 
#> ------checking: ASV84
#> 
#> ------relative cooccurence: 14
#>  which is sufficient!4
#> 
#> ------min avg abundance: 103.6254
#>  which is OK!4
#> 
#> SETTING ASV831 to be an ERROR of ASV8
#> 4
#> 
#> ------checking: ASV184
#> 
#> ------checking: ASV1704
#> 
#> ------checking: ASV6664
#> 
#> ------checking: ASV944
#> 
#> ------checking: ASV2094
#> 
#> ------checking: ASV2614
#> 
#> ------checking: ASV3484
#> 
#> ------checking: ASV4934
#>   |                                                          |============                                      |  25%
#> 
#> ####processing: ASV1078 #####4
#> 
#> ---hits: ASV1467
#> ---hits: ASV178
#> ---hits: ASV1131
#> ---hits: ASV1010
#> ---hits: ASV19
#> ---hits: ASV333
#> ---hits: ASV594
#> ---hits: ASV986
#> ---hits: ASV731
#> ---hits: ASV8584
#> 
#> ---potential parent: ASV19
#> ---potential parent: ASV178
#> ---potential parent: ASV3334
#> 
#> ------checking: ASV194
#> 
#> ------relative cooccurence: 14
#>  which is sufficient!4
#> 
#> ------min avg abundance: 44
#>  which is OK!4
#> 
#> SETTING ASV1078 to be an ERROR of ASV19
#> 4
#> 
#> ------checking: ASV1784
#> 
#> ------checking: ASV3334
#> 
#> ####processing: ASV1230 #####4
#> 
#> ---hits: ASV1303
#> ---hits: ASV746
#> ---hits: ASV1332
#> ---hits: ASV378
#> ---hits: ASV1458
#> ---hits: ASV1562
#> ---hits: ASV650
#> ---hits: ASV1030
#> ---hits: ASV693
#> ---hits: ASV14214
#> 
#> ---potential parent: ASV378
#> ---potential parent: ASV746
#> ---potential parent: ASV14214
#> 
#> ------checking: ASV3784
#> 
#> ------relative cooccurence: 04
#> 
#> ------checking: ASV7464
#> 
#> ------relative cooccurence: 0.3333333333333334
#> 
#> ------checking: ASV14214
#> 
#> ------relative cooccurence: 04
#> 
#> No parent found!
#> 4
#>   |                                                          |=============                                     |  25%
#> 
#> ####processing: ASV1004 #####4
#> 
#> ---hits: 4
#> 
#> ---potential parent: 4
#> 
#> No parent found!
#> 4
#> 
#> ####processing: ASV586 #####4
#> 
#> ---hits: ASV15764
#> 
#> ---potential parent: 4
#> 
#> No parent found!
#> 4
#>   |                                                          |=============                                     |  26%
#> 
#> ####processing: ASV859 #####4
#> 
#> ---hits: ASV1427
#> ---hits: ASV420
#> ---hits: ASV724
#> ---hits: ASV662
#> ---hits: ASV1052
#> ---hits: ASV505
#> ---hits: ASV509
#> ---hits: ASV790
#> ---hits: ASV1128
#> ---hits: ASV11984
#> 
#> ---potential parent: ASV662
#> ---potential parent: ASV505
#> ---potential parent: ASV509
#> ---potential parent: ASV7244
#> 
#> ------checking: ASV6624
#> 
#> ------relative cooccurence: 0.6666666666666674
#> 
#> ------checking: ASV5054
#> 
#> ------relative cooccurence: 14
#>  which is sufficient!4
#> 
#> ------min avg abundance: 14
#> 
#> ------checking: ASV5094
#> 
#> ------relative cooccurence: 14
#>  which is sufficient!4
#> 
#> ------min avg abundance: 3.54
#>  which is OK!4
#> 
#> SETTING ASV859 to be an ERROR of ASV509
#> 4
#> 
#> ------checking: ASV7244
#> 
#> ####processing: ASV1278 #####4
#> 
#> ---hits: ASV1359
#> ---hits: ASV1074
#> ---hits: ASV1224
#> ---hits: ASV906
#> ---hits: ASV1007
#> ---hits: ASV828
#> ---hits: ASV1192
#> ---hits: ASV950
#> ---hits: ASV756
#> ---hits: ASV13194
#> 
#> ---potential parent: ASV756
#> ---potential parent: ASV1192
#> ---potential parent: ASV1074
#> ---potential parent: ASV906
#> ---potential parent: ASV1007
#> ---potential parent: ASV13594
#> 
#> ------checking: ASV7564
#> 
#> ------relative cooccurence: 0.6666666666666674
#> 
#> ------checking: ASV11924
#> 
#> ------relative cooccurence: 0.6666666666666674
#> 
#> ------checking: ASV10744
#> 
#> ------relative cooccurence: 0.3333333333333334
#> 
#> ------checking: ASV9064
#> 
#> ------relative cooccurence: 04
#> 
#> ------checking: ASV10074
#> 
#> ------relative cooccurence: 0.3333333333333334
#> 
#> ------checking: ASV13594
#> 
#> ------relative cooccurence: 0.3333333333333334
#> 
#> No parent found!
#> 4
#> 
#> ####processing: ASV1356 #####4
#> 
#> ---hits: ASV635
#> ---hits: ASV12794
#> 
#> ---potential parent: ASV635
#> ---potential parent: ASV12794
#> 
#> ------checking: ASV6354
#> 
#> ------relative cooccurence: 14
#>  which is sufficient!4
#> 
#> ------min avg abundance: 14
#> 
#> ------checking: ASV12794
#> 
#> ------relative cooccurence: 14
#>  which is sufficient!4
#> 
#> ------min avg abundance: 24
#>  which is OK!4
#> 
#> SETTING ASV1356 to be an ERROR of ASV1279
#> 4
#> 
#> ####processing: ASV801 #####4
#> 
#> ---hits: ASV1461
#> ---hits: ASV1155
#> ---hits: ASV1632
#> ---hits: ASV13404
#> 
#> ---potential parent: 4
#> 
#> No parent found!
#> 4
#>   |                                                          |=============                                     |  27%
#> 
#> ####processing: ASV64 #####4
#> 
#> ---hits: ASV168
#> ---hits: ASV249
#> ---hits: ASV580
#> ---hits: ASV2034
#> 
#> ---potential parent: 4
#> 
#> No parent found!
#> 4
#> 
#> ####processing: ASV144 #####4
#> 
#> ---hits: ASV6164
#> 
#> ---potential parent: 4
#> 
#> No parent found!
#> 4
#>   |                                                          |==============                                    |  27%
#> 
#> ####processing: ASV131 #####4
#> 
#> ---hits: 4
#> 
#> ---potential parent: 4
#> 
#> No parent found!
#> 4
#> 
#> ####processing: ASV202 #####4
#> 
#> ---hits: 4
#> 
#> ---potential parent: 4
#> 
#> No parent found!
#> 4
#> 
#> ####processing: ASV107 #####4
#> 
#> ---hits: ASV12
#> ---hits: ASV67
#> ---hits: ASV1542
#> ---hits: ASV624
#> ---hits: ASV1262
#> ---hits: ASV847
#> ---hits: ASV1054
#> 
#> ---potential parent: ASV124
#> 
#> ------checking: ASV124
#> 
#> ------relative cooccurence: 14
#>  which is sufficient!4
#> 
#> ------min avg abundance: 0.024
#> 
#> No parent found!
#> 4
#>   |                                                          |==============                                    |  28%
#> 
#> ####processing: ASV322 #####4
#> 
#> ---hits: ASV65
#> ---hits: ASV3654
#> 
#> ---potential parent: ASV65
#> ---potential parent: ASV3654
#> 
#> ------checking: ASV654
#> 
#> ------relative cooccurence: 14
#>  which is sufficient!4
#> 
#> ------min avg abundance: 0.1707317073170734
#> 
#> ------checking: ASV3654
#> 
#> ------relative cooccurence: 14
#>  which is sufficient!4
#> 
#> ------min avg abundance: 0.1707317073170734
#> 
#> No parent found!
#> 4
#> 
#> ####processing: ASV444 #####4
#> 
#> ---hits: ASV9984
#> 
#> ---potential parent: 4
#> 
#> No parent found!
#> 4
#> 
#> ####processing: ASV313 #####4
#> 
#> ---hits: ASV462
#> ---hits: ASV534
#> ---hits: ASV89
#> ---hits: ASV1101
#> ---hits: ASV1557
#> ---hits: ASV1459
#> ---hits: ASV1058
#> ---hits: ASV1548
#> ---hits: ASV926
#> ---hits: ASV10594
#> 
#> ---potential parent: ASV89
#> ---potential parent: ASV1058
#> ---potential parent: ASV926
#> ---potential parent: ASV1459
#> ---potential parent: ASV1059
#> ---potential parent: ASV462
#> ---potential parent: ASV534
#> ---potential parent: ASV1101
#> ---potential parent: ASV15574
#> 
#> ------checking: ASV894
#> 
#> ------relative cooccurence: 14
#>  which is sufficient!4
#> 
#> ------min avg abundance: 2.124324324324324
#>  which is OK!4
#> 
#> SETTING ASV313 to be an ERROR of ASV89
#> 4
#> 
#> ------checking: ASV10584
#> 
#> ------checking: ASV9264
#> 
#> ------checking: ASV14594
#> 
#> ------checking: ASV10594
#> 
#> ------checking: ASV4624
#> 
#> ------checking: ASV5344
#> 
#> ------checking: ASV11014
#> 
#> ------checking: ASV15574
#> 
#> ####processing: ASV569 #####4
#> 
#> ---hits: ASV546
#> ---hits: ASV727
#> ---hits: ASV1233
#> ---hits: ASV1396
#> ---hits: ASV13114
#> 
#> ---potential parent: ASV727
#> ---potential parent: ASV546
#> ---potential parent: ASV13964
#> 
#> ------checking: ASV7274
#> 
#> ------relative cooccurence: 04
#> 
#> ------checking: ASV5464
#> 
#> ------relative cooccurence: 04
#> 
#> ------checking: ASV13964
#> 
#> ------relative cooccurence: 0.54
#> 
#> No parent found!
#> 4
#>   |                                                          |==============                                    |  29%
#> 
#> ####processing: ASV672 #####4
#> 
#> ---hits: ASV41
#> ---hits: ASV251
#> ---hits: ASV10924
#> 
#> ---potential parent: ASV41
#> ---potential parent: ASV2514
#> 
#> ------checking: ASV414
#> 
#> ------relative cooccurence: 14
#>  which is sufficient!4
#> 
#> ------min avg abundance: 0.01880877742946714
#> 
#> ------checking: ASV2514
#> 
#> ------relative cooccurence: 0.54
#> 
#> No parent found!
#> 4
#> 
#> ####processing: ASV462 #####4
#> 
#> ---hits: ASV89
#> ---hits: ASV313
#> ---hits: ASV534
#> ---hits: ASV1101
#> ---hits: ASV1557
#> ---hits: ASV1459
#> ---hits: ASV1058
#> ---hits: ASV1548
#> ---hits: ASV926
#> ---hits: ASV10594
#> 
#> ---potential parent: ASV89
#> ---potential parent: ASV1058
#> ---potential parent: ASV926
#> ---potential parent: ASV1459
#> ---potential parent: ASV1059
#> ---potential parent: ASV313
#> ---potential parent: ASV534
#> ---potential parent: ASV1101
#> ---potential parent: ASV15574
#> 
#> ------checking: ASV894
#> 
#> ------relative cooccurence: 14
#>  which is sufficient!4
#> 
#> ------min avg abundance: 3.661490683229814
#>  which is OK!4
#> 
#> SETTING ASV462 to be an ERROR of ASV89
#> 4
#> 
#> ------checking: ASV10584
#> 
#> ------checking: ASV9264
#> 
#> ------checking: ASV14594
#> 
#> ------checking: ASV10594
#> 
#> ------checking: ASV3134
#> 
#> ------checking: ASV5344
#> 
#> ------checking: ASV11014
#> 
#> ------checking: ASV15574
#>   |                                                          |===============                                   |  29%
#> 
#> ####processing: ASV223 #####4
#> 
#> ---hits: 4
#> 
#> ---potential parent: 4
#> 
#> No parent found!
#> 4
#> 
#> ####processing: ASV594 #####4
#> 
#> ---hits: ASV178
#> ---hits: ASV19
#> ---hits: ASV333
#> ---hits: ASV1131
#> ---hits: ASV1467
#> ---hits: ASV1010
#> ---hits: ASV858
#> ---hits: ASV1078
#> ---hits: ASV986
#> ---hits: ASV9794
#> 
#> ---potential parent: ASV19
#> ---potential parent: ASV178
#> ---potential parent: ASV333
#> ---potential parent: ASV1078
#> ---potential parent: ASV11314
#> 
#> ------checking: ASV194
#> 
#> ------relative cooccurence: 14
#>  which is sufficient!4
#> 
#> ------min avg abundance: 10.18150684931514
#>  which is OK!4
#> 
#> SETTING ASV594 to be an ERROR of ASV19
#> 4
#> 
#> ------checking: ASV1784
#> 
#> ------checking: ASV3334
#> 
#> ------checking: ASV10784
#> 
#> ------checking: ASV11314
#>   |                                                          |===============                                   |  30%
#> 
#> ####processing: ASV534 #####4
#> 
#> ---hits: ASV313
#> ---hits: ASV89
#> ---hits: ASV462
#> ---hits: ASV1101
#> ---hits: ASV1557
#> ---hits: ASV1058
#> ---hits: ASV1459
#> ---hits: ASV1548
#> ---hits: ASV926
#> ---hits: ASV10594
#> 
#> ---potential parent: ASV89
#> ---potential parent: ASV1058
#> ---potential parent: ASV926
#> ---potential parent: ASV1459
#> ---potential parent: ASV1059
#> ---potential parent: ASV313
#> ---potential parent: ASV462
#> ---potential parent: ASV1101
#> ---potential parent: ASV15574
#> 
#> ------checking: ASV894
#> 
#> ------relative cooccurence: 14
#>  which is sufficient!4
#> 
#> ------min avg abundance: 4.334558823529414
#>  which is OK!4
#> 
#> SETTING ASV534 to be an ERROR of ASV89
#> 4
#> 
#> ------checking: ASV10584
#> 
#> ------checking: ASV9264
#> 
#> ------checking: ASV14594
#> 
#> ------checking: ASV10594
#> 
#> ------checking: ASV3134
#> 
#> ------checking: ASV4624
#> 
#> ------checking: ASV11014
#> 
#> ------checking: ASV15574
#> 
#> ####processing: ASV950 #####4
#> 
#> ---hits: ASV906
#> ---hits: ASV1278
#> ---hits: ASV1074
#> ---hits: ASV828
#> ---hits: ASV1007
#> ---hits: ASV1359
#> ---hits: ASV1192
#> ---hits: ASV1319
#> ---hits: ASV1224
#> ---hits: ASV7564
#> 
#> ---potential parent: ASV756
#> ---potential parent: ASV1192
#> ---potential parent: ASV1074
#> ---potential parent: ASV906
#> ---potential parent: ASV1007
#> ---potential parent: ASV1359
#> ---potential parent: ASV1278
#> ---potential parent: ASV13194
#> 
#> ------checking: ASV7564
#> 
#> ------relative cooccurence: 0.54
#> 
#> ------checking: ASV11924
#> 
#> ------relative cooccurence: 14
#>  which is sufficient!4
#> 
#> ------min avg abundance: 0.02298850574712644
#> 
#> ------checking: ASV10744
#> 
#> ------relative cooccurence: 14
#>  which is sufficient!4
#> 
#> ------min avg abundance: 0.04022988505747134
#> 
#> ------checking: ASV9064
#> 
#> ------relative cooccurence: 14
#>  which is sufficient!4
#> 
#> ------min avg abundance: 0.01149425287356324
#> 
#> ------checking: ASV10074
#> 
#> ------relative cooccurence: 0.54
#> 
#> ------checking: ASV13594
#> 
#> ------relative cooccurence: 04
#> 
#> ------checking: ASV12784
#> 
#> ------relative cooccurence: 04
#> 
#> ------checking: ASV13194
#> 
#> ------relative cooccurence: 0.54
#> 
#> No parent found!
#> 4
#> 
#> ####processing: ASV428 #####4
#> 
#> ---hits: ASV756
#> ---hits: ASV552
#> ---hits: ASV904
#> ---hits: ASV1468
#> ---hits: ASV1035
#> ---hits: ASV1278
#> ---hits: ASV1074
#> ---hits: ASV906
#> ---hits: ASV1007
#> ---hits: ASV13594
#> 
#> ---potential parent: ASV756
#> ---potential parent: ASV1074
#> ---potential parent: ASV906
#> ---potential parent: ASV1007
#> ---potential parent: ASV1359
#> ---potential parent: ASV904
#> ---potential parent: ASV12784
#> 
#> ------checking: ASV7564
#> 
#> ------relative cooccurence: 04
#> 
#> ------checking: ASV10744
#> 
#> ------relative cooccurence: 0.54
#> 
#> ------checking: ASV9064
#> 
#> ------relative cooccurence: 04
#> 
#> ------checking: ASV10074
#> 
#> ------relative cooccurence: 0.54
#> 
#> ------checking: ASV13594
#> 
#> ------relative cooccurence: 0.54
#> 
#> ------checking: ASV9044
#> 
#> ------relative cooccurence: 04
#> 
#> ------checking: ASV12784
#> 
#> ------relative cooccurence: 0.54
#> 
#> No parent found!
#> 4
#> 
#> ####processing: ASV915 #####4
#> 
#> ---hits: ASV201
#> ---hits: ASV43
#> ---hits: ASV542
#> ---hits: ASV1387
#> ---hits: ASV507
#> ---hits: ASV999
#> ---hits: ASV1117
#> ---hits: ASV7664
#> 
#> ---potential parent: ASV43
#> ---potential parent: ASV201
#> ---potential parent: ASV507
#> ---potential parent: ASV542
#> ---potential parent: ASV1387
#> ---potential parent: ASV7664
#> 
#> ------checking: ASV434
#> 
#> ------relative cooccurence: 14
#>  which is sufficient!4
#> 
#> ------min avg abundance: 7.754
#>  which is OK!4
#> 
#> SETTING ASV915 to be an ERROR of ASV43
#> 4
#> 
#> ------checking: ASV2014
#> 
#> ------checking: ASV5074
#> 
#> ------checking: ASV5424
#> 
#> ------checking: ASV13874
#> 
#> ------checking: ASV7664
#>   |                                                          |===============                                   |  31%
#> 
#> ####processing: ASV346 #####4
#> 
#> ---hits: ASV592
#> ---hits: ASV1261
#> ---hits: ASV618
#> ---hits: ASV953
#> ---hits: ASV906
#> ---hits: ASV38
#> ---hits: ASV466
#> ---hits: ASV989
#> ---hits: ASV1388
#> ---hits: ASV14684
#> 
#> ---potential parent: ASV38
#> ---potential parent: ASV618
#> ---potential parent: ASV953
#> ---potential parent: ASV989
#> ---potential parent: ASV906
#> ---potential parent: ASV592
#> ---potential parent: ASV1261
#> ---potential parent: ASV4664
#> 
#> ------checking: ASV384
#> 
#> ------relative cooccurence: 0.54
#> 
#> ------checking: ASV6184
#> 
#> ------relative cooccurence: 14
#>  which is sufficient!4
#> 
#> ------min avg abundance: 1.070967741935484
#>  which is OK!4
#> 
#> SETTING ASV346 to be an ERROR of ASV618
#> 4
#> 
#> ------checking: ASV9534
#> 
#> ------checking: ASV9894
#> 
#> ------checking: ASV9064
#> 
#> ------checking: ASV5924
#> 
#> ------checking: ASV12614
#> 
#> ------checking: ASV4664
#> 
#> ####processing: ASV1146 #####4
#> 
#> ---hits: ASV1609
#> ---hits: ASV4544
#> 
#> ---potential parent: 4
#> 
#> No parent found!
#> 4
#>   |                                                          |================                                  |  31%
#> 
#> ####processing: ASV942 #####4
#> 
#> ---hits: 4
#> 
#> ---potential parent: 4
#> 
#> No parent found!
#> 4
#> 
#> ####processing: ASV1164 #####4
#> 
#> ---hits: 4
#> 
#> ---potential parent: 4
#> 
#> No parent found!
#> 4
#>   |                                                          |================                                  |  32%
#> 
#> ####processing: ASV1268 #####4
#> 
#> ---hits: ASV244
#> ---hits: ASV139
#> ---hits: ASV766
#> ---hits: ASV744
#> ---hits: ASV461
#> ---hits: ASV13874
#> 
#> ---potential parent: ASV139
#> ---potential parent: ASV244
#> ---potential parent: ASV1387
#> ---potential parent: ASV766
#> ---potential parent: ASV7444
#> 
#> ------checking: ASV1394
#> 
#> ------relative cooccurence: 14
#>  which is sufficient!4
#> 
#> ------min avg abundance: 21.22471910112364
#>  which is OK!4
#> 
#> SETTING ASV1268 to be an ERROR of ASV139
#> 4
#> 
#> ------checking: ASV2444
#> 
#> ------checking: ASV13874
#> 
#> ------checking: ASV7664
#> 
#> ------checking: ASV7444
#> 
#> ####processing: ASV1152 #####4
#> 
#> ---hits: ASV839
#> ---hits: ASV272
#> ---hits: ASV339
#> ---hits: ASV81
#> ---hits: ASV383
#> ---hits: ASV959
#> ---hits: ASV1263
#> ---hits: ASV11444
#> 
#> ---potential parent: ASV81
#> ---potential parent: ASV3394
#> 
#> ------checking: ASV814
#> 
#> ------relative cooccurence: 04
#> 
#> ------checking: ASV3394
#> 
#> ------relative cooccurence: 04
#> 
#> No parent found!
#> 4
#> 
#> ####processing: ASV1387 #####4
#> 
#> ---hits: ASV201
#> ---hits: ASV43
#> ---hits: ASV507
#> ---hits: ASV542
#> ---hits: ASV915
#> ---hits: ASV999
#> ---hits: ASV1268
#> ---hits: ASV11174
#> 
#> ---potential parent: ASV43
#> ---potential parent: ASV201
#> ---potential parent: ASV507
#> ---potential parent: ASV542
#> ---potential parent: ASV915
#> ---potential parent: ASV12684
#> 
#> ------checking: ASV434
#> 
#> ------relative cooccurence: 14
#>  which is sufficient!4
#> 
#> ------min avg abundance: 314
#>  which is OK!4
#> 
#> SETTING ASV1387 to be an ERROR of ASV43
#> 4
#> 
#> ------checking: ASV2014
#> 
#> ------checking: ASV5074
#> 
#> ------checking: ASV5424
#> 
#> ------checking: ASV9154
#> 
#> ------checking: ASV12684
#> 
#> ####processing: ASV592 #####4
#> 
#> ---hits: ASV346
#> ---hits: ASV1261
#> ---hits: ASV953
#> ---hits: ASV618
#> ---hits: ASV906
#> ---hits: ASV38
#> ---hits: ASV466
#> ---hits: ASV989
#> ---hits: ASV1468
#> ---hits: ASV13884
#> 
#> ---potential parent: ASV38
#> ---potential parent: ASV618
#> ---potential parent: ASV953
#> ---potential parent: ASV989
#> ---potential parent: ASV906
#> ---potential parent: ASV346
#> ---potential parent: ASV1261
#> ---potential parent: ASV4664
#> 
#> ------checking: ASV384
#> 
#> ------relative cooccurence: 0.54
#> 
#> ------checking: ASV6184
#> 
#> ------relative cooccurence: 14
#>  which is sufficient!4
#> 
#> ------min avg abundance: 2.243243243243244
#>  which is OK!4
#> 
#> SETTING ASV592 to be an ERROR of ASV618
#> 4
#> 
#> ------checking: ASV9534
#> 
#> ------checking: ASV9894
#> 
#> ------checking: ASV9064
#> 
#> ------checking: ASV3464
#> 
#> ------checking: ASV12614
#> 
#> ------checking: ASV4664
#>   |                                                          |================                                  |  33%
#> 
#> ####processing: ASV737 #####4
#> 
#> ---hits: ASV1369
#> ---hits: ASV715
#> ---hits: ASV981
#> ---hits: ASV963
#> ---hits: ASV172
#> ---hits: ASV313
#> ---hits: ASV4624
#> 
#> ---potential parent: ASV172
#> ---potential parent: ASV1369
#> ---potential parent: ASV313
#> ---potential parent: ASV462
#> ---potential parent: ASV963
#> ---potential parent: ASV7154
#> 
#> ------checking: ASV1724
#> 
#> ------relative cooccurence: 0.54
#> 
#> ------checking: ASV13694
#> 
#> ------relative cooccurence: 14
#>  which is sufficient!4
#> 
#> ------min avg abundance: 0.244
#> 
#> ------checking: ASV3134
#> 
#> ------relative cooccurence: 0.54
#> 
#> ------checking: ASV4624
#> 
#> ------relative cooccurence: 0.54
#> 
#> ------checking: ASV9634
#> 
#> ------relative cooccurence: 0.54
#> 
#> ------checking: ASV7154
#> 
#> ------relative cooccurence: 0.54
#> 
#> No parent found!
#> 4
#> 
#> ####processing: ASV1458 #####4
#> 
#> ---hits: ASV1030
#> ---hits: ASV1303
#> ---hits: ASV1332
#> ---hits: ASV693
#> ---hits: ASV1562
#> ---hits: ASV772
#> ---hits: ASV860
#> ---hits: ASV1230
#> ---hits: ASV504
#> ---hits: ASV14094
#> 
#> ---potential parent: ASV1409
#> ---potential parent: ASV1230
#> ---potential parent: ASV772
#> ---potential parent: ASV693
#> ---potential parent: ASV13324
#> 
#> ------checking: ASV14094
#> 
#> ------relative cooccurence: 0.54
#> 
#> ------checking: ASV12304
#> 
#> ------relative cooccurence: 04
#> 
#> ------checking: ASV7724
#> 
#> ------relative cooccurence: 0.54
#> 
#> ------checking: ASV6934
#> 
#> ------relative cooccurence: 0.54
#> 
#> ------checking: ASV13324
#> 
#> ------relative cooccurence: 04
#> 
#> No parent found!
#> 4
#>   |                                                          |=================                                 |  33%
#> 
#> ####processing: ASV772 #####4
#> 
#> ---hits: ASV1409
#> ---hits: ASV504
#> ---hits: ASV860
#> ---hits: ASV832
#> ---hits: ASV693
#> ---hits: ASV1458
#> ---hits: ASV1200
#> ---hits: ASV1332
#> ---hits: ASV1303
#> ---hits: ASV10304
#> 
#> ---potential parent: ASV1409
#> ---potential parent: ASV1458
#> ---potential parent: ASV693
#> ---potential parent: ASV13324
#> 
#> ------checking: ASV14094
#> 
#> ------relative cooccurence: 0.54
#> 
#> ------checking: ASV14584
#> 
#> ------relative cooccurence: 0.54
#> 
#> ------checking: ASV6934
#> 
#> ------relative cooccurence: 0.54
#> 
#> ------checking: ASV13324
#> 
#> ------relative cooccurence: 0.54
#> 
#> No parent found!
#> 4
#> 
#> ####processing: ASV523 #####4
#> 
#> ---hits: ASV9874
#> 
#> ---potential parent: 4
#> 
#> No parent found!
#> 4
#>   |                                                          |=================                                 |  34%
#> 
#> ####processing: ASV766 #####4
#> 
#> ---hits: ASV244
#> ---hits: ASV139
#> ---hits: ASV1268
#> ---hits: ASV744
#> ---hits: ASV461
#> ---hits: ASV9154
#> 
#> ---potential parent: ASV139
#> ---potential parent: ASV244
#> ---potential parent: ASV915
#> ---potential parent: ASV1268
#> ---potential parent: ASV7444
#> 
#> ------checking: ASV1394
#> 
#> ------relative cooccurence: 14
#>  which is sufficient!4
#> 
#> ------min avg abundance: 4.333333333333334
#>  which is OK!4
#> 
#> SETTING ASV766 to be an ERROR of ASV139
#> 4
#> 
#> ------checking: ASV2444
#> 
#> ------checking: ASV9154
#> 
#> ------checking: ASV12684
#> 
#> ------checking: ASV7444
#> 
#> ####processing: ASV55 #####4
#> 
#> ---hits: ASV500
#> ---hits: ASV254
#> ---hits: ASV15264
#> 
#> ---potential parent: ASV5004
#> 
#> ------checking: ASV5004
#> 
#> ------relative cooccurence: 04
#> 
#> No parent found!
#> 4
#> 
#> ####processing: ASV1236 #####4
#> 
#> ---hits: ASV1300
#> ---hits: ASV15884
#> 
#> ---potential parent: ASV13004
#> 
#> ------checking: ASV13004
#> 
#> ------relative cooccurence: 0.54
#> 
#> No parent found!
#> 4
#> 
#> ####processing: ASV749 #####4
#> 
#> ---hits: ASV662
#> ---hits: ASV1052
#> ---hits: ASV724
#> ---hits: ASV420
#> ---hits: ASV505
#> ---hits: ASV509
#> ---hits: ASV790
#> ---hits: ASV1427
#> ---hits: ASV1410
#> ---hits: ASV8594
#> 
#> ---potential parent: ASV662
#> ---potential parent: ASV505
#> ---potential parent: ASV509
#> ---potential parent: ASV724
#> ---potential parent: ASV859
#> ---potential parent: ASV1052
#> ---potential parent: ASV420
#> ---potential parent: ASV14274
#> 
#> ------checking: ASV6624
#> 
#> ------relative cooccurence: 0.54
#> 
#> ------checking: ASV5054
#> 
#> ------relative cooccurence: 0.54
#> 
#> ------checking: ASV5094
#> 
#> ------relative cooccurence: 0.54
#> 
#> ------checking: ASV7244
#> 
#> ------relative cooccurence: 0.54
#> 
#> ------checking: ASV8594
#> 
#> ------relative cooccurence: 0.54
#> 
#> ------checking: ASV10524
#> 
#> ------relative cooccurence: 04
#> 
#> ------checking: ASV4204
#> 
#> ------relative cooccurence: 04
#> 
#> ------checking: ASV14274
#> 
#> ------relative cooccurence: 04
#> 
#> No parent found!
#> 4
#>   |                                                          |=================                                 |  35%
#> 
#> ####processing: ASV854 #####4
#> 
#> ---hits: ASV255
#> ---hits: ASV570
#> ---hits: ASV14594
#> 
#> ---potential parent: ASV1459
#> ---potential parent: ASV2554
#> 
#> ------checking: ASV14594
#> 
#> ------relative cooccurence: 04
#> 
#> ------checking: ASV2554
#> 
#> ------relative cooccurence: 0.54
#> 
#> No parent found!
#> 4
#> 
#> ####processing: ASV743 #####4
#> 
#> ---hits: 4
#> 
#> ---potential parent: 4
#> 
#> No parent found!
#> 4
#>   |                                                          |==================                                |  35%
#> 
#> ####processing: ASV975 #####4
#> 
#> ---hits: ASV38
#> ---hits: ASV466
#> ---hits: ASV989
#> ---hits: ASV816
#> ---hits: ASV853
#> ---hits: ASV1526
#> ---hits: ASV1107
#> ---hits: ASV254
#> ---hits: ASV500
#> ---hits: ASV5924
#> 
#> ---potential parent: ASV38
#> ---potential parent: ASV989
#> ---potential parent: ASV816
#> ---potential parent: ASV500
#> ---potential parent: ASV592
#> ---potential parent: ASV4664
#> 
#> ------checking: ASV384
#> 
#> ------relative cooccurence: 14
#>  which is sufficient!4
#> 
#> ------min avg abundance: 7.659090909090914
#>  which is OK!4
#> 
#> SETTING ASV975 to be an ERROR of ASV38
#> 4
#> 
#> ------checking: ASV9894
#> 
#> ------checking: ASV8164
#> 
#> ------checking: ASV5004
#> 
#> ------checking: ASV5924
#> 
#> ------checking: ASV4664
#> 
#> ####processing: ASV111 #####4
#> 
#> ---hits: ASV2534
#> 
#> ---potential parent: 4
#> 
#> No parent found!
#> 4
#>   |                                                          |==================                                |  36%
#> 
#> ####processing: ASV1532 #####4
#> 
#> ---hits: ASV60
#> ---hits: ASV113
#> ---hits: ASV415
#> ---hits: ASV984
#> ---hits: ASV796
#> ---hits: ASV1380
#> ---hits: ASV1020
#> ---hits: ASV1011
#> ---hits: ASV1036
#> ---hits: ASV10664
#> 
#> ---potential parent: ASV113
#> ---potential parent: ASV60
#> ---potential parent: ASV10364
#> 
#> ------checking: ASV1134
#> 
#> ------relative cooccurence: 0.54
#> 
#> ------checking: ASV604
#> 
#> ------relative cooccurence: 0.54
#> 
#> ------checking: ASV10364
#> 
#> ------relative cooccurence: 0.54
#> 
#> No parent found!
#> 4
#> 
#> ####processing: ASV1493 #####4
#> 
#> ---hits: ASV1022
#> ---hits: ASV1058
#> ---hits: ASV1059
#> ---hits: ASV717
#> ---hits: ASV926
#> ---hits: ASV344
#> ---hits: ASV722
#> ---hits: ASV14594
#> 
#> ---potential parent: ASV717
#> ---potential parent: ASV1022
#> ---potential parent: ASV1058
#> ---potential parent: ASV344
#> ---potential parent: ASV926
#> ---potential parent: ASV1459
#> ---potential parent: ASV10594
#> 
#> ------checking: ASV7174
#> 
#> ------relative cooccurence: 14
#>  which is sufficient!4
#> 
#> ------min avg abundance: 0.5405405405405414
#> 
#> ------checking: ASV10224
#> 
#> ------relative cooccurence: 14
#>  which is sufficient!4
#> 
#> ------min avg abundance: 0.6216216216216224
#> 
#> ------checking: ASV10584
#> 
#> ------relative cooccurence: 0.54
#> 
#> ------checking: ASV3444
#> 
#> ------relative cooccurence: 14
#>  which is sufficient!4
#> 
#> ------min avg abundance: 2.972972972972974
#>  which is OK!4
#> 
#> SETTING ASV1493 to be an ERROR of ASV344
#> 4
#> 
#> ------checking: ASV9264
#> 
#> ------checking: ASV14594
#> 
#> ------checking: ASV10594
#> 
#> ####processing: ASV1131 #####4
#> 
#> ---hits: ASV178
#> ---hits: ASV19
#> ---hits: ASV594
#> ---hits: ASV333
#> ---hits: ASV1467
#> ---hits: ASV1010
#> ---hits: ASV1078
#> ---hits: ASV858
#> ---hits: ASV986
#> ---hits: ASV9794
#> 
#> ---potential parent: ASV19
#> ---potential parent: ASV178
#> ---potential parent: ASV333
#> ---potential parent: ASV1078
#> ---potential parent: ASV5944
#> 
#> ------checking: ASV194
#> 
#> ------relative cooccurence: 14
#>  which is sufficient!4
#> 
#> ------min avg abundance: 514
#>  which is OK!4
#> 
#> SETTING ASV1131 to be an ERROR of ASV19
#> 4
#> 
#> ------checking: ASV1784
#> 
#> ------checking: ASV3334
#> 
#> ------checking: ASV10784
#> 
#> ------checking: ASV5944
#> 
#> ####processing: ASV1101 #####4
#> 
#> ---hits: ASV89
#> ---hits: ASV462
#> ---hits: ASV534
#> ---hits: ASV313
#> ---hits: ASV1557
#> ---hits: ASV1058
#> ---hits: ASV1459
#> ---hits: ASV15484
#> 
#> ---potential parent: ASV89
#> ---potential parent: ASV1058
#> ---potential parent: ASV1459
#> ---potential parent: ASV313
#> ---potential parent: ASV462
#> ---potential parent: ASV534
#> ---potential parent: ASV15574
#> 
#> ------checking: ASV894
#> 
#> ------relative cooccurence: 14
#>  which is sufficient!4
#> 
#> ------min avg abundance: 39.34
#>  which is OK!4
#> 
#> SETTING ASV1101 to be an ERROR of ASV89
#> 4
#> 
#> ------checking: ASV10584
#> 
#> ------checking: ASV14594
#> 
#> ------checking: ASV3134
#> 
#> ------checking: ASV4624
#> 
#> ------checking: ASV5344
#> 
#> ------checking: ASV15574
#> 
#> ####processing: ASV474 #####4
#> 
#> ---hits: ASV170
#> ---hits: ASV18
#> ---hits: ASV831
#> ---hits: ASV8
#> ---hits: ASV94
#> ---hits: ASV209
#> ---hits: ASV666
#> ---hits: ASV261
#> ---hits: ASV493
#> ---hits: ASV3484
#> 
#> ---potential parent: ASV8
#> ---potential parent: ASV18
#> ---potential parent: ASV170
#> ---potential parent: ASV666
#> ---potential parent: ASV94
#> ---potential parent: ASV209
#> ---potential parent: ASV261
#> ---potential parent: ASV348
#> ---potential parent: ASV493
#> ---potential parent: ASV8314
#> 
#> ------checking: ASV84
#> 
#> ------relative cooccurence: 14
#>  which is sufficient!4
#> 
#> ------min avg abundance: 72.254
#>  which is OK!4
#> 
#> SETTING ASV474 to be an ERROR of ASV8
#> 4
#> 
#> ------checking: ASV184
#> 
#> ------checking: ASV1704
#> 
#> ------checking: ASV6664
#> 
#> ------checking: ASV944
#> 
#> ------checking: ASV2094
#> 
#> ------checking: ASV2614
#> 
#> ------checking: ASV3484
#> 
#> ------checking: ASV4934
#> 
#> ------checking: ASV8314
#>   |                                                          |==================                                |  37%
#> 
#> ####processing: ASV744 #####4
#> 
#> ---hits: ASV766
#> ---hits: ASV244
#> ---hits: ASV139
#> ---hits: ASV1268
#> ---hits: ASV4614
#> 
#> ---potential parent: ASV139
#> ---potential parent: ASV244
#> ---potential parent: ASV1268
#> ---potential parent: ASV7664
#> 
#> ------checking: ASV1394
#> 
#> ------relative cooccurence: 14
#>  which is sufficient!4
#> 
#> ------min avg abundance: 1.368421052631584
#>  which is OK!4
#> 
#> SETTING ASV744 to be an ERROR of ASV139
#> 4
#> 
#> ------checking: ASV2444
#> 
#> ------checking: ASV12684
#> 
#> ------checking: ASV7664
#> 
#> ####processing: ASV1396 #####4
#> 
#> ---hits: ASV1311
#> ---hits: ASV546
#> ---hits: ASV1233
#> ---hits: ASV727
#> ---hits: ASV927
#> ---hits: ASV569
#> ---hits: ASV6044
#> 
#> ---potential parent: ASV727
#> ---potential parent: ASV546
#> ---potential parent: ASV569
#> ---potential parent: ASV9274
#> 
#> ------checking: ASV7274
#> 
#> ------relative cooccurence: 04
#> 
#> ------checking: ASV5464
#> 
#> ------relative cooccurence: 04
#> 
#> ------checking: ASV5694
#> 
#> ------relative cooccurence: 0.54
#> 
#> ------checking: ASV9274
#> 
#> ------relative cooccurence: 04
#> 
#> No parent found!
#> 4
#>   |                                                          |===================                               |  37%
#> 
#> ####processing: ASV758 #####4
#> 
#> ---hits: ASV505
#> ---hits: ASV509
#> ---hits: ASV724
#> ---hits: ASV420
#> ---hits: ASV790
#> ---hits: ASV662
#> ---hits: ASV1052
#> ---hits: ASV884
#> ---hits: ASV1410
#> ---hits: ASV14274
#> 
#> ---potential parent: ASV662
#> ---potential parent: ASV505
#> ---potential parent: ASV509
#> ---potential parent: ASV724
#> ---potential parent: ASV1052
#> ---potential parent: ASV420
#> ---potential parent: ASV884
#> ---potential parent: ASV14274
#> 
#> ------checking: ASV6624
#> 
#> ------relative cooccurence: 14
#>  which is sufficient!4
#> 
#> ------min avg abundance: 0.1666666666666674
#> 
#> ------checking: ASV5054
#> 
#> ------relative cooccurence: 14
#>  which is sufficient!4
#> 
#> ------min avg abundance: 0.2222222222222224
#> 
#> ------checking: ASV5094
#> 
#> ------relative cooccurence: 14
#>  which is sufficient!4
#> 
#> ------min avg abundance: 0.8888888888888894
#> 
#> ------checking: ASV7244
#> 
#> ------relative cooccurence: 14
#>  which is sufficient!4
#> 
#> ------min avg abundance: 0.3333333333333334
#> 
#> ------checking: ASV10524
#> 
#> ------relative cooccurence: 14
#>  which is sufficient!4
#> 
#> ------min avg abundance: 0.3333333333333334
#> 
#> ------checking: ASV4204
#> 
#> ------relative cooccurence: 14
#>  which is sufficient!4
#> 
#> ------min avg abundance: 0.2222222222222224
#> 
#> ------checking: ASV8844
#> 
#> ------relative cooccurence: 14
#>  which is sufficient!4
#> 
#> ------min avg abundance: 0.3333333333333334
#> 
#> ------checking: ASV14274
#> 
#> ------relative cooccurence: 14
#>  which is sufficient!4
#> 
#> ------min avg abundance: 0.1666666666666674
#> 
#> No parent found!
#> 4
#> 
#> ####processing: ASV784 #####4
#> 
#> ---hits: ASV1315
#> ---hits: ASV1421
#> ---hits: ASV13034
#> 
#> ---potential parent: ASV1421
#> ---potential parent: ASV13154
#> 
#> ------checking: ASV14214
#> 
#> ------relative cooccurence: 04
#> 
#> ------checking: ASV13154
#> 
#> ------relative cooccurence: 14
#>  which is sufficient!4
#> 
#> ------min avg abundance: 0.24
#> 
#> No parent found!
#> 4
#>   |                                                          |===================                               |  38%
#> 
#> ####processing: ASV297 #####4
#> 
#> ---hits: 4
#> 
#> ---potential parent: 4
#> 
#> No parent found!
#> 4
#> 
#> ####processing: ASV91 #####4
#> 
#> ---hits: 4
#> 
#> ---potential parent: 4
#> 
#> No parent found!
#> 4
#> 
#> ####processing: ASV199 #####4
#> 
#> ---hits: ASV543
#> ---hits: ASV685
#> ---hits: ASV692
#> ---hits: ASV1037
#> ---hits: ASV3104
#> 
#> ---potential parent: ASV310
#> ---potential parent: ASV543
#> ---potential parent: ASV6854
#> 
#> ------checking: ASV3104
#> 
#> ------relative cooccurence: 14
#>  which is sufficient!4
#> 
#> ------min avg abundance: 14
#> 
#> ------checking: ASV5434
#> 
#> ------relative cooccurence: 14
#>  which is sufficient!4
#> 
#> ------min avg abundance: 0.4166666666666674
#> 
#> ------checking: ASV6854
#> 
#> ------relative cooccurence: 04
#> 
#> No parent found!
#> 4
#> 
#> ####processing: ASV643 #####4
#> 
#> ---hits: ASV5004
#> 
#> ---potential parent: ASV5004
#> 
#> ------checking: ASV5004
#> 
#> ------relative cooccurence: 04
#> 
#> No parent found!
#> 4
#>   |                                                          |===================                               |  39%
#> 
#> ####processing: ASV787 #####4
#> 
#> ---hits: ASV8754
#> 
#> ---potential parent: 4
#> 
#> No parent found!
#> 4
#> 
#> ####processing: ASV1261 #####4
#> 
#> ---hits: ASV592
#> ---hits: ASV346
#> ---hits: ASV953
#> ---hits: ASV618
#> ---hits: ASV756
#> ---hits: ASV1074
#> ---hits: ASV906
#> ---hits: ASV904
#> ---hits: ASV1468
#> ---hits: ASV12784
#> 
#> ---potential parent: ASV756
#> ---potential parent: ASV618
#> ---potential parent: ASV953
#> ---potential parent: ASV1074
#> ---potential parent: ASV906
#> ---potential parent: ASV904
#> ---potential parent: ASV1278
#> ---potential parent: ASV346
#> ---potential parent: ASV5924
#> 
#> ------checking: ASV7564
#> 
#> ------relative cooccurence: 0.54
#> 
#> ------checking: ASV6184
#> 
#> ------relative cooccurence: 14
#>  which is sufficient!4
#> 
#> ------min avg abundance: 114
#>  which is OK!4
#> 
#> SETTING ASV1261 to be an ERROR of ASV618
#> 4
#> 
#> ------checking: ASV9534
#> 
#> ------checking: ASV10744
#> 
#> ------checking: ASV9064
#> 
#> ------checking: ASV9044
#> 
#> ------checking: ASV12784
#> 
#> ------checking: ASV3464
#> 
#> ------checking: ASV5924
#>   |                                                          |====================                              |  39%
#> 
#> ####processing: ASV1434 #####4
#> 
#> ---hits: 4
#> 
#> ---potential parent: 4
#> 
#> No parent found!
#> 4
#> 
#> ####processing: ASV357 #####4
#> 
#> ---hits: ASV6264
#> 
#> ---potential parent: 4
#> 
#> No parent found!
#> 4
#>   |                                                          |====================                              |  40%
#> 
#> ####processing: ASV927 #####4
#> 
#> ---hits: ASV1311
#> ---hits: ASV1396
#> ---hits: ASV12334
#> 
#> ---potential parent: ASV13964
#> 
#> ------checking: ASV13964
#> 
#> ------relative cooccurence: 04
#> 
#> No parent found!
#> 4
#> 
#> ####processing: ASV248 #####4
#> 
#> ---hits: ASV404
#> ---hits: ASV377
#> ---hits: ASV488
#> ---hits: ASV840
#> ---hits: ASV440
#> ---hits: ASV5454
#> 
#> ---potential parent: 4
#> 
#> No parent found!
#> 4
#> 
#> ####processing: ASV963 #####4
#> 
#> ---hits: ASV172
#> ---hits: ASV981
#> ---hits: ASV715
#> ---hits: ASV737
#> ---hits: ASV13694
#> 
#> ---potential parent: ASV172
#> ---potential parent: ASV1369
#> ---potential parent: ASV737
#> ---potential parent: ASV7154
#> 
#> ------checking: ASV1724
#> 
#> ------relative cooccurence: 14
#>  which is sufficient!4
#> 
#> ------min avg abundance: 54
#>  which is OK!4
#> 
#> SETTING ASV963 to be an ERROR of ASV172
#> 4
#> 
#> ------checking: ASV13694
#> 
#> ------checking: ASV7374
#> 
#> ------checking: ASV7154
#> 
#> ####processing: ASV1212 #####4
#> 
#> ---hits: 4
#> 
#> ---potential parent: 4
#> 
#> No parent found!
#> 4
#>   |                                                          |====================                              |  41%
#> 
#> ####processing: ASV1321 #####4
#> 
#> ---hits: 4
#> 
#> ---potential parent: 4
#> 
#> No parent found!
#> 4
#> 
#> ####processing: ASV880 #####4
#> 
#> ---hits: ASV159
#> ---hits: ASV52
#> ---hits: ASV1276
#> ---hits: ASV10704
#> 
#> ---potential parent: ASV52
#> ---potential parent: ASV1594
#> 
#> ------checking: ASV524
#> 
#> ------relative cooccurence: 14
#>  which is sufficient!4
#> 
#> ------min avg abundance: 10.71428571428574
#>  which is OK!4
#> 
#> SETTING ASV880 to be an ERROR of ASV52
#> 4
#> 
#> ------checking: ASV1594
#>   |                                                          |=====================                             |  41%
#> 
#> ####processing: ASV1052 #####4
#> 
#> ---hits: ASV724
#> ---hits: ASV420
#> ---hits: ASV505
#> ---hits: ASV509
#> ---hits: ASV790
#> ---hits: ASV662
#> ---hits: ASV1410
#> ---hits: ASV1427
#> ---hits: ASV1198
#> ---hits: ASV8594
#> 
#> ---potential parent: ASV662
#> ---potential parent: ASV505
#> ---potential parent: ASV509
#> ---potential parent: ASV724
#> ---potential parent: ASV859
#> ---potential parent: ASV420
#> ---potential parent: ASV1198
#> ---potential parent: ASV14274
#> 
#> ------checking: ASV6624
#> 
#> ------relative cooccurence: 14
#>  which is sufficient!4
#> 
#> ------min avg abundance: 0.1666666666666674
#> 
#> ------checking: ASV5054
#> 
#> ------relative cooccurence: 14
#>  which is sufficient!4
#> 
#> ------min avg abundance: 0.6666666666666674
#> 
#> ------checking: ASV5094
#> 
#> ------relative cooccurence: 14
#>  which is sufficient!4
#> 
#> ------min avg abundance: 1.166666666666674
#>  which is OK!4
#> 
#> SETTING ASV1052 to be an ERROR of ASV509
#> 4
#> 
#> ------checking: ASV7244
#> 
#> ------checking: ASV8594
#> 
#> ------checking: ASV4204
#> 
#> ------checking: ASV11984
#> 
#> ------checking: ASV14274
#> 
#> ####processing: ASV1155 #####4
#> 
#> ---hits: ASV1461
#> ---hits: ASV801
#> ---hits: ASV16324
#> 
#> ---potential parent: ASV8014
#> 
#> ------checking: ASV8014
#> 
#> ------relative cooccurence: 14
#>  which is sufficient!4
#> 
#> ------min avg abundance: 0.1428571428571434
#> 
#> No parent found!
#> 4
#>   |                                                          |=====================                             |  42%
#> 
#> ####processing: ASV685 #####4
#> 
#> ---hits: ASV199
#> ---hits: ASV310
#> ---hits: ASV543
#> ---hits: ASV692
#> ---hits: ASV10374
#> 
#> ---potential parent: ASV310
#> ---potential parent: ASV543
#> ---potential parent: ASV1994
#> 
#> ------checking: ASV3104
#> 
#> ------relative cooccurence: 04
#> 
#> ------checking: ASV5434
#> 
#> ------relative cooccurence: 0.54
#> 
#> ------checking: ASV1994
#> 
#> ------relative cooccurence: 04
#> 
#> No parent found!
#> 4
#> 
#> ####processing: ASV1128 #####4
#> 
#> ---hits: ASV662
#> ---hits: ASV724
#> ---hits: ASV859
#> ---hits: ASV420
#> ---hits: ASV505
#> ---hits: ASV509
#> ---hits: ASV790
#> ---hits: ASV1198
#> ---hits: ASV1427
#> ---hits: ASV10524
#> 
#> ---potential parent: ASV662
#> ---potential parent: ASV505
#> ---potential parent: ASV509
#> ---potential parent: ASV724
#> ---potential parent: ASV859
#> ---potential parent: ASV1052
#> ---potential parent: ASV420
#> ---potential parent: ASV1198
#> ---potential parent: ASV14274
#> 
#> ------checking: ASV6624
#> 
#> ------relative cooccurence: 14
#>  which is sufficient!4
#> 
#> ------min avg abundance: 0.24
#> 
#> ------checking: ASV5054
#> 
#> ------relative cooccurence: 14
#>  which is sufficient!4
#> 
#> ------min avg abundance: 14
#> 
#> ------checking: ASV5094
#> 
#> ------relative cooccurence: 14
#>  which is sufficient!4
#> 
#> ------min avg abundance: 1.44
#>  which is OK!4
#> 
#> SETTING ASV1128 to be an ERROR of ASV509
#> 4
#> 
#> ------checking: ASV7244
#> 
#> ------checking: ASV8594
#> 
#> ------checking: ASV10524
#> 
#> ------checking: ASV4204
#> 
#> ------checking: ASV11984
#> 
#> ------checking: ASV14274
#> 
#> ####processing: ASV81 #####4
#> 
#> ---hits: ASV272
#> ---hits: ASV339
#> ---hits: ASV959
#> ---hits: ASV383
#> ---hits: ASV839
#> ---hits: ASV1152
#> ---hits: ASV1144
#> ---hits: ASV12634
#> 
#> ---potential parent: ASV1152
#> ---potential parent: ASV3394
#> 
#> ------checking: ASV11524
#> 
#> ------relative cooccurence: 04
#> 
#> ------checking: ASV3394
#> 
#> ------relative cooccurence: 14
#>  which is sufficient!4
#> 
#> ------min avg abundance: 0.54
#> 
#> No parent found!
#> 4
#> 
#> ####processing: ASV339 #####4
#> 
#> ---hits: ASV272
#> ---hits: ASV81
#> ---hits: ASV383
#> ---hits: ASV959
#> ---hits: ASV839
#> ---hits: ASV1152
#> ---hits: ASV1263
#> ---hits: ASV11444
#> 
#> ---potential parent: ASV1152
#> ---potential parent: ASV814
#> 
#> ------checking: ASV11524
#> 
#> ------relative cooccurence: 04
#> 
#> ------checking: ASV814
#> 
#> ------relative cooccurence: 14
#>  which is sufficient!4
#> 
#> ------min avg abundance: 0.84
#> 
#> No parent found!
#> 4
#>   |                                                          |=====================                             |  43%
#> 
#> ####processing: ASV358 #####4
#> 
#> ---hits: ASV405
#> ---hits: ASV69
#> ---hits: ASV930
#> ---hits: ASV914
#> ---hits: ASV833
#> ---hits: ASV221
#> ---hits: ASV1420
#> ---hits: ASV673
#> ---hits: ASV1218
#> ---hits: ASV4314
#> 
#> ---potential parent: ASV69
#> ---potential parent: ASV12184
#> 
#> ------checking: ASV694
#> 
#> ------relative cooccurence: 14
#>  which is sufficient!4
#> 
#> ------min avg abundance: 1.54
#>  which is OK!4
#> 
#> SETTING ASV358 to be an ERROR of ASV69
#> 4
#> 
#> ------checking: ASV12184
#> 
#> ####processing: ASV420 #####4
#> 
#> ---hits: ASV1427
#> ---hits: ASV505
#> ---hits: ASV509
#> ---hits: ASV790
#> ---hits: ASV1410
#> ---hits: ASV662
#> ---hits: ASV1052
#> ---hits: ASV724
#> ---hits: ASV1198
#> ---hits: ASV8594
#> 
#> ---potential parent: ASV662
#> ---potential parent: ASV505
#> ---potential parent: ASV509
#> ---potential parent: ASV724
#> ---potential parent: ASV859
#> ---potential parent: ASV1052
#> ---potential parent: ASV1198
#> ---potential parent: ASV14274
#> 
#> ------checking: ASV6624
#> 
#> ------relative cooccurence: 14
#>  which is sufficient!4
#> 
#> ------min avg abundance: 0.254
#> 
#> ------checking: ASV5054
#> 
#> ------relative cooccurence: 14
#>  which is sufficient!4
#> 
#> ------min avg abundance: 14
#> 
#> ------checking: ASV5094
#> 
#> ------relative cooccurence: 14
#>  which is sufficient!4
#> 
#> ------min avg abundance: 1.754
#>  which is OK!4
#> 
#> SETTING ASV420 to be an ERROR of ASV509
#> 4
#> 
#> ------checking: ASV7244
#> 
#> ------checking: ASV8594
#> 
#> ------checking: ASV10524
#> 
#> ------checking: ASV11984
#> 
#> ------checking: ASV14274
#>   |                                                          |======================                            |  43%
#> 
#> ####processing: ASV693 #####4
#> 
#> ---hits: ASV1409
#> ---hits: ASV772
#> ---hits: ASV504
#> ---hits: ASV860
#> ---hits: ASV832
#> ---hits: ASV1458
#> ---hits: ASV1303
#> ---hits: ASV1200
#> ---hits: ASV1332
#> ---hits: ASV10304
#> 
#> ---potential parent: ASV1409
#> ---potential parent: ASV1458
#> ---potential parent: ASV772
#> ---potential parent: ASV13324
#> 
#> ------checking: ASV14094
#> 
#> ------relative cooccurence: 14
#>  which is sufficient!4
#> 
#> ------min avg abundance: 2.84
#>  which is OK!4
#> 
#> SETTING ASV693 to be an ERROR of ASV1409
#> 4
#> 
#> ------checking: ASV14584
#> 
#> ------checking: ASV7724
#> 
#> ------checking: ASV13324
#> 
#> ####processing: ASV715 #####4
#> 
#> ---hits: ASV981
#> ---hits: ASV737
#> ---hits: ASV1369
#> ---hits: ASV963
#> ---hits: ASV1724
#> 
#> ---potential parent: ASV172
#> ---potential parent: ASV1369
#> ---potential parent: ASV737
#> ---potential parent: ASV9634
#> 
#> ------checking: ASV1724
#> 
#> ------relative cooccurence: 0.54
#> 
#> ------checking: ASV13694
#> 
#> ------relative cooccurence: 0.54
#> 
#> ------checking: ASV7374
#> 
#> ------relative cooccurence: 0.54
#> 
#> ------checking: ASV9634
#> 
#> ------relative cooccurence: 0.54
#> 
#> No parent found!
#> 4
#>   |                                                          |======================                            |  44%
#> 
#> ####processing: ASV884 #####4
#> 
#> ---hits: ASV641
#> ---hits: ASV790
#> ---hits: ASV420
#> ---hits: ASV505
#> ---hits: ASV509
#> ---hits: ASV662
#> ---hits: ASV1052
#> ---hits: ASV1410
#> ---hits: ASV1427
#> ---hits: ASV7244
#> 
#> ---potential parent: ASV662
#> ---potential parent: ASV505
#> ---potential parent: ASV509
#> ---potential parent: ASV724
#> ---potential parent: ASV1052
#> ---potential parent: ASV420
#> ---potential parent: ASV641
#> ---potential parent: ASV14274
#> 
#> ------checking: ASV6624
#> 
#> ------relative cooccurence: 14
#>  which is sufficient!4
#> 
#> ------min avg abundance: 0.3333333333333334
#> 
#> ------checking: ASV5054
#> 
#> ------relative cooccurence: 14
#>  which is sufficient!4
#> 
#> ------min avg abundance: 0.6666666666666674
#> 
#> ------checking: ASV5094
#> 
#> ------relative cooccurence: 14
#>  which is sufficient!4
#> 
#> ------min avg abundance: 2.333333333333334
#>  which is OK!4
#> 
#> SETTING ASV884 to be an ERROR of ASV509
#> 4
#> 
#> ------checking: ASV7244
#> 
#> ------checking: ASV10524
#> 
#> ------checking: ASV4204
#> 
#> ------checking: ASV6414
#> 
#> ------checking: ASV14274
#> 
#> ####processing: ASV973 #####4
#> 
#> ---hits: ASV1468
#> ---hits: ASV1388
#> ---hits: ASV1653
#> ---hits: ASV45
#> ---hits: ASV137
#> ---hits: ASV996
#> ---hits: ASV1067
#> ---hits: ASV711
#> ---hits: ASV760
#> ---hits: ASV12614
#> 
#> ---potential parent: ASV12614
#> 
#> ------checking: ASV12614
#> 
#> ------relative cooccurence: 04
#> 
#> No parent found!
#> 4
#> 
#> ####processing: ASV1218 #####4
#> 
#> ---hits: ASV673
#> ---hits: ASV1420
#> ---hits: ASV563
#> ---hits: ASV431
#> ---hits: ASV242
#> ---hits: ASV358
#> ---hits: ASV405
#> ---hits: ASV930
#> ---hits: ASV694
#> 
#> ---potential parent: ASV563
#> ---potential parent: ASV69
#> ---potential parent: ASV3584
#> 
#> ------checking: ASV5634
#> 
#> ------relative cooccurence: 04
#> 
#> ------checking: ASV694
#> 
#> ------relative cooccurence: 04
#> 
#> ------checking: ASV3584
#> 
#> ------relative cooccurence: 04
#> 
#> No parent found!
#> 4
#> 
#> ####processing: ASV1319 #####4
#> 
#> ---hits: ASV828
#> ---hits: ASV906
#> ---hits: ASV1074
#> ---hits: ASV1278
#> ---hits: ASV1192
#> ---hits: ASV1007
#> ---hits: ASV1359
#> ---hits: ASV950
#> ---hits: ASV756
#> ---hits: ASV12244
#> 
#> ---potential parent: ASV756
#> ---potential parent: ASV1192
#> ---potential parent: ASV1074
#> ---potential parent: ASV906
#> ---potential parent: ASV1007
#> ---potential parent: ASV1359
#> ---potential parent: ASV1278
#> ---potential parent: ASV9504
#> 
#> ------checking: ASV7564
#> 
#> ------relative cooccurence: 0.54
#> 
#> ------checking: ASV11924
#> 
#> ------relative cooccurence: 0.54
#> 
#> ------checking: ASV10744
#> 
#> ------relative cooccurence: 14
#>  which is sufficient!4
#> 
#> ------min avg abundance: 1.44
#>  which is OK!4
#> 
#> SETTING ASV1319 to be an ERROR of ASV1074
#> 4
#> 
#> ------checking: ASV9064
#> 
#> ------checking: ASV10074
#> 
#> ------checking: ASV13594
#> 
#> ------checking: ASV12784
#> 
#> ------checking: ASV9504
#>   |                                                          |======================                            |  45%
#> 
#> ####processing: ASV1557 #####4
#> 
#> ---hits: ASV89
#> ---hits: ASV1101
#> ---hits: ASV462
#> ---hits: ASV534
#> ---hits: ASV3134
#> 
#> ---potential parent: ASV89
#> ---potential parent: ASV313
#> ---potential parent: ASV462
#> ---potential parent: ASV534
#> ---potential parent: ASV11014
#> 
#> ------checking: ASV894
#> 
#> ------relative cooccurence: 14
#>  which is sufficient!4
#> 
#> ------min avg abundance: 24
#>  which is OK!4
#> 
#> SETTING ASV1557 to be an ERROR of ASV89
#> 4
#> 
#> ------checking: ASV3134
#> 
#> ------checking: ASV4624
#> 
#> ------checking: ASV5344
#> 
#> ------checking: ASV11014
#> 
#> ####processing: ASV1315 #####4
#> 
#> ---hits: ASV784
#> ---hits: ASV1421
#> ---hits: ASV1303
#> ---hits: ASV12304
#> 
#> ---potential parent: ASV1421
#> ---potential parent: ASV1230
#> ---potential parent: ASV7844
#> 
#> ------checking: ASV14214
#> 
#> ------relative cooccurence: 04
#> 
#> ------checking: ASV12304
#> 
#> ------relative cooccurence: 04
#> 
#> ------checking: ASV7844
#> 
#> ------relative cooccurence: 14
#>  which is sufficient!4
#> 
#> ------min avg abundance: 1.666666666666674
#>  which is OK!4
#> 
#> SETTING ASV1315 to be an ERROR of ASV784
#> 4
#>   |                                                          |=======================                           |  45%
#> 
#> ####processing: ASV576 #####4
#> 
#> ---hits: ASV2
#> ---hits: ASV962
#> ---hits: ASV1027
#> ---hits: ASV1082
#> ---hits: ASV1247
#> ---hits: ASV11024
#> 
#> ---potential parent: ASV2
#> ---potential parent: ASV12474
#> 
#> ------checking: ASV24
#> 
#> ------relative cooccurence: 14
#>  which is sufficient!4
#> 
#> ------min avg abundance: 93.66666666666674
#>  which is OK!4
#> 
#> SETTING ASV576 to be an ERROR of ASV2
#> 4
#> 
#> ------checking: ASV12474
#> 
#> ####processing: ASV641 #####4
#> 
#> ---hits: ASV790
#> ---hits: ASV884
#> ---hits: ASV420
#> ---hits: ASV505
#> ---hits: ASV509
#> ---hits: ASV662
#> ---hits: ASV1052
#> ---hits: ASV1410
#> ---hits: ASV1427
#> ---hits: ASV7244
#> 
#> ---potential parent: ASV662
#> ---potential parent: ASV505
#> ---potential parent: ASV509
#> ---potential parent: ASV724
#> ---potential parent: ASV1052
#> ---potential parent: ASV420
#> ---potential parent: ASV884
#> ---potential parent: ASV14274
#> 
#> ------checking: ASV6624
#> 
#> ------relative cooccurence: 14
#>  which is sufficient!4
#> 
#> ------min avg abundance: 0.3333333333333334
#> 
#> ------checking: ASV5054
#> 
#> ------relative cooccurence: 14
#>  which is sufficient!4
#> 
#> ------min avg abundance: 24
#>  which is OK!4
#> 
#> SETTING ASV641 to be an ERROR of ASV505
#> 4
#> 
#> ------checking: ASV5094
#> 
#> ------checking: ASV7244
#> 
#> ------checking: ASV10524
#> 
#> ------checking: ASV4204
#> 
#> ------checking: ASV8844
#> 
#> ------checking: ASV14274
#> 
#> ####processing: ASV1198 #####4
#> 
#> ---hits: ASV662
#> ---hits: ASV724
#> ---hits: ASV420
#> ---hits: ASV505
#> ---hits: ASV509
#> ---hits: ASV790
#> ---hits: ASV1052
#> ---hits: ASV1410
#> ---hits: ASV1427
#> ---hits: ASV8594
#> 
#> ---potential parent: ASV662
#> ---potential parent: ASV505
#> ---potential parent: ASV509
#> ---potential parent: ASV724
#> ---potential parent: ASV859
#> ---potential parent: ASV1052
#> ---potential parent: ASV420
#> ---potential parent: ASV14274
#> 
#> ------checking: ASV6624
#> 
#> ------relative cooccurence: 0.54
#> 
#> ------checking: ASV5054
#> 
#> ------relative cooccurence: 14
#>  which is sufficient!4
#> 
#> ------min avg abundance: 14
#> 
#> ------checking: ASV5094
#> 
#> ------relative cooccurence: 14
#>  which is sufficient!4
#> 
#> ------min avg abundance: 3.54
#>  which is OK!4
#> 
#> SETTING ASV1198 to be an ERROR of ASV509
#> 4
#> 
#> ------checking: ASV7244
#> 
#> ------checking: ASV8594
#> 
#> ------checking: ASV10524
#> 
#> ------checking: ASV4204
#> 
#> ------checking: ASV14274
#>   |                                                          |=======================                           |  46%
#> 
#> ####processing: ASV29 #####4
#> 
#> ---hits: 4
#> 
#> ---potential parent: 4
#> 
#> No parent found!
#> 4
#> 
#> ####processing: ASV286 #####4
#> 
#> ---hits: 4
#> 
#> ---potential parent: 4
#> 
#> No parent found!
#> 4
#> 
#> ####processing: ASV464 #####4
#> 
#> ---hits: ASV13144
#> 
#> ---potential parent: 4
#> 
#> No parent found!
#> 4
#> 
#> ####processing: ASV466 #####4
#> 
#> ---hits: ASV38
#> ---hits: ASV989
#> ---hits: ASV853
#> ---hits: ASV816
#> ---hits: ASV975
#> ---hits: ASV1526
#> ---hits: ASV254
#> ---hits: ASV500
#> ---hits: ASV592
#> ---hits: ASV11074
#> 
#> ---potential parent: ASV38
#> ---potential parent: ASV989
#> ---potential parent: ASV816
#> ---potential parent: ASV500
#> ---potential parent: ASV592
#> ---potential parent: ASV9754
#> 
#> ------checking: ASV384
#> 
#> ------relative cooccurence: 14
#>  which is sufficient!4
#> 
#> ------min avg abundance: 218.54
#>  which is OK!4
#> 
#> SETTING ASV466 to be an ERROR of ASV38
#> 4
#> 
#> ------checking: ASV9894
#> 
#> ------checking: ASV8164
#> 
#> ------checking: ASV5004
#> 
#> ------checking: ASV5924
#> 
#> ------checking: ASV9754
#>   |                                                          |=======================                           |  47%
#> 
#> ####processing: ASV1036 #####4
#> 
#> ---hits: ASV1020
#> ---hits: ASV1603
#> ---hits: ASV60
#> ---hits: ASV796
#> ---hits: ASV113
#> ---hits: ASV415
#> ---hits: ASV984
#> ---hits: ASV1380
#> ---hits: ASV1011
#> ---hits: ASV10664
#> 
#> ---potential parent: ASV113
#> ---potential parent: ASV604
#> 
#> ------checking: ASV1134
#> 
#> ------relative cooccurence: 14
#>  which is sufficient!4
#> 
#> ------min avg abundance: 24
#>  which is OK!4
#> 
#> SETTING ASV1036 to be an ERROR of ASV113
#> 4
#> 
#> ------checking: ASV604
#> 
#> ####processing: ASV1247 #####4
#> 
#> ---hits: ASV962
#> ---hits: ASV1027
#> ---hits: ASV1082
#> ---hits: ASV576
#> ---hits: ASV24
#> 
#> ---potential parent: ASV2
#> ---potential parent: ASV5764
#> 
#> ------checking: ASV24
#> 
#> ------relative cooccurence: 14
#>  which is sufficient!4
#> 
#> ------min avg abundance: 1304
#>  which is OK!4
#> 
#> SETTING ASV1247 to be an ERROR of ASV2
#> 4
#> 
#> ------checking: ASV5764
#>   |                                                          |========================                          |  47%
#> 
#> ####processing: ASV1300 #####4
#> 
#> ---hits: ASV1236
#> ---hits: ASV1031
#> ---hits: ASV1239
#> ---hits: ASV892
#> ---hits: ASV1510
#> ---hits: ASV1624
#> ---hits: ASV1644
#> ---hits: ASV15884
#> 
#> ---potential parent: ASV1624
#> ---potential parent: ASV1644
#> ---potential parent: ASV12364
#> 
#> ------checking: ASV16244
#> 
#> ------relative cooccurence: 14
#>  which is sufficient!4
#> 
#> ------min avg abundance: 0.54
#> 
#> ------checking: ASV16444
#> 
#> ------relative cooccurence: 0.54
#> 
#> ------checking: ASV12364
#> 
#> ------relative cooccurence: 0.54
#> 
#> No parent found!
#> 4
#> 
#> ####processing: ASV1427 #####4
#> 
#> ---hits: ASV420
#> ---hits: ASV859
#> ---hits: ASV505
#> ---hits: ASV509
#> ---hits: ASV790
#> ---hits: ASV1410
#> ---hits: ASV662
#> ---hits: ASV1052
#> ---hits: ASV724
#> ---hits: ASV11984
#> 
#> ---potential parent: ASV662
#> ---potential parent: ASV505
#> ---potential parent: ASV509
#> ---potential parent: ASV724
#> ---potential parent: ASV859
#> ---potential parent: ASV1052
#> ---potential parent: ASV420
#> ---potential parent: ASV11984
#> 
#> ------checking: ASV6624
#> 
#> ------relative cooccurence: 14
#>  which is sufficient!4
#> 
#> ------min avg abundance: 14
#> 
#> ------checking: ASV5054
#> 
#> ------relative cooccurence: 14
#>  which is sufficient!4
#> 
#> ------min avg abundance: 14
#> 
#> ------checking: ASV5094
#> 
#> ------relative cooccurence: 14
#>  which is sufficient!4
#> 
#> ------min avg abundance: 44
#>  which is OK!4
#> 
#> SETTING ASV1427 to be an ERROR of ASV509
#> 4
#> 
#> ------checking: ASV7244
#> 
#> ------checking: ASV8594
#> 
#> ------checking: ASV10524
#> 
#> ------checking: ASV4204
#> 
#> ------checking: ASV11984
#>   |                                                          |========================                          |  48%
#> 
#> ####processing: ASV1690 #####4
#> 
#> ---hits: 4
#> 
#> ---potential parent: 4
#> 
#> No parent found!
#> 4
#> 
#> ####processing: ASV834 #####4
#> 
#> ---hits: ASV355
#> ---hits: ASV1514
#> 
#> ---potential parent: ASV1514
#> 
#> ------checking: ASV1514
#> 
#> ------relative cooccurence: 0.54
#> 
#> No parent found!
#> 4
#> 
#> ####processing: ASV1176 #####4
#> 
#> ---hits: ASV443
#> ---hits: ASV1111
#> ---hits: ASV368
#> ---hits: ASV292
#> ---hits: ASV566
#> ---hits: ASV961
#> ---hits: ASV984
#> 
#> ---potential parent: ASV292
#> ---potential parent: ASV443
#> ---potential parent: ASV98
#> ---potential parent: ASV368
#> ---potential parent: ASV5664
#> 
#> ------checking: ASV2924
#> 
#> ------relative cooccurence: 14
#>  which is sufficient!4
#> 
#> ------min avg abundance: 104
#>  which is OK!4
#> 
#> SETTING ASV1176 to be an ERROR of ASV292
#> 4
#> 
#> ------checking: ASV4434
#> 
#> ------checking: ASV984
#> 
#> ------checking: ASV3684
#> 
#> ------checking: ASV5664
#> 
#> ####processing: ASV1223 #####4
#> 
#> ---hits: ASV491
#> ---hits: ASV2104
#> 
#> ---potential parent: 4
#> 
#> No parent found!
#> 4
#>   |                                                          |========================                          |  49%
#> 
#> ####processing: ASV1332 #####4
#> 
#> ---hits: ASV1303
#> ---hits: ASV746
#> ---hits: ASV378
#> ---hits: ASV1562
#> ---hits: ASV1200
#> ---hits: ASV1458
#> ---hits: ASV650
#> ---hits: ASV1230
#> ---hits: ASV860
#> ---hits: ASV14094
#> 
#> ---potential parent: ASV378
#> ---potential parent: ASV746
#> ---potential parent: ASV1409
#> ---potential parent: ASV1230
#> ---potential parent: ASV14584
#> 
#> ------checking: ASV3784
#> 
#> ------relative cooccurence: 0.54
#> 
#> ------checking: ASV7464
#> 
#> ------relative cooccurence: 0.54
#> 
#> ------checking: ASV14094
#> 
#> ------relative cooccurence: 04
#> 
#> ------checking: ASV12304
#> 
#> ------relative cooccurence: 0.54
#> 
#> ------checking: ASV14584
#> 
#> ------relative cooccurence: 04
#> 
#> No parent found!
#> 4
#> 
#> ####processing: ASV1365 #####4
#> 
#> ---hits: ASV348
#> ---hits: ASV94
#> ---hits: ASV209
#> ---hits: ASV261
#> ---hits: ASV493
#> ---hits: ASV1182
#> ---hits: ASV170
#> ---hits: ASV26
#> ---hits: ASV474
#> ---hits: ASV184
#> 
#> ---potential parent: ASV18
#> ---potential parent: ASV170
#> ---potential parent: ASV94
#> ---potential parent: ASV209
#> ---potential parent: ASV261
#> ---potential parent: ASV348
#> ---potential parent: ASV493
#> ---potential parent: ASV4744
#> 
#> ------checking: ASV184
#> 
#> ------relative cooccurence: 0.54
#> 
#> ------checking: ASV1704
#> 
#> ------relative cooccurence: 0.54
#> 
#> ------checking: ASV944
#> 
#> ------relative cooccurence: 0.54
#> 
#> ------checking: ASV2094
#> 
#> ------relative cooccurence: 0.54
#> 
#> ------checking: ASV2614
#> 
#> ------relative cooccurence: 0.54
#> 
#> ------checking: ASV3484
#> 
#> ------relative cooccurence: 0.54
#> 
#> ------checking: ASV4934
#> 
#> ------relative cooccurence: 0.54
#> 
#> ------checking: ASV4744
#> 
#> ------relative cooccurence: 0.54
#> 
#> No parent found!
#> 4
#>   |                                                          |=========================                         |  49%
#> 
#> ####processing: ASV78 #####4
#> 
#> ---hits: ASV33
#> ---hits: ASV1313
#> ---hits: ASV941
#> ---hits: ASV16254
#> 
#> ---potential parent: ASV1313
#> ---potential parent: ASV33
#> ---potential parent: ASV1625
#> ---potential parent: ASV9414
#> 
#> ------checking: ASV13134
#> 
#> ------relative cooccurence: 14
#>  which is sufficient!4
#> 
#> ------min avg abundance: 0.01165695253955044
#> 
#> ------checking: ASV334
#> 
#> ------relative cooccurence: 14
#>  which is sufficient!4
#> 
#> ------min avg abundance: 0.003164029975020824
#> 
#> ------checking: ASV16254
#> 
#> ------relative cooccurence: 14
#>  which is sufficient!4
#> 
#> ------min avg abundance: 0.0004995836802664454
#> 
#> ------checking: ASV9414
#> 
#> ------relative cooccurence: 14
#>  which is sufficient!4
#> 
#> ------min avg abundance: 0.0001665278934221484
#> 
#> No parent found!
#> 4
#> 
#> ####processing: ASV101 #####4
#> 
#> ---hits: 4
#> 
#> ---potential parent: 4
#> 
#> No parent found!
#> 4
#>   |                                                          |=========================                         |  50%
#> 
#> ####processing: ASV168 #####4
#> 
#> ---hits: ASV64
#> ---hits: ASV249
#> ---hits: ASV580
#> ---hits: ASV2034
#> 
#> ---potential parent: ASV64
#> ---potential parent: ASV580
#> ---potential parent: ASV203
#> ---potential parent: ASV2494
#> 
#> ------checking: ASV644
#> 
#> ------relative cooccurence: 14
#>  which is sufficient!4
#> 
#> ------min avg abundance: 2.845672575599584
#>  which is OK!4
#> 
#> SETTING ASV168 to be an ERROR of ASV64
#> 4
#> 
#> ------checking: ASV5804
#> 
#> ------checking: ASV2034
#> 
#> ------checking: ASV2494
#> 
#> ####processing: ASV221 #####4
#> 
#> ---hits: ASV833
#> ---hits: ASV69
#> ---hits: ASV405
#> ---hits: ASV358
#> ---hits: ASV914
#> ---hits: ASV930
#> ---hits: ASV14204
#> 
#> ---potential parent: ASV69
#> ---potential parent: ASV358
#> ---potential parent: ASV833
#> ---potential parent: ASV1420
#> ---potential parent: ASV930
#> ---potential parent: ASV914
#> ---potential parent: ASV4054
#> 
#> ------checking: ASV694
#> 
#> ------relative cooccurence: 14
#>  which is sufficient!4
#> 
#> ------min avg abundance: 0.001556824078879094
#> 
#> ------checking: ASV3584
#> 
#> ------relative cooccurence: 14
#>  which is sufficient!4
#> 
#> ------min avg abundance: 0.001037882719252724
#> 
#> ------checking: ASV8334
#> 
#> ------relative cooccurence: 14
#>  which is sufficient!4
#> 
#> ------min avg abundance: 0.1286974571873384
#> 
#> ------checking: ASV14204
#> 
#> ------relative cooccurence: 04
#> 
#> ------checking: ASV9304
#> 
#> ------relative cooccurence: 04
#> 
#> ------checking: ASV9144
#> 
#> ------relative cooccurence: 04
#> 
#> ------checking: ASV4054
#> 
#> ------relative cooccurence: 14
#>  which is sufficient!4
#> 
#> ------min avg abundance: 0.0005189413596263624
#> 
#> No parent found!
#> 4
#> 
#> ####processing: ASV242 #####4
#> 
#> ---hits: ASV1166
#> ---hits: ASV1218
#> ---hits: ASV1420
#> ---hits: ASV3584
#> 
#> ---potential parent: ASV358
#> ---potential parent: ASV1218
#> ---potential parent: ASV1166
#> ---potential parent: ASV14204
#> 
#> ------checking: ASV3584
#> 
#> ------relative cooccurence: 04
#> 
#> ------checking: ASV12184
#> 
#> ------relative cooccurence: 14
#>  which is sufficient!4
#> 
#> ------min avg abundance: 0.0005910165484633574
#> 
#> ------checking: ASV11664
#> 
#> ------relative cooccurence: 04
#> 
#> ------checking: ASV14204
#> 
#> ------relative cooccurence: 14
#>  which is sufficient!4
#> 
#> ------min avg abundance: 0.0413711583924354
#> 
#> No parent found!
#> 4
#> 
#> ####processing: ASV253 #####4
#> 
#> ---hits: ASV1114
#> 
#> ---potential parent: ASV1114
#> 
#> ------checking: ASV1114
#> 
#> ------relative cooccurence: 14
#>  which is sufficient!4
#> 
#> ------min avg abundance: 0.02321204516938524
#> 
#> No parent found!
#> 4
#>   |                                                          |=========================                         |  51%
#> 
#> ####processing: ASV59 #####4
#> 
#> ---hits: 4
#> 
#> ---potential parent: 4
#> 
#> No parent found!
#> 4
#> 
#> ####processing: ASV264 #####4
#> 
#> ---hits: 4
#> 
#> ---potential parent: 4
#> 
#> No parent found!
#> 4
#>   |                                                          |==========================                        |  51%
#> 
#> ####processing: ASV45 #####4
#> 
#> ---hits: ASV137
#> ---hits: ASV711
#> ---hits: ASV996
#> ---hits: ASV760
#> ---hits: ASV1067
#> ---hits: ASV1167
#> ---hits: ASV1484
#> ---hits: ASV940
#> ---hits: ASV1468
#> ---hits: ASV9734
#> 
#> ---potential parent: ASV973
#> ---potential parent: ASV711
#> ---potential parent: ASV1067
#> ---potential parent: ASV137
#> ---potential parent: ASV1468
#> ---potential parent: ASV996
#> ---potential parent: ASV1167
#> ---potential parent: ASV760
#> ---potential parent: ASV940
#> ---potential parent: ASV14844
#> 
#> ------checking: ASV9734
#> 
#> ------relative cooccurence: 04
#> 
#> ------checking: ASV7114
#> 
#> ------relative cooccurence: 14
#>  which is sufficient!4
#> 
#> ------min avg abundance: 0.02511961722488044
#> 
#> ------checking: ASV10674
#> 
#> ------relative cooccurence: 14
#>  which is sufficient!4
#> 
#> ------min avg abundance: 0.01435406698564594
#> 
#> ------checking: ASV1374
#> 
#> ------relative cooccurence: 04
#> 
#> ------checking: ASV14684
#> 
#> ------relative cooccurence: 04
#> 
#> ------checking: ASV9964
#> 
#> ------relative cooccurence: 14
#>  which is sufficient!4
#> 
#> ------min avg abundance: 0.004784688995215314
#> 
#> ------checking: ASV11674
#> 
#> ------relative cooccurence: 04
#> 
#> ------checking: ASV7604
#> 
#> ------relative cooccurence: 14
#>  which is sufficient!4
#> 
#> ------min avg abundance: 0.002392344497607664
#> 
#> ------checking: ASV9404
#> 
#> ------relative cooccurence: 14
#>  which is sufficient!4
#> 
#> ------min avg abundance: 0.001196172248803834
#> 
#> ------checking: ASV14844
#> 
#> ------relative cooccurence: 14
#>  which is sufficient!4
#> 
#> ------min avg abundance: 0.001196172248803834
#> 
#> No parent found!
#> 4
#> 
#> ####processing: ASV454 #####4
#> 
#> ---hits: ASV1609
#> ---hits: ASV11464
#> 
#> ---potential parent: ASV1146
#> ---potential parent: ASV16094
#> 
#> ------checking: ASV11464
#> 
#> ------relative cooccurence: 04
#> 
#> ------checking: ASV16094
#> 
#> ------relative cooccurence: 04
#> 
#> No parent found!
#> 4
#>   |                                                          |==========================                        |  52%
#> 
#> ####processing: ASV489 #####4
#> 
#> ---hits: ASV762
#> ---hits: ASV12344
#> 
#> ---potential parent: ASV762
#> ---potential parent: ASV12344
#> 
#> ------checking: ASV7624
#> 
#> ------relative cooccurence: 14
#>  which is sufficient!4
#> 
#> ------min avg abundance: 0.5025728987993144
#> 
#> ------checking: ASV12344
#> 
#> ------relative cooccurence: 04
#> 
#> No parent found!
#> 4
#> 
#> ####processing: ASV524 #####4
#> 
#> ---hits: 4
#> 
#> ---potential parent: 4
#> 
#> No parent found!
#> 4
#> 
#> ####processing: ASV552 #####4
#> 
#> ---hits: ASV428
#> ---hits: ASV1468
#> ---hits: ASV1653
#> ---hits: ASV1388
#> ---hits: ASV756
#> ---hits: ASV973
#> ---hits: ASV9534
#> 
#> ---potential parent: ASV756
#> ---potential parent: ASV953
#> ---potential parent: ASV428
#> ---potential parent: ASV973
#> ---potential parent: ASV1468
#> ---potential parent: ASV1388
#> ---potential parent: ASV16534
#> 
#> ------checking: ASV7564
#> 
#> ------relative cooccurence: 04
#> 
#> ------checking: ASV9534
#> 
#> ------relative cooccurence: 14
#>  which is sufficient!4
#> 
#> ------min avg abundance: 0.01431492842535794
#> 
#> ------checking: ASV4284
#> 
#> ------relative cooccurence: 04
#> 
#> ------checking: ASV9734
#> 
#> ------relative cooccurence: 04
#> 
#> ------checking: ASV14684
#> 
#> ------relative cooccurence: 04
#> 
#> ------checking: ASV13884
#> 
#> ------relative cooccurence: 04
#> 
#> ------checking: ASV16534
#> 
#> ------relative cooccurence: 04
#> 
#> No parent found!
#> 4
#> 
#> ####processing: ASV533 #####4
#> 
#> ---hits: 4
#> 
#> ---potential parent: 4
#> 
#> No parent found!
#> 4
#>   |                                                          |==========================                        |  53%
#> 
#> ####processing: ASV619 #####4
#> 
#> ---hits: 4
#> 
#> ---potential parent: 4
#> 
#> No parent found!
#> 4
#> 
#> ####processing: ASV504 #####4
#> 
#> ---hits: ASV860
#> ---hits: ASV832
#> ---hits: ASV772
#> ---hits: ASV1409
#> ---hits: ASV693
#> ---hits: ASV1200
#> ---hits: ASV1332
#> ---hits: ASV1458
#> ---hits: ASV650
#> ---hits: ASV13034
#> 
#> ---potential parent: ASV1409
#> ---potential parent: ASV1458
#> ---potential parent: ASV772
#> ---potential parent: ASV693
#> ---potential parent: ASV1332
#> ---potential parent: ASV832
#> ---potential parent: ASV1200
#> ---potential parent: ASV650
#> ---potential parent: ASV860
#> ---potential parent: ASV13034
#> 
#> ------checking: ASV14094
#> 
#> ------relative cooccurence: 04
#> 
#> ------checking: ASV14584
#> 
#> ------relative cooccurence: 04
#> 
#> ------checking: ASV7724
#> 
#> ------relative cooccurence: 04
#> 
#> ------checking: ASV6934
#> 
#> ------relative cooccurence: 04
#> 
#> ------checking: ASV13324
#> 
#> ------relative cooccurence: 04
#> 
#> ------checking: ASV8324
#> 
#> ------relative cooccurence: 14
#>  which is sufficient!4
#> 
#> ------min avg abundance: 0.3217821782178224
#> 
#> ------checking: ASV12004
#> 
#> ------relative cooccurence: 04
#> 
#> ------checking: ASV6504
#> 
#> ------relative cooccurence: 04
#> 
#> ------checking: ASV8604
#> 
#> ------relative cooccurence: 14
#>  which is sufficient!4
#> 
#> ------min avg abundance: 0.002475247524752484
#> 
#> ------checking: ASV13034
#> 
#> ------relative cooccurence: 04
#> 
#> No parent found!
#> 4
#>   |                                                          |===========================                       |  53%
#> 
#> ####processing: ASV26 #####4
#> 
#> ---hits: ASV1182
#> ---hits: ASV348
#> ---hits: ASV1365
#> ---hits: ASV94
#> ---hits: ASV209
#> ---hits: ASV261
#> ---hits: ASV474
#> ---hits: ASV170
#> ---hits: ASV493
#> ---hits: ASV184
#> 
#> ---potential parent: ASV18
#> ---potential parent: ASV170
#> ---potential parent: ASV94
#> ---potential parent: ASV209
#> ---potential parent: ASV261
#> ---potential parent: ASV348
#> ---potential parent: ASV493
#> ---potential parent: ASV474
#> ---potential parent: ASV1365
#> ---potential parent: ASV11824
#> 
#> ------checking: ASV184
#> 
#> ------relative cooccurence: 14
#>  which is sufficient!4
#> 
#> ------min avg abundance: 0.08547008547008554
#> 
#> ------checking: ASV1704
#> 
#> ------relative cooccurence: 14
#>  which is sufficient!4
#> 
#> ------min avg abundance: 0.2962962962962964
#> 
#> ------checking: ASV944
#> 
#> ------relative cooccurence: 14
#>  which is sufficient!4
#> 
#> ------min avg abundance: 0.4358974358974364
#> 
#> ------checking: ASV2094
#> 
#> ------relative cooccurence: 14
#>  which is sufficient!4
#> 
#> ------min avg abundance: 0.1481481481481484
#> 
#> ------checking: ASV2614
#> 
#> ------relative cooccurence: 14
#>  which is sufficient!4
#> 
#> ------min avg abundance: 0.07122507122507124
#> 
#> ------checking: ASV3484
#> 
#> ------relative cooccurence: 14
#>  which is sufficient!4
#> 
#> ------min avg abundance: 0.06837606837606844
#> 
#> ------checking: ASV4934
#> 
#> ------relative cooccurence: 14
#>  which is sufficient!4
#> 
#> ------min avg abundance: 0.01994301994301994
#> 
#> ------checking: ASV4744
#> 
#> ------relative cooccurence: 14
#>  which is sufficient!4
#> 
#> ------min avg abundance: 0.02279202279202284
#> 
#> ------checking: ASV13654
#> 
#> ------relative cooccurence: 04
#> 
#> ------checking: ASV11824
#> 
#> ------relative cooccurence: 14
#>  which is sufficient!4
#> 
#> ------min avg abundance: 0.3247863247863254
#> 
#> No parent found!
#> 4
#> 
#> ####processing: ASV692 #####4
#> 
#> ---hits: ASV543
#> ---hits: ASV199
#> ---hits: ASV310
#> ---hits: ASV685
#> ---hits: ASV10374
#> 
#> ---potential parent: ASV310
#> ---potential parent: ASV543
#> ---potential parent: ASV199
#> ---potential parent: ASV685
#> ---potential parent: ASV10374
#> 
#> ------checking: ASV3104
#> 
#> ------relative cooccurence: 14
#>  which is sufficient!4
#> 
#> ------min avg abundance: 0.03468208092485554
#> 
#> ------checking: ASV5434
#> 
#> ------relative cooccurence: 14
#>  which is sufficient!4
#> 
#> ------min avg abundance: 0.01445086705202314
#> 
#> ------checking: ASV1994
#> 
#> ------relative cooccurence: 14
#>  which is sufficient!4
#> 
#> ------min avg abundance: 0.03468208092485554
#> 
#> ------checking: ASV6854
#> 
#> ------relative cooccurence: 04
#> 
#> ------checking: ASV10374
#> 
#> ------relative cooccurence: 14
#>  which is sufficient!4
#> 
#> ------min avg abundance: 0.06647398843930644
#> 
#> No parent found!
#> 4
#>   |                                                          |===========================                       |  54%
#> 
#> ####processing: ASV47 #####4
#> 
#> ---hits: 4
#> 
#> ---potential parent: 4
#> 
#> No parent found!
#> 4
#> 
#> ####processing: ASV762 #####4
#> 
#> ---hits: ASV489
#> ---hits: ASV12344
#> 
#> ---potential parent: ASV489
#> ---potential parent: ASV12344
#> 
#> ------checking: ASV4894
#> 
#> ------relative cooccurence: 14
#>  which is sufficient!4
#> 
#> ------min avg abundance: 1.989761092150174
#>  which is OK!4
#> 
#> SETTING ASV762 to be an ERROR of ASV489
#> 4
#> 
#> ------checking: ASV12344
#> 
#> ####processing: ASV461 #####4
#> 
#> ---hits: ASV139
#> ---hits: ASV244
#> ---hits: ASV1268
#> ---hits: ASV766
#> ---hits: ASV744
#> ---hits: ASV434
#> 
#> ---potential parent: ASV43
#> ---potential parent: ASV139
#> ---potential parent: ASV244
#> ---potential parent: ASV1268
#> ---potential parent: ASV766
#> ---potential parent: ASV7444
#> 
#> ------checking: ASV434
#> 
#> ------relative cooccurence: 14
#>  which is sufficient!4
#> 
#> ------min avg abundance: 0.120155038759694
#> 
#> ------checking: ASV1394
#> 
#> ------relative cooccurence: 14
#>  which is sufficient!4
#> 
#> ------min avg abundance: 7.321705426356594
#>  which is OK!4
#> 
#> SETTING ASV461 to be an ERROR of ASV139
#> 4
#> 
#> ------checking: ASV2444
#> 
#> ------checking: ASV12684
#> 
#> ------checking: ASV7664
#> 
#> ------checking: ASV7444
#> 
#> ####processing: ASV833 #####4
#> 
#> ---hits: ASV221
#> ---hits: ASV358
#> ---hits: ASV405
#> ---hits: ASV69
#> ---hits: ASV930
#> ---hits: ASV914
#> ---hits: ASV14204
#> 
#> ---potential parent: ASV69
#> ---potential parent: ASV358
#> ---potential parent: ASV221
#> ---potential parent: ASV1420
#> ---potential parent: ASV930
#> ---potential parent: ASV914
#> ---potential parent: ASV4054
#> 
#> ------checking: ASV694
#> 
#> ------relative cooccurence: 14
#>  which is sufficient!4
#> 
#> ------min avg abundance: 0.01209677419354844
#> 
#> ------checking: ASV3584
#> 
#> ------relative cooccurence: 14
#>  which is sufficient!4
#> 
#> ------min avg abundance: 0.008064516129032264
#> 
#> ------checking: ASV2214
#> 
#> ------relative cooccurence: 14
#>  which is sufficient!4
#> 
#> ------min avg abundance: 7.770161290322584
#>  which is OK!4
#> 
#> SETTING ASV833 to be an ERROR of ASV221
#> 4
#> 
#> ------checking: ASV14204
#> 
#> ------checking: ASV9304
#> 
#> ------checking: ASV9144
#> 
#> ------checking: ASV4054
#>   |                                                          |===========================                       |  55%
#> 
#> ####processing: ASV338 #####4
#> 
#> ---hits: 4
#> 
#> ---potential parent: 4
#> 
#> No parent found!
#> 4
#> 
#> ####processing: ASV898 #####4
#> 
#> ---hits: ASV13204
#> 
#> ---potential parent: ASV13204
#> 
#> ------checking: ASV13204
#> 
#> ------relative cooccurence: 14
#>  which is sufficient!4
#> 
#> ------min avg abundance: 0.452736318407964
#> 
#> No parent found!
#> 4
#> 
#> ####processing: ASV946 #####4
#> 
#> ---hits: ASV1069
#> ---hits: ASV16324
#> 
#> ---potential parent: ASV1632
#> ---potential parent: ASV10694
#> 
#> ------checking: ASV16324
#> 
#> ------relative cooccurence: 04
#> 
#> ------checking: ASV10694
#> 
#> ------relative cooccurence: 04
#> 
#> No parent found!
#> 4
#>   |                                                          |============================                      |  55%
#> 
#> ####processing: ASV987 #####4
#> 
#> ---hits: ASV5234
#> 
#> ---potential parent: ASV5234
#> 
#> ------checking: ASV5234
#> 
#> ------relative cooccurence: 14
#>  which is sufficient!4
#> 
#> ------min avg abundance: 0.3295454545454554
#> 
#> No parent found!
#> 4
#> 
#> ####processing: ASV998 #####4
#> 
#> ---hits: ASV4444
#> 
#> ---potential parent: ASV4444
#> 
#> ------checking: ASV4444
#> 
#> ------relative cooccurence: 14
#>  which is sufficient!4
#> 
#> ------min avg abundance: 3.635838150289024
#>  which is OK!4
#> 
#> SETTING ASV998 to be an ERROR of ASV444
#> 4
#>   |                                                          |============================                      |  56%
#> 
#> ####processing: ASV254 #####4
#> 
#> ---hits: ASV500
#> ---hits: ASV816
#> ---hits: ASV989
#> ---hits: ASV55
#> ---hits: ASV1526
#> ---hits: ASV38
#> ---hits: ASV116
#> ---hits: ASV466
#> ---hits: ASV1107
#> ---hits: ASV8534
#> 
#> ---potential parent: ASV38
#> ---potential parent: ASV989
#> ---potential parent: ASV816
#> ---potential parent: ASV500
#> ---potential parent: ASV55
#> ---potential parent: ASV466
#> ---potential parent: ASV1526
#> ---potential parent: ASV853
#> ---potential parent: ASV116
#> ---potential parent: ASV11074
#> 
#> ------checking: ASV384
#> 
#> ------relative cooccurence: 14
#>  which is sufficient!4
#> 
#> ------min avg abundance: 2.601190476190484
#>  which is OK!4
#> 
#> SETTING ASV254 to be an ERROR of ASV38
#> 4
#> 
#> ------checking: ASV9894
#> 
#> ------checking: ASV8164
#> 
#> ------checking: ASV5004
#> 
#> ------checking: ASV554
#> 
#> ------checking: ASV4664
#> 
#> ------checking: ASV15264
#> 
#> ------checking: ASV8534
#> 
#> ------checking: ASV1164
#> 
#> ------checking: ASV11074
#> 
#> ####processing: ASV637 #####4
#> 
#> ---hits: ASV1015
#> ---hits: ASV9034
#> 
#> ---potential parent: ASV903
#> ---potential parent: ASV10154
#> 
#> ------checking: ASV9034
#> 
#> ------relative cooccurence: 14
#>  which is sufficient!4
#> 
#> ------min avg abundance: 0.03012048192771084
#> 
#> ------checking: ASV10154
#> 
#> ------relative cooccurence: 14
#>  which is sufficient!4
#> 
#> ------min avg abundance: 0.006024096385542174
#> 
#> No parent found!
#> 4
#> 
#> ####processing: ASV1070 #####4
#> 
#> ---hits: ASV52
#> ---hits: ASV159
#> ---hits: ASV8804
#> 
#> ---potential parent: ASV52
#> ---potential parent: ASV159
#> ---potential parent: ASV8804
#> 
#> ------checking: ASV524
#> 
#> ------relative cooccurence: 04
#> 
#> ------checking: ASV1594
#> 
#> ------relative cooccurence: 04
#> 
#> ------checking: ASV8804
#> 
#> ------relative cooccurence: 04
#> 
#> No parent found!
#> 4
#> 
#> ####processing: ASV1092 #####4
#> 
#> ---hits: ASV251
#> ---hits: ASV41
#> ---hits: ASV6724
#> 
#> ---potential parent: ASV41
#> ---potential parent: ASV251
#> ---potential parent: ASV6724
#> 
#> ------checking: ASV414
#> 
#> ------relative cooccurence: 14
#>  which is sufficient!4
#> 
#> ------min avg abundance: 0.03424657534246584
#> 
#> ------checking: ASV2514
#> 
#> ------relative cooccurence: 14
#>  which is sufficient!4
#> 
#> ------min avg abundance: 4.260273972602744
#>  which is OK!4
#> 
#> SETTING ASV1092 to be an ERROR of ASV251
#> 4
#> 
#> ------checking: ASV6724
#>   |                                                          |============================                      |  57%
#> 
#> ####processing: ASV1102 #####4
#> 
#> ---hits: ASV2
#> ---hits: ASV5764
#> 
#> ---potential parent: ASV2
#> ---potential parent: ASV5764
#> 
#> ------checking: ASV24
#> 
#> ------relative cooccurence: 14
#>  which is sufficient!4
#> 
#> ------min avg abundance: 0.04861111111111114
#> 
#> ------checking: ASV5764
#> 
#> ------relative cooccurence: 04
#> 
#> No parent found!
#> 4
#> 
#> ####processing: ASV49 #####4
#> 
#> ---hits: 4
#> 
#> ---potential parent: 4
#> 
#> No parent found!
#> 4
#>   |                                                          |=============================                     |  57%
#> 
#> ####processing: ASV1117 #####4
#> 
#> ---hits: ASV999
#> ---hits: ASV43
#> ---hits: ASV507
#> ---hits: ASV201
#> ---hits: ASV542
#> ---hits: ASV1387
#> ---hits: ASV915
#> ---hits: ASV14194
#> 
#> ---potential parent: ASV43
#> ---potential parent: ASV201
#> ---potential parent: ASV507
#> ---potential parent: ASV542
#> ---potential parent: ASV915
#> ---potential parent: ASV1387
#> ---potential parent: ASV999
#> ---potential parent: ASV14194
#> 
#> ------checking: ASV434
#> 
#> ------relative cooccurence: 14
#>  which is sufficient!4
#> 
#> ------min avg abundance: 78.87142857142864
#>  which is OK!4
#> 
#> SETTING ASV1117 to be an ERROR of ASV43
#> 4
#> 
#> ------checking: ASV2014
#> 
#> ------checking: ASV5074
#> 
#> ------checking: ASV5424
#> 
#> ------checking: ASV9154
#> 
#> ------checking: ASV13874
#> 
#> ------checking: ASV9994
#> 
#> ------checking: ASV14194
#> 
#> ####processing: ASV1118 #####4
#> 
#> ---hits: 4
#> 
#> ---potential parent: 4
#> 
#> No parent found!
#> 4
#>   |                                                          |=============================                     |  58%
#> 
#> ####processing: ASV1132 #####4
#> 
#> ---hits: ASV8474
#> 
#> ---potential parent: ASV8474
#> 
#> ------checking: ASV8474
#> 
#> ------relative cooccurence: 14
#>  which is sufficient!4
#> 
#> ------min avg abundance: 0.3382352941176474
#> 
#> No parent found!
#> 4
#> 
#> ####processing: ASV964 #####4
#> 
#> ---hits: 4
#> 
#> ---potential parent: 4
#> 
#> No parent found!
#> 4
#> 
#> ####processing: ASV832 #####4
#> 
#> ---hits: ASV504
#> ---hits: ASV860
#> ---hits: ASV772
#> ---hits: ASV1409
#> ---hits: ASV693
#> ---hits: ASV1200
#> ---hits: ASV1332
#> ---hits: ASV1458
#> ---hits: ASV1303
#> ---hits: ASV3784
#> 
#> ---potential parent: ASV378
#> ---potential parent: ASV1409
#> ---potential parent: ASV1458
#> ---potential parent: ASV772
#> ---potential parent: ASV693
#> ---potential parent: ASV1332
#> ---potential parent: ASV504
#> ---potential parent: ASV1200
#> ---potential parent: ASV860
#> ---potential parent: ASV13034
#> 
#> ------checking: ASV3784
#> 
#> ------relative cooccurence: 04
#> 
#> ------checking: ASV14094
#> 
#> ------relative cooccurence: 04
#> 
#> ------checking: ASV14584
#> 
#> ------relative cooccurence: 04
#> 
#> ------checking: ASV7724
#> 
#> ------relative cooccurence: 04
#> 
#> ------checking: ASV6934
#> 
#> ------relative cooccurence: 04
#> 
#> ------checking: ASV13324
#> 
#> ------relative cooccurence: 04
#> 
#> ------checking: ASV5044
#> 
#> ------relative cooccurence: 14
#>  which is sufficient!4
#> 
#> ------min avg abundance: 3.107692307692314
#>  which is OK!4
#> 
#> SETTING ASV832 to be an ERROR of ASV504
#> 4
#> 
#> ------checking: ASV12004
#> 
#> ------checking: ASV8604
#> 
#> ------checking: ASV13034
#> 
#> ####processing: ASV1166 #####4
#> 
#> ---hits: ASV2424
#> 
#> ---potential parent: ASV2424
#> 
#> ------checking: ASV2424
#> 
#> ------relative cooccurence: 04
#> 
#> No parent found!
#> 4
#>   |                                                          |=============================                     |  59%
#> 
#> ####processing: ASV1182 #####4
#> 
#> ---hits: ASV26
#> ---hits: ASV348
#> ---hits: ASV1365
#> ---hits: ASV94
#> ---hits: ASV474
#> ---hits: ASV209
#> ---hits: ASV261
#> ---hits: ASV170
#> ---hits: ASV493
#> ---hits: ASV184
#> 
#> ---potential parent: ASV18
#> ---potential parent: ASV170
#> ---potential parent: ASV94
#> ---potential parent: ASV209
#> ---potential parent: ASV261
#> ---potential parent: ASV348
#> ---potential parent: ASV493
#> ---potential parent: ASV474
#> ---potential parent: ASV1365
#> ---potential parent: ASV264
#> 
#> ------checking: ASV184
#> 
#> ------relative cooccurence: 14
#>  which is sufficient!4
#> 
#> ------min avg abundance: 0.2631578947368424
#> 
#> ------checking: ASV1704
#> 
#> ------relative cooccurence: 14
#>  which is sufficient!4
#> 
#> ------min avg abundance: 0.9122807017543864
#> 
#> ------checking: ASV944
#> 
#> ------relative cooccurence: 14
#>  which is sufficient!4
#> 
#> ------min avg abundance: 1.342105263157894
#>  which is OK!4
#> 
#> SETTING ASV1182 to be an ERROR of ASV8
#> 4
#> 
#> ------checking: ASV2094
#> 
#> ------checking: ASV2614
#> 
#> ------checking: ASV3484
#> 
#> ------checking: ASV4934
#> 
#> ------checking: ASV4744
#> 
#> ------checking: ASV13654
#> 
#> ------checking: ASV264
#> 
#> ####processing: ASV1314 #####4
#> 
#> ---hits: ASV4644
#> 
#> ---potential parent: ASV4644
#> 
#> ------checking: ASV4644
#> 
#> ------relative cooccurence: 14
#>  which is sufficient!4
#> 
#> ------min avg abundance: 0.02061855670103094
#> 
#> No parent found!
#> 4
#>   |                                                          |==============================                    |  59%
#> 
#> ####processing: ASV969 #####4
#> 
#> ---hits: ASV1014
#> ---hits: ASV13864
#> 
#> ---potential parent: ASV1386
#> ---potential parent: ASV10144
#> 
#> ------checking: ASV13864
#> 
#> ------relative cooccurence: 14
#>  which is sufficient!4
#> 
#> ------min avg abundance: 0.3152173913043484
#> 
#> ------checking: ASV10144
#> 
#> ------relative cooccurence: 14
#>  which is sufficient!4
#> 
#> ------min avg abundance: 0.03260869565217394
#> 
#> No parent found!
#> 4
#> 
#> ####processing: ASV1200 #####4
#> 
#> ---hits: ASV1332
#> ---hits: ASV1562
#> ---hits: ASV1303
#> ---hits: ASV693
#> ---hits: ASV860
#> ---hits: ASV504
#> ---hits: ASV772
#> ---hits: ASV1409
#> ---hits: ASV1458
#> ---hits: ASV7464
#> 
#> ---potential parent: ASV746
#> ---potential parent: ASV1409
#> ---potential parent: ASV1458
#> ---potential parent: ASV772
#> ---potential parent: ASV693
#> ---potential parent: ASV1332
#> ---potential parent: ASV504
#> ---potential parent: ASV860
#> ---potential parent: ASV1303
#> ---potential parent: ASV15624
#> 
#> ------checking: ASV7464
#> 
#> ------relative cooccurence: 14
#>  which is sufficient!4
#> 
#> ------min avg abundance: 0.02173913043478264
#> 
#> ------checking: ASV14094
#> 
#> ------relative cooccurence: 14
#>  which is sufficient!4
#> 
#> ------min avg abundance: 0.1521739130434784
#> 
#> ------checking: ASV14584
#> 
#> ------relative cooccurence: 04
#> 
#> ------checking: ASV7724
#> 
#> ------relative cooccurence: 04
#> 
#> ------checking: ASV6934
#> 
#> ------relative cooccurence: 14
#>  which is sufficient!4
#> 
#> ------min avg abundance: 0.05434782608695654
#> 
#> ------checking: ASV13324
#> 
#> ------relative cooccurence: 04
#> 
#> ------checking: ASV5044
#> 
#> ------relative cooccurence: 04
#> 
#> ------checking: ASV8604
#> 
#> ------relative cooccurence: 04
#> 
#> ------checking: ASV13034
#> 
#> ------relative cooccurence: 04
#> 
#> ------checking: ASV15624
#> 
#> ------relative cooccurence: 04
#> 
#> No parent found!
#> 4
#>   |                                                          |==============================                    |  60%
#> 
#> ####processing: ASV624 #####4
#> 
#> ---hits: ASV12
#> ---hits: ASV67
#> ---hits: ASV107
#> ---hits: ASV1542
#> ---hits: ASV1262
#> ---hits: ASV847
#> ---hits: ASV1054
#> 
#> ---potential parent: ASV12
#> ---potential parent: ASV107
#> ---potential parent: ASV105
#> ---potential parent: ASV847
#> ---potential parent: ASV1542
#> ---potential parent: ASV1262
#> ---potential parent: ASV674
#> 
#> ------checking: ASV124
#> 
#> ------relative cooccurence: 14
#>  which is sufficient!4
#> 
#> ------min avg abundance: 35.0109890109894
#>  which is OK!4
#> 
#> SETTING ASV624 to be an ERROR of ASV12
#> 4
#> 
#> ------checking: ASV1074
#> 
#> ------checking: ASV1054
#> 
#> ------checking: ASV8474
#> 
#> ------checking: ASV15424
#> 
#> ------checking: ASV12624
#> 
#> ------checking: ASV674
#> 
#> ####processing: ASV1320 #####4
#> 
#> ---hits: ASV8984
#> 
#> ---potential parent: ASV8984
#> 
#> ------checking: ASV8984
#> 
#> ------relative cooccurence: 14
#>  which is sufficient!4
#> 
#> ------min avg abundance: 2.208791208791214
#>  which is OK!4
#> 
#> SETTING ASV1320 to be an ERROR of ASV898
#> 4
#> 
#> ####processing: ASV999 #####4
#> 
#> ---hits: ASV507
#> ---hits: ASV542
#> ---hits: ASV43
#> ---hits: ASV201
#> ---hits: ASV1387
#> ---hits: ASV915
#> ---hits: ASV1117
#> ---hits: ASV139
#> ---hits: ASV14194
#> 
#> ---potential parent: ASV43
#> ---potential parent: ASV139
#> ---potential parent: ASV201
#> ---potential parent: ASV507
#> ---potential parent: ASV542
#> ---potential parent: ASV915
#> ---potential parent: ASV1387
#> ---potential parent: ASV1117
#> ---potential parent: ASV14194
#> 
#> ------checking: ASV434
#> 
#> ------relative cooccurence: 14
#>  which is sufficient!4
#> 
#> ------min avg abundance: 0.3690476190476194
#> 
#> ------checking: ASV1394
#> 
#> ------relative cooccurence: 14
#>  which is sufficient!4
#> 
#> ------min avg abundance: 22.48809523809524
#>  which is OK!4
#> 
#> SETTING ASV999 to be an ERROR of ASV139
#> 4
#> 
#> ------checking: ASV2014
#> 
#> ------checking: ASV5074
#> 
#> ------checking: ASV5424
#> 
#> ------checking: ASV9154
#> 
#> ------checking: ASV13874
#> 
#> ------checking: ASV11174
#> 
#> ------checking: ASV14194
#> 
#> ####processing: ASV1485 #####4
#> 
#> ---hits: 4
#> 
#> ---potential parent: 4
#> 
#> No parent found!
#> 4
#>   |                                                          |==============================                    |  61%
#> 
#> ####processing: ASV1461 #####4
#> 
#> ---hits: ASV1155
#> ---hits: ASV801
#> ---hits: ASV16324
#> 
#> ---potential parent: ASV801
#> ---potential parent: ASV1155
#> ---potential parent: ASV16324
#> 
#> ------checking: ASV8014
#> 
#> ------relative cooccurence: 04
#> 
#> ------checking: ASV11554
#> 
#> ------relative cooccurence: 04
#> 
#> ------checking: ASV16324
#> 
#> ------relative cooccurence: 04
#> 
#> No parent found!
#> 4
#> 
#> ####processing: ASV1313 #####4
#> 
#> ---hits: ASV78
#> ---hits: ASV33
#> ---hits: ASV941
#> ---hits: ASV16254
#> 
#> ---potential parent: ASV78
#> ---potential parent: ASV33
#> ---potential parent: ASV1625
#> ---potential parent: ASV9414
#> 
#> ------checking: ASV784
#> 
#> ------relative cooccurence: 14
#>  which is sufficient!4
#> 
#> ------min avg abundance: 85.78571428571434
#>  which is OK!4
#> 
#> SETTING ASV1313 to be an ERROR of ASV78
#> 4
#> 
#> ------checking: ASV334
#> 
#> ------checking: ASV16254
#> 
#> ------checking: ASV9414
#>   |                                                          |===============================                   |  61%
#> 
#> ####processing: ASV1420 #####4
#> 
#> ---hits: ASV431
#> ---hits: ASV563
#> ---hits: ASV673
#> ---hits: ASV1218
#> ---hits: ASV833
#> ---hits: ASV221
#> ---hits: ASV358
#> ---hits: ASV930
#> ---hits: ASV405
#> ---hits: ASV694
#> 
#> ---potential parent: ASV563
#> ---potential parent: ASV69
#> ---potential parent: ASV358
#> ---potential parent: ASV1218
#> ---potential parent: ASV221
#> ---potential parent: ASV833
#> ---potential parent: ASV431
#> ---potential parent: ASV673
#> ---potential parent: ASV930
#> ---potential parent: ASV4054
#> 
#> ------checking: ASV5634
#> 
#> ------relative cooccurence: 04
#> 
#> ------checking: ASV694
#> 
#> ------relative cooccurence: 04
#> 
#> ------checking: ASV3584
#> 
#> ------relative cooccurence: 04
#> 
#> ------checking: ASV12184
#> 
#> ------relative cooccurence: 14
#>  which is sufficient!4
#> 
#> ------min avg abundance: 0.01428571428571434
#> 
#> ------checking: ASV2214
#> 
#> ------relative cooccurence: 04
#> 
#> ------checking: ASV8334
#> 
#> ------relative cooccurence: 04
#> 
#> ------checking: ASV4314
#> 
#> ------relative cooccurence: 04
#> 
#> ------checking: ASV6734
#> 
#> ------relative cooccurence: 04
#> 
#> ------checking: ASV9304
#> 
#> ------relative cooccurence: 04
#> 
#> ------checking: ASV4054
#> 
#> ------relative cooccurence: 04
#> 
#> No parent found!
#> 4
#> 
#> ####processing: ASV105 #####4
#> 
#> ---hits: ASV12
#> ---hits: ASV67
#> ---hits: ASV107
#> ---hits: ASV1542
#> ---hits: ASV624
#> ---hits: ASV847
#> ---hits: ASV12624
#> 
#> ---potential parent: ASV12
#> ---potential parent: ASV107
#> ---potential parent: ASV624
#> ---potential parent: ASV847
#> ---potential parent: ASV1542
#> ---potential parent: ASV1262
#> ---potential parent: ASV674
#> 
#> ------checking: ASV124
#> 
#> ------relative cooccurence: 14
#>  which is sufficient!4
#> 
#> ------min avg abundance: 0.01470588235294124
#> 
#> ------checking: ASV1074
#> 
#> ------relative cooccurence: 04
#> 
#> ------checking: ASV6244
#> 
#> ------relative cooccurence: 04
#> 
#> ------checking: ASV8474
#> 
#> ------relative cooccurence: 04
#> 
#> ------checking: ASV15424
#> 
#> ------relative cooccurence: 04
#> 
#> ------checking: ASV12624
#> 
#> ------relative cooccurence: 04
#> 
#> ------checking: ASV674
#> 
#> ------relative cooccurence: 04
#> 
#> No parent found!
#> 4
#>   |                                                          |===============================                   |  62%
#> 
#> ####processing: ASV431 #####4
#> 
#> ---hits: ASV563
#> ---hits: ASV1420
#> ---hits: ASV673
#> ---hits: ASV1218
#> ---hits: ASV358
#> ---hits: ASV405
#> ---hits: ASV694
#> 
#> ---potential parent: ASV563
#> ---potential parent: ASV69
#> ---potential parent: ASV358
#> ---potential parent: ASV1218
#> ---potential parent: ASV1420
#> ---potential parent: ASV673
#> ---potential parent: ASV4054
#> 
#> ------checking: ASV5634
#> 
#> ------relative cooccurence: 14
#>  which is sufficient!4
#> 
#> ------min avg abundance: 0.03076923076923084
#> 
#> ------checking: ASV694
#> 
#> ------relative cooccurence: 04
#> 
#> ------checking: ASV3584
#> 
#> ------relative cooccurence: 04
#> 
#> ------checking: ASV12184
#> 
#> ------relative cooccurence: 04
#> 
#> ------checking: ASV14204
#> 
#> ------relative cooccurence: 04
#> 
#> ------checking: ASV6734
#> 
#> ------relative cooccurence: 04
#> 
#> ------checking: ASV4054
#> 
#> ------relative cooccurence: 04
#> 
#> No parent found!
#> 4
#> 
#> ####processing: ASV377 #####4
#> 
#> ---hits: ASV404
#> ---hits: ASV248
#> ---hits: ASV840
#> ---hits: ASV440
#> ---hits: ASV545
#> ---hits: ASV4884
#> 
#> ---potential parent: ASV248
#> ---potential parent: ASV404
#> ---potential parent: ASV545
#> ---potential parent: ASV440
#> ---potential parent: ASV488
#> ---potential parent: ASV8404
#> 
#> ------checking: ASV2484
#> 
#> ------relative cooccurence: 04
#> 
#> ------checking: ASV4044
#> 
#> ------relative cooccurence: 04
#> 
#> ------checking: ASV5454
#> 
#> ------relative cooccurence: 04
#> 
#> ------checking: ASV4404
#> 
#> ------relative cooccurence: 04
#> 
#> ------checking: ASV4884
#> 
#> ------relative cooccurence: 04
#> 
#> ------checking: ASV8404
#> 
#> ------relative cooccurence: 04
#> 
#> No parent found!
#> 4
#> 
#> ####processing: ASV1561 #####4
#> 
#> ---hits: 4
#> 
#> ---potential parent: 4
#> 
#> No parent found!
#> 4
#> 
#> ####processing: ASV477 #####4
#> 
#> ---hits: 4
#> 
#> ---potential parent: 4
#> 
#> No parent found!
#> 4
#>   |                                                          |===============================                   |  63%
#> 
#> ####processing: ASV1537 #####4
#> 
#> ---hits: 4
#> 
#> ---potential parent: 4
#> 
#> No parent found!
#> 4
#> 
#> ####processing: ASV847 #####4
#> 
#> ---hits: ASV12
#> ---hits: ASV67
#> ---hits: ASV107
#> ---hits: ASV1542
#> ---hits: ASV624
#> ---hits: ASV105
#> ---hits: ASV1262
#> ---hits: ASV11324
#> 
#> ---potential parent: ASV12
#> ---potential parent: ASV107
#> ---potential parent: ASV1132
#> ---potential parent: ASV624
#> ---potential parent: ASV105
#> ---potential parent: ASV1542
#> ---potential parent: ASV1262
#> ---potential parent: ASV674
#> 
#> ------checking: ASV124
#> 
#> ------relative cooccurence: 14
#>  which is sufficient!4
#> 
#> ------min avg abundance: 0.02173913043478264
#> 
#> ------checking: ASV1074
#> 
#> ------relative cooccurence: 14
#>  which is sufficient!4
#> 
#> ------min avg abundance: 1.086956521739134
#>  which is OK!4
#> 
#> SETTING ASV847 to be an ERROR of ASV107
#> 4
#> 
#> ------checking: ASV11324
#> 
#> ------checking: ASV6244
#> 
#> ------checking: ASV1054
#> 
#> ------checking: ASV15424
#> 
#> ------checking: ASV12624
#> 
#> ------checking: ASV674
#>   |                                                          |================================                  |  63%
#> 
#> ####processing: ASV722 #####4
#> 
#> ---hits: ASV344
#> ---hits: ASV926
#> ---hits: ASV1058
#> ---hits: ASV1059
#> ---hits: ASV717
#> ---hits: ASV1493
#> ---hits: ASV1022
#> ---hits: ASV1459
#> ---hits: ASV1548
#> ---hits: ASV894
#> 
#> ---potential parent: ASV717
#> ---potential parent: ASV89
#> ---potential parent: ASV1022
#> ---potential parent: ASV1058
#> ---potential parent: ASV344
#> ---potential parent: ASV926
#> ---potential parent: ASV1459
#> ---potential parent: ASV1059
#> ---potential parent: ASV1493
#> ---potential parent: ASV15484
#> 
#> ------checking: ASV7174
#> 
#> ------relative cooccurence: 14
#>  which is sufficient!4
#> 
#> ------min avg abundance: 0.08888888888888894
#> 
#> ------checking: ASV894
#> 
#> ------relative cooccurence: 04
#> 
#> ------checking: ASV10224
#> 
#> ------relative cooccurence: 14
#>  which is sufficient!4
#> 
#> ------min avg abundance: 0.04444444444444444
#> 
#> ------checking: ASV10584
#> 
#> ------relative cooccurence: 04
#> 
#> ------checking: ASV3444
#> 
#> ------relative cooccurence: 14
#>  which is sufficient!4
#> 
#> ------min avg abundance: 7.533333333333334
#>  which is OK!4
#> 
#> SETTING ASV722 to be an ERROR of ASV344
#> 4
#> 
#> ------checking: ASV9264
#> 
#> ------checking: ASV14594
#> 
#> ------checking: ASV10594
#> 
#> ------checking: ASV14934
#> 
#> ------checking: ASV15484
#> 
#> ####processing: ASV404 #####4
#> 
#> ---hits: ASV248
#> ---hits: ASV377
#> ---hits: ASV840
#> ---hits: ASV545
#> ---hits: ASV488
#> ---hits: ASV4404
#> 
#> ---potential parent: ASV248
#> ---potential parent: ASV377
#> ---potential parent: ASV545
#> ---potential parent: ASV440
#> ---potential parent: ASV488
#> ---potential parent: ASV8404
#> 
#> ------checking: ASV2484
#> 
#> ------relative cooccurence: 14
#>  which is sufficient!4
#> 
#> ------min avg abundance: 0.02564102564102564
#> 
#> ------checking: ASV3774
#> 
#> ------relative cooccurence: 04
#> 
#> ------checking: ASV5454
#> 
#> ------relative cooccurence: 04
#> 
#> ------checking: ASV4404
#> 
#> ------relative cooccurence: 04
#> 
#> ------checking: ASV4884
#> 
#> ------relative cooccurence: 14
#>  which is sufficient!4
#> 
#> ------min avg abundance: 0.05128205128205134
#> 
#> ------checking: ASV8404
#> 
#> ------relative cooccurence: 14
#>  which is sufficient!4
#> 
#> ------min avg abundance: 0.02564102564102564
#> 
#> No parent found!
#> 4
#>   |                                                          |================================                  |  64%
#> 
#> ####processing: ASV1510 #####4
#> 
#> ---hits: ASV892
#> ---hits: ASV1239
#> ---hits: ASV1031
#> ---hits: ASV13004
#> 
#> ---potential parent: ASV1300
#> ---potential parent: ASV892
#> ---potential parent: ASV1031
#> ---potential parent: ASV12394
#> 
#> ------checking: ASV13004
#> 
#> ------relative cooccurence: 04
#> 
#> ------checking: ASV8924
#> 
#> ------relative cooccurence: 04
#> 
#> ------checking: ASV10314
#> 
#> ------relative cooccurence: 04
#> 
#> ------checking: ASV12394
#> 
#> ------relative cooccurence: 04
#> 
#> No parent found!
#> 4
#> 
#> ####processing: ASV1526 #####4
#> 
#> ---hits: ASV1107
#> ---hits: ASV116
#> ---hits: ASV989
#> ---hits: ASV816
#> ---hits: ASV38
#> ---hits: ASV466
#> ---hits: ASV254
#> ---hits: ASV500
#> ---hits: ASV853
#> ---hits: ASV554
#> 
#> ---potential parent: ASV38
#> ---potential parent: ASV989
#> ---potential parent: ASV816
#> ---potential parent: ASV500
#> ---potential parent: ASV55
#> ---potential parent: ASV466
#> ---potential parent: ASV254
#> ---potential parent: ASV853
#> ---potential parent: ASV116
#> ---potential parent: ASV11074
#> 
#> ------checking: ASV384
#> 
#> ------relative cooccurence: 14
#>  which is sufficient!4
#> 
#> ------min avg abundance: 0.3076923076923084
#> 
#> ------checking: ASV9894
#> 
#> ------relative cooccurence: 04
#> 
#> ------checking: ASV8164
#> 
#> ------relative cooccurence: 04
#> 
#> ------checking: ASV5004
#> 
#> ------relative cooccurence: 04
#> 
#> ------checking: ASV554
#> 
#> ------relative cooccurence: 14
#>  which is sufficient!4
#> 
#> ------min avg abundance: 1.205128205128214
#>  which is OK!4
#> 
#> SETTING ASV1526 to be an ERROR of ASV55
#> 4
#> 
#> ------checking: ASV4664
#> 
#> ------checking: ASV2544
#> 
#> ------checking: ASV8534
#> 
#> ------checking: ASV1164
#> 
#> ------checking: ASV11074
#> 
#> ####processing: ASV1379 #####4
#> 
#> ---hits: ASV5174
#> 
#> ---potential parent: ASV5174
#> 
#> ------checking: ASV5174
#> 
#> ------relative cooccurence: 04
#> 
#> No parent found!
#> 4
#> 
#> ####processing: ASV367 #####4
#> 
#> ---hits: ASV284
#> 
#> ---potential parent: ASV284
#> 
#> ------checking: ASV284
#> 
#> ------relative cooccurence: 14
#>  which is sufficient!4
#> 
#> ------min avg abundance: 2.166666666666674
#>  which is OK!4
#> 
#> SETTING ASV367 to be an ERROR of ASV28
#> 4
#> 
#> ####processing: ASV1642 #####4
#> 
#> ---hits: ASV1984
#> 
#> ---potential parent: ASV1984
#> 
#> ------checking: ASV1984
#> 
#> ------relative cooccurence: 04
#> 
#> No parent found!
#> 4
#>   |                                                          |================================                  |  65%
#> 
#> ####processing: ASV1233 #####4
#> 
#> ---hits: ASV1396
#> ---hits: ASV569
#> ---hits: ASV546
#> ---hits: ASV727
#> ---hits: ASV1311
#> ---hits: ASV1628
#> ---hits: ASV9274
#> 
#> ---potential parent: ASV727
#> ---potential parent: ASV546
#> ---potential parent: ASV569
#> ---potential parent: ASV1396
#> ---potential parent: ASV927
#> ---potential parent: ASV1311
#> ---potential parent: ASV16284
#> 
#> ------checking: ASV7274
#> 
#> ------relative cooccurence: 14
#>  which is sufficient!4
#> 
#> ------min avg abundance: 0.02941176470588244
#> 
#> ------checking: ASV5464
#> 
#> ------relative cooccurence: 14
#>  which is sufficient!4
#> 
#> ------min avg abundance: 0.2941176470588244
#> 
#> ------checking: ASV5694
#> 
#> ------relative cooccurence: 04
#> 
#> ------checking: ASV13964
#> 
#> ------relative cooccurence: 04
#> 
#> ------checking: ASV9274
#> 
#> ------relative cooccurence: 14
#>  which is sufficient!4
#> 
#> ------min avg abundance: 0.1176470588235294
#> 
#> ------checking: ASV13114
#> 
#> ------relative cooccurence: 04
#> 
#> ------checking: ASV16284
#> 
#> ------relative cooccurence: 04
#> 
#> No parent found!
#> 4
#> 
#> ####processing: ASV1386 #####4
#> 
#> ---hits: ASV1014
#> ---hits: ASV9694
#> 
#> ---potential parent: ASV969
#> ---potential parent: ASV10144
#> 
#> ------checking: ASV9694
#> 
#> ------relative cooccurence: 14
#>  which is sufficient!4
#> 
#> ------min avg abundance: 3.172413793103454
#>  which is OK!4
#> 
#> SETTING ASV1386 to be an ERROR of ASV969
#> 4
#> 
#> ------checking: ASV10144
#>   |                                                          |=================================                 |  65%
#> 
#> ####processing: ASV1234 #####4
#> 
#> ---hits: ASV489
#> ---hits: ASV7624
#> 
#> ---potential parent: ASV489
#> ---potential parent: ASV7624
#> 
#> ------checking: ASV4894
#> 
#> ------relative cooccurence: 04
#> 
#> ------checking: ASV7624
#> 
#> ------relative cooccurence: 04
#> 
#> No parent found!
#> 4
#> 
#> ####processing: ASV1340 #####4
#> 
#> ---hits: ASV1632
#> ---hits: ASV8014
#> 
#> ---potential parent: ASV801
#> ---potential parent: ASV16324
#> 
#> ------checking: ASV8014
#> 
#> ------relative cooccurence: 04
#> 
#> ------checking: ASV16324
#> 
#> ------relative cooccurence: 04
#> 
#> No parent found!
#> 4
#>   |                                                          |=================================                 |  66%
#> 
#> ####processing: ASV1035 #####4
#> 
#> ---hits: ASV904
#> ---hits: ASV756
#> ---hits: ASV1007
#> ---hits: ASV1359
#> ---hits: ASV1278
#> ---hits: ASV1074
#> ---hits: ASV906
#> ---hits: ASV1224
#> ---hits: ASV828
#> ---hits: ASV11924
#> 
#> ---potential parent: ASV756
#> ---potential parent: ASV1192
#> ---potential parent: ASV1074
#> ---potential parent: ASV906
#> ---potential parent: ASV1007
#> ---potential parent: ASV1359
#> ---potential parent: ASV904
#> ---potential parent: ASV1278
#> ---potential parent: ASV1224
#> ---potential parent: ASV8284
#> 
#> ------checking: ASV7564
#> 
#> ------relative cooccurence: 14
#>  which is sufficient!4
#> 
#> ------min avg abundance: 0.164
#> 
#> ------checking: ASV11924
#> 
#> ------relative cooccurence: 04
#> 
#> ------checking: ASV10744
#> 
#> ------relative cooccurence: 14
#>  which is sufficient!4
#> 
#> ------min avg abundance: 0.084
#> 
#> ------checking: ASV9064
#> 
#> ------relative cooccurence: 14
#>  which is sufficient!4
#> 
#> ------min avg abundance: 0.24
#> 
#> ------checking: ASV10074
#> 
#> ------relative cooccurence: 14
#>  which is sufficient!4
#> 
#> ------min avg abundance: 2.724
#>  which is OK!4
#> 
#> SETTING ASV1035 to be an ERROR of ASV1007
#> 4
#> 
#> ------checking: ASV13594
#> 
#> ------checking: ASV9044
#> 
#> ------checking: ASV12784
#> 
#> ------checking: ASV12244
#> 
#> ------checking: ASV8284
#> 
#> ####processing: ASV1037 #####4
#> 
#> ---hits: ASV199
#> ---hits: ASV543
#> ---hits: ASV685
#> ---hits: ASV692
#> ---hits: ASV3104
#> 
#> ---potential parent: ASV310
#> ---potential parent: ASV543
#> ---potential parent: ASV199
#> ---potential parent: ASV685
#> ---potential parent: ASV6924
#> 
#> ------checking: ASV3104
#> 
#> ------relative cooccurence: 14
#>  which is sufficient!4
#> 
#> ------min avg abundance: 0.5217391304347834
#> 
#> ------checking: ASV5434
#> 
#> ------relative cooccurence: 14
#>  which is sufficient!4
#> 
#> ------min avg abundance: 0.2173913043478264
#> 
#> ------checking: ASV1994
#> 
#> ------relative cooccurence: 14
#>  which is sufficient!4
#> 
#> ------min avg abundance: 0.5217391304347834
#> 
#> ------checking: ASV6854
#> 
#> ------relative cooccurence: 04
#> 
#> ------checking: ASV6924
#> 
#> ------relative cooccurence: 14
#>  which is sufficient!4
#> 
#> ------min avg abundance: 15.04347826086964
#>  which is OK!4
#> 
#> SETTING ASV1037 to be an ERROR of ASV692
#> 4
#> 
#> ####processing: ASV48 #####4
#> 
#> ---hits: 4
#> 
#> ---potential parent: 4
#> 
#> No parent found!
#> 4
#> 
#> ####processing: ASV711 #####4
#> 
#> ---hits: ASV45
#> ---hits: ASV137
#> ---hits: ASV996
#> ---hits: ASV760
#> ---hits: ASV1067
#> ---hits: ASV1484
#> ---hits: ASV1167
#> ---hits: ASV940
#> ---hits: ASV1468
#> ---hits: ASV9734
#> 
#> ---potential parent: ASV973
#> ---potential parent: ASV45
#> ---potential parent: ASV1067
#> ---potential parent: ASV137
#> ---potential parent: ASV1468
#> ---potential parent: ASV996
#> ---potential parent: ASV1167
#> ---potential parent: ASV760
#> ---potential parent: ASV940
#> ---potential parent: ASV14844
#> 
#> ------checking: ASV9734
#> 
#> ------relative cooccurence: 04
#> 
#> ------checking: ASV454
#> 
#> ------relative cooccurence: 14
#>  which is sufficient!4
#> 
#> ------min avg abundance: 39.80952380952384
#>  which is OK!4
#> 
#> SETTING ASV711 to be an ERROR of ASV45
#> 4
#> 
#> ------checking: ASV10674
#> 
#> ------checking: ASV1374
#> 
#> ------checking: ASV14684
#> 
#> ------checking: ASV9964
#> 
#> ------checking: ASV11674
#> 
#> ------checking: ASV7604
#> 
#> ------checking: ASV9404
#> 
#> ------checking: ASV14844
#>   |                                                          |=================================                 |  67%
#> 
#> ####processing: ASV1609 #####4
#> 
#> ---hits: ASV1146
#> ---hits: ASV4544
#> 
#> ---potential parent: ASV1146
#> ---potential parent: ASV4544
#> 
#> ------checking: ASV11464
#> 
#> ------relative cooccurence: 04
#> 
#> ------checking: ASV4544
#> 
#> ------relative cooccurence: 04
#> 
#> No parent found!
#> 4
#> 
#> ####processing: ASV33 #####4
#> 
#> ---hits: ASV78
#> ---hits: ASV1313
#> ---hits: ASV941
#> ---hits: ASV16254
#> 
#> ---potential parent: ASV78
#> ---potential parent: ASV1313
#> ---potential parent: ASV1625
#> ---potential parent: ASV9414
#> 
#> ------checking: ASV784
#> 
#> ------relative cooccurence: 14
#>  which is sufficient!4
#> 
#> ------min avg abundance: 316.0526315789474
#>  which is OK!4
#> 
#> SETTING ASV33 to be an ERROR of ASV78
#> 4
#> 
#> ------checking: ASV13134
#> 
#> ------checking: ASV16254
#> 
#> ------checking: ASV9414
#>   |                                                          |==================================                |  67%
#> 
#> ####processing: ASV1683 #####4
#> 
#> ---hits: 4
#> 
#> ---potential parent: 4
#> 
#> No parent found!
#> 4
#> 
#> ####processing: ASV673 #####4
#> 
#> ---hits: ASV1218
#> ---hits: ASV1420
#> ---hits: ASV563
#> ---hits: ASV431
#> ---hits: ASV358
#> ---hits: ASV405
#> ---hits: ASV9304
#> 
#> ---potential parent: ASV563
#> ---potential parent: ASV358
#> ---potential parent: ASV1218
#> ---potential parent: ASV1420
#> ---potential parent: ASV431
#> ---potential parent: ASV930
#> ---potential parent: ASV4054
#> 
#> ------checking: ASV5634
#> 
#> ------relative cooccurence: 04
#> 
#> ------checking: ASV3584
#> 
#> ------relative cooccurence: 04
#> 
#> ------checking: ASV12184
#> 
#> ------relative cooccurence: 04
#> 
#> ------checking: ASV14204
#> 
#> ------relative cooccurence: 04
#> 
#> ------checking: ASV4314
#> 
#> ------relative cooccurence: 04
#> 
#> ------checking: ASV9304
#> 
#> ------relative cooccurence: 04
#> 
#> ------checking: ASV4054
#> 
#> ------relative cooccurence: 04
#> 
#> No parent found!
#> 4
#>   |                                                          |==================================                |  68%
#> 
#> ####processing: ASV1548 #####4
#> 
#> ---hits: ASV89
#> ---hits: ASV462
#> ---hits: ASV534
#> ---hits: ASV1459
#> ---hits: ASV313
#> ---hits: ASV722
#> ---hits: ASV1058
#> ---hits: ASV344
#> ---hits: ASV717
#> ---hits: ASV11014
#> 
#> ---potential parent: ASV717
#> ---potential parent: ASV89
#> ---potential parent: ASV1058
#> ---potential parent: ASV344
#> ---potential parent: ASV1459
#> ---potential parent: ASV313
#> ---potential parent: ASV462
#> ---potential parent: ASV534
#> ---potential parent: ASV1101
#> ---potential parent: ASV7224
#> 
#> ------checking: ASV7174
#> 
#> ------relative cooccurence: 04
#> 
#> ------checking: ASV894
#> 
#> ------relative cooccurence: 04
#> 
#> ------checking: ASV10584
#> 
#> ------relative cooccurence: 04
#> 
#> ------checking: ASV3444
#> 
#> ------relative cooccurence: 04
#> 
#> ------checking: ASV14594
#> 
#> ------relative cooccurence: 04
#> 
#> ------checking: ASV3134
#> 
#> ------relative cooccurence: 04
#> 
#> ------checking: ASV4624
#> 
#> ------relative cooccurence: 04
#> 
#> ------checking: ASV5344
#> 
#> ------relative cooccurence: 04
#> 
#> ------checking: ASV11014
#> 
#> ------relative cooccurence: 04
#> 
#> ------checking: ASV7224
#> 
#> ------relative cooccurence: 04
#> 
#> No parent found!
#> 4
#> 
#> ####processing: ASV1460 #####4
#> 
#> ---hits: ASV13414
#> 
#> ---potential parent: ASV13414
#> 
#> ------checking: ASV13414
#> 
#> ------relative cooccurence: 14
#>  which is sufficient!4
#> 
#> ------min avg abundance: 0.07142857142857144
#> 
#> No parent found!
#> 4
#> 
#> ####processing: ASV545 #####4
#> 
#> ---hits: ASV404
#> ---hits: ASV248
#> ---hits: ASV377
#> ---hits: ASV840
#> ---hits: ASV440
#> ---hits: ASV4884
#> 
#> ---potential parent: ASV248
#> ---potential parent: ASV377
#> ---potential parent: ASV404
#> ---potential parent: ASV440
#> ---potential parent: ASV488
#> ---potential parent: ASV8404
#> 
#> ------checking: ASV2484
#> 
#> ------relative cooccurence: 04
#> 
#> ------checking: ASV3774
#> 
#> ------relative cooccurence: 04
#> 
#> ------checking: ASV4044
#> 
#> ------relative cooccurence: 04
#> 
#> ------checking: ASV4404
#> 
#> ------relative cooccurence: 14
#>  which is sufficient!4
#> 
#> ------min avg abundance: 0.7692307692307694
#> 
#> ------checking: ASV4884
#> 
#> ------relative cooccurence: 04
#> 
#> ------checking: ASV8404
#> 
#> ------relative cooccurence: 04
#> 
#> No parent found!
#> 4
#> 
#> ####processing: ASV1410 #####4
#> 
#> ---hits: ASV420
#> ---hits: ASV505
#> ---hits: ASV509
#> ---hits: ASV790
#> ---hits: ASV1427
#> ---hits: ASV662
#> ---hits: ASV1052
#> ---hits: ASV724
#> ---hits: ASV1198
#> ---hits: ASV8594
#> 
#> ---potential parent: ASV662
#> ---potential parent: ASV505
#> ---potential parent: ASV509
#> ---potential parent: ASV724
#> ---potential parent: ASV859
#> ---potential parent: ASV1052
#> ---potential parent: ASV420
#> ---potential parent: ASV1198
#> ---potential parent: ASV1427
#> ---potential parent: ASV7904
#> 
#> ------checking: ASV6624
#> 
#> ------relative cooccurence: 14
#>  which is sufficient!4
#> 
#> ------min avg abundance: 0.07692307692307694
#> 
#> ------checking: ASV5054
#> 
#> ------relative cooccurence: 14
#>  which is sufficient!4
#> 
#> ------min avg abundance: 3.076923076923084
#>  which is OK!4
#> 
#> SETTING ASV1410 to be an ERROR of ASV505
#> 4
#> 
#> ------checking: ASV5094
#> 
#> ------checking: ASV7244
#> 
#> ------checking: ASV8594
#> 
#> ------checking: ASV10524
#> 
#> ------checking: ASV4204
#> 
#> ------checking: ASV11984
#> 
#> ------checking: ASV14274
#> 
#> ------checking: ASV7904
#>   |                                                          |==================================                |  69%
#> 
#> ####processing: ASV1462 #####4
#> 
#> ---hits: 4
#> 
#> ---potential parent: 4
#> 
#> No parent found!
#> 4
#> 
#> ####processing: ASV853 #####4
#> 
#> ---hits: ASV38
#> ---hits: ASV466
#> ---hits: ASV989
#> ---hits: ASV816
#> ---hits: ASV975
#> ---hits: ASV1526
#> ---hits: ASV1107
#> ---hits: ASV254
#> ---hits: ASV500
#> ---hits: ASV5924
#> 
#> ---potential parent: ASV38
#> ---potential parent: ASV989
#> ---potential parent: ASV816
#> ---potential parent: ASV500
#> ---potential parent: ASV592
#> ---potential parent: ASV975
#> ---potential parent: ASV466
#> ---potential parent: ASV254
#> ---potential parent: ASV1526
#> ---potential parent: ASV11074
#> 
#> ------checking: ASV384
#> 
#> ------relative cooccurence: 14
#>  which is sufficient!4
#> 
#> ------min avg abundance: 28.08333333333334
#>  which is OK!4
#> 
#> SETTING ASV853 to be an ERROR of ASV38
#> 4
#> 
#> ------checking: ASV9894
#> 
#> ------checking: ASV8164
#> 
#> ------checking: ASV5004
#> 
#> ------checking: ASV5924
#> 
#> ------checking: ASV9754
#> 
#> ------checking: ASV4664
#> 
#> ------checking: ASV2544
#> 
#> ------checking: ASV15264
#> 
#> ------checking: ASV11074
#>   |                                                          |===================================               |  69%
#> 
#> ####processing: ASV1067 #####4
#> 
#> ---hits: ASV45
#> ---hits: ASV137
#> ---hits: ASV996
#> ---hits: ASV1167
#> ---hits: ASV711
#> ---hits: ASV760
#> ---hits: ASV940
#> ---hits: ASV1484
#> ---hits: ASV1468
#> ---hits: ASV9734
#> 
#> ---potential parent: ASV973
#> ---potential parent: ASV45
#> ---potential parent: ASV711
#> ---potential parent: ASV137
#> ---potential parent: ASV1468
#> ---potential parent: ASV996
#> ---potential parent: ASV1167
#> ---potential parent: ASV760
#> ---potential parent: ASV940
#> ---potential parent: ASV14844
#> 
#> ------checking: ASV9734
#> 
#> ------relative cooccurence: 04
#> 
#> ------checking: ASV454
#> 
#> ------relative cooccurence: 14
#>  which is sufficient!4
#> 
#> ------min avg abundance: 69.66666666666674
#>  which is OK!4
#> 
#> SETTING ASV1067 to be an ERROR of ASV45
#> 4
#> 
#> ------checking: ASV7114
#> 
#> ------checking: ASV1374
#> 
#> ------checking: ASV14684
#> 
#> ------checking: ASV9964
#> 
#> ------checking: ASV11674
#> 
#> ------checking: ASV7604
#> 
#> ------checking: ASV9404
#> 
#> ------checking: ASV14844
#> 
#> ####processing: ASV1311 #####4
#> 
#> ---hits: ASV927
#> ---hits: ASV1396
#> ---hits: ASV604
#> ---hits: ASV891
#> ---hits: ASV1233
#> ---hits: ASV546
#> ---hits: ASV569
#> ---hits: ASV7274
#> 
#> ---potential parent: ASV727
#> ---potential parent: ASV546
#> ---potential parent: ASV569
#> ---potential parent: ASV1396
#> ---potential parent: ASV927
#> ---potential parent: ASV1233
#> ---potential parent: ASV604
#> ---potential parent: ASV8914
#> 
#> ------checking: ASV7274
#> 
#> ------relative cooccurence: 04
#> 
#> ------checking: ASV5464
#> 
#> ------relative cooccurence: 04
#> 
#> ------checking: ASV5694
#> 
#> ------relative cooccurence: 04
#> 
#> ------checking: ASV13964
#> 
#> ------relative cooccurence: 04
#> 
#> ------checking: ASV9274
#> 
#> ------relative cooccurence: 04
#> 
#> ------checking: ASV12334
#> 
#> ------relative cooccurence: 04
#> 
#> ------checking: ASV6044
#> 
#> ------relative cooccurence: 14
#>  which is sufficient!4
#> 
#> ------min avg abundance: 0.1666666666666674
#> 
#> ------checking: ASV8914
#> 
#> ------relative cooccurence: 04
#> 
#> No parent found!
#> 4
#>   |                                                          |===================================               |  70%
#> 
#> ####processing: ASV1552 #####4
#> 
#> ---hits: ASV526
#> ---hits: ASV11334
#> 
#> ---potential parent: ASV526
#> ---potential parent: ASV11334
#> 
#> ------checking: ASV5264
#> 
#> ------relative cooccurence: 14
#>  which is sufficient!4
#> 
#> ------min avg abundance: 0.754
#> 
#> ------checking: ASV11334
#> 
#> ------relative cooccurence: 14
#>  which is sufficient!4
#> 
#> ------min avg abundance: 0.1666666666666674
#> 
#> No parent found!
#> 4
#> 
#> ####processing: ASV1588 #####4
#> 
#> ---hits: ASV1236
#> ---hits: ASV13004
#> 
#> ---potential parent: ASV1236
#> ---potential parent: ASV13004
#> 
#> ------checking: ASV12364
#> 
#> ------relative cooccurence: 04
#> 
#> ------checking: ASV13004
#> 
#> ------relative cooccurence: 04
#> 
#> No parent found!
#> 4
#> 
#> ####processing: ASV359 #####4
#> 
#> ---hits: ASV15504
#> 
#> ---potential parent: ASV15504
#> 
#> ------checking: ASV15504
#> 
#> ------relative cooccurence: 04
#> 
#> No parent found!
#> 4
#> 
#> ####processing: ASV440 #####4
#> 
#> ---hits: ASV248
#> ---hits: ASV377
#> ---hits: ASV404
#> ---hits: ASV545
#> ---hits: ASV488
#> ---hits: ASV8404
#> 
#> ---potential parent: ASV248
#> ---potential parent: ASV377
#> ---potential parent: ASV404
#> ---potential parent: ASV545
#> ---potential parent: ASV488
#> ---potential parent: ASV8404
#> 
#> ------checking: ASV2484
#> 
#> ------relative cooccurence: 04
#> 
#> ------checking: ASV3774
#> 
#> ------relative cooccurence: 04
#> 
#> ------checking: ASV4044
#> 
#> ------relative cooccurence: 04
#> 
#> ------checking: ASV5454
#> 
#> ------relative cooccurence: 14
#>  which is sufficient!4
#> 
#> ------min avg abundance: 1.34
#>  which is OK!4
#> 
#> SETTING ASV440 to be an ERROR of ASV545
#> 4
#> 
#> ------checking: ASV4884
#> 
#> ------checking: ASV8404
#>   |                                                          |===================================               |  71%
#> 
#> ####processing: ASV580 #####4
#> 
#> ---hits: ASV64
#> ---hits: ASV168
#> ---hits: ASV203
#> ---hits: ASV2494
#> 
#> ---potential parent: ASV64
#> ---potential parent: ASV168
#> ---potential parent: ASV203
#> ---potential parent: ASV2494
#> 
#> ------checking: ASV644
#> 
#> ------relative cooccurence: 14
#>  which is sufficient!4
#> 
#> ------min avg abundance: 818.74
#>  which is OK!4
#> 
#> SETTING ASV580 to be an ERROR of ASV64
#> 4
#> 
#> ------checking: ASV1684
#> 
#> ------checking: ASV2034
#> 
#> ------checking: ASV2494
#> 
#> ####processing: ASV858 #####4
#> 
#> ---hits: ASV19
#> ---hits: ASV979
#> ---hits: ASV333
#> ---hits: ASV736
#> ---hits: ASV178
#> ---hits: ASV594
#> ---hits: ASV24
#> ---hits: ASV1131
#> ---hits: ASV1010
#> ---hits: ASV2664
#> 
#> ---potential parent: ASV19
#> ---potential parent: ASV178
#> ---potential parent: ASV24
#> ---potential parent: ASV266
#> ---potential parent: ASV333
#> ---potential parent: ASV594
#> ---potential parent: ASV1131
#> ---potential parent: ASV1010
#> ---potential parent: ASV736
#> ---potential parent: ASV9794
#> 
#> ------checking: ASV194
#> 
#> ------relative cooccurence: 14
#>  which is sufficient!4
#> 
#> ------min avg abundance: 297.34
#>  which is OK!4
#> 
#> SETTING ASV858 to be an ERROR of ASV19
#> 4
#> 
#> ------checking: ASV1784
#> 
#> ------checking: ASV244
#> 
#> ------checking: ASV2664
#> 
#> ------checking: ASV3334
#> 
#> ------checking: ASV5944
#> 
#> ------checking: ASV11314
#> 
#> ------checking: ASV10104
#> 
#> ------checking: ASV7364
#> 
#> ------checking: ASV9794
#>   |                                                          |====================================              |  71%
#> 
#> ####processing: ASV1632 #####4
#> 
#> ---hits: ASV1340
#> ---hits: ASV801
#> ---hits: ASV1461
#> ---hits: ASV1155
#> ---hits: ASV9464
#> 
#> ---potential parent: ASV801
#> ---potential parent: ASV1155
#> ---potential parent: ASV946
#> ---potential parent: ASV1461
#> ---potential parent: ASV13404
#> 
#> ------checking: ASV8014
#> 
#> ------relative cooccurence: 04
#> 
#> ------checking: ASV11554
#> 
#> ------relative cooccurence: 04
#> 
#> ------checking: ASV9464
#> 
#> ------relative cooccurence: 04
#> 
#> ------checking: ASV14614
#> 
#> ------relative cooccurence: 04
#> 
#> ------checking: ASV13404
#> 
#> ------relative cooccurence: 04
#> 
#> No parent found!
#> 4
#> 
#> ####processing: ASV526 #####4
#> 
#> ---hits: ASV1133
#> ---hits: ASV15524
#> 
#> ---potential parent: ASV1552
#> ---potential parent: ASV11334
#> 
#> ------checking: ASV15524
#> 
#> ------relative cooccurence: 14
#>  which is sufficient!4
#> 
#> ------min avg abundance: 1.333333333333334
#>  which is OK!4
#> 
#> SETTING ASV526 to be an ERROR of ASV1552
#> 4
#> 
#> ------checking: ASV11334
#>   |                                                          |====================================              |  72%
#> 
#> ####processing: ASV1179 #####4
#> 
#> ---hits: 4
#> 
#> ---potential parent: 4
#> 
#> No parent found!
#> 4
#> 
#> ####processing: ASV1428 #####4
#> 
#> ---hits: 4
#> 
#> ---potential parent: 4
#> 
#> No parent found!
#> 4
#> 
#> ####processing: ASV355 #####4
#> 
#> ---hits: ASV151
#> ---hits: ASV8344
#> 
#> ---potential parent: ASV151
#> ---potential parent: ASV8344
#> 
#> ------checking: ASV1514
#> 
#> ------relative cooccurence: 14
#>  which is sufficient!4
#> 
#> ------min avg abundance: 1.6254
#>  which is OK!4
#> 
#> SETTING ASV355 to be an ERROR of ASV151
#> 4
#> 
#> ------checking: ASV8344
#> 
#> ####processing: ASV1010 #####4
#> 
#> ---hits: ASV19
#> ---hits: ASV333
#> ---hits: ASV986
#> ---hits: ASV178
#> ---hits: ASV594
#> ---hits: ASV1131
#> ---hits: ASV1078
#> ---hits: ASV1467
#> ---hits: ASV858
#> ---hits: ASV244
#> 
#> ---potential parent: ASV19
#> ---potential parent: ASV178
#> ---potential parent: ASV24
#> ---potential parent: ASV333
#> ---potential parent: ASV1078
#> ---potential parent: ASV594
#> ---potential parent: ASV1131
#> ---potential parent: ASV858
#> ---potential parent: ASV1467
#> ---potential parent: ASV9864
#> 
#> ------checking: ASV194
#> 
#> ------relative cooccurence: 14
#>  which is sufficient!4
#> 
#> ------min avg abundance: 371.6254
#>  which is OK!4
#> 
#> SETTING ASV1010 to be an ERROR of ASV19
#> 4
#> 
#> ------checking: ASV1784
#> 
#> ------checking: ASV244
#> 
#> ------checking: ASV3334
#> 
#> ------checking: ASV10784
#> 
#> ------checking: ASV5944
#> 
#> ------checking: ASV11314
#> 
#> ------checking: ASV8584
#> 
#> ------checking: ASV14674
#> 
#> ------checking: ASV9864
#>   |                                                          |====================================              |  73%
#> 
#> ####processing: ASV1542 #####4
#> 
#> ---hits: ASV12
#> ---hits: ASV67
#> ---hits: ASV107
#> ---hits: ASV624
#> ---hits: ASV1262
#> ---hits: ASV847
#> ---hits: ASV1054
#> 
#> ---potential parent: ASV12
#> ---potential parent: ASV107
#> ---potential parent: ASV624
#> ---potential parent: ASV105
#> ---potential parent: ASV847
#> ---potential parent: ASV1262
#> ---potential parent: ASV674
#> 
#> ------checking: ASV124
#> 
#> ------relative cooccurence: 14
#>  which is sufficient!4
#> 
#> ------min avg abundance: 398.254
#>  which is OK!4
#> 
#> SETTING ASV1542 to be an ERROR of ASV12
#> 4
#> 
#> ------checking: ASV1074
#> 
#> ------checking: ASV6244
#> 
#> ------checking: ASV1054
#> 
#> ------checking: ASV8474
#> 
#> ------checking: ASV12624
#> 
#> ------checking: ASV674
#> 
#> ####processing: ASV1577 #####4
#> 
#> ---hits: ASV12424
#> 
#> ---potential parent: ASV12424
#> 
#> ------checking: ASV12424
#> 
#> ------relative cooccurence: 04
#> 
#> No parent found!
#> 4
#> 
#> ####processing: ASV1262 #####4
#> 
#> ---hits: ASV624
#> ---hits: ASV67
#> ---hits: ASV12
#> ---hits: ASV107
#> ---hits: ASV1542
#> ---hits: ASV847
#> ---hits: ASV1054
#> 
#> ---potential parent: ASV12
#> ---potential parent: ASV107
#> ---potential parent: ASV624
#> ---potential parent: ASV105
#> ---potential parent: ASV847
#> ---potential parent: ASV1542
#> ---potential parent: ASV674
#> 
#> ------checking: ASV124
#> 
#> ------relative cooccurence: 14
#>  which is sufficient!4
#> 
#> ------min avg abundance: 455.1428571428574
#>  which is OK!4
#> 
#> SETTING ASV1262 to be an ERROR of ASV12
#> 4
#> 
#> ------checking: ASV1074
#> 
#> ------checking: ASV6244
#> 
#> ------checking: ASV1054
#> 
#> ------checking: ASV8474
#> 
#> ------checking: ASV15424
#> 
#> ------checking: ASV674
#>   |                                                          |=====================================             |  73%
#> 
#> ####processing: ASV27 #####4
#> 
#> ---hits: 4
#> 
#> ---potential parent: 4
#> 
#> No parent found!
#> 4
#> 
#> ####processing: ASV137 #####4
#> 
#> ---hits: ASV45
#> ---hits: ASV760
#> ---hits: ASV711
#> ---hits: ASV996
#> ---hits: ASV1067
#> ---hits: ASV1167
#> ---hits: ASV1484
#> ---hits: ASV940
#> ---hits: ASV1468
#> ---hits: ASV9734
#> 
#> ---potential parent: ASV973
#> ---potential parent: ASV45
#> ---potential parent: ASV711
#> ---potential parent: ASV1067
#> ---potential parent: ASV1468
#> ---potential parent: ASV996
#> ---potential parent: ASV1167
#> ---potential parent: ASV760
#> ---potential parent: ASV940
#> ---potential parent: ASV14844
#> 
#> ------checking: ASV9734
#> 
#> ------relative cooccurence: 04
#> 
#> ------checking: ASV454
#> 
#> ------relative cooccurence: 04
#> 
#> ------checking: ASV7114
#> 
#> ------relative cooccurence: 04
#> 
#> ------checking: ASV10674
#> 
#> ------relative cooccurence: 04
#> 
#> ------checking: ASV14684
#> 
#> ------relative cooccurence: 04
#> 
#> ------checking: ASV9964
#> 
#> ------relative cooccurence: 04
#> 
#> ------checking: ASV11674
#> 
#> ------relative cooccurence: 04
#> 
#> ------checking: ASV7604
#> 
#> ------relative cooccurence: 04
#> 
#> ------checking: ASV9404
#> 
#> ------relative cooccurence: 04
#> 
#> ------checking: ASV14844
#> 
#> ------relative cooccurence: 04
#> 
#> No parent found!
#> 4
#>   |                                                          |=====================================             |  74%
#> 
#> ####processing: ASV930 #####4
#> 
#> ---hits: ASV914
#> ---hits: ASV358
#> ---hits: ASV405
#> ---hits: ASV69
#> ---hits: ASV833
#> ---hits: ASV221
#> ---hits: ASV1420
#> ---hits: ASV1218
#> ---hits: ASV6734
#> 
#> ---potential parent: ASV69
#> ---potential parent: ASV358
#> ---potential parent: ASV1218
#> ---potential parent: ASV221
#> ---potential parent: ASV833
#> ---potential parent: ASV1420
#> ---potential parent: ASV673
#> ---potential parent: ASV914
#> ---potential parent: ASV4054
#> 
#> ------checking: ASV694
#> 
#> ------relative cooccurence: 04
#> 
#> ------checking: ASV3584
#> 
#> ------relative cooccurence: 04
#> 
#> ------checking: ASV12184
#> 
#> ------relative cooccurence: 04
#> 
#> ------checking: ASV2214
#> 
#> ------relative cooccurence: 04
#> 
#> ------checking: ASV8334
#> 
#> ------relative cooccurence: 04
#> 
#> ------checking: ASV14204
#> 
#> ------relative cooccurence: 04
#> 
#> ------checking: ASV6734
#> 
#> ------relative cooccurence: 04
#> 
#> ------checking: ASV9144
#> 
#> ------relative cooccurence: 04
#> 
#> ------checking: ASV4054
#> 
#> ------relative cooccurence: 04
#> 
#> No parent found!
#> 4
#> 
#> ####processing: ASV962 #####4
#> 
#> ---hits: ASV1027
#> ---hits: ASV1247
#> ---hits: ASV1082
#> ---hits: ASV576
#> ---hits: ASV24
#> 
#> ---potential parent: ASV2
#> ---potential parent: ASV576
#> ---potential parent: ASV1247
#> ---potential parent: ASV1082
#> ---potential parent: ASV10274
#> 
#> ------checking: ASV24
#> 
#> ------relative cooccurence: 14
#>  which is sufficient!4
#> 
#> ------min avg abundance: 46.83333333333334
#>  which is OK!4
#> 
#> SETTING ASV962 to be an ERROR of ASV2
#> 4
#> 
#> ------checking: ASV5764
#> 
#> ------checking: ASV12474
#> 
#> ------checking: ASV10824
#> 
#> ------checking: ASV10274
#> 
#> ####processing: ASV1242 #####4
#> 
#> ---hits: ASV15774
#> 
#> ---potential parent: ASV15774
#> 
#> ------checking: ASV15774
#> 
#> ------relative cooccurence: 04
#> 
#> No parent found!
#> 4
#> 
#> ####processing: ASV1628 #####4
#> 
#> ---hits: ASV1109
#> ---hits: ASV12334
#> 
#> ---potential parent: ASV1233
#> ---potential parent: ASV11094
#> 
#> ------checking: ASV12334
#> 
#> ------relative cooccurence: 04
#> 
#> ------checking: ASV11094
#> 
#> ------relative cooccurence: 04
#> 
#> No parent found!
#> 4
#>   |                                                          |=====================================             |  75%
#> 
#> ####processing: ASV1650 #####4
#> 
#> ---hits: 4
#> 
#> ---potential parent: 4
#> 
#> No parent found!
#> 4
#> 
#> ####processing: ASV903 #####4
#> 
#> ---hits: ASV637
#> ---hits: ASV10154
#> 
#> ---potential parent: ASV637
#> ---potential parent: ASV10154
#> 
#> ------checking: ASV6374
#> 
#> ------relative cooccurence: 14
#>  which is sufficient!4
#> 
#> ------min avg abundance: 33.24
#>  which is OK!4
#> 
#> SETTING ASV903 to be an ERROR of ASV637
#> 4
#> 
#> ------checking: ASV10154
#>   |                                                          |======================================            |  75%
#> 
#> ####processing: ASV984 #####4
#> 
#> ---hits: ASV796
#> ---hits: ASV60
#> ---hits: ASV113
#> ---hits: ASV1020
#> ---hits: ASV1380
#> ---hits: ASV415
#> ---hits: ASV1036
#> ---hits: ASV1011
#> ---hits: ASV1603
#> ---hits: ASV10664
#> 
#> ---potential parent: ASV113
#> ---potential parent: ASV60
#> ---potential parent: ASV1036
#> ---potential parent: ASV1011
#> ---potential parent: ASV1603
#> ---potential parent: ASV415
#> ---potential parent: ASV796
#> ---potential parent: ASV1020
#> ---potential parent: ASV1066
#> ---potential parent: ASV13804
#> 
#> ------checking: ASV1134
#> 
#> ------relative cooccurence: 14
#>  which is sufficient!4
#> 
#> ------min avg abundance: 16.64
#>  which is OK!4
#> 
#> SETTING ASV984 to be an ERROR of ASV113
#> 4
#> 
#> ------checking: ASV604
#> 
#> ------checking: ASV10364
#> 
#> ------checking: ASV10114
#> 
#> ------checking: ASV16034
#> 
#> ------checking: ASV4154
#> 
#> ------checking: ASV7964
#> 
#> ------checking: ASV10204
#> 
#> ------checking: ASV10664
#> 
#> ------checking: ASV13804
#> 
#> ####processing: ASV1011 #####4
#> 
#> ---hits: ASV1066
#> ---hits: ASV60
#> ---hits: ASV113
#> ---hits: ASV984
#> ---hits: ASV415
#> ---hits: ASV796
#> ---hits: ASV1020
#> ---hits: ASV1380
#> ---hits: ASV1036
#> ---hits: ASV15324
#> 
#> ---potential parent: ASV113
#> ---potential parent: ASV60
#> ---potential parent: ASV1532
#> ---potential parent: ASV1036
#> ---potential parent: ASV984
#> ---potential parent: ASV415
#> ---potential parent: ASV796
#> ---potential parent: ASV1020
#> ---potential parent: ASV1066
#> ---potential parent: ASV13804
#> 
#> ------checking: ASV1134
#> 
#> ------relative cooccurence: 04
#> 
#> ------checking: ASV604
#> 
#> ------relative cooccurence: 04
#> 
#> ------checking: ASV15324
#> 
#> ------relative cooccurence: 14
#>  which is sufficient!4
#> 
#> ------min avg abundance: 7.24
#>  which is OK!4
#> 
#> SETTING ASV1011 to be an ERROR of ASV1532
#> 4
#> 
#> ------checking: ASV10364
#> 
#> ------checking: ASV9844
#> 
#> ------checking: ASV4154
#> 
#> ------checking: ASV7964
#> 
#> ------checking: ASV10204
#> 
#> ------checking: ASV10664
#> 
#> ------checking: ASV13804
#>   |                                                          |======================================            |  76%
#> 
#> ####processing: ASV1224 #####4
#> 
#> ---hits: ASV1278
#> ---hits: ASV1359
#> ---hits: ASV1074
#> ---hits: ASV1007
#> ---hits: ASV906
#> ---hits: ASV1192
#> ---hits: ASV828
#> ---hits: ASV950
#> ---hits: ASV756
#> ---hits: ASV13194
#> 
#> ---potential parent: ASV756
#> ---potential parent: ASV1192
#> ---potential parent: ASV1074
#> ---potential parent: ASV906
#> ---potential parent: ASV1007
#> ---potential parent: ASV1359
#> ---potential parent: ASV1278
#> ---potential parent: ASV950
#> ---potential parent: ASV1319
#> ---potential parent: ASV8284
#> 
#> ------checking: ASV7564
#> 
#> ------relative cooccurence: 04
#> 
#> ------checking: ASV11924
#> 
#> ------relative cooccurence: 14
#>  which is sufficient!4
#> 
#> ------min avg abundance: 14
#> 
#> ------checking: ASV10744
#> 
#> ------relative cooccurence: 14
#>  which is sufficient!4
#> 
#> ------min avg abundance: 0.24
#> 
#> ------checking: ASV9064
#> 
#> ------relative cooccurence: 04
#> 
#> ------checking: ASV10074
#> 
#> ------relative cooccurence: 14
#>  which is sufficient!4
#> 
#> ------min avg abundance: 0.84
#> 
#> ------checking: ASV13594
#> 
#> ------relative cooccurence: 14
#>  which is sufficient!4
#> 
#> ------min avg abundance: 1.84
#>  which is OK!4
#> 
#> SETTING ASV1224 to be an ERROR of ASV1359
#> 4
#> 
#> ------checking: ASV12784
#> 
#> ------checking: ASV9504
#> 
#> ------checking: ASV13194
#> 
#> ------checking: ASV8284
#> 
#> ####processing: ASV1276 #####4
#> 
#> ---hits: ASV159
#> ---hits: ASV52
#> ---hits: ASV8804
#> 
#> ---potential parent: ASV52
#> ---potential parent: ASV159
#> ---potential parent: ASV8804
#> 
#> ------checking: ASV524
#> 
#> ------relative cooccurence: 14
#>  which is sufficient!4
#> 
#> ------min avg abundance: 154
#>  which is OK!4
#> 
#> SETTING ASV1276 to be an ERROR of ASV52
#> 4
#> 
#> ------checking: ASV1594
#> 
#> ------checking: ASV8804
#> 
#> ####processing: ASV1468 #####4
#> 
#> ---hits: ASV973
#> ---hits: ASV1388
#> ---hits: ASV1653
#> ---hits: ASV45
#> ---hits: ASV137
#> ---hits: ASV996
#> ---hits: ASV1067
#> ---hits: ASV711
#> ---hits: ASV760
#> ---hits: ASV12614
#> 
#> ---potential parent: ASV1261
#> ---potential parent: ASV973
#> ---potential parent: ASV45
#> ---potential parent: ASV711
#> ---potential parent: ASV1067
#> ---potential parent: ASV137
#> ---potential parent: ASV996
#> ---potential parent: ASV760
#> ---potential parent: ASV1388
#> ---potential parent: ASV16534
#> 
#> ------checking: ASV12614
#> 
#> ------relative cooccurence: 04
#> 
#> ------checking: ASV9734
#> 
#> ------relative cooccurence: 14
#>  which is sufficient!4
#> 
#> ------min avg abundance: 14
#> 
#> ------checking: ASV454
#> 
#> ------relative cooccurence: 04
#> 
#> ------checking: ASV7114
#> 
#> ------relative cooccurence: 04
#> 
#> ------checking: ASV10674
#> 
#> ------relative cooccurence: 04
#> 
#> ------checking: ASV1374
#> 
#> ------relative cooccurence: 04
#> 
#> ------checking: ASV9964
#> 
#> ------relative cooccurence: 04
#> 
#> ------checking: ASV7604
#> 
#> ------relative cooccurence: 04
#> 
#> ------checking: ASV13884
#> 
#> ------relative cooccurence: 04
#> 
#> ------checking: ASV16534
#> 
#> ------relative cooccurence: 04
#> 
#> No parent found!
#> 4
#> 
#> ####processing: ASV1574 #####4
#> 
#> ---hits: 4
#> 
#> ---potential parent: 4
#> 
#> No parent found!
#> 4
#>   |                                                          |======================================            |  77%
#> 
#> ####processing: ASV1630 #####4
#> 
#> ---hits: 4
#> 
#> ---potential parent: 4
#> 
#> No parent found!
#> 4
#> 
#> ####processing: ASV210 #####4
#> 
#> ---hits: ASV491
#> ---hits: ASV12234
#> 
#> ---potential parent: ASV1223
#> ---potential parent: ASV4914
#> 
#> ------checking: ASV12234
#> 
#> ------relative cooccurence: 14
#>  which is sufficient!4
#> 
#> ------min avg abundance: 0.254
#> 
#> ------checking: ASV4914
#> 
#> ------relative cooccurence: 14
#>  which is sufficient!4
#> 
#> ------min avg abundance: 0.754
#> 
#> No parent found!
#> 4
#>   |                                                          |=======================================           |  77%
#> 
#> ####processing: ASV383 #####4
#> 
#> ---hits: ASV272
#> ---hits: ASV339
#> ---hits: ASV81
#> ---hits: ASV959
#> ---hits: ASV839
#> ---hits: ASV1152
#> ---hits: ASV1263
#> ---hits: ASV11444
#> 
#> ---potential parent: ASV1152
#> ---potential parent: ASV81
#> ---potential parent: ASV339
#> ---potential parent: ASV1144
#> ---potential parent: ASV839
#> ---potential parent: ASV272
#> ---potential parent: ASV959
#> ---potential parent: ASV12634
#> 
#> ------checking: ASV11524
#> 
#> ------relative cooccurence: 04
#> 
#> ------checking: ASV814
#> 
#> ------relative cooccurence: 04
#> 
#> ------checking: ASV3394
#> 
#> ------relative cooccurence: 04
#> 
#> ------checking: ASV11444
#> 
#> ------relative cooccurence: 04
#> 
#> ------checking: ASV8394
#> 
#> ------relative cooccurence: 04
#> 
#> ------checking: ASV2724
#> 
#> ------relative cooccurence: 04
#> 
#> ------checking: ASV9594
#> 
#> ------relative cooccurence: 04
#> 
#> ------checking: ASV12634
#> 
#> ------relative cooccurence: 04
#> 
#> No parent found!
#> 4
#> 
#> ####processing: ASV753 #####4
#> 
#> ---hits: ASV171
#> ---hits: ASV7484
#> 
#> ---potential parent: ASV748
#> ---potential parent: ASV1714
#> 
#> ------checking: ASV7484
#> 
#> ------relative cooccurence: 14
#>  which is sufficient!4
#> 
#> ------min avg abundance: 0.54
#> 
#> ------checking: ASV1714
#> 
#> ------relative cooccurence: 04
#> 
#> No parent found!
#> 4
#>   |                                                          |=======================================           |  78%
#> 
#> ####processing: ASV996 #####4
#> 
#> ---hits: ASV45
#> ---hits: ASV137
#> ---hits: ASV1067
#> ---hits: ASV711
#> ---hits: ASV760
#> ---hits: ASV1167
#> ---hits: ASV940
#> ---hits: ASV1484
#> ---hits: ASV1468
#> ---hits: ASV9734
#> 
#> ---potential parent: ASV973
#> ---potential parent: ASV45
#> ---potential parent: ASV711
#> ---potential parent: ASV1067
#> ---potential parent: ASV137
#> ---potential parent: ASV1468
#> ---potential parent: ASV1167
#> ---potential parent: ASV760
#> ---potential parent: ASV940
#> ---potential parent: ASV14844
#> 
#> ------checking: ASV9734
#> 
#> ------relative cooccurence: 04
#> 
#> ------checking: ASV454
#> 
#> ------relative cooccurence: 14
#>  which is sufficient!4
#> 
#> ------min avg abundance: 2094
#>  which is OK!4
#> 
#> SETTING ASV996 to be an ERROR of ASV45
#> 4
#> 
#> ------checking: ASV7114
#> 
#> ------checking: ASV10674
#> 
#> ------checking: ASV1374
#> 
#> ------checking: ASV14684
#> 
#> ------checking: ASV11674
#> 
#> ------checking: ASV7604
#> 
#> ------checking: ASV9404
#> 
#> ------checking: ASV14844
#> 
#> ####processing: ASV1069 #####4
#> 
#> ---hits: ASV9464
#> 
#> ---potential parent: ASV9464
#> 
#> ------checking: ASV9464
#> 
#> ------relative cooccurence: 04
#> 
#> No parent found!
#> 4
#> 
#> ####processing: ASV1082 #####4
#> 
#> ---hits: ASV1027
#> ---hits: ASV962
#> ---hits: ASV2
#> ---hits: ASV1247
#> ---hits: ASV5764
#> 
#> ---potential parent: ASV2
#> ---potential parent: ASV576
#> ---potential parent: ASV1247
#> ---potential parent: ASV962
#> ---potential parent: ASV10274
#> 
#> ------checking: ASV24
#> 
#> ------relative cooccurence: 14
#>  which is sufficient!4
#> 
#> ------min avg abundance: 70.254
#>  which is OK!4
#> 
#> SETTING ASV1082 to be an ERROR of ASV2
#> 4
#> 
#> ------checking: ASV5764
#> 
#> ------checking: ASV12474
#> 
#> ------checking: ASV9624
#> 
#> ------checking: ASV10274
#> 
#> ####processing: ASV1109 #####4
#> 
#> ---hits: ASV1628
#> ---hits: ASV12334
#> 
#> ---potential parent: ASV1233
#> ---potential parent: ASV16284
#> 
#> ------checking: ASV12334
#> 
#> ------relative cooccurence: 14
#>  which is sufficient!4
#> 
#> ------min avg abundance: 8.54
#>  which is OK!4
#> 
#> SETTING ASV1109 to be an ERROR of ASV1233
#> 4
#> 
#> ------checking: ASV16284
#>   |                                                          |=======================================           |  79%
#> 
#> ####processing: ASV1144 #####4
#> 
#> ---hits: ASV1263
#> ---hits: ASV839
#> ---hits: ASV81
#> ---hits: ASV272
#> ---hits: ASV959
#> ---hits: ASV339
#> ---hits: ASV383
#> ---hits: ASV11524
#> 
#> ---potential parent: ASV1152
#> ---potential parent: ASV81
#> ---potential parent: ASV339
#> ---potential parent: ASV383
#> ---potential parent: ASV839
#> ---potential parent: ASV272
#> ---potential parent: ASV959
#> ---potential parent: ASV12634
#> 
#> ------checking: ASV11524
#> 
#> ------relative cooccurence: 04
#> 
#> ------checking: ASV814
#> 
#> ------relative cooccurence: 04
#> 
#> ------checking: ASV3394
#> 
#> ------relative cooccurence: 04
#> 
#> ------checking: ASV3834
#> 
#> ------relative cooccurence: 04
#> 
#> ------checking: ASV8394
#> 
#> ------relative cooccurence: 04
#> 
#> ------checking: ASV2724
#> 
#> ------relative cooccurence: 04
#> 
#> ------checking: ASV9594
#> 
#> ------relative cooccurence: 04
#> 
#> ------checking: ASV12634
#> 
#> ------relative cooccurence: 14
#>  which is sufficient!4
#> 
#> ------min avg abundance: 0.254
#> 
#> No parent found!
#> 4
#> 
#> ####processing: ASV1419 #####4
#> 
#> ---hits: ASV43
#> ---hits: ASV999
#> ---hits: ASV11174
#> 
#> ---potential parent: ASV43
#> ---potential parent: ASV1117
#> ---potential parent: ASV9994
#> 
#> ------checking: ASV434
#> 
#> ------relative cooccurence: 04
#> 
#> ------checking: ASV11174
#> 
#> ------relative cooccurence: 04
#> 
#> ------checking: ASV9994
#> 
#> ------relative cooccurence: 04
#> 
#> No parent found!
#> 4
#>   |                                                          |========================================          |  79%
#> 
#> ####processing: ASV1467 #####4
#> 
#> ---hits: ASV1078
#> ---hits: ASV178
#> ---hits: ASV1131
#> ---hits: ASV19
#> ---hits: ASV594
#> ---hits: ASV333
#> ---hits: ASV1010
#> ---hits: ASV858
#> ---hits: ASV986
#> ---hits: ASV7314
#> 
#> ---potential parent: ASV19
#> ---potential parent: ASV178
#> ---potential parent: ASV333
#> ---potential parent: ASV1078
#> ---potential parent: ASV594
#> ---potential parent: ASV1131
#> ---potential parent: ASV858
#> ---potential parent: ASV1010
#> ---potential parent: ASV986
#> ---potential parent: ASV7314
#> 
#> ------checking: ASV194
#> 
#> ------relative cooccurence: 14
#>  which is sufficient!4
#> 
#> ------min avg abundance: 743.254
#>  which is OK!4
#> 
#> SETTING ASV1467 to be an ERROR of ASV19
#> 4
#> 
#> ------checking: ASV1784
#> 
#> ------checking: ASV3334
#> 
#> ------checking: ASV10784
#> 
#> ------checking: ASV5944
#> 
#> ------checking: ASV11314
#> 
#> ------checking: ASV8584
#> 
#> ------checking: ASV10104
#> 
#> ------checking: ASV9864
#> 
#> ------checking: ASV7314
#> 
#> ####processing: ASV42 #####4
#> 
#> ---hits: ASV254
#> 
#> ---potential parent: ASV254
#> 
#> ------checking: ASV254
#> 
#> ------relative cooccurence: 04
#> 
#> No parent found!
#> 4
#>   |                                                          |========================================          |  80%
#> 
#> ####processing: ASV104 #####4
#> 
#> ---hits: 4
#> 
#> ---potential parent: 4
#> 
#> No parent found!
#> 4
#> 
#> ####processing: ASV203 #####4
#> 
#> ---hits: ASV580
#> ---hits: ASV64
#> ---hits: ASV168
#> ---hits: ASV2494
#> 
#> ---potential parent: ASV64
#> ---potential parent: ASV168
#> ---potential parent: ASV580
#> ---potential parent: ASV2494
#> 
#> ------checking: ASV644
#> 
#> ------relative cooccurence: 14
#>  which is sufficient!4
#> 
#> ------min avg abundance: 27294
#>  which is OK!4
#> 
#> SETTING ASV203 to be an ERROR of ASV64
#> 4
#> 
#> ------checking: ASV1684
#> 
#> ------checking: ASV5804
#> 
#> ------checking: ASV2494
#> 
#> ####processing: ASV491 #####4
#> 
#> ---hits: ASV210
#> ---hits: ASV12234
#> 
#> ---potential parent: ASV1223
#> ---potential parent: ASV2104
#> 
#> ------checking: ASV12234
#> 
#> ------relative cooccurence: 14
#>  which is sufficient!4
#> 
#> ------min avg abundance: 0.3333333333333334
#> 
#> ------checking: ASV2104
#> 
#> ------relative cooccurence: 14
#>  which is sufficient!4
#> 
#> ------min avg abundance: 1.333333333333334
#>  which is OK!4
#> 
#> SETTING ASV491 to be an ERROR of ASV210
#> 4
#> 
#> ####processing: ASV517 #####4
#> 
#> ---hits: ASV13794
#> 
#> ---potential parent: ASV13794
#> 
#> ------checking: ASV13794
#> 
#> ------relative cooccurence: 04
#> 
#> No parent found!
#> 4
#>   |                                                          |========================================          |  81%
#> 
#> ####processing: ASV598 #####4
#> 
#> ---hits: 4
#> 
#> ---potential parent: 4
#> 
#> No parent found!
#> 4
#> 
#> ####processing: ASV626 #####4
#> 
#> ---hits: ASV3574
#> 
#> ---potential parent: ASV3574
#> 
#> ------checking: ASV3574
#> 
#> ------relative cooccurence: 14
#>  which is sufficient!4
#> 
#> ------min avg abundance: 3.666666666666674
#>  which is OK!4
#> 
#> SETTING ASV626 to be an ERROR of ASV357
#> 4
#>   |                                                          |=========================================         |  81%
#> 
#> ####processing: ASV650 #####4
#> 
#> ---hits: ASV1332
#> ---hits: ASV1303
#> ---hits: ASV1230
#> ---hits: ASV378
#> ---hits: ASV504
#> ---hits: ASV860
#> ---hits: ASV772
#> ---hits: ASV746
#> ---hits: ASV1409
#> ---hits: ASV8324
#> 
#> ---potential parent: ASV378
#> ---potential parent: ASV746
#> ---potential parent: ASV1409
#> ---potential parent: ASV1230
#> ---potential parent: ASV772
#> ---potential parent: ASV1332
#> ---potential parent: ASV504
#> ---potential parent: ASV832
#> ---potential parent: ASV860
#> ---potential parent: ASV13034
#> 
#> ------checking: ASV3784
#> 
#> ------relative cooccurence: 14
#>  which is sufficient!4
#> 
#> ------min avg abundance: 10.33333333333334
#>  which is OK!4
#> 
#> SETTING ASV650 to be an ERROR of ASV378
#> 4
#> 
#> ------checking: ASV7464
#> 
#> ------checking: ASV14094
#> 
#> ------checking: ASV12304
#> 
#> ------checking: ASV7724
#> 
#> ------checking: ASV13324
#> 
#> ------checking: ASV5044
#> 
#> ------checking: ASV8324
#> 
#> ------checking: ASV8604
#> 
#> ------checking: ASV13034
#> 
#> ####processing: ASV790 #####4
#> 
#> ---hits: ASV420
#> ---hits: ASV505
#> ---hits: ASV509
#> ---hits: ASV641
#> ---hits: ASV662
#> ---hits: ASV1052
#> ---hits: ASV1410
#> ---hits: ASV1427
#> ---hits: ASV724
#> ---hits: ASV11984
#> 
#> ---potential parent: ASV662
#> ---potential parent: ASV505
#> ---potential parent: ASV509
#> ---potential parent: ASV724
#> ---potential parent: ASV1052
#> ---potential parent: ASV420
#> ---potential parent: ASV641
#> ---potential parent: ASV1198
#> ---potential parent: ASV1427
#> ---potential parent: ASV14104
#> 
#> ------checking: ASV6624
#> 
#> ------relative cooccurence: 14
#>  which is sufficient!4
#> 
#> ------min avg abundance: 0.3333333333333334
#> 
#> ------checking: ASV5054
#> 
#> ------relative cooccurence: 14
#>  which is sufficient!4
#> 
#> ------min avg abundance: 13.33333333333334
#>  which is OK!4
#> 
#> SETTING ASV790 to be an ERROR of ASV505
#> 4
#> 
#> ------checking: ASV5094
#> 
#> ------checking: ASV7244
#> 
#> ------checking: ASV10524
#> 
#> ------checking: ASV4204
#> 
#> ------checking: ASV6414
#> 
#> ------checking: ASV11984
#> 
#> ------checking: ASV14274
#> 
#> ------checking: ASV14104
#>   |                                                          |=========================================         |  82%
#> 
#> ####processing: ASV839 #####4
#> 
#> ---hits: ASV81
#> ---hits: ASV272
#> ---hits: ASV339
#> ---hits: ASV959
#> ---hits: ASV383
#> ---hits: ASV1152
#> ---hits: ASV1144
#> ---hits: ASV12634
#> 
#> ---potential parent: ASV1152
#> ---potential parent: ASV81
#> ---potential parent: ASV339
#> ---potential parent: ASV383
#> ---potential parent: ASV1144
#> ---potential parent: ASV272
#> ---potential parent: ASV959
#> ---potential parent: ASV12634
#> 
#> ------checking: ASV11524
#> 
#> ------relative cooccurence: 04
#> 
#> ------checking: ASV814
#> 
#> ------relative cooccurence: 14
#>  which is sufficient!4
#> 
#> ------min avg abundance: 0.6666666666666674
#> 
#> ------checking: ASV3394
#> 
#> ------relative cooccurence: 14
#>  which is sufficient!4
#> 
#> ------min avg abundance: 0.3333333333333334
#> 
#> ------checking: ASV3834
#> 
#> ------relative cooccurence: 04
#> 
#> ------checking: ASV11444
#> 
#> ------relative cooccurence: 04
#> 
#> ------checking: ASV2724
#> 
#> ------relative cooccurence: 14
#>  which is sufficient!4
#> 
#> ------min avg abundance: 0.6666666666666674
#> 
#> ------checking: ASV9594
#> 
#> ------relative cooccurence: 04
#> 
#> ------checking: ASV12634
#> 
#> ------relative cooccurence: 04
#> 
#> No parent found!
#> 4
#> 
#> ####processing: ASV1014 #####4
#> 
#> ---hits: ASV1386
#> ---hits: ASV9694
#> 
#> ---potential parent: ASV969
#> ---potential parent: ASV13864
#> 
#> ------checking: ASV9694
#> 
#> ------relative cooccurence: 14
#>  which is sufficient!4
#> 
#> ------min avg abundance: 30.66666666666674
#>  which is OK!4
#> 
#> SETTING ASV1014 to be an ERROR of ASV969
#> 4
#> 
#> ------checking: ASV13864
#> 
#> ####processing: ASV1017 #####4
#> 
#> ---hits: 4
#> 
#> ---potential parent: 4
#> 
#> No parent found!
#> 4
#> 
#> ####processing: ASV1030 #####4
#> 
#> ---hits: ASV1303
#> ---hits: ASV1458
#> ---hits: ASV1562
#> ---hits: ASV746
#> ---hits: ASV1332
#> ---hits: ASV693
#> ---hits: ASV1230
#> ---hits: ASV1409
#> ---hits: ASV378
#> ---hits: ASV7724
#> 
#> ---potential parent: ASV378
#> ---potential parent: ASV746
#> ---potential parent: ASV1409
#> ---potential parent: ASV1230
#> ---potential parent: ASV1458
#> ---potential parent: ASV772
#> ---potential parent: ASV693
#> ---potential parent: ASV1332
#> ---potential parent: ASV1303
#> ---potential parent: ASV15624
#> 
#> ------checking: ASV3784
#> 
#> ------relative cooccurence: 04
#> 
#> ------checking: ASV7464
#> 
#> ------relative cooccurence: 04
#> 
#> ------checking: ASV14094
#> 
#> ------relative cooccurence: 04
#> 
#> ------checking: ASV12304
#> 
#> ------relative cooccurence: 04
#> 
#> ------checking: ASV14584
#> 
#> ------relative cooccurence: 04
#> 
#> ------checking: ASV7724
#> 
#> ------relative cooccurence: 04
#> 
#> ------checking: ASV6934
#> 
#> ------relative cooccurence: 04
#> 
#> ------checking: ASV13324
#> 
#> ------relative cooccurence: 04
#> 
#> ------checking: ASV13034
#> 
#> ------relative cooccurence: 04
#> 
#> ------checking: ASV15624
#> 
#> ------relative cooccurence: 04
#> 
#> No parent found!
#> 4
#> 
#> ####processing: ASV1085 #####4
#> 
#> ---hits: 4
#> 
#> ---potential parent: 4
#> 
#> No parent found!
#> 4
#>   |                                                          |=========================================         |  83%
#> 
#> ####processing: ASV1167 #####4
#> 
#> ---hits: ASV45
#> ---hits: ASV137
#> ---hits: ASV1067
#> ---hits: ASV996
#> ---hits: ASV711
#> ---hits: ASV760
#> ---hits: ASV940
#> ---hits: ASV1484
#> ---hits: ASV1468
#> ---hits: ASV9734
#> 
#> ---potential parent: ASV973
#> ---potential parent: ASV45
#> ---potential parent: ASV711
#> ---potential parent: ASV1067
#> ---potential parent: ASV137
#> ---potential parent: ASV1468
#> ---potential parent: ASV996
#> ---potential parent: ASV760
#> ---potential parent: ASV940
#> ---potential parent: ASV14844
#> 
#> ------checking: ASV9734
#> 
#> ------relative cooccurence: 04
#> 
#> ------checking: ASV454
#> 
#> ------relative cooccurence: 04
#> 
#> ------checking: ASV7114
#> 
#> ------relative cooccurence: 04
#> 
#> ------checking: ASV10674
#> 
#> ------relative cooccurence: 04
#> 
#> ------checking: ASV1374
#> 
#> ------relative cooccurence: 04
#> 
#> ------checking: ASV14684
#> 
#> ------relative cooccurence: 04
#> 
#> ------checking: ASV9964
#> 
#> ------relative cooccurence: 04
#> 
#> ------checking: ASV7604
#> 
#> ------relative cooccurence: 04
#> 
#> ------checking: ASV9404
#> 
#> ------relative cooccurence: 04
#> 
#> ------checking: ASV14844
#> 
#> ------relative cooccurence: 04
#> 
#> No parent found!
#> 4
#> 
#> ####processing: ASV1625 #####4
#> 
#> ---hits: ASV33
#> ---hits: ASV78
#> ---hits: ASV941
#> ---hits: ASV13134
#> 
#> ---potential parent: ASV78
#> ---potential parent: ASV1313
#> ---potential parent: ASV33
#> ---potential parent: ASV9414
#> 
#> ------checking: ASV784
#> 
#> ------relative cooccurence: 14
#>  which is sufficient!4
#> 
#> ------min avg abundance: 2001.666666666674
#>  which is OK!4
#> 
#> SETTING ASV1625 to be an ERROR of ASV78
#> 4
#> 
#> ------checking: ASV13134
#> 
#> ------checking: ASV334
#> 
#> ------checking: ASV9414
#>   |                                                          |==========================================        |  83%
#> 
#> ####processing: ASV25 #####4
#> 
#> ---hits: ASV424
#> 
#> ---potential parent: ASV424
#> 
#> ------checking: ASV424
#> 
#> ------relative cooccurence: 04
#> 
#> No parent found!
#> 4
#> 
#> ####processing: ASV249 #####4
#> 
#> ---hits: ASV64
#> ---hits: ASV168
#> ---hits: ASV580
#> ---hits: ASV2034
#> 
#> ---potential parent: ASV64
#> ---potential parent: ASV168
#> ---potential parent: ASV580
#> ---potential parent: ASV2034
#> 
#> ------checking: ASV644
#> 
#> ------relative cooccurence: 14
#>  which is sufficient!4
#> 
#> ------min avg abundance: 4093.54
#>  which is OK!4
#> 
#> SETTING ASV249 to be an ERROR of ASV64
#> 4
#> 
#> ------checking: ASV1684
#> 
#> ------checking: ASV5804
#> 
#> ------checking: ASV2034
#>   |                                                          |==========================================        |  84%
#> 
#> ####processing: ASV272 #####4
#> 
#> ---hits: ASV81
#> ---hits: ASV339
#> ---hits: ASV383
#> ---hits: ASV959
#> ---hits: ASV839
#> ---hits: ASV1152
#> ---hits: ASV1263
#> ---hits: ASV11444
#> 
#> ---potential parent: ASV1152
#> ---potential parent: ASV81
#> ---potential parent: ASV339
#> ---potential parent: ASV383
#> ---potential parent: ASV1144
#> ---potential parent: ASV839
#> ---potential parent: ASV959
#> ---potential parent: ASV12634
#> 
#> ------checking: ASV11524
#> 
#> ------relative cooccurence: 04
#> 
#> ------checking: ASV814
#> 
#> ------relative cooccurence: 14
#>  which is sufficient!4
#> 
#> ------min avg abundance: 14
#> 
#> ------checking: ASV3394
#> 
#> ------relative cooccurence: 14
#>  which is sufficient!4
#> 
#> ------min avg abundance: 0.54
#> 
#> ------checking: ASV3834
#> 
#> ------relative cooccurence: 04
#> 
#> ------checking: ASV11444
#> 
#> ------relative cooccurence: 04
#> 
#> ------checking: ASV8394
#> 
#> ------relative cooccurence: 14
#>  which is sufficient!4
#> 
#> ------min avg abundance: 1.54
#>  which is OK!4
#> 
#> SETTING ASV272 to be an ERROR of ASV839
#> 4
#> 
#> ------checking: ASV9594
#> 
#> ------checking: ASV12634
#> 
#> ####processing: ASV320 #####4
#> 
#> ---hits: 4
#> 
#> ---potential parent: 4
#> 
#> No parent found!
#> 4
#> 
#> ####processing: ASV334 #####4
#> 
#> ---hits: 4
#> 
#> ---potential parent: 4
#> 
#> No parent found!
#> 4
#> 
#> ####processing: ASV488 #####4
#> 
#> ---hits: ASV248
#> ---hits: ASV840
#> ---hits: ASV404
#> ---hits: ASV377
#> ---hits: ASV440
#> ---hits: ASV5454
#> 
#> ---potential parent: ASV248
#> ---potential parent: ASV377
#> ---potential parent: ASV404
#> ---potential parent: ASV545
#> ---potential parent: ASV440
#> ---potential parent: ASV8404
#> 
#> ------checking: ASV2484
#> 
#> ------relative cooccurence: 14
#>  which is sufficient!4
#> 
#> ------min avg abundance: 0.54
#> 
#> ------checking: ASV3774
#> 
#> ------relative cooccurence: 04
#> 
#> ------checking: ASV4044
#> 
#> ------relative cooccurence: 14
#>  which is sufficient!4
#> 
#> ------min avg abundance: 19.54
#>  which is OK!4
#> 
#> SETTING ASV488 to be an ERROR of ASV404
#> 4
#> 
#> ------checking: ASV5454
#> 
#> ------checking: ASV4404
#> 
#> ------checking: ASV8404
#>   |                                                          |==========================================        |  85%
#> 
#> ####processing: ASV515 #####4
#> 
#> ---hits: 4
#> 
#> ---potential parent: 4
#> 
#> No parent found!
#> 4
#> 
#> ####processing: ASV549 #####4
#> 
#> ---hits: 4
#> 
#> ---potential parent: 4
#> 
#> No parent found!
#> 4
#>   |                                                          |===========================================       |  85%
#> 
#> ####processing: ASV570 #####4
#> 
#> ---hits: ASV854
#> ---hits: ASV2554
#> 
#> ---potential parent: ASV255
#> ---potential parent: ASV8544
#> 
#> ------checking: ASV2554
#> 
#> ------relative cooccurence: 04
#> 
#> ------checking: ASV8544
#> 
#> ------relative cooccurence: 04
#> 
#> No parent found!
#> 4
#> 
#> ####processing: ASV604 #####4
#> 
#> ---hits: ASV891
#> ---hits: ASV1311
#> ---hits: ASV13964
#> 
#> ---potential parent: ASV1396
#> ---potential parent: ASV1311
#> ---potential parent: ASV8914
#> 
#> ------checking: ASV13964
#> 
#> ------relative cooccurence: 04
#> 
#> ------checking: ASV13114
#> 
#> ------relative cooccurence: 14
#>  which is sufficient!4
#> 
#> ------min avg abundance: 64
#>  which is OK!4
#> 
#> SETTING ASV604 to be an ERROR of ASV1311
#> 4
#> 
#> ------checking: ASV8914
#>   |                                                          |===========================================       |  86%
#> 
#> ####processing: ASV616 #####4
#> 
#> ---hits: ASV1444
#> 
#> ---potential parent: ASV1444
#> 
#> ------checking: ASV1444
#> 
#> ------relative cooccurence: 14
#>  which is sufficient!4
#> 
#> ------min avg abundance: 14624
#>  which is OK!4
#> 
#> SETTING ASV616 to be an ERROR of ASV144
#> 4
#> 
#> ####processing: ASV736 #####4
#> 
#> ---hits: ASV979
#> ---hits: ASV24
#> ---hits: ASV266
#> ---hits: ASV858
#> ---hits: ASV731
#> ---hits: ASV986
#> ---hits: ASV19
#> ---hits: ASV333
#> ---hits: ASV1010
#> ---hits: ASV1784
#> 
#> ---potential parent: ASV19
#> ---potential parent: ASV178
#> ---potential parent: ASV24
#> ---potential parent: ASV266
#> ---potential parent: ASV333
#> ---potential parent: ASV858
#> ---potential parent: ASV1010
#> ---potential parent: ASV979
#> ---potential parent: ASV986
#> ---potential parent: ASV7314
#> 
#> ------checking: ASV194
#> 
#> ------relative cooccurence: 14
#>  which is sufficient!4
#> 
#> ------min avg abundance: 1486.54
#>  which is OK!4
#> 
#> SETTING ASV736 to be an ERROR of ASV19
#> 4
#> 
#> ------checking: ASV1784
#> 
#> ------checking: ASV244
#> 
#> ------checking: ASV2664
#> 
#> ------checking: ASV3334
#> 
#> ------checking: ASV8584
#> 
#> ------checking: ASV10104
#> 
#> ------checking: ASV9794
#> 
#> ------checking: ASV9864
#> 
#> ------checking: ASV7314
#> 
#> ####processing: ASV748 #####4
#> 
#> ---hits: ASV171
#> ---hits: ASV7534
#> 
#> ---potential parent: ASV753
#> ---potential parent: ASV1714
#> 
#> ------checking: ASV7534
#> 
#> ------relative cooccurence: 14
#>  which is sufficient!4
#> 
#> ------min avg abundance: 24
#>  which is OK!4
#> 
#> SETTING ASV748 to be an ERROR of ASV753
#> 4
#> 
#> ------checking: ASV1714
#> 
#> ####processing: ASV760 #####4
#> 
#> ---hits: ASV137
#> ---hits: ASV45
#> ---hits: ASV711
#> ---hits: ASV996
#> ---hits: ASV1067
#> ---hits: ASV1167
#> ---hits: ASV1484
#> ---hits: ASV940
#> ---hits: ASV1468
#> ---hits: ASV9734
#> 
#> ---potential parent: ASV973
#> ---potential parent: ASV45
#> ---potential parent: ASV711
#> ---potential parent: ASV1067
#> ---potential parent: ASV137
#> ---potential parent: ASV1468
#> ---potential parent: ASV996
#> ---potential parent: ASV1167
#> ---potential parent: ASV940
#> ---potential parent: ASV14844
#> 
#> ------checking: ASV9734
#> 
#> ------relative cooccurence: 04
#> 
#> ------checking: ASV454
#> 
#> ------relative cooccurence: 14
#>  which is sufficient!4
#> 
#> ------min avg abundance: 4184
#>  which is OK!4
#> 
#> SETTING ASV760 to be an ERROR of ASV45
#> 4
#> 
#> ------checking: ASV7114
#> 
#> ------checking: ASV10674
#> 
#> ------checking: ASV1374
#> 
#> ------checking: ASV14684
#> 
#> ------checking: ASV9964
#> 
#> ------checking: ASV11674
#> 
#> ------checking: ASV9404
#> 
#> ------checking: ASV14844
#>   |                                                          |===========================================       |  87%
#> 
#> ####processing: ASV823 #####4
#> 
#> ---hits: 4
#> 
#> ---potential parent: 4
#> 
#> No parent found!
#> 4
#> 
#> ####processing: ASV828 #####4
#> 
#> ---hits: ASV1319
#> ---hits: ASV906
#> ---hits: ASV1278
#> ---hits: ASV1074
#> ---hits: ASV1007
#> ---hits: ASV1359
#> ---hits: ASV1192
#> ---hits: ASV950
#> ---hits: ASV1224
#> ---hits: ASV7564
#> 
#> ---potential parent: ASV756
#> ---potential parent: ASV1192
#> ---potential parent: ASV1074
#> ---potential parent: ASV906
#> ---potential parent: ASV1007
#> ---potential parent: ASV1359
#> ---potential parent: ASV1278
#> ---potential parent: ASV950
#> ---potential parent: ASV1319
#> ---potential parent: ASV12244
#> 
#> ------checking: ASV7564
#> 
#> ------relative cooccurence: 04
#> 
#> ------checking: ASV11924
#> 
#> ------relative cooccurence: 14
#>  which is sufficient!4
#> 
#> ------min avg abundance: 24
#>  which is OK!4
#> 
#> SETTING ASV828 to be an ERROR of ASV1192
#> 4
#> 
#> ------checking: ASV10744
#> 
#> ------checking: ASV9064
#> 
#> ------checking: ASV10074
#> 
#> ------checking: ASV13594
#> 
#> ------checking: ASV12784
#> 
#> ------checking: ASV9504
#> 
#> ------checking: ASV13194
#> 
#> ------checking: ASV12244
#>   |                                                          |============================================      |  87%
#> 
#> ####processing: ASV829 #####4
#> 
#> ---hits: 4
#> 
#> ---potential parent: 4
#> 
#> No parent found!
#> 4
#> 
#> ####processing: ASV891 #####4
#> 
#> ---hits: ASV604
#> ---hits: ASV13114
#> 
#> ---potential parent: ASV1311
#> ---potential parent: ASV6044
#> 
#> ------checking: ASV13114
#> 
#> ------relative cooccurence: 04
#> 
#> ------checking: ASV6044
#> 
#> ------relative cooccurence: 04
#> 
#> No parent found!
#> 4
#>   |                                                          |============================================      |  88%
#> 
#> ####processing: ASV914 #####4
#> 
#> ---hits: ASV930
#> ---hits: ASV69
#> ---hits: ASV405
#> ---hits: ASV358
#> ---hits: ASV221
#> ---hits: ASV833
#> ---hits: ASV1420
#> ---hits: ASV12184
#> 
#> ---potential parent: ASV69
#> ---potential parent: ASV358
#> ---potential parent: ASV1218
#> ---potential parent: ASV221
#> ---potential parent: ASV833
#> ---potential parent: ASV1420
#> ---potential parent: ASV930
#> ---potential parent: ASV4054
#> 
#> ------checking: ASV694
#> 
#> ------relative cooccurence: 14
#>  which is sufficient!4
#> 
#> ------min avg abundance: 14
#> 
#> ------checking: ASV3584
#> 
#> ------relative cooccurence: 04
#> 
#> ------checking: ASV12184
#> 
#> ------relative cooccurence: 04
#> 
#> ------checking: ASV2214
#> 
#> ------relative cooccurence: 04
#> 
#> ------checking: ASV8334
#> 
#> ------relative cooccurence: 04
#> 
#> ------checking: ASV14204
#> 
#> ------relative cooccurence: 04
#> 
#> ------checking: ASV9304
#> 
#> ------relative cooccurence: 04
#> 
#> ------checking: ASV4054
#> 
#> ------relative cooccurence: 04
#> 
#> No parent found!
#> 4
#> 
#> ####processing: ASV961 #####4
#> 
#> ---hits: ASV98
#> ---hits: ASV292
#> ---hits: ASV368
#> ---hits: ASV566
#> ---hits: ASV443
#> ---hits: ASV1111
#> ---hits: ASV11764
#> 
#> ---potential parent: ASV292
#> ---potential parent: ASV443
#> ---potential parent: ASV98
#> ---potential parent: ASV368
#> ---potential parent: ASV566
#> ---potential parent: ASV1176
#> ---potential parent: ASV11114
#> 
#> ------checking: ASV2924
#> 
#> ------relative cooccurence: 14
#>  which is sufficient!4
#> 
#> ------min avg abundance: 224
#>  which is OK!4
#> 
#> SETTING ASV961 to be an ERROR of ASV292
#> 4
#> 
#> ------checking: ASV4434
#> 
#> ------checking: ASV984
#> 
#> ------checking: ASV3684
#> 
#> ------checking: ASV5664
#> 
#> ------checking: ASV11764
#> 
#> ------checking: ASV11114
#> 
#> ####processing: ASV979 #####4
#> 
#> ---hits: ASV736
#> ---hits: ASV24
#> ---hits: ASV858
#> ---hits: ASV266
#> ---hits: ASV986
#> ---hits: ASV731
#> ---hits: ASV19
#> ---hits: ASV333
#> ---hits: ASV178
#> ---hits: ASV10104
#> 
#> ---potential parent: ASV19
#> ---potential parent: ASV178
#> ---potential parent: ASV24
#> ---potential parent: ASV266
#> ---potential parent: ASV333
#> ---potential parent: ASV858
#> ---potential parent: ASV1010
#> ---potential parent: ASV736
#> ---potential parent: ASV986
#> ---potential parent: ASV7314
#> 
#> ------checking: ASV194
#> 
#> ------relative cooccurence: 14
#>  which is sufficient!4
#> 
#> ------min avg abundance: 1486.54
#>  which is OK!4
#> 
#> SETTING ASV979 to be an ERROR of ASV19
#> 4
#> 
#> ------checking: ASV1784
#> 
#> ------checking: ASV244
#> 
#> ------checking: ASV2664
#> 
#> ------checking: ASV3334
#> 
#> ------checking: ASV8584
#> 
#> ------checking: ASV10104
#> 
#> ------checking: ASV7364
#> 
#> ------checking: ASV9864
#> 
#> ------checking: ASV7314
#> 
#> ####processing: ASV986 #####4
#> 
#> ---hits: ASV1010
#> ---hits: ASV24
#> ---hits: ASV736
#> ---hits: ASV266
#> ---hits: ASV979
#> ---hits: ASV19
#> ---hits: ASV333
#> ---hits: ASV731
#> ---hits: ASV178
#> ---hits: ASV5944
#> 
#> ---potential parent: ASV19
#> ---potential parent: ASV178
#> ---potential parent: ASV24
#> ---potential parent: ASV266
#> ---potential parent: ASV333
#> ---potential parent: ASV594
#> ---potential parent: ASV1010
#> ---potential parent: ASV736
#> ---potential parent: ASV979
#> ---potential parent: ASV7314
#> 
#> ------checking: ASV194
#> 
#> ------relative cooccurence: 14
#>  which is sufficient!4
#> 
#> ------min avg abundance: 1486.54
#>  which is OK!4
#> 
#> SETTING ASV986 to be an ERROR of ASV19
#> 4
#> 
#> ------checking: ASV1784
#> 
#> ------checking: ASV244
#> 
#> ------checking: ASV2664
#> 
#> ------checking: ASV3334
#> 
#> ------checking: ASV5944
#> 
#> ------checking: ASV10104
#> 
#> ------checking: ASV7364
#> 
#> ------checking: ASV9794
#> 
#> ------checking: ASV7314
#>   |                                                          |============================================      |  89%
#> 
#> ####processing: ASV1133 #####4
#> 
#> ---hits: ASV526
#> ---hits: ASV15524
#> 
#> ---potential parent: ASV1552
#> ---potential parent: ASV5264
#> 
#> ------checking: ASV15524
#> 
#> ------relative cooccurence: 14
#>  which is sufficient!4
#> 
#> ------min avg abundance: 64
#>  which is OK!4
#> 
#> SETTING ASV1133 to be an ERROR of ASV1552
#> 4
#> 
#> ------checking: ASV5264
#> 
#> ####processing: ASV1270 #####4
#> 
#> ---hits: 4
#> 
#> ---potential parent: 4
#> 
#> No parent found!
#> 4
#>   |                                                          |=============================================     |  89%
#> 
#> ####processing: ASV1422 #####4
#> 
#> ---hits: 4
#> 
#> ---potential parent: 4
#> 
#> No parent found!
#> 4
#> 
#> ####processing: ASV1523 #####4
#> 
#> ---hits: 4
#> 
#> ---potential parent: 4
#> 
#> No parent found!
#> 4
#>   |                                                          |=============================================     |  90%
#> 
#> ####processing: ASV1576 #####4
#> 
#> ---hits: ASV5864
#> 
#> ---potential parent: ASV5864
#> 
#> ------checking: ASV5864
#> 
#> ------relative cooccurence: 14
#>  which is sufficient!4
#> 
#> ------min avg abundance: 0.54
#> 
#> No parent found!
#> 4
#> 
#> ####processing: ASV1603 #####4
#> 
#> ---hits: ASV1036
#> ---hits: ASV1020
#> ---hits: ASV60
#> ---hits: ASV796
#> ---hits: ASV984
#> ---hits: ASV113
#> ---hits: ASV415
#> ---hits: ASV1380
#> ---hits: ASV1011
#> ---hits: ASV10664
#> 
#> ---potential parent: ASV113
#> ---potential parent: ASV60
#> ---potential parent: ASV1036
#> ---potential parent: ASV984
#> ---potential parent: ASV1011
#> ---potential parent: ASV415
#> ---potential parent: ASV796
#> ---potential parent: ASV1020
#> ---potential parent: ASV1066
#> ---potential parent: ASV13804
#> 
#> ------checking: ASV1134
#> 
#> ------relative cooccurence: 14
#>  which is sufficient!4
#> 
#> ------min avg abundance: 24
#>  which is OK!4
#> 
#> SETTING ASV1603 to be an ERROR of ASV113
#> 4
#> 
#> ------checking: ASV604
#> 
#> ------checking: ASV10364
#> 
#> ------checking: ASV9844
#> 
#> ------checking: ASV10114
#> 
#> ------checking: ASV4154
#> 
#> ------checking: ASV7964
#> 
#> ------checking: ASV10204
#> 
#> ------checking: ASV10664
#> 
#> ------checking: ASV13804
#> 
#> ####processing: ASV1726 #####4
#> 
#> ---hits: 4
#> 
#> ---potential parent: 4
#> 
#> No parent found!
#> 4
#> 
#> ####processing: ASV67 #####4
#> 
#> ---hits: ASV12
#> ---hits: ASV107
#> ---hits: ASV1542
#> ---hits: ASV624
#> ---hits: ASV1262
#> ---hits: ASV847
#> ---hits: ASV1054
#> 
#> ---potential parent: ASV12
#> ---potential parent: ASV107
#> ---potential parent: ASV624
#> ---potential parent: ASV105
#> ---potential parent: ASV847
#> ---potential parent: ASV1542
#> ---potential parent: ASV12624
#> 
#> ------checking: ASV124
#> 
#> ------relative cooccurence: 14
#>  which is sufficient!4
#> 
#> ------min avg abundance: 31864
#>  which is OK!4
#> 
#> SETTING ASV67 to be an ERROR of ASV12
#> 4
#> 
#> ------checking: ASV1074
#> 
#> ------checking: ASV6244
#> 
#> ------checking: ASV1054
#> 
#> ------checking: ASV8474
#> 
#> ------checking: ASV15424
#> 
#> ------checking: ASV12624
#>   |                                                          |=============================================     |  91%
#> 
#> ####processing: ASV87 #####4
#> 
#> ---hits: 4
#> 
#> ---potential parent: 4
#> 
#> No parent found!
#> 4
#> 
#> ####processing: ASV116 #####4
#> 
#> ---hits: ASV1107
#> ---hits: ASV1526
#> ---hits: ASV254
#> ---hits: ASV500
#> ---hits: ASV989
#> ---hits: ASV38
#> ---hits: ASV816
#> ---hits: ASV466
#> ---hits: ASV8534
#> 
#> ---potential parent: ASV38
#> ---potential parent: ASV989
#> ---potential parent: ASV816
#> ---potential parent: ASV500
#> ---potential parent: ASV466
#> ---potential parent: ASV254
#> ---potential parent: ASV1526
#> ---potential parent: ASV853
#> ---potential parent: ASV11074
#> 
#> ------checking: ASV384
#> 
#> ------relative cooccurence: 14
#>  which is sufficient!4
#> 
#> ------min avg abundance: 94
#>  which is OK!4
#> 
#> SETTING ASV116 to be an ERROR of ASV38
#> 4
#> 
#> ------checking: ASV9894
#> 
#> ------checking: ASV8164
#> 
#> ------checking: ASV5004
#> 
#> ------checking: ASV4664
#> 
#> ------checking: ASV2544
#> 
#> ------checking: ASV15264
#> 
#> ------checking: ASV8534
#> 
#> ------checking: ASV11074
#> 
#> ####processing: ASV171 #####4
#> 
#> ---hits: ASV753
#> ---hits: ASV7484
#> 
#> ---potential parent: ASV753
#> ---potential parent: ASV7484
#> 
#> ------checking: ASV7534
#> 
#> ------relative cooccurence: 04
#> 
#> ------checking: ASV7484
#> 
#> ------relative cooccurence: 04
#> 
#> No parent found!
#> 4
#>   |                                                          |==============================================    |  91%
#> 
#> ####processing: ASV193 #####4
#> 
#> ---hits: ASV6314
#> 
#> ---potential parent: ASV6314
#> 
#> ------checking: ASV6314
#> 
#> ------relative cooccurence: 04
#> 
#> No parent found!
#> 4
#> 
#> ####processing: ASV198 #####4
#> 
#> ---hits: ASV16424
#> 
#> ---potential parent: ASV16424
#> 
#> ------checking: ASV16424
#> 
#> ------relative cooccurence: 04
#> 
#> No parent found!
#> 4
#>   |                                                          |==============================================    |  92%
#> 
#> ####processing: ASV314 #####4
#> 
#> ---hits: 4
#> 
#> ---potential parent: 4
#> 
#> No parent found!
#> 4
#> 
#> ####processing: ASV405 #####4
#> 
#> ---hits: ASV69
#> ---hits: ASV358
#> ---hits: ASV914
#> ---hits: ASV930
#> ---hits: ASV221
#> ---hits: ASV833
#> ---hits: ASV1420
#> ---hits: ASV673
#> ---hits: ASV1218
#> ---hits: ASV4314
#> 
#> ---potential parent: ASV69
#> ---potential parent: ASV358
#> ---potential parent: ASV1218
#> ---potential parent: ASV221
#> ---potential parent: ASV833
#> ---potential parent: ASV1420
#> ---potential parent: ASV431
#> ---potential parent: ASV673
#> ---potential parent: ASV930
#> ---potential parent: ASV9144
#> 
#> ------checking: ASV694
#> 
#> ------relative cooccurence: 14
#>  which is sufficient!4
#> 
#> ------min avg abundance: 34
#>  which is OK!4
#> 
#> SETTING ASV405 to be an ERROR of ASV69
#> 4
#> 
#> ------checking: ASV3584
#> 
#> ------checking: ASV12184
#> 
#> ------checking: ASV2214
#> 
#> ------checking: ASV8334
#> 
#> ------checking: ASV14204
#> 
#> ------checking: ASV4314
#> 
#> ------checking: ASV6734
#> 
#> ------checking: ASV9304
#> 
#> ------checking: ASV9144
#> 
#> ####processing: ASV415 #####4
#> 
#> ---hits: ASV113
#> ---hits: ASV60
#> ---hits: ASV1380
#> ---hits: ASV1020
#> ---hits: ASV984
#> ---hits: ASV796
#> ---hits: ASV1036
#> ---hits: ASV1066
#> ---hits: ASV1011
#> ---hits: ASV15324
#> 
#> ---potential parent: ASV113
#> ---potential parent: ASV60
#> ---potential parent: ASV1532
#> ---potential parent: ASV1036
#> ---potential parent: ASV984
#> ---potential parent: ASV1011
#> ---potential parent: ASV796
#> ---potential parent: ASV1020
#> ---potential parent: ASV1066
#> ---potential parent: ASV13804
#> 
#> ------checking: ASV1134
#> 
#> ------relative cooccurence: 14
#>  which is sufficient!4
#> 
#> ------min avg abundance: 834
#>  which is OK!4
#> 
#> SETTING ASV415 to be an ERROR of ASV113
#> 4
#> 
#> ------checking: ASV604
#> 
#> ------checking: ASV15324
#> 
#> ------checking: ASV10364
#> 
#> ------checking: ASV9844
#> 
#> ------checking: ASV10114
#> 
#> ------checking: ASV7964
#> 
#> ------checking: ASV10204
#> 
#> ------checking: ASV10664
#> 
#> ------checking: ASV13804
#> 
#> ####processing: ASV731 #####4
#> 
#> ---hits: ASV266
#> ---hits: ASV24
#> ---hits: ASV736
#> ---hits: ASV979
#> ---hits: ASV986
#> ---hits: ASV1078
#> ---hits: ASV1467
#> ---hits: ASV858
#> ---hits: ASV1010
#> ---hits: ASV1784
#> 
#> ---potential parent: ASV178
#> ---potential parent: ASV24
#> ---potential parent: ASV266
#> ---potential parent: ASV1078
#> ---potential parent: ASV858
#> ---potential parent: ASV1010
#> ---potential parent: ASV1467
#> ---potential parent: ASV736
#> ---potential parent: ASV979
#> ---potential parent: ASV9864
#> 
#> ------checking: ASV1784
#> 
#> ------relative cooccurence: 04
#> 
#> ------checking: ASV244
#> 
#> ------relative cooccurence: 14
#>  which is sufficient!4
#> 
#> ------min avg abundance: 14
#> 
#> ------checking: ASV2664
#> 
#> ------relative cooccurence: 14
#>  which is sufficient!4
#> 
#> ------min avg abundance: 14
#> 
#> ------checking: ASV10784
#> 
#> ------relative cooccurence: 04
#> 
#> ------checking: ASV8584
#> 
#> ------relative cooccurence: 04
#> 
#> ------checking: ASV10104
#> 
#> ------relative cooccurence: 04
#> 
#> ------checking: ASV14674
#> 
#> ------relative cooccurence: 04
#> 
#> ------checking: ASV7364
#> 
#> ------relative cooccurence: 04
#> 
#> ------checking: ASV9794
#> 
#> ------relative cooccurence: 04
#> 
#> ------checking: ASV9864
#> 
#> ------relative cooccurence: 04
#> 
#> No parent found!
#> 4
#>   |                                                          |==============================================    |  93%
#> 
#> ####processing: ASV796 #####4
#> 
#> ---hits: ASV984
#> ---hits: ASV60
#> ---hits: ASV113
#> ---hits: ASV1020
#> ---hits: ASV1036
#> ---hits: ASV415
#> ---hits: ASV1380
#> ---hits: ASV1011
#> ---hits: ASV1603
#> ---hits: ASV10664
#> 
#> ---potential parent: ASV113
#> ---potential parent: ASV60
#> ---potential parent: ASV1036
#> ---potential parent: ASV984
#> ---potential parent: ASV1011
#> ---potential parent: ASV1603
#> ---potential parent: ASV415
#> ---potential parent: ASV1020
#> ---potential parent: ASV1066
#> ---potential parent: ASV13804
#> 
#> ------checking: ASV1134
#> 
#> ------relative cooccurence: 04
#> 
#> ------checking: ASV604
#> 
#> ------relative cooccurence: 04
#> 
#> ------checking: ASV10364
#> 
#> ------relative cooccurence: 04
#> 
#> ------checking: ASV9844
#> 
#> ------relative cooccurence: 04
#> 
#> ------checking: ASV10114
#> 
#> ------relative cooccurence: 14
#>  which is sufficient!4
#> 
#> ------min avg abundance: 54
#>  which is OK!4
#> 
#> SETTING ASV796 to be an ERROR of ASV1532
#> 4
#> 
#> ------checking: ASV16034
#> 
#> ------checking: ASV4154
#> 
#> ------checking: ASV10204
#> 
#> ------checking: ASV10664
#> 
#> ------checking: ASV13804
#> 
#> ####processing: ASV840 #####4
#> 
#> ---hits: ASV404
#> ---hits: ASV248
#> ---hits: ASV377
#> ---hits: ASV488
#> ---hits: ASV545
#> ---hits: ASV4404
#> 
#> ---potential parent: ASV248
#> ---potential parent: ASV377
#> ---potential parent: ASV404
#> ---potential parent: ASV545
#> ---potential parent: ASV440
#> ---potential parent: ASV4884
#> 
#> ------checking: ASV2484
#> 
#> ------relative cooccurence: 14
#>  which is sufficient!4
#> 
#> ------min avg abundance: 14
#> 
#> ------checking: ASV3774
#> 
#> ------relative cooccurence: 04
#> 
#> ------checking: ASV4044
#> 
#> ------relative cooccurence: 14
#>  which is sufficient!4
#> 
#> ------min avg abundance: 394
#>  which is OK!4
#> 
#> SETTING ASV840 to be an ERROR of ASV404
#> 4
#> 
#> ------checking: ASV5454
#> 
#> ------checking: ASV4404
#> 
#> ------checking: ASV4884
#>   |                                                          |===============================================   |  93%
#> 
#> ####processing: ASV860 #####4
#> 
#> ---hits: ASV504
#> ---hits: ASV832
#> ---hits: ASV772
#> ---hits: ASV1409
#> ---hits: ASV693
#> ---hits: ASV1200
#> ---hits: ASV1332
#> ---hits: ASV1458
#> ---hits: ASV1303
#> ---hits: ASV6504
#> 
#> ---potential parent: ASV1409
#> ---potential parent: ASV1458
#> ---potential parent: ASV772
#> ---potential parent: ASV693
#> ---potential parent: ASV1332
#> ---potential parent: ASV504
#> ---potential parent: ASV832
#> ---potential parent: ASV1200
#> ---potential parent: ASV650
#> ---potential parent: ASV13034
#> 
#> ------checking: ASV14094
#> 
#> ------relative cooccurence: 04
#> 
#> ------checking: ASV14584
#> 
#> ------relative cooccurence: 04
#> 
#> ------checking: ASV7724
#> 
#> ------relative cooccurence: 04
#> 
#> ------checking: ASV6934
#> 
#> ------relative cooccurence: 04
#> 
#> ------checking: ASV13324
#> 
#> ------relative cooccurence: 04
#> 
#> ------checking: ASV5044
#> 
#> ------relative cooccurence: 14
#>  which is sufficient!4
#> 
#> ------min avg abundance: 4044
#>  which is OK!4
#> 
#> SETTING ASV860 to be an ERROR of ASV504
#> 4
#> 
#> ------checking: ASV8324
#> 
#> ------checking: ASV12004
#> 
#> ------checking: ASV6504
#> 
#> ------checking: ASV13034
#> 
#> ####processing: ASV875 #####4
#> 
#> ---hits: ASV7874
#> 
#> ---potential parent: ASV7874
#> 
#> ------checking: ASV7874
#> 
#> ------relative cooccurence: 14
#>  which is sufficient!4
#> 
#> ------min avg abundance: 74
#>  which is OK!4
#> 
#> SETTING ASV875 to be an ERROR of ASV787
#> 4
#>   |                                                          |===============================================   |  94%
#> 
#> ####processing: ASV892 #####4
#> 
#> ---hits: ASV1031
#> ---hits: ASV1239
#> ---hits: ASV1510
#> ---hits: ASV13004
#> 
#> ---potential parent: ASV1300
#> ---potential parent: ASV1510
#> ---potential parent: ASV1031
#> ---potential parent: ASV12394
#> 
#> ------checking: ASV13004
#> 
#> ------relative cooccurence: 04
#> 
#> ------checking: ASV15104
#> 
#> ------relative cooccurence: 04
#> 
#> ------checking: ASV10314
#> 
#> ------relative cooccurence: 04
#> 
#> ------checking: ASV12394
#> 
#> ------relative cooccurence: 14
#>  which is sufficient!4
#> 
#> ------min avg abundance: 14
#> 
#> No parent found!
#> 4
#> 
#> ####processing: ASV940 #####4
#> 
#> ---hits: ASV996
#> ---hits: ASV45
#> ---hits: ASV137
#> ---hits: ASV1067
#> ---hits: ASV711
#> ---hits: ASV760
#> ---hits: ASV1167
#> ---hits: ASV1484
#> ---hits: ASV14684
#> 
#> ---potential parent: ASV45
#> ---potential parent: ASV711
#> ---potential parent: ASV1067
#> ---potential parent: ASV137
#> ---potential parent: ASV1468
#> ---potential parent: ASV996
#> ---potential parent: ASV1167
#> ---potential parent: ASV760
#> ---potential parent: ASV14844
#> 
#> ------checking: ASV454
#> 
#> ------relative cooccurence: 14
#>  which is sufficient!4
#> 
#> ------min avg abundance: 8364
#>  which is OK!4
#> 
#> SETTING ASV940 to be an ERROR of ASV45
#> 4
#> 
#> ------checking: ASV7114
#> 
#> ------checking: ASV10674
#> 
#> ------checking: ASV1374
#> 
#> ------checking: ASV14684
#> 
#> ------checking: ASV9964
#> 
#> ------checking: ASV11674
#> 
#> ------checking: ASV7604
#> 
#> ------checking: ASV14844
#> 
#> ####processing: ASV941 #####4
#> 
#> ---hits: ASV33
#> ---hits: ASV78
#> ---hits: ASV1625
#> ---hits: ASV13134
#> 
#> ---potential parent: ASV78
#> ---potential parent: ASV1313
#> ---potential parent: ASV33
#> ---potential parent: ASV16254
#> 
#> ------checking: ASV784
#> 
#> ------relative cooccurence: 14
#>  which is sufficient!4
#> 
#> ------min avg abundance: 60054
#>  which is OK!4
#> 
#> SETTING ASV941 to be an ERROR of ASV78
#> 4
#> 
#> ------checking: ASV13134
#> 
#> ------checking: ASV334
#> 
#> ------checking: ASV16254
#> 
#> ####processing: ASV959 #####4
#> 
#> ---hits: ASV81
#> ---hits: ASV272
#> ---hits: ASV339
#> ---hits: ASV383
#> ---hits: ASV839
#> ---hits: ASV1152
#> ---hits: ASV1144
#> ---hits: ASV12634
#> 
#> ---potential parent: ASV1152
#> ---potential parent: ASV81
#> ---potential parent: ASV339
#> ---potential parent: ASV383
#> ---potential parent: ASV1144
#> ---potential parent: ASV839
#> ---potential parent: ASV272
#> ---potential parent: ASV12634
#> 
#> ------checking: ASV11524
#> 
#> ------relative cooccurence: 04
#> 
#> ------checking: ASV814
#> 
#> ------relative cooccurence: 14
#>  which is sufficient!4
#> 
#> ------min avg abundance: 44
#>  which is OK!4
#> 
#> SETTING ASV959 to be an ERROR of ASV81
#> 4
#> 
#> ------checking: ASV3394
#> 
#> ------checking: ASV3834
#> 
#> ------checking: ASV11444
#> 
#> ------checking: ASV8394
#> 
#> ------checking: ASV2724
#> 
#> ------checking: ASV12634
#>   |                                                          |===============================================   |  95%
#> 
#> ####processing: ASV981 #####4
#> 
#> ---hits: ASV715
#> ---hits: ASV963
#> ---hits: ASV737
#> ---hits: ASV1369
#> ---hits: ASV1724
#> 
#> ---potential parent: ASV172
#> ---potential parent: ASV1369
#> ---potential parent: ASV737
#> ---potential parent: ASV963
#> ---potential parent: ASV7154
#> 
#> ------checking: ASV1724
#> 
#> ------relative cooccurence: 14
#>  which is sufficient!4
#> 
#> ------min avg abundance: 2554
#>  which is OK!4
#> 
#> SETTING ASV981 to be an ERROR of ASV172
#> 4
#> 
#> ------checking: ASV13694
#> 
#> ------checking: ASV7374
#> 
#> ------checking: ASV9634
#> 
#> ------checking: ASV7154
#> 
#> ####processing: ASV988 #####4
#> 
#> ---hits: 4
#> 
#> ---potential parent: 4
#> 
#> No parent found!
#> 4
#>   |                                                          |================================================  |  95%
#> 
#> ####processing: ASV1015 #####4
#> 
#> ---hits: ASV637
#> ---hits: ASV9034
#> 
#> ---potential parent: ASV637
#> ---potential parent: ASV9034
#> 
#> ------checking: ASV6374
#> 
#> ------relative cooccurence: 14
#>  which is sufficient!4
#> 
#> ------min avg abundance: 1664
#>  which is OK!4
#> 
#> SETTING ASV1015 to be an ERROR of ASV637
#> 4
#> 
#> ------checking: ASV9034
#> 
#> ####processing: ASV1020 #####4
#> 
#> ---hits: ASV1036
#> ---hits: ASV60
#> ---hits: ASV984
#> ---hits: ASV796
#> ---hits: ASV1603
#> ---hits: ASV113
#> ---hits: ASV415
#> ---hits: ASV1380
#> ---hits: ASV1011
#> ---hits: ASV10664
#> 
#> ---potential parent: ASV113
#> ---potential parent: ASV60
#> ---potential parent: ASV1036
#> ---potential parent: ASV984
#> ---potential parent: ASV1011
#> ---potential parent: ASV1603
#> ---potential parent: ASV415
#> ---potential parent: ASV796
#> ---potential parent: ASV1066
#> ---potential parent: ASV13804
#> 
#> ------checking: ASV1134
#> 
#> ------relative cooccurence: 14
#>  which is sufficient!4
#> 
#> ------min avg abundance: 44
#>  which is OK!4
#> 
#> SETTING ASV1020 to be an ERROR of ASV113
#> 4
#> 
#> ------checking: ASV604
#> 
#> ------checking: ASV10364
#> 
#> ------checking: ASV9844
#> 
#> ------checking: ASV10114
#> 
#> ------checking: ASV16034
#> 
#> ------checking: ASV4154
#> 
#> ------checking: ASV7964
#> 
#> ------checking: ASV10664
#> 
#> ------checking: ASV13804
#>   |                                                          |================================================  |  96%
#> 
#> ####processing: ASV1027 #####4
#> 
#> ---hits: ASV1082
#> ---hits: ASV962
#> ---hits: ASV2
#> ---hits: ASV1247
#> ---hits: ASV5764
#> 
#> ---potential parent: ASV2
#> ---potential parent: ASV576
#> ---potential parent: ASV1247
#> ---potential parent: ASV962
#> ---potential parent: ASV10824
#> 
#> ------checking: ASV24
#> 
#> ------relative cooccurence: 14
#>  which is sufficient!4
#> 
#> ------min avg abundance: 2814
#>  which is OK!4
#> 
#> SETTING ASV1027 to be an ERROR of ASV2
#> 4
#> 
#> ------checking: ASV5764
#> 
#> ------checking: ASV12474
#> 
#> ------checking: ASV9624
#> 
#> ------checking: ASV10824
#> 
#> ####processing: ASV1031 #####4
#> 
#> ---hits: ASV892
#> ---hits: ASV1239
#> ---hits: ASV1510
#> ---hits: ASV13004
#> 
#> ---potential parent: ASV1300
#> ---potential parent: ASV1510
#> ---potential parent: ASV892
#> ---potential parent: ASV12394
#> 
#> ------checking: ASV13004
#> 
#> ------relative cooccurence: 14
#>  which is sufficient!4
#> 
#> ------min avg abundance: 24
#>  which is OK!4
#> 
#> SETTING ASV1031 to be an ERROR of ASV1300
#> 4
#> 
#> ------checking: ASV15104
#> 
#> ------checking: ASV8924
#> 
#> ------checking: ASV12394
#> 
#> ####processing: ASV1066 #####4
#> 
#> ---hits: ASV1011
#> ---hits: ASV113
#> ---hits: ASV415
#> ---hits: ASV60
#> ---hits: ASV1380
#> ---hits: ASV984
#> ---hits: ASV796
#> ---hits: ASV1020
#> ---hits: ASV1036
#> ---hits: ASV15324
#> 
#> ---potential parent: ASV113
#> ---potential parent: ASV60
#> ---potential parent: ASV1532
#> ---potential parent: ASV1036
#> ---potential parent: ASV984
#> ---potential parent: ASV1011
#> ---potential parent: ASV415
#> ---potential parent: ASV796
#> ---potential parent: ASV1020
#> ---potential parent: ASV13804
#> 
#> ------checking: ASV1134
#> 
#> ------relative cooccurence: 14
#>  which is sufficient!4
#> 
#> ------min avg abundance: 44
#>  which is OK!4
#> 
#> SETTING ASV1066 to be an ERROR of ASV113
#> 4
#> 
#> ------checking: ASV604
#> 
#> ------checking: ASV15324
#> 
#> ------checking: ASV10364
#> 
#> ------checking: ASV9844
#> 
#> ------checking: ASV10114
#> 
#> ------checking: ASV4154
#> 
#> ------checking: ASV7964
#> 
#> ------checking: ASV10204
#> 
#> ------checking: ASV13804
#> 
#> ####processing: ASV1107 #####4
#> 
#> ---hits: ASV116
#> ---hits: ASV1526
#> ---hits: ASV500
#> ---hits: ASV989
#> ---hits: ASV816
#> ---hits: ASV254
#> ---hits: ASV38
#> ---hits: ASV4664
#> 
#> ---potential parent: ASV38
#> ---potential parent: ASV989
#> ---potential parent: ASV816
#> ---potential parent: ASV500
#> ---potential parent: ASV466
#> ---potential parent: ASV254
#> ---potential parent: ASV1526
#> ---potential parent: ASV1164
#> 
#> ------checking: ASV384
#> 
#> ------relative cooccurence: 14
#>  which is sufficient!4
#> 
#> ------min avg abundance: 4374
#>  which is OK!4
#> 
#> SETTING ASV1107 to be an ERROR of ASV38
#> 4
#> 
#> ------checking: ASV9894
#> 
#> ------checking: ASV8164
#> 
#> ------checking: ASV5004
#> 
#> ------checking: ASV4664
#> 
#> ------checking: ASV2544
#> 
#> ------checking: ASV15264
#> 
#> ------checking: ASV1164
#>   |                                                          |================================================  |  97%
#> 
#> ####processing: ASV1111 #####4
#> 
#> ---hits: ASV443
#> ---hits: ASV1176
#> ---hits: ASV368
#> ---hits: ASV292
#> ---hits: ASV961
#> ---hits: ASV566
#> ---hits: ASV984
#> 
#> ---potential parent: ASV292
#> ---potential parent: ASV443
#> ---potential parent: ASV98
#> ---potential parent: ASV368
#> ---potential parent: ASV566
#> ---potential parent: ASV1176
#> ---potential parent: ASV9614
#> 
#> ------checking: ASV2924
#> 
#> ------relative cooccurence: 14
#>  which is sufficient!4
#> 
#> ------min avg abundance: 14
#> 
#> ------checking: ASV4434
#> 
#> ------relative cooccurence: 14
#>  which is sufficient!4
#> 
#> ------min avg abundance: 24
#>  which is OK!4
#> 
#> SETTING ASV1111 to be an ERROR of ASV443
#> 4
#> 
#> ------checking: ASV984
#> 
#> ------checking: ASV3684
#> 
#> ------checking: ASV5664
#> 
#> ------checking: ASV11764
#> 
#> ------checking: ASV9614
#> 
#> ####processing: ASV1239 #####4
#> 
#> ---hits: ASV1031
#> ---hits: ASV892
#> ---hits: ASV1510
#> ---hits: ASV13004
#> 
#> ---potential parent: ASV1300
#> ---potential parent: ASV1510
#> ---potential parent: ASV892
#> ---potential parent: ASV10314
#> 
#> ------checking: ASV13004
#> 
#> ------relative cooccurence: 04
#> 
#> ------checking: ASV15104
#> 
#> ------relative cooccurence: 04
#> 
#> ------checking: ASV8924
#> 
#> ------relative cooccurence: 14
#>  which is sufficient!4
#> 
#> ------min avg abundance: 14
#> 
#> ------checking: ASV10314
#> 
#> ------relative cooccurence: 04
#> 
#> No parent found!
#> 4
#>   |                                                          |================================================= |  97%
#> 
#> ####processing: ASV1263 #####4
#> 
#> ---hits: ASV1144
#> ---hits: ASV272
#> ---hits: ASV339
#> ---hits: ASV839
#> ---hits: ASV81
#> ---hits: ASV383
#> ---hits: ASV959
#> ---hits: ASV11524
#> 
#> ---potential parent: ASV1152
#> ---potential parent: ASV81
#> ---potential parent: ASV339
#> ---potential parent: ASV383
#> ---potential parent: ASV1144
#> ---potential parent: ASV839
#> ---potential parent: ASV272
#> ---potential parent: ASV9594
#> 
#> ------checking: ASV11524
#> 
#> ------relative cooccurence: 04
#> 
#> ------checking: ASV814
#> 
#> ------relative cooccurence: 04
#> 
#> ------checking: ASV3394
#> 
#> ------relative cooccurence: 04
#> 
#> ------checking: ASV3834
#> 
#> ------relative cooccurence: 04
#> 
#> ------checking: ASV11444
#> 
#> ------relative cooccurence: 14
#>  which is sufficient!4
#> 
#> ------min avg abundance: 44
#>  which is OK!4
#> 
#> SETTING ASV1263 to be an ERROR of ASV1144
#> 4
#> 
#> ------checking: ASV8394
#> 
#> ------checking: ASV2724
#> 
#> ------checking: ASV9594
#> 
#> ####processing: ASV1302 #####4
#> 
#> ---hits: 4
#> 
#> ---potential parent: 4
#> 
#> No parent found!
#> 4
#>   |                                                          |================================================= |  98%
#> 
#> ####processing: ASV1303 #####4
#> 
#> ---hits: ASV378
#> ---hits: ASV1332
#> ---hits: ASV746
#> ---hits: ASV1562
#> ---hits: ASV1030
#> ---hits: ASV1230
#> ---hits: ASV1458
#> ---hits: ASV1200
#> ---hits: ASV693
#> ---hits: ASV14214
#> 
#> ---potential parent: ASV378
#> ---potential parent: ASV746
#> ---potential parent: ASV1421
#> ---potential parent: ASV1230
#> ---potential parent: ASV1458
#> ---potential parent: ASV693
#> ---potential parent: ASV1332
#> ---potential parent: ASV1200
#> ---potential parent: ASV1030
#> ---potential parent: ASV15624
#> 
#> ------checking: ASV3784
#> 
#> ------relative cooccurence: 04
#> 
#> ------checking: ASV7464
#> 
#> ------relative cooccurence: 14
#>  which is sufficient!4
#> 
#> ------min avg abundance: 1564
#>  which is OK!4
#> 
#> SETTING ASV1303 to be an ERROR of ASV746
#> 4
#> 
#> ------checking: ASV14214
#> 
#> ------checking: ASV12304
#> 
#> ------checking: ASV14584
#> 
#> ------checking: ASV6934
#> 
#> ------checking: ASV13324
#> 
#> ------checking: ASV12004
#> 
#> ------checking: ASV10304
#> 
#> ------checking: ASV15624
#> 
#> ####processing: ASV1341 #####4
#> 
#> ---hits: ASV14604
#> 
#> ---potential parent: ASV14604
#> 
#> ------checking: ASV14604
#> 
#> ------relative cooccurence: 14
#>  which is sufficient!4
#> 
#> ------min avg abundance: 144
#>  which is OK!4
#> 
#> SETTING ASV1341 to be an ERROR of ASV1460
#> 4
#> 
#> ####processing: ASV1380 #####4
#> 
#> ---hits: ASV113
#> ---hits: ASV415
#> ---hits: ASV60
#> ---hits: ASV984
#> ---hits: ASV796
#> ---hits: ASV1066
#> ---hits: ASV1020
#> ---hits: ASV1036
#> ---hits: ASV1011
#> ---hits: ASV16034
#> 
#> ---potential parent: ASV113
#> ---potential parent: ASV60
#> ---potential parent: ASV1036
#> ---potential parent: ASV984
#> ---potential parent: ASV1011
#> ---potential parent: ASV1603
#> ---potential parent: ASV415
#> ---potential parent: ASV796
#> ---potential parent: ASV1020
#> ---potential parent: ASV10664
#> 
#> ------checking: ASV1134
#> 
#> ------relative cooccurence: 14
#>  which is sufficient!4
#> 
#> ------min avg abundance: 114
#>  which is OK!4
#> 
#> SETTING ASV1380 to be an ERROR of ASV113
#> 4
#> 
#> ------checking: ASV604
#> 
#> ------checking: ASV10364
#> 
#> ------checking: ASV9844
#> 
#> ------checking: ASV10114
#> 
#> ------checking: ASV16034
#> 
#> ------checking: ASV4154
#> 
#> ------checking: ASV7964
#> 
#> ------checking: ASV10204
#> 
#> ------checking: ASV10664
#> 
#> ####processing: ASV1388 #####4
#> 
#> ---hits: ASV1653
#> ---hits: ASV973
#> ---hits: ASV1468
#> ---hits: ASV552
#> ---hits: ASV592
#> ---hits: ASV1261
#> ---hits: ASV346
#> ---hits: ASV45
#> ---hits: ASV137
#> ---hits: ASV9964
#> 
#> ---potential parent: ASV346
#> ---potential parent: ASV592
#> ---potential parent: ASV1261
#> ---potential parent: ASV973
#> ---potential parent: ASV45
#> ---potential parent: ASV552
#> ---potential parent: ASV137
#> ---potential parent: ASV1468
#> ---potential parent: ASV996
#> ---potential parent: ASV16534
#> 
#> ------checking: ASV3464
#> 
#> ------relative cooccurence: 14
#>  which is sufficient!4
#> 
#> ------min avg abundance: 1554
#>  which is OK!4
#> 
#> SETTING ASV1388 to be an ERROR of ASV618
#> 4
#> 
#> ------checking: ASV5924
#> 
#> ------checking: ASV12614
#> 
#> ------checking: ASV9734
#> 
#> ------checking: ASV454
#> 
#> ------checking: ASV5524
#> 
#> ------checking: ASV1374
#> 
#> ------checking: ASV14684
#> 
#> ------checking: ASV9964
#> 
#> ------checking: ASV16534
#>   |                                                          |================================================= |  99%
#> 
#> ####processing: ASV1484 #####4
#> 
#> ---hits: ASV711
#> ---hits: ASV45
#> ---hits: ASV137
#> ---hits: ASV996
#> ---hits: ASV760
#> ---hits: ASV1067
#> ---hits: ASV1167
#> ---hits: ASV940
#> ---hits: ASV14684
#> 
#> ---potential parent: ASV45
#> ---potential parent: ASV711
#> ---potential parent: ASV1067
#> ---potential parent: ASV137
#> ---potential parent: ASV1468
#> ---potential parent: ASV996
#> ---potential parent: ASV1167
#> ---potential parent: ASV760
#> ---potential parent: ASV9404
#> 
#> ------checking: ASV454
#> 
#> ------relative cooccurence: 14
#>  which is sufficient!4
#> 
#> ------min avg abundance: 8364
#>  which is OK!4
#> 
#> SETTING ASV1484 to be an ERROR of ASV45
#> 4
#> 
#> ------checking: ASV7114
#> 
#> ------checking: ASV10674
#> 
#> ------checking: ASV1374
#> 
#> ------checking: ASV14684
#> 
#> ------checking: ASV9964
#> 
#> ------checking: ASV11674
#> 
#> ------checking: ASV7604
#> 
#> ------checking: ASV9404
#> 
#> ####processing: ASV1550 #####4
#> 
#> ---hits: ASV3594
#> 
#> ---potential parent: ASV3594
#> 
#> ------checking: ASV3594
#> 
#> ------relative cooccurence: 04
#> 
#> No parent found!
#> 4
#>   |                                                          |==================================================|  99%
#> 
#> ####processing: ASV1562 #####4
#> 
#> ---hits: ASV1303
#> ---hits: ASV1332
#> ---hits: ASV1200
#> ---hits: ASV1458
#> ---hits: ASV378
#> ---hits: ASV1030
#> ---hits: ASV746
#> ---hits: ASV1230
#> ---hits: ASV1409
#> ---hits: ASV6934
#> 
#> ---potential parent: ASV378
#> ---potential parent: ASV746
#> ---potential parent: ASV1409
#> ---potential parent: ASV1230
#> ---potential parent: ASV1458
#> ---potential parent: ASV693
#> ---potential parent: ASV1332
#> ---potential parent: ASV1200
#> ---potential parent: ASV1030
#> ---potential parent: ASV13034
#> 
#> ------checking: ASV3784
#> 
#> ------relative cooccurence: 04
#> 
#> ------checking: ASV7464
#> 
#> ------relative cooccurence: 04
#> 
#> ------checking: ASV14094
#> 
#> ------relative cooccurence: 04
#> 
#> ------checking: ASV12304
#> 
#> ------relative cooccurence: 04
#> 
#> ------checking: ASV14584
#> 
#> ------relative cooccurence: 04
#> 
#> ------checking: ASV6934
#> 
#> ------relative cooccurence: 04
#> 
#> ------checking: ASV13324
#> 
#> ------relative cooccurence: 04
#> 
#> ------checking: ASV12004
#> 
#> ------relative cooccurence: 04
#> 
#> ------checking: ASV10304
#> 
#> ------relative cooccurence: 04
#> 
#> ------checking: ASV13034
#> 
#> ------relative cooccurence: 04
#> 
#> No parent found!
#> 4
#> 
#> ####processing: ASV1635 #####4
#> 
#> ---hits: 4
#> 
#> ---potential parent: 4
#> 
#> No parent found!
#> 4
#>   |                                                          |==================================================| 100%
#> 
#> ####processing: ASV1653 #####4
#> 
#> ---hits: ASV1388
#> ---hits: ASV973
#> ---hits: ASV1468
#> ---hits: ASV1261
#> ---hits: ASV552
#> ---hits: ASV592
#> ---hits: ASV346
#> ---hits: ASV1167
#> ---hits: ASV45
#> ---hits: ASV9964
#> 
#> ---potential parent: ASV346
#> ---potential parent: ASV592
#> ---potential parent: ASV1261
#> ---potential parent: ASV973
#> ---potential parent: ASV45
#> ---potential parent: ASV552
#> ---potential parent: ASV1468
#> ---potential parent: ASV996
#> ---potential parent: ASV1167
#> ---potential parent: ASV13884
#> 
#> ------checking: ASV3464
#> 
#> ------relative cooccurence: 04
#> 
#> ------checking: ASV5924
#> 
#> ------relative cooccurence: 04
#> 
#> ------checking: ASV12614
#> 
#> ------relative cooccurence: 04
#> 
#> ------checking: ASV9734
#> 
#> ------relative cooccurence: 14
#>  which is sufficient!4
#> 
#> ------min avg abundance: 14
#> 
#> ------checking: ASV454
#> 
#> ------relative cooccurence: 04
#> 
#> ------checking: ASV5524
#> 
#> ------relative cooccurence: 04
#> 
#> ------checking: ASV14684
#> 
#> ------relative cooccurence: 04
#> 
#> ------checking: ASV9964
#> 
#> ------relative cooccurence: 04
#> 
#> ------checking: ASV11674
#> 
#> ------relative cooccurence: 04
#> 
#> ------checking: ASV13884
#> 
#> ------relative cooccurence: 04
#> 
#> No parent found!
#> 4
#> 
#> ####processing: ASV1674 #####4
#> 
#> ---hits: 4
#> 
#> ---potential parent: 4
#> 
#> No parent found!
#> 4
#> 
#> ####processing: ASV1712 #####4
#> 
#> ---hits: 4
#> 
#> ---potential parent: 4
#> 
#> No parent found!
#> 4
#> $new_physeq
#> phyloseq-class experiment-level object
#> otu_table()   OTU Table:         [ 264 taxa and 20 samples ]
#> sample_data() Sample Data:       [ 20 samples by 7 sample variables ]
#> tax_table()   Taxonomy Table:    [ 264 taxa by 12 taxonomic ranks ]
#> refseq()      DNAStringSet:      [ 264 reference sequences ]
#> 
#> $discrepancy_vector
#>             Domain             Phylum              Class              Order 
#>          1.0000000          0.8299320          0.5850340          0.5578231 
#>             Family              Genus            Species       Trophic.Mode 
#>          0.5034014          0.4829932          0.4149660          0.6598639 
#>              Guild              Trait Confidence.Ranking      Genus_species 
#>          0.6054422          0.7823129          0.7142857          0.4149660 
#> 
#> $res_lulu
#> $res_lulu$curated_table
#>         A10.005.B_S188_MERGED.fastq.gz A10.005.H_S189_MERGED.fastq.gz
#> ASV1004                              1                              0
#> ASV1007                              0                              0
#> ASV101                               0                              0
#> ASV1017                              0                              0
#> ASV1022                              0                              0
#> ASV1030                              0                              0
#> ASV104                               0                              0
#> ASV105                               0                              0
#> ASV1059                              0                              0
#> ASV1069                              0                              0
#> ASV107                               0                              0
#> ASV1070                              0                              0
#> ASV1074                              0                              0
#> ASV1085                              0                              0
#> ASV1102                              0                              0
#> ASV1108                              0                              0
#> ASV111                               0                              7
#> ASV1118                              0                              0
#> ASV113                               0                              5
#> ASV1132                              0                              0
#> ASV1144                              0                              0
#> ASV1146                              0                              0
#> ASV1152                              0                             83
#> ASV1155                              0                              0
#> ASV1164                              0                              0
#> ASV1166                              0                              0
#> ASV1167                              0                              0
#> ASV1179                              0                              0
#> ASV1192                              5                              0
#> ASV12                                0                              1
#> ASV1200                              0                              0
#> ASV1212                              0                              0
#> ASV1218                              0                              0
#> ASV1223                              0                              0
#> ASV1230                              0                              3
#> ASV1233                              0                              0
#> ASV1234                              0                              0
#> ASV1236                             11                              0
#> ASV1239                              0                              0
#> ASV1242                              6                              0
#> ASV1253                              0                              0
#> ASV1270                              0                              0
#> ASV1278                              1                              0
#> ASV1279                              0                              0
#> ASV1300                              0                              0
#> ASV1302                              0                              0
#> ASV131                               0                              0
#> ASV1311                              0                              0
#> ASV1314                              0                              0
#> ASV1321                              0                              0
#> ASV1332                              0                              0
#> ASV1340                              0                              0
#> ASV1359                              0                              0
#> ASV1365                              0                              1
#> ASV1369                              0                              0
#> ASV137                               0                              0
#> ASV1370                              0                              0
#> ASV1379                              0                              0
#> ASV139                               0                              0
#> ASV1396                              0                              0
#> ASV1409                              0                              0
#> ASV1419                              0                              0
#> ASV1420                              0                              0
#> ASV1421                              0                              0
#> ASV1422                              0                              0
#> ASV1428                              0                              0
#> ASV1434                              0                              0
#> ASV144                               0                              0
#> ASV1458                              0                              0
#> ASV1459                              0                              0
#> ASV1460                              0                              0
#> ASV1461                              0                              0
#> ASV1462                              0                              0
#> ASV1468                              0                              0
#> ASV1483                              0                              1
#> ASV1485                              0                              0
#> ASV151                               0                              0
#> ASV1510                              0                              0
#> ASV1523                              0                              0
#> ASV1532                              0                              0
#> ASV1537                              0                              0
#> ASV1548                              0                              0
#> ASV1550                              0                              0
#> ASV1552                              0                              0
#> ASV1561                              0                              0
#> ASV1562                              0                              0
#> ASV1574                              0                              0
#> ASV1576                              0                              0
#> ASV1577                              0                              0
#> ASV1588                              0                              0
#> ASV159                               0                              0
#> ASV1609                              0                              0
#> ASV162                               0                              0
#> ASV1624                              0                              0
#> ASV1628                              0                              0
#> ASV1630                              0                              0
#> ASV1632                              0                              0
#> ASV1635                              0                              0
#> ASV1642                              0                              0
#> ASV1644                             11                              0
#> ASV1650                              0                              0
#> ASV1653                              0                              0
#> ASV1674                              0                              0
#> ASV1683                              0                              0
#> ASV1690                              0                              0
#> ASV170                             104                              0
#> ASV171                               0                              0
#> ASV1712                              0                              0
#> ASV172                               0                              0
#> ASV1726                              0                              0
#> ASV175                               0                              0
#> ASV18                               30                              0
#> ASV19                                0                              0
#> ASV193                               0                              0
#> ASV198                               0                              0
#> ASV199                               0                              0
#> ASV2                                 1                             23
#> ASV201                               0                              0
#> ASV202                               0                              0
#> ASV210                               0                              0
#> ASV219                               0                              0
#> ASV221                               0                              0
#> ASV223                               0                              0
#> ASV24                                0                              0
#> ASV242                               0                              0
#> ASV244                               0                              0
#> ASV248                               0                              0
#> ASV25                                0                              0
#> ASV251                               3                              0
#> ASV253                               0                              0
#> ASV255                               0                              0
#> ASV26                              351                              0
#> ASV264                               0                              0
#> ASV266                               0                              0
#> ASV27                                0                              0
#> ASV28                                1                              0
#> ASV286                               0                              0
#> ASV29                                0                              0
#> ASV292                               0                              0
#> ASV297                               0                              0
#> ASV310                               0                              0
#> ASV314                               0                              0
#> ASV320                               0                              0
#> ASV322                               0                              0
#> ASV329                               0                              0
#> ASV334                               0                              0
#> ASV338                               0                              0
#> ASV339                               0                              0
#> ASV344                               0                              0
#> ASV357                               0                              0
#> ASV359                               0                              0
#> ASV365                               0                              0
#> ASV368                               0                              0
#> ASV377                               0                              0
#> ASV378                               0                              0
#> ASV38                                0                             10
#> ASV383                               0                              0
#> ASV404                               0                              0
#> ASV41                                6                              0
#> ASV42                                0                              0
#> ASV428                               0                              0
#> ASV43                                0                              0
#> ASV431                               0                              0
#> ASV443                               0                              0
#> ASV444                               0                              0
#> ASV45                                0                              0
#> ASV454                               0                              0
#> ASV46                                0                              0
#> ASV464                               0                              0
#> ASV47                                0                              0
#> ASV477                               0                              0
#> ASV48                                0                              0
#> ASV483                               0                              0
#> ASV489                               0                              0
#> ASV49                                0                              0
#> ASV500                               0                             41
#> ASV504                               0                              0
#> ASV505                               0                              0
#> ASV507                               0                              0
#> ASV509                               0                              0
#> ASV515                               0                              0
#> ASV517                               0                              0
#> ASV52                                0                              0
#> ASV523                               0                              0
#> ASV524                               0                              0
#> ASV533                               0                              0
#> ASV543                               0                              0
#> ASV545                               0                              0
#> ASV546                               0                              0
#> ASV549                               0                              0
#> ASV55                                0                              0
#> ASV552                               0                            489
#> ASV559                               0                              0
#> ASV563                              13                              0
#> ASV569                               0                              0
#> ASV570                               0                              2
#> ASV577                               0                              0
#> ASV586                               0                              0
#> ASV59                                0                              0
#> ASV598                               0                              0
#> ASV60                                0                              0
#> ASV602                               0                              0
#> ASV618                               0                             26
#> ASV619                               0                              0
#> ASV631                               0                              0
#> ASV635                               0                              0
#> ASV637                               0                              0
#> ASV64                                0                              0
#> ASV643                               0                              0
#> ASV65                                0                              0
#> ASV662                               0                              0
#> ASV672                             319                              0
#> ASV673                               0                              0
#> ASV685                               0                              0
#> ASV69                                0                              0
#> ASV692                               0                              0
#> ASV694                               0                              0
#> ASV715                               0                              0
#> ASV717                               0                              0
#> ASV727                               0                              0
#> ASV731                               0                              0
#> ASV737                               0                              0
#> ASV741                               0                              0
#> ASV743                               0                              0
#> ASV746                               0                              0
#> ASV749                               0                              0
#> ASV753                               0                              0
#> ASV756                              49                              0
#> ASV758                               0                              0
#> ASV772                               0                              0
#> ASV78                                0                              0
#> ASV784                               0                              0
#> ASV787                               0                              0
#> ASV8                               996                              9
#> ASV801                               0                              0
#> ASV81                                0                              0
#> ASV817                               0                              0
#> ASV823                               0                              0
#> ASV829                               2                              0
#> ASV834                               0                              0
#> ASV839                               0                              0
#> ASV854                               0                              0
#> ASV87                                0                              0
#> ASV89                                0                              0
#> ASV891                               0                              0
#> ASV892                               0                              0
#> ASV898                               0                              0
#> ASV899                               0                              0
#> ASV904                               0                              0
#> ASV906                               0                              0
#> ASV91                                0                              0
#> ASV914                               0                              0
#> ASV926                               0                              0
#> ASV927                               0                              0
#> ASV930                               0                              0
#> ASV942                               0                              0
#> ASV946                             197                              0
#> ASV950                               0                              0
#> ASV953                               0                              7
#> ASV964                               0                              0
#> ASV969                               0                            124
#> ASV973                               0                              0
#> ASV987                               0                              0
#> ASV988                               0                              1
#>         A10.005.M_S190_MERGED.fastq.gz A12.007_S191_MERGED.fastq.gz
#> ASV1004                              0                            0
#> ASV1007                              0                           93
#> ASV101                               0                            0
#> ASV1017                              0                            0
#> ASV1022                              0                            0
#> ASV1030                              0                            0
#> ASV104                               0                            0
#> ASV105                               0                            0
#> ASV1059                              0                            0
#> ASV1069                              0                            0
#> ASV107                               0                            0
#> ASV1070                              0                            0
#> ASV1074                              1                            3
#> ASV1085                              0                            3
#> ASV1102                              0                            0
#> ASV1108                              1                            1
#> ASV111                               0                            0
#> ASV1118                              0                            0
#> ASV113                               0                            0
#> ASV1132                              0                            0
#> ASV1144                              0                            0
#> ASV1146                              0                            0
#> ASV1152                              0                            1
#> ASV1155                              0                            7
#> ASV1164                              0                            0
#> ASV1166                              0                            0
#> ASV1167                              0                            0
#> ASV1179                              0                            0
#> ASV1192                              0                            0
#> ASV12                                0                         3293
#> ASV1200                              0                            0
#> ASV1212                              0                            0
#> ASV1218                              0                            0
#> ASV1223                              0                            0
#> ASV1230                              0                            0
#> ASV1233                              0                            0
#> ASV1234                              0                            0
#> ASV1236                              0                            0
#> ASV1239                              0                            0
#> ASV1242                              0                            0
#> ASV1253                              0                            0
#> ASV1270                              0                            0
#> ASV1278                              0                            0
#> ASV1279                              0                            0
#> ASV1300                              0                            0
#> ASV1302                              0                            0
#> ASV131                               0                            0
#> ASV1311                              0                            0
#> ASV1314                              0                            0
#> ASV1321                              0                            0
#> ASV1332                              0                            0
#> ASV1340                              0                            0
#> ASV1359                              0                           46
#> ASV1365                              0                            0
#> ASV1369                              0                            0
#> ASV137                               0                            0
#> ASV1370                              0                            0
#> ASV1379                              0                            0
#> ASV139                               0                            0
#> ASV1396                              0                            0
#> ASV1409                              0                            5
#> ASV1419                              0                            0
#> ASV1420                              0                            0
#> ASV1421                              0                            0
#> ASV1422                              0                            0
#> ASV1428                              0                            0
#> ASV1434                              0                           12
#> ASV144                               0                            0
#> ASV1458                             73                            3
#> ASV1459                              0                            0
#> ASV1460                              0                            0
#> ASV1461                              0                            0
#> ASV1462                              0                            0
#> ASV1468                              0                            0
#> ASV1483                              0                            1
#> ASV1485                              0                            0
#> ASV151                               0                            0
#> ASV1510                              0                            0
#> ASV1523                              0                            0
#> ASV1532                              0                            0
#> ASV1537                              0                           53
#> ASV1548                              0                            0
#> ASV1550                              0                            0
#> ASV1552                              0                            0
#> ASV1561                              0                            0
#> ASV1562                              0                            0
#> ASV1574                              0                            5
#> ASV1576                              0                            0
#> ASV1577                              0                            0
#> ASV1588                              0                            0
#> ASV159                               0                            0
#> ASV1609                              0                            0
#> ASV162                               0                            0
#> ASV1624                              0                            2
#> ASV1628                              0                            0
#> ASV1630                              0                            0
#> ASV1632                              0                            0
#> ASV1635                              0                            0
#> ASV1642                              0                            0
#> ASV1644                              0                            0
#> ASV1650                              0                            0
#> ASV1653                              0                            0
#> ASV1674                              0                            0
#> ASV1683                              0                            0
#> ASV1690                              0                            1
#> ASV170                              22                            4
#> ASV171                               0                            0
#> ASV1712                              0                            0
#> ASV172                               0                            0
#> ASV1726                              0                            0
#> ASV175                               0                            0
#> ASV18                                3                           29
#> ASV19                                0                         3774
#> ASV193                               0                            0
#> ASV198                               0                            0
#> ASV199                               0                            0
#> ASV2                                 1                          296
#> ASV201                               0                            0
#> ASV202                               0                         2283
#> ASV210                               0                            0
#> ASV219                               0                            0
#> ASV221                               0                            0
#> ASV223                               0                            0
#> ASV24                                0                            9
#> ASV242                               0                            0
#> ASV244                               0                            0
#> ASV248                               0                            9
#> ASV25                                0                            0
#> ASV251                             768                            0
#> ASV253                               0                            0
#> ASV255                               0                            0
#> ASV26                                0                            0
#> ASV264                               0                            0
#> ASV266                               0                            1
#> ASV27                                0                            0
#> ASV28                                0                            0
#> ASV286                               0                            0
#> ASV29                                0                            0
#> ASV292                               0                            0
#> ASV297                               0                            0
#> ASV310                               0                            0
#> ASV314                               0                            1
#> ASV320                               0                            2
#> ASV322                               0                            0
#> ASV329                               0                            0
#> ASV334                               0                            0
#> ASV338                               0                            0
#> ASV339                               0                            0
#> ASV344                               0                            0
#> ASV357                               0                            0
#> ASV359                               0                            0
#> ASV365                               0                            0
#> ASV368                               0                            0
#> ASV377                               0                            0
#> ASV378                               0                            0
#> ASV38                                0                          404
#> ASV383                               0                            0
#> ASV404                               0                            0
#> ASV41                                5                           22
#> ASV42                                0                            0
#> ASV428                               0                            0
#> ASV43                                0                            0
#> ASV431                               0                            0
#> ASV443                               0                            0
#> ASV444                               0                            0
#> ASV45                                0                            0
#> ASV454                               0                            0
#> ASV46                                0                            0
#> ASV464                               0                            0
#> ASV47                                0                            0
#> ASV477                               0                            0
#> ASV48                                0                            0
#> ASV483                               0                          134
#> ASV489                               0                            0
#> ASV49                                0                            0
#> ASV500                             435                            0
#> ASV504                               0                            0
#> ASV505                               0                            0
#> ASV507                               0                            0
#> ASV509                               0                            0
#> ASV515                               0                            0
#> ASV517                               0                            0
#> ASV52                                0                            0
#> ASV523                               0                            0
#> ASV524                               0                            0
#> ASV533                               0                            0
#> ASV543                               0                            0
#> ASV545                               0                            0
#> ASV546                               0                            0
#> ASV549                               0                            0
#> ASV55                                0                            0
#> ASV552                               0                            0
#> ASV559                               0                            0
#> ASV563                               0                           17
#> ASV569                               0                            0
#> ASV570                               0                            0
#> ASV577                               0                          220
#> ASV586                               0                            0
#> ASV59                                0                            0
#> ASV598                               0                            0
#> ASV60                                0                            0
#> ASV602                               0                            2
#> ASV618                             131                           56
#> ASV619                               0                            0
#> ASV631                              63                            0
#> ASV635                               0                            0
#> ASV637                               0                          172
#> ASV64                                0                            0
#> ASV643                               0                           10
#> ASV65                                0                            0
#> ASV662                               0                            0
#> ASV672                               0                           23
#> ASV673                               0                            0
#> ASV685                               0                            0
#> ASV69                                0                            0
#> ASV692                               0                            0
#> ASV694                               0                            0
#> ASV715                               0                            0
#> ASV717                               0                            0
#> ASV727                               0                            0
#> ASV731                               0                            0
#> ASV737                               0                            0
#> ASV741                               0                            1
#> ASV743                               0                            0
#> ASV746                               0                            8
#> ASV749                               0                            0
#> ASV753                               0                            0
#> ASV756                               0                            4
#> ASV758                               0                            0
#> ASV772                               0                           71
#> ASV78                                0                            0
#> ASV784                               0                            0
#> ASV787                               0                            0
#> ASV8                               463                          120
#> ASV801                               0                            1
#> ASV81                                0                            0
#> ASV817                               0                            0
#> ASV823                               0                            0
#> ASV829                               0                            0
#> ASV834                               0                            0
#> ASV839                               0                            0
#> ASV854                               0                            0
#> ASV87                                0                            0
#> ASV89                                0                            0
#> ASV891                               0                            0
#> ASV892                               0                            0
#> ASV898                               0                            0
#> ASV899                               0                            0
#> ASV904                               0                            1
#> ASV906                               2                            5
#> ASV91                                0                            0
#> ASV914                               0                            0
#> ASV926                               0                            0
#> ASV927                               0                            0
#> ASV930                               0                            0
#> ASV942                             128                            0
#> ASV946                               0                            0
#> ASV950                               0                            0
#> ASV953                              53                           26
#> ASV964                               0                            0
#> ASV969                               0                            0
#> ASV973                               0                            0
#> ASV987                               0                            0
#> ASV988                               0                            0
#>         A12.007.B_S2_MERGED.fastq.gz A15.004_S3_MERGED.fastq.gz
#> ASV1004                            7                          0
#> ASV1007                            0                          0
#> ASV101                             0                          0
#> ASV1017                            0                          0
#> ASV1022                            0                          0
#> ASV1030                            0                          0
#> ASV104                             0                          0
#> ASV105                             0                          0
#> ASV1059                            0                          0
#> ASV1069                            0                          0
#> ASV107                             0                          0
#> ASV1070                            0                          0
#> ASV1074                            0                          0
#> ASV1085                            0                          0
#> ASV1102                            0                          0
#> ASV1108                            0                          0
#> ASV111                             0                          0
#> ASV1118                            0                          0
#> ASV113                             0                          9
#> ASV1132                            0                          0
#> ASV1144                            5                          0
#> ASV1146                            0                          0
#> ASV1152                            0                          0
#> ASV1155                            0                          0
#> ASV1164                            0                          0
#> ASV1166                          129                          0
#> ASV1167                            0                          0
#> ASV1179                            0                          0
#> ASV1192                            0                          0
#> ASV12                              0                          0
#> ASV1200                            0                          0
#> ASV1212                            0                          0
#> ASV1218                            0                          0
#> ASV1223                            1                          0
#> ASV1230                            0                          0
#> ASV1233                            0                          0
#> ASV1234                            0                          0
#> ASV1236                           47                          0
#> ASV1239                            0                          1
#> ASV1242                            0                          0
#> ASV1253                            0                          0
#> ASV1270                            0                          0
#> ASV1278                            0                          0
#> ASV1279                            0                          0
#> ASV1300                            1                          0
#> ASV1302                            0                          0
#> ASV131                             0                          0
#> ASV1311                            0                         14
#> ASV1314                            0                          0
#> ASV1321                            0                          0
#> ASV1332                            0                          0
#> ASV1340                            0                          0
#> ASV1359                            0                          0
#> ASV1365                            0                          0
#> ASV1369                            0                          0
#> ASV137                             0                          0
#> ASV1370                            0                          0
#> ASV1379                            0                          0
#> ASV139                             0                          0
#> ASV1396                            0                          0
#> ASV1409                            0                          0
#> ASV1419                            0                          4
#> ASV1420                            0                          0
#> ASV1421                            0                          0
#> ASV1422                            0                          0
#> ASV1428                            0                          0
#> ASV1434                            0                          0
#> ASV144                             0                          0
#> ASV1458                            0                          0
#> ASV1459                            0                          0
#> ASV1460                            0                          0
#> ASV1461                            0                          0
#> ASV1462                            0                          0
#> ASV1468                            0                          0
#> ASV1483                            0                          0
#> ASV1485                            0                          0
#> ASV151                             0                          0
#> ASV1510                            0                          0
#> ASV1523                            0                          0
#> ASV1532                           42                          0
#> ASV1537                            0                          0
#> ASV1548                            0                          0
#> ASV1550                            0                          0
#> ASV1552                            0                          0
#> ASV1561                            0                          0
#> ASV1562                            0                          0
#> ASV1574                            0                          0
#> ASV1576                            0                          0
#> ASV1577                            0                          0
#> ASV1588                            0                          0
#> ASV159                             0                          0
#> ASV1609                            0                          0
#> ASV162                             0                          0
#> ASV1624                           37                          0
#> ASV1628                            0                          6
#> ASV1630                            0                          5
#> ASV1632                            0                          0
#> ASV1635                            0                          0
#> ASV1642                            0                          0
#> ASV1644                            1                          0
#> ASV1650                            0                          0
#> ASV1653                            0                          0
#> ASV1674                            0                          0
#> ASV1683                            0                          0
#> ASV1690                            0                          0
#> ASV170                             0                          0
#> ASV171                             0                          0
#> ASV1712                            0                          0
#> ASV172                             0                          0
#> ASV1726                            0                          0
#> ASV175                             0                          0
#> ASV18                           1040                          0
#> ASV19                              0                          0
#> ASV193                             0                          0
#> ASV198                             0                          1
#> ASV199                             0                          0
#> ASV2                             263                          0
#> ASV201                             0                          0
#> ASV202                             0                          0
#> ASV210                             7                          0
#> ASV219                             0                          0
#> ASV221                             0                          0
#> ASV223                             0                          0
#> ASV24                              0                          1
#> ASV242                             0                          0
#> ASV244                             0                          0
#> ASV248                             0                          0
#> ASV25                              0                          0
#> ASV251                             0                          0
#> ASV253                             0                          0
#> ASV255                             0                          0
#> ASV26                              0                          0
#> ASV264                             0                          0
#> ASV266                             0                          0
#> ASV27                              0                          0
#> ASV28                              0                          0
#> ASV286                             0                          0
#> ASV29                              0                          0
#> ASV292                             0                          0
#> ASV297                             0                          0
#> ASV310                             0                          0
#> ASV314                             0                          0
#> ASV320                             0                          0
#> ASV322                             0                          0
#> ASV329                             0                          1
#> ASV334                             0                          0
#> ASV338                             0                          0
#> ASV339                             0                          0
#> ASV344                             0                          0
#> ASV357                             0                          0
#> ASV359                             0                          0
#> ASV365                             0                          0
#> ASV368                             0                          0
#> ASV377                             0                          0
#> ASV378                             0                          0
#> ASV38                              0                          0
#> ASV383                             0                          0
#> ASV404                             0                          0
#> ASV41                              0                          0
#> ASV42                              0                          3
#> ASV428                             0                          0
#> ASV43                              0                          0
#> ASV431                             0                          0
#> ASV443                             0                          0
#> ASV444                             0                          0
#> ASV45                              0                          0
#> ASV454                             0                          0
#> ASV46                              0                          0
#> ASV464                             0                          0
#> ASV47                              0                          0
#> ASV477                             0                          0
#> ASV48                              0                          0
#> ASV483                             0                          0
#> ASV489                             0                          0
#> ASV49                              0                          0
#> ASV500                             0                          0
#> ASV504                             0                          0
#> ASV505                             0                         59
#> ASV507                             0                          0
#> ASV509                             0                         32
#> ASV515                             0                          0
#> ASV517                             0                          0
#> ASV52                              0                          0
#> ASV523                             0                          2
#> ASV524                             0                          0
#> ASV533                             0                          0
#> ASV543                             0                          0
#> ASV545                             0                          0
#> ASV546                             0                          0
#> ASV549                             0                          0
#> ASV55                              0                          0
#> ASV552                             0                          0
#> ASV559                             0                          0
#> ASV563                             0                          0
#> ASV569                             0                          0
#> ASV570                             0                          0
#> ASV577                             0                          0
#> ASV586                             0                          1
#> ASV59                              0                          0
#> ASV598                             0                          0
#> ASV60                              0                         21
#> ASV602                             0                          0
#> ASV618                             0                          0
#> ASV619                             0                          0
#> ASV631                             0                          0
#> ASV635                             0                          0
#> ASV637                             0                          0
#> ASV64                              0                          0
#> ASV643                             0                          0
#> ASV65                              0                          0
#> ASV662                             6                          1
#> ASV672                             0                          0
#> ASV673                             0                          0
#> ASV685                             0                          0
#> ASV69                              0                          2
#> ASV692                             0                          0
#> ASV694                             0                          0
#> ASV715                             0                          0
#> ASV717                             0                          0
#> ASV727                             0                          0
#> ASV731                             0                          0
#> ASV737                             0                          0
#> ASV741                             0                          0
#> ASV743                             0                          0
#> ASV746                             0                          0
#> ASV749                             7                          0
#> ASV753                             0                          6
#> ASV756                             0                          0
#> ASV758                             0                          6
#> ASV772                             0                          0
#> ASV78                              0                          0
#> ASV784                             0                          0
#> ASV787                             0                          0
#> ASV8                             841                          2
#> ASV801                             0                          0
#> ASV81                              0                          0
#> ASV817                             0                          0
#> ASV823                             0                          0
#> ASV829                             0                          0
#> ASV834                             0                          0
#> ASV839                             0                          0
#> ASV854                             0                          0
#> ASV87                              0                          0
#> ASV89                              0                          0
#> ASV891                             0                          0
#> ASV892                             0                          1
#> ASV898                             0                          0
#> ASV899                             0                          0
#> ASV904                             0                          0
#> ASV906                             0                          0
#> ASV91                              0                         12
#> ASV914                             0                          2
#> ASV926                             0                          0
#> ASV927                             0                          0
#> ASV930                             6                          0
#> ASV942                             0                          0
#> ASV946                             0                          0
#> ASV950                             0                          0
#> ASV953                             0                          0
#> ASV964                             0                          0
#> ASV969                             0                          0
#> ASV973                             0                          0
#> ASV987                             0                          0
#> ASV988                             0                          0
#>         A8.005_S4_MERGED.fastq.gz A9.012_S5_MERGED.fastq.gz
#> ASV1004                         0                         0
#> ASV1007                         0                         0
#> ASV101                          0                         0
#> ASV1017                         0                         0
#> ASV1022                         0                         0
#> ASV1030                         0                         3
#> ASV104                          0                         0
#> ASV105                          0                         0
#> ASV1059                         0                         0
#> ASV1069                         0                         0
#> ASV107                          0                         0
#> ASV1070                       153                         0
#> ASV1074                         0                         0
#> ASV1085                         0                         0
#> ASV1102                         0                         0
#> ASV1108                         0                         0
#> ASV111                          0                         0
#> ASV1118                         0                         0
#> ASV113                          0                         0
#> ASV1132                         0                         0
#> ASV1144                         0                         0
#> ASV1146                         0                         0
#> ASV1152                         0                         0
#> ASV1155                         0                         0
#> ASV1164                         0                         0
#> ASV1166                         0                         0
#> ASV1167                         0                         0
#> ASV1179                         0                         0
#> ASV1192                         0                         0
#> ASV12                           1                         0
#> ASV1200                         0                         0
#> ASV1212                         0                         0
#> ASV1218                         0                         1
#> ASV1223                         0                         0
#> ASV1230                         1                         0
#> ASV1233                         0                         0
#> ASV1234                         0                         0
#> ASV1236                         0                         0
#> ASV1239                         0                         0
#> ASV1242                         0                         0
#> ASV1253                         0                         0
#> ASV1270                         0                         0
#> ASV1278                         0                         0
#> ASV1279                         0                         0
#> ASV1300                         0                         0
#> ASV1302                         0                         0
#> ASV131                          0                         0
#> ASV1311                         0                         0
#> ASV1314                         0                         0
#> ASV1321                         0                         0
#> ASV1332                         0                         0
#> ASV1340                         0                         0
#> ASV1359                         0                         0
#> ASV1365                         1                         0
#> ASV1369                         0                         0
#> ASV137                          0                         0
#> ASV1370                         0                         0
#> ASV1379                         0                         0
#> ASV139                          0                         0
#> ASV1396                         7                         0
#> ASV1409                         0                         0
#> ASV1419                         0                         0
#> ASV1420                         0                        70
#> ASV1421                         0                         0
#> ASV1422                         0                         0
#> ASV1428                         0                         0
#> ASV1434                         0                         0
#> ASV144                          0                         0
#> ASV1458                         0                         0
#> ASV1459                         0                         0
#> ASV1460                         0                         0
#> ASV1461                         0                         0
#> ASV1462                         0                         0
#> ASV1468                         0                         0
#> ASV1483                         0                         0
#> ASV1485                         0                         0
#> ASV151                          0                         0
#> ASV1510                         0                        39
#> ASV1523                         0                         0
#> ASV1532                         0                         0
#> ASV1537                         0                         0
#> ASV1548                         0                         0
#> ASV1550                         0                         0
#> ASV1552                         0                         0
#> ASV1561                         0                         0
#> ASV1562                         0                         0
#> ASV1574                         0                         0
#> ASV1576                         0                         2
#> ASV1577                         0                         0
#> ASV1588                         0                         0
#> ASV159                          0                         0
#> ASV1609                         0                         0
#> ASV162                          0                         0
#> ASV1624                         0                         0
#> ASV1628                         0                         0
#> ASV1630                         0                         0
#> ASV1632                         0                         0
#> ASV1635                         0                         0
#> ASV1642                         0                        36
#> ASV1644                         0                         0
#> ASV1650                         0                         0
#> ASV1653                         0                         0
#> ASV1674                         0                         0
#> ASV1683                         0                         0
#> ASV1690                         0                         0
#> ASV170                        218                         0
#> ASV171                          0                         0
#> ASV1712                         0                         0
#> ASV172                          0                         0
#> ASV1726                         0                         0
#> ASV175                          0                         0
#> ASV18                          42                         0
#> ASV19                           0                         0
#> ASV193                          0                         0
#> ASV198                          0                         0
#> ASV199                          0                         0
#> ASV2                            5                         1
#> ASV201                          0                         0
#> ASV202                          0                         0
#> ASV210                          0                         0
#> ASV219                          0                         0
#> ASV221                          0                         0
#> ASV223                          0                         0
#> ASV24                           0                         0
#> ASV242                          0                      1692
#> ASV244                          0                         0
#> ASV248                          0                         0
#> ASV25                           0                         0
#> ASV251                          0                         0
#> ASV253                          0                         0
#> ASV255                          0                         0
#> ASV26                           0                         0
#> ASV264                          0                         0
#> ASV266                          0                         0
#> ASV27                           0                         0
#> ASV28                           0                         0
#> ASV286                          0                         0
#> ASV29                           0                         0
#> ASV292                          0                         0
#> ASV297                          0                         0
#> ASV310                          0                         0
#> ASV314                          0                         0
#> ASV320                          0                         0
#> ASV322                          0                         0
#> ASV329                          5                         0
#> ASV334                          0                         0
#> ASV338                          0                         0
#> ASV339                          0                         0
#> ASV344                          0                         0
#> ASV357                          0                         0
#> ASV359                         11                         0
#> ASV365                          0                         0
#> ASV368                          0                         0
#> ASV377                          0                         0
#> ASV378                          0                         0
#> ASV38                           0                         0
#> ASV383                          0                         0
#> ASV404                          0                         0
#> ASV41                           0                         0
#> ASV42                           0                         0
#> ASV428                          0                         0
#> ASV43                           0                         0
#> ASV431                         65                         0
#> ASV443                          0                         0
#> ASV444                          0                         0
#> ASV45                           0                         0
#> ASV454                          0                         0
#> ASV46                           0                         0
#> ASV464                          0                         0
#> ASV47                           0                         0
#> ASV477                         57                         0
#> ASV48                           0                         0
#> ASV483                          0                         0
#> ASV489                          0                         0
#> ASV49                           0                         0
#> ASV500                          0                         1
#> ASV504                        535                         0
#> ASV505                          0                         3
#> ASV507                          0                         0
#> ASV509                          0                        26
#> ASV515                          0                         2
#> ASV517                          0                         0
#> ASV52                           0                         0
#> ASV523                          0                         0
#> ASV524                          0                         0
#> ASV533                          0                         0
#> ASV543                          0                         0
#> ASV545                          0                         0
#> ASV546                          0                         0
#> ASV549                          0                         0
#> ASV55                           0                         0
#> ASV552                          0                         0
#> ASV559                          0                         0
#> ASV563                          2                         0
#> ASV569                          0                         0
#> ASV570                          0                         0
#> ASV577                          0                         0
#> ASV586                          0                         1
#> ASV59                           0                         0
#> ASV598                          0                         0
#> ASV60                           0                         0
#> ASV602                          0                         0
#> ASV618                          0                         0
#> ASV619                          0                         0
#> ASV631                         28                        21
#> ASV635                          0                         0
#> ASV637                          0                         0
#> ASV64                         116                         0
#> ASV643                          0                         0
#> ASV65                           0                         0
#> ASV662                          0                         4
#> ASV672                          0                         0
#> ASV673                          0                         0
#> ASV685                          0                         0
#> ASV69                           0                         0
#> ASV692                          0                         0
#> ASV694                          0                         0
#> ASV715                          3                         0
#> ASV717                          0                         0
#> ASV727                          0                         0
#> ASV731                          0                         0
#> ASV737                          0                         0
#> ASV741                          0                         0
#> ASV743                          0                         0
#> ASV746                          0                         0
#> ASV749                          0                         0
#> ASV753                          0                         0
#> ASV756                          2                         0
#> ASV758                          0                         9
#> ASV772                          0                         0
#> ASV78                           0                         0
#> ASV784                          0                         0
#> ASV787                          0                         8
#> ASV8                         2214                         0
#> ASV801                          0                         0
#> ASV81                           0                         0
#> ASV817                          0                         0
#> ASV823                          0                         0
#> ASV829                          0                         0
#> ASV834                          0                         0
#> ASV839                          0                         0
#> ASV854                          0                         0
#> ASV87                           0                         0
#> ASV89                        2358                         0
#> ASV891                          0                         0
#> ASV892                          0                         0
#> ASV898                          0                       292
#> ASV899                          0                         0
#> ASV904                          1                         0
#> ASV906                          0                         0
#> ASV91                           0                         0
#> ASV914                          0                         0
#> ASV926                          0                         0
#> ASV927                          0                         0
#> ASV930                          0                         0
#> ASV942                          0                         0
#> ASV946                          0                         0
#> ASV950                          0                         0
#> ASV953                          0                         0
#> ASV964                          0                         0
#> ASV969                          0                         0
#> ASV973                          0                         0
#> ASV987                          0                         0
#> ASV988                          0                         0
#>         AB29.ABMX.H_S6_MERGED.fastq.gz AC27.013_S7_MERGED.fastq.gz
#> ASV1004                              0                           0
#> ASV1007                              0                           0
#> ASV101                               0                           0
#> ASV1017                              0                           0
#> ASV1022                              0                           0
#> ASV1030                              0                           0
#> ASV104                               0                           0
#> ASV105                              68                           0
#> ASV1059                              0                           0
#> ASV1069                              0                           0
#> ASV107                               0                          96
#> ASV1070                              0                           0
#> ASV1074                              0                           0
#> ASV1085                              0                           0
#> ASV1102                              0                           0
#> ASV1108                              0                           0
#> ASV111                               0                           0
#> ASV1118                              0                           0
#> ASV113                               0                          14
#> ASV1132                              0                         136
#> ASV1144                              0                           0
#> ASV1146                            132                           0
#> ASV1152                              0                           0
#> ASV1155                              0                           0
#> ASV1164                              0                           0
#> ASV1166                              0                           0
#> ASV1167                              0                           0
#> ASV1179                              9                           0
#> ASV1192                              0                           0
#> ASV12                                1                           1
#> ASV1200                              0                           0
#> ASV1212                              0                           0
#> ASV1218                              0                           0
#> ASV1223                              0                           0
#> ASV1230                              0                           0
#> ASV1233                              0                           0
#> ASV1234                              0                           0
#> ASV1236                              0                           0
#> ASV1239                              0                           0
#> ASV1242                              0                           0
#> ASV1253                              0                           0
#> ASV1270                              0                           0
#> ASV1278                              0                           0
#> ASV1279                              0                           0
#> ASV1300                              0                           0
#> ASV1302                              0                           0
#> ASV131                               0                           0
#> ASV1311                              0                           0
#> ASV1314                              0                           0
#> ASV1321                              0                           0
#> ASV1332                              0                           0
#> ASV1340                             26                           0
#> ASV1359                              0                           0
#> ASV1365                              0                           0
#> ASV1369                              0                           0
#> ASV137                               0                           0
#> ASV1370                              0                           0
#> ASV1379                              0                           0
#> ASV139                               0                           0
#> ASV1396                              0                           0
#> ASV1409                              0                           0
#> ASV1419                              0                           0
#> ASV1420                              0                           0
#> ASV1421                              3                           0
#> ASV1422                              0                           0
#> ASV1428                              9                           0
#> ASV1434                              0                           0
#> ASV144                               0                           0
#> ASV1458                              0                           0
#> ASV1459                              0                           0
#> ASV1460                              0                           0
#> ASV1461                              0                           0
#> ASV1462                              0                           0
#> ASV1468                              0                           0
#> ASV1483                              0                           0
#> ASV1485                              0                           0
#> ASV151                              21                           0
#> ASV1510                              0                           0
#> ASV1523                              0                           0
#> ASV1532                              0                           0
#> ASV1537                              0                           0
#> ASV1548                              0                           0
#> ASV1550                              0                           0
#> ASV1552                              0                           0
#> ASV1561                              0                           0
#> ASV1562                              0                           0
#> ASV1574                              0                           0
#> ASV1576                              0                           0
#> ASV1577                              0                           8
#> ASV1588                              0                           0
#> ASV159                               0                           0
#> ASV1609                              0                           0
#> ASV162                               2                           0
#> ASV1624                              0                           0
#> ASV1628                              0                           0
#> ASV1630                              0                           0
#> ASV1632                              0                           0
#> ASV1635                              0                           0
#> ASV1642                              0                           0
#> ASV1644                              0                           0
#> ASV1650                              0                           0
#> ASV1653                              0                           0
#> ASV1674                              0                           1
#> ASV1683                              0                           0
#> ASV1690                              0                           0
#> ASV170                               0                           0
#> ASV171                               0                           0
#> ASV1712                              0                           0
#> ASV172                               0                           0
#> ASV1726                              0                           0
#> ASV175                            1295                           0
#> ASV18                                0                           0
#> ASV19                                0                           0
#> ASV193                               0                           0
#> ASV198                               0                           0
#> ASV199                               0                           0
#> ASV2                                 5                           5
#> ASV201                               0                           0
#> ASV202                               0                           0
#> ASV210                               0                           0
#> ASV219                            1961                           0
#> ASV221                               0                           0
#> ASV223                               0                           0
#> ASV24                                0                           0
#> ASV242                               0                           0
#> ASV244                               0                           0
#> ASV248                               0                           0
#> ASV25                                0                           0
#> ASV251                               0                           0
#> ASV253                               0                           0
#> ASV255                               0                           0
#> ASV26                                0                           0
#> ASV264                               0                           0
#> ASV266                               0                           0
#> ASV27                                0                           6
#> ASV28                                0                           0
#> ASV286                               0                           0
#> ASV29                                0                           0
#> ASV292                               0                           0
#> ASV297                               0                           1
#> ASV310                              12                           0
#> ASV314                               0                           0
#> ASV320                               0                           0
#> ASV322                             957                           0
#> ASV329                               0                           4
#> ASV334                               0                           0
#> ASV338                               0                           0
#> ASV339                               0                           0
#> ASV344                               0                           0
#> ASV357                               0                           0
#> ASV359                               0                           0
#> ASV365                             705                           0
#> ASV368                               2                           0
#> ASV377                               0                           0
#> ASV378                               0                           0
#> ASV38                                0                           0
#> ASV383                               0                           0
#> ASV404                               0                           0
#> ASV41                                0                           0
#> ASV42                                0                           0
#> ASV428                              26                           0
#> ASV43                                0                           0
#> ASV431                               0                           0
#> ASV443                               0                           0
#> ASV444                               0                           0
#> ASV45                              877                           0
#> ASV454                               0                           0
#> ASV46                             1230                           0
#> ASV464                               0                           0
#> ASV47                                0                           0
#> ASV477                               0                           0
#> ASV48                                0                           0
#> ASV483                               0                           0
#> ASV489                               0                           0
#> ASV49                                0                           0
#> ASV500                               0                           0
#> ASV504                               0                           0
#> ASV505                               0                           2
#> ASV507                               0                           0
#> ASV509                               0                          40
#> ASV515                               0                           0
#> ASV517                               0                           0
#> ASV52                                0                           0
#> ASV523                               0                          58
#> ASV524                             531                           0
#> ASV533                               0                           0
#> ASV543                               0                           0
#> ASV545                               0                           0
#> ASV546                               8                           0
#> ASV549                               0                           0
#> ASV55                                0                           0
#> ASV552                               0                           0
#> ASV559                               0                           0
#> ASV563                               0                           0
#> ASV569                               0                           0
#> ASV570                               0                           0
#> ASV577                             444                           0
#> ASV586                               0                           0
#> ASV59                                0                           0
#> ASV598                               0                           0
#> ASV60                                0                           0
#> ASV602                               0                           0
#> ASV618                               0                           0
#> ASV619                               0                           0
#> ASV631                               0                           0
#> ASV635                               0                           0
#> ASV637                               0                           0
#> ASV64                                0                       11079
#> ASV643                               0                           0
#> ASV65                             8505                           0
#> ASV662                               0                           0
#> ASV672                               0                           0
#> ASV673                               0                          15
#> ASV685                               0                           0
#> ASV69                                0                           0
#> ASV692                               0                           0
#> ASV694                             101                           0
#> ASV715                               0                           0
#> ASV717                               0                           0
#> ASV727                               4                           0
#> ASV731                               0                           0
#> ASV737                               0                           0
#> ASV741                               0                           0
#> ASV743                               0                           0
#> ASV746                               0                           0
#> ASV749                               0                          50
#> ASV753                               0                           0
#> ASV756                               0                           0
#> ASV758                               0                           0
#> ASV772                               0                           0
#> ASV78                                0                        6098
#> ASV784                               0                           0
#> ASV787                               0                           0
#> ASV8                                 1                           3
#> ASV801                               0                           0
#> ASV81                                0                           0
#> ASV817                              13                           0
#> ASV823                               0                           0
#> ASV829                               0                           0
#> ASV834                               0                           0
#> ASV839                               0                           0
#> ASV854                               0                           0
#> ASV87                                0                           0
#> ASV89                                0                           0
#> ASV891                               0                           0
#> ASV892                               0                           0
#> ASV898                               0                           0
#> ASV899                             210                           0
#> ASV904                               0                           0
#> ASV906                               0                           0
#> ASV91                                0                           0
#> ASV914                               0                           0
#> ASV926                               0                           0
#> ASV927                               0                           0
#> ASV930                               0                           0
#> ASV942                               0                           0
#> ASV946                               0                           0
#> ASV950                               0                           0
#> ASV953                               0                           0
#> ASV964                               0                           0
#> ASV969                               0                           0
#> ASV973                               0                           0
#> ASV987                               0                         176
#> ASV988                               0                           0
#>         AC29033_S8_MERGED.fastq.gz AD26.005.B_S9_MERGED.fastq.gz
#> ASV1004                          1                             0
#> ASV1007                          0                             0
#> ASV101                           0                             0
#> ASV1017                          0                             0
#> ASV1022                         14                             0
#> ASV1030                          0                             0
#> ASV104                           0                             0
#> ASV105                           0                             0
#> ASV1059                          0                             0
#> ASV1069                          0                             0
#> ASV107                           0                             0
#> ASV1070                          0                             0
#> ASV1074                          0                             0
#> ASV1085                          0                             0
#> ASV1102                          0                             0
#> ASV1108                          0                             0
#> ASV111                           0                             0
#> ASV1118                          0                             0
#> ASV113                          90                             0
#> ASV1132                          0                             0
#> ASV1144                          0                             0
#> ASV1146                          0                             0
#> ASV1152                          0                             0
#> ASV1155                          0                             0
#> ASV1164                          0                            10
#> ASV1166                          0                             0
#> ASV1167                          0                             3
#> ASV1179                          0                             0
#> ASV1192                          3                             0
#> ASV12                            0                             0
#> ASV1200                          0                             0
#> ASV1212                          0                             0
#> ASV1218                          0                             5
#> ASV1223                          1                             0
#> ASV1230                          0                             0
#> ASV1233                         38                             0
#> ASV1234                          0                             0
#> ASV1236                          0                             0
#> ASV1239                          0                             0
#> ASV1242                          0                             0
#> ASV1253                          3                             0
#> ASV1270                          0                             0
#> ASV1278                          0                             0
#> ASV1279                          0                            24
#> ASV1300                          3                             0
#> ASV1302                          0                             1
#> ASV131                           0                             0
#> ASV1311                          0                             0
#> ASV1314                          0                             0
#> ASV1321                          0                             3
#> ASV1332                          1                             0
#> ASV1340                          0                             0
#> ASV1359                          0                             1
#> ASV1365                          0                             0
#> ASV1369                          6                             1
#> ASV137                           0                             0
#> ASV1370                          0                             0
#> ASV1379                          0                            38
#> ASV139                           0                           142
#> ASV1396                          0                            12
#> ASV1409                          0                             0
#> ASV1419                          0                             0
#> ASV1420                          0                             0
#> ASV1421                          1                             0
#> ASV1422                          2                             0
#> ASV1428                          0                             0
#> ASV1434                          0                             1
#> ASV144                           0                             0
#> ASV1458                          0                             0
#> ASV1459                          0                            48
#> ASV1460                          0                             0
#> ASV1461                          0                             0
#> ASV1462                          0                             0
#> ASV1468                          0                             0
#> ASV1483                          0                             1
#> ASV1485                          0                            73
#> ASV151                           0                             1
#> ASV1510                          0                             0
#> ASV1523                          2                             0
#> ASV1532                          8                             0
#> ASV1537                          0                             0
#> ASV1548                          0                             0
#> ASV1550                          1                             0
#> ASV1552                          0                             0
#> ASV1561                          0                             0
#> ASV1562                          0                             0
#> ASV1574                          0                             0
#> ASV1576                          0                             0
#> ASV1577                          0                             0
#> ASV1588                          0                             0
#> ASV159                           0                            32
#> ASV1609                          0                             0
#> ASV162                           0                             0
#> ASV1624                          1                             0
#> ASV1628                          0                             0
#> ASV1630                          0                             0
#> ASV1632                          0                            10
#> ASV1635                          0                             0
#> ASV1642                          0                             0
#> ASV1644                          0                             0
#> ASV1650                          0                             0
#> ASV1653                          0                             0
#> ASV1674                          0                             0
#> ASV1683                          0                             0
#> ASV1690                          0                             0
#> ASV170                           0                             0
#> ASV171                           0                             0
#> ASV1712                          0                             0
#> ASV172                           0                             0
#> ASV1726                          0                             2
#> ASV175                           0                            73
#> ASV18                            0                             0
#> ASV19                            0                             1
#> ASV193                           1                             0
#> ASV198                           0                             0
#> ASV199                          12                             0
#> ASV2                            19                             8
#> ASV201                           0                          1784
#> ASV202                           0                             5
#> ASV210                           0                             0
#> ASV219                           0                             0
#> ASV221                        2175                             0
#> ASV223                           0                             0
#> ASV24                            0                             0
#> ASV242                           0                             0
#> ASV244                           0                            13
#> ASV248                           0                             0
#> ASV25                            0                             0
#> ASV251                           0                             0
#> ASV253                           0                             0
#> ASV255                           0                             0
#> ASV26                            0                             0
#> ASV264                           0                             0
#> ASV266                           1                             0
#> ASV27                            0                             0
#> ASV28                          114                             0
#> ASV286                           2                             0
#> ASV29                            0                             0
#> ASV292                           1                             0
#> ASV297                          13                             0
#> ASV310                          12                             2
#> ASV314                           0                             0
#> ASV320                           0                             0
#> ASV322                           0                             0
#> ASV329                          67                             0
#> ASV334                           0                             0
#> ASV338                           0                             0
#> ASV339                           0                             0
#> ASV344                           2                             0
#> ASV357                           0                             0
#> ASV359                           0                             0
#> ASV365                           0                             0
#> ASV368                           0                             0
#> ASV377                           0                             0
#> ASV378                           1                             0
#> ASV38                            0                           617
#> ASV383                           4                             0
#> ASV404                           0                             0
#> ASV41                            1                             0
#> ASV42                            0                             0
#> ASV428                           0                             0
#> ASV43                            0                         11426
#> ASV431                           0                             0
#> ASV443                           5                             0
#> ASV444                         802                             0
#> ASV45                            0                             0
#> ASV454                           0                             0
#> ASV46                            0                             1
#> ASV464                           1                             0
#> ASV47                            0                             0
#> ASV477                           0                             0
#> ASV48                            0                             0
#> ASV483                           0                             0
#> ASV489                           0                             0
#> ASV49                            0                             0
#> ASV500                           0                             0
#> ASV504                           0                             0
#> ASV505                           0                             0
#> ASV507                           0                             8
#> ASV509                           0                             0
#> ASV515                           0                             0
#> ASV517                           0                             0
#> ASV52                            0                            86
#> ASV523                           0                             0
#> ASV524                           0                             0
#> ASV533                           0                           434
#> ASV543                           5                             0
#> ASV545                           0                             0
#> ASV546                          10                             0
#> ASV549                           0                             2
#> ASV55                            0                            11
#> ASV552                           0                             0
#> ASV559                           3                             0
#> ASV563                           0                             0
#> ASV569                           0                             4
#> ASV570                           0                             0
#> ASV577                           0                             0
#> ASV586                           0                             0
#> ASV59                            0                             0
#> ASV598                           0                             0
#> ASV60                           53                             0
#> ASV602                           0                             0
#> ASV618                           0                             5
#> ASV619                           0                             0
#> ASV631                           0                             0
#> ASV635                           1                             9
#> ASV637                           0                             0
#> ASV64                            0                             0
#> ASV643                           0                             0
#> ASV65                            1                             0
#> ASV662                           2                             0
#> ASV672                           0                             0
#> ASV673                           0                             0
#> ASV685                           0                             0
#> ASV69                            6                             0
#> ASV692                         369                             0
#> ASV694                           0                            16
#> ASV715                           0                             0
#> ASV717                           3                             1
#> ASV727                           1                             0
#> ASV731                           0                             0
#> ASV737                           2                             0
#> ASV741                           0                             0
#> ASV743                           0                            25
#> ASV746                           0                             0
#> ASV749                           0                             0
#> ASV753                           0                             0
#> ASV756                           0                             1
#> ASV758                           0                             0
#> ASV772                           0                             0
#> ASV78                            0                             0
#> ASV784                           0                             0
#> ASV787                           0                             0
#> ASV8                             1                             0
#> ASV801                           0                             0
#> ASV81                            0                             0
#> ASV817                           0                             2
#> ASV823                           0                             0
#> ASV829                           0                             0
#> ASV834                           0                             0
#> ASV839                           0                             0
#> ASV854                           0                             0
#> ASV87                            0                             1
#> ASV89                          281                             0
#> ASV891                           0                             0
#> ASV892                           0                             0
#> ASV898                           0                             0
#> ASV899                           0                             0
#> ASV904                           0                             0
#> ASV906                           0                             0
#> ASV91                            0                             0
#> ASV914                           0                             0
#> ASV926                           3                             0
#> ASV927                           4                             0
#> ASV930                           0                             0
#> ASV942                           0                             2
#> ASV946                           0                             0
#> ASV950                           0                             0
#> ASV953                           0                             8
#> ASV964                           0                             0
#> ASV969                           0                             0
#> ASV973                           0                             0
#> ASV987                           0                             0
#> ASV988                           0                             0
#>         AD26.005.H_S10_MERGED.fastq.gz AD26.005.M_S11_MERGED.fastq.gz
#> ASV1004                              0                              0
#> ASV1007                              0                              0
#> ASV101                            4669                              0
#> ASV1017                              0                              0
#> ASV1022                             23                              2
#> ASV1030                              0                              0
#> ASV104                               0                              0
#> ASV105                               0                              0
#> ASV1059                             38                              2
#> ASV1069                              0                              4
#> ASV107                               0                           1632
#> ASV1070                              0                              0
#> ASV1074                              1                              0
#> ASV1085                              0                              0
#> ASV1102                              0                              0
#> ASV1108                              0                              0
#> ASV111                               0                              0
#> ASV1118                              0                            140
#> ASV113                               0                              0
#> ASV1132                              0                              0
#> ASV1144                              0                              0
#> ASV1146                              0                              0
#> ASV1152                              0                              0
#> ASV1155                              0                              1
#> ASV1164                            117                              0
#> ASV1166                              0                              0
#> ASV1167                              0                              0
#> ASV1179                              0                              0
#> ASV1192                              3                              0
#> ASV12                                0                           1338
#> ASV1200                              0                              0
#> ASV1212                              1                              0
#> ASV1218                              0                              0
#> ASV1223                              0                              0
#> ASV1230                              0                              0
#> ASV1233                              0                              0
#> ASV1234                              0                              0
#> ASV1236                              0                              0
#> ASV1239                              0                              0
#> ASV1242                              0                              0
#> ASV1253                              2                              0
#> ASV1270                              0                              0
#> ASV1278                              0                              2
#> ASV1279                              0                              0
#> ASV1300                              0                              0
#> ASV1302                              0                              0
#> ASV131                               0                              0
#> ASV1311                              0                              0
#> ASV1314                              0                             97
#> ASV1321                              7                              0
#> ASV1332                              0                              0
#> ASV1340                              0                              0
#> ASV1359                              0                              0
#> ASV1365                              0                              0
#> ASV1369                              0                              0
#> ASV137                               6                              0
#> ASV1370                              3                              3
#> ASV1379                              0                              0
#> ASV139                               0                              0
#> ASV1396                              0                              0
#> ASV1409                              0                              0
#> ASV1419                              0                              0
#> ASV1420                              0                              0
#> ASV1421                              8                              7
#> ASV1422                              0                              0
#> ASV1428                              0                              0
#> ASV1434                              0                              0
#> ASV144                               0                            672
#> ASV1458                              0                              0
#> ASV1459                             16                              0
#> ASV1460                              0                              0
#> ASV1461                              0                              0
#> ASV1462                              0                              0
#> ASV1468                              0                              0
#> ASV1483                              2                              0
#> ASV1485                              0                              0
#> ASV151                               0                              0
#> ASV1510                              0                              0
#> ASV1523                              0                              0
#> ASV1532                              0                              0
#> ASV1537                              0                              0
#> ASV1548                              0                              0
#> ASV1550                              0                              0
#> ASV1552                              0                              0
#> ASV1561                              0                              0
#> ASV1562                              0                              0
#> ASV1574                              0                              0
#> ASV1576                              0                              0
#> ASV1577                              0                              0
#> ASV1588                              0                             12
#> ASV159                               0                             14
#> ASV1609                              0                              0
#> ASV162                               0                              0
#> ASV1624                              0                              0
#> ASV1628                              0                              0
#> ASV1630                              0                              0
#> ASV1632                              0                              0
#> ASV1635                              1                              0
#> ASV1642                              0                              0
#> ASV1644                              0                             17
#> ASV1650                              0                              0
#> ASV1653                              0                              0
#> ASV1674                              0                              0
#> ASV1683                              0                             18
#> ASV1690                              0                              0
#> ASV170                               2                              0
#> ASV171                               0                              0
#> ASV1712                              0                              0
#> ASV172                               3                             21
#> ASV1726                              0                              0
#> ASV175                               0                              0
#> ASV18                                0                              0
#> ASV19                                0                              0
#> ASV193                               0                              0
#> ASV198                               0                              0
#> ASV199                               0                              0
#> ASV2                                 7                              5
#> ASV201                               0                              0
#> ASV202                               0                              0
#> ASV210                               0                              0
#> ASV219                               0                              0
#> ASV221                               0                              0
#> ASV223                               0                              0
#> ASV24                                0                              0
#> ASV242                               0                              0
#> ASV244                               0                              0
#> ASV248                               0                              0
#> ASV25                                0                              0
#> ASV251                               0                              0
#> ASV253                               0                              0
#> ASV255                               0                              0
#> ASV26                                0                              0
#> ASV264                               0                              0
#> ASV266                               0                              0
#> ASV27                                0                              0
#> ASV28                                0                              0
#> ASV286                               0                              0
#> ASV29                                0                              0
#> ASV292                               0                             47
#> ASV297                               0                              0
#> ASV310                               0                              0
#> ASV314                               0                              0
#> ASV320                               0                              0
#> ASV322                               0                              0
#> ASV329                               0                              0
#> ASV334                               2                              0
#> ASV338                               0                              0
#> ASV339                               0                              0
#> ASV344                             147                            385
#> ASV357                               0                              0
#> ASV359                               0                              0
#> ASV365                               0                              0
#> ASV368                               0                             28
#> ASV377                               0                             60
#> ASV378                               2                             10
#> ASV38                                0                             12
#> ASV383                               0                              0
#> ASV404                               0                              0
#> ASV41                                0                              0
#> ASV42                                0                              0
#> ASV428                               0                              0
#> ASV43                                0                              0
#> ASV431                               0                              0
#> ASV443                               0                             72
#> ASV444                               0                              0
#> ASV45                                0                              0
#> ASV454                               0                            646
#> ASV46                                0                              0
#> ASV464                               0                              2
#> ASV47                                0                              0
#> ASV477                               0                              0
#> ASV48                                0                              0
#> ASV483                               0                              0
#> ASV489                               0                              0
#> ASV49                                0                              0
#> ASV500                               0                              0
#> ASV504                               0                              0
#> ASV505                               0                              0
#> ASV507                               0                              0
#> ASV509                               0                              0
#> ASV515                               0                              0
#> ASV517                               0                              0
#> ASV52                                0                              6
#> ASV523                               0                              0
#> ASV524                               0                              0
#> ASV533                               0                              0
#> ASV543                               0                              0
#> ASV545                              23                              0
#> ASV546                               0                              0
#> ASV549                               0                              0
#> ASV55                                0                             86
#> ASV552                               0                              0
#> ASV559                               1                              0
#> ASV563                               0                              0
#> ASV569                               0                              0
#> ASV570                               0                              0
#> ASV577                               0                              0
#> ASV586                               6                              0
#> ASV59                                0                              0
#> ASV598                               3                              0
#> ASV60                                9                              0
#> ASV602                               0                             36
#> ASV618                             408                              0
#> ASV619                               0                              0
#> ASV631                               0                              0
#> ASV635                               0                              0
#> ASV637                               0                              0
#> ASV64                                0                              0
#> ASV643                               0                              0
#> ASV65                                0                              0
#> ASV662                               0                              0
#> ASV672                               0                              0
#> ASV673                               0                              0
#> ASV685                               0                              0
#> ASV69                                0                              0
#> ASV692                               0                              0
#> ASV694                               0                              0
#> ASV715                               0                              0
#> ASV717                              20                              4
#> ASV727                               0                              0
#> ASV731                               0                              0
#> ASV737                               0                              0
#> ASV741                               0                              0
#> ASV743                               0                             22
#> ASV746                               0                              0
#> ASV749                               0                              0
#> ASV753                               0                              0
#> ASV756                               9                              4
#> ASV758                               0                              0
#> ASV772                               0                              0
#> ASV78                                0                              0
#> ASV784                               0                              0
#> ASV787                               6                              0
#> ASV8                                 0                              0
#> ASV801                               0                              2
#> ASV81                                0                              0
#> ASV817                               0                             78
#> ASV823                               0                              0
#> ASV829                               0                              0
#> ASV834                               0                              0
#> ASV839                               0                              0
#> ASV854                               0                             13
#> ASV87                                0                              0
#> ASV89                              131                              0
#> ASV891                               0                              0
#> ASV892                               0                              0
#> ASV898                               0                              0
#> ASV899                               0                              0
#> ASV904                               0                             33
#> ASV906                              24                              0
#> ASV91                                0                              0
#> ASV914                               0                              0
#> ASV926                             104                             58
#> ASV927                               0                              0
#> ASV930                               0                              0
#> ASV942                               0                              0
#> ASV946                               0                              0
#> ASV950                               1                              0
#> ASV953                              77                              0
#> ASV964                               0                              0
#> ASV969                               0                              0
#> ASV973                               0                              0
#> ASV987                               0                              0
#> ASV988                               0                              0
#>         AD30.ABMX.M_S12_MERGED.fastq.gz AD32.007.M_S13_MERGED.fastq.gz
#> ASV1004                               0                              0
#> ASV1007                               0                              4
#> ASV101                                0                              0
#> ASV1017                               0                              0
#> ASV1022                               2                              0
#> ASV1030                               0                              0
#> ASV104                                0                              0
#> ASV105                                0                              0
#> ASV1059                               0                              0
#> ASV1069                               0                              0
#> ASV107                                0                              0
#> ASV1070                               0                              0
#> ASV1074                               0                              1
#> ASV1085                               0                              0
#> ASV1102                               0                              0
#> ASV1108                               0                              1
#> ASV111                                0                              0
#> ASV1118                               0                              0
#> ASV113                                2                              0
#> ASV1132                               0                              0
#> ASV1144                               0                              0
#> ASV1146                               0                              0
#> ASV1152                               0                              0
#> ASV1155                               0                              0
#> ASV1164                               0                              0
#> ASV1166                               0                              0
#> ASV1167                               0                              0
#> ASV1179                               0                              0
#> ASV1192                               0                              5
#> ASV12                                 0                              0
#> ASV1200                               0                             92
#> ASV1212                               0                              9
#> ASV1218                               0                              0
#> ASV1223                               0                              0
#> ASV1230                               0                              0
#> ASV1233                               0                              0
#> ASV1234                               0                              0
#> ASV1236                               0                              0
#> ASV1239                               0                              0
#> ASV1242                               0                              0
#> ASV1253                               8                              0
#> ASV1270                               2                              0
#> ASV1278                               0                              2
#> ASV1279                               0                              0
#> ASV1300                               0                              0
#> ASV1302                               0                              0
#> ASV131                                0                              0
#> ASV1311                               0                              0
#> ASV1314                               0                              0
#> ASV1321                               0                              0
#> ASV1332                               0                              0
#> ASV1340                               0                              0
#> ASV1359                               0                             14
#> ASV1365                               0                              0
#> ASV1369                               0                              0
#> ASV137                                0                              0
#> ASV1370                               0                              0
#> ASV1379                               0                              0
#> ASV139                                0                             51
#> ASV1396                               0                              0
#> ASV1409                               0                             19
#> ASV1419                               0                              0
#> ASV1420                               0                              0
#> ASV1421                               0                              0
#> ASV1422                               0                              0
#> ASV1428                               0                              0
#> ASV1434                               0                              0
#> ASV144                                0                              0
#> ASV1458                               0                              0
#> ASV1459                               0                              0
#> ASV1460                               0                              0
#> ASV1461                               0                              0
#> ASV1462                               0                              0
#> ASV1468                               0                              0
#> ASV1483                               0                              0
#> ASV1485                               0                              0
#> ASV151                                1                              0
#> ASV1510                               0                              0
#> ASV1523                               0                              0
#> ASV1532                               0                              0
#> ASV1537                               0                              0
#> ASV1548                               0                              0
#> ASV1550                               0                              0
#> ASV1552                               0                              0
#> ASV1561                              58                              0
#> ASV1562                               1                              0
#> ASV1574                               0                              0
#> ASV1576                               0                              0
#> ASV1577                               0                              0
#> ASV1588                               0                              0
#> ASV159                                0                            114
#> ASV1609                               0                              0
#> ASV162                              985                              0
#> ASV1624                               0                              0
#> ASV1628                               0                              0
#> ASV1630                               0                              0
#> ASV1632                               0                              0
#> ASV1635                               0                              0
#> ASV1642                               0                              0
#> ASV1644                               0                              0
#> ASV1650                               0                              0
#> ASV1653                               0                              0
#> ASV1674                               0                              0
#> ASV1683                               0                              0
#> ASV1690                               0                              0
#> ASV170                                0                              0
#> ASV171                                0                              0
#> ASV1712                               0                              0
#> ASV172                                6                              0
#> ASV1726                               0                              0
#> ASV175                              165                              0
#> ASV18                                 0                              0
#> ASV19                                 0                             65
#> ASV193                                0                              0
#> ASV198                                0                              0
#> ASV199                                0                              1
#> ASV2                                 20                              5
#> ASV201                                1                              0
#> ASV202                                0                              0
#> ASV210                                0                              0
#> ASV219                                0                              0
#> ASV221                                0                              0
#> ASV223                                0                              0
#> ASV24                                 0                              0
#> ASV242                                0                              0
#> ASV244                                0                             31
#> ASV248                                0                              1
#> ASV25                                 0                              0
#> ASV251                                0                            379
#> ASV253                                0                              0
#> ASV255                                2                             12
#> ASV26                                 0                              0
#> ASV264                                0                              0
#> ASV266                                0                              0
#> ASV27                                 0                              0
#> ASV28                                 0                              1
#> ASV286                                1                              0
#> ASV29                                 2                              0
#> ASV292                                0                              0
#> ASV297                                0                              0
#> ASV310                                0                             47
#> ASV314                                0                              0
#> ASV320                                0                              0
#> ASV322                                0                              0
#> ASV329                                0                              0
#> ASV334                                0                              0
#> ASV338                                0                              0
#> ASV339                                1                              0
#> ASV344                                0                              0
#> ASV357                                1                             14
#> ASV359                                0                              0
#> ASV365                               14                              0
#> ASV368                                0                              0
#> ASV377                                0                              0
#> ASV378                                0                              1
#> ASV38                                 0                              0
#> ASV383                                0                              0
#> ASV404                                0                             42
#> ASV41                                 0                          10572
#> ASV42                                 0                              0
#> ASV428                                0                            146
#> ASV43                                 1                              1
#> ASV431                                0                              0
#> ASV443                                0                              0
#> ASV444                                0                              0
#> ASV45                                 0                              0
#> ASV454                                0                              0
#> ASV46                               589                              0
#> ASV464                                0                              0
#> ASV47                                 0                              0
#> ASV477                                0                              0
#> ASV48                                 0                              0
#> ASV483                                0                              0
#> ASV489                                0                            876
#> ASV49                                 0                            141
#> ASV500                                0                              0
#> ASV504                                0                              0
#> ASV505                                0                              0
#> ASV507                                0                            773
#> ASV509                                0                              0
#> ASV515                                0                              0
#> ASV517                                0                              0
#> ASV52                                 0                             87
#> ASV523                                0                              0
#> ASV524                                0                              0
#> ASV533                                0                              0
#> ASV543                                0                             78
#> ASV545                                0                              0
#> ASV546                                8                              0
#> ASV549                                0                              0
#> ASV55                                 0                              0
#> ASV552                                0                              0
#> ASV559                                0                              0
#> ASV563                                0                              0
#> ASV569                                0                              0
#> ASV570                                0                              0
#> ASV577                                0                             19
#> ASV586                                0                              0
#> ASV59                                 0                              0
#> ASV598                                0                              0
#> ASV60                                 0                              0
#> ASV602                                0                              0
#> ASV618                                0                              0
#> ASV619                                0                              0
#> ASV631                                0                              0
#> ASV635                                1                              0
#> ASV637                                0                              0
#> ASV64                                 0                              0
#> ASV643                                0                              0
#> ASV65                                18                              0
#> ASV662                                0                              0
#> ASV672                                0                              0
#> ASV673                                0                              0
#> ASV685                                6                              0
#> ASV69                                 0                              0
#> ASV692                                0                              0
#> ASV694                                5                              0
#> ASV715                                0                              0
#> ASV717                                3                              0
#> ASV727                                5                              0
#> ASV731                                0                              0
#> ASV737                                0                              0
#> ASV741                                0                              1
#> ASV743                                0                              0
#> ASV746                                0                              2
#> ASV749                                0                              0
#> ASV753                                0                              0
#> ASV756                                0                              0
#> ASV758                                0                              0
#> ASV772                                0                              0
#> ASV78                                 0                              0
#> ASV784                                0                              8
#> ASV787                                0                              0
#> ASV8                                  1                              7
#> ASV801                                0                              0
#> ASV81                                 2                              0
#> ASV817                                0                              0
#> ASV823                                0                              0
#> ASV829                                0                              0
#> ASV834                                1                              0
#> ASV839                                5                              0
#> ASV854                                0                             35
#> ASV87                                 0                              0
#> ASV89                                13                              0
#> ASV891                                0                              0
#> ASV892                                0                              0
#> ASV898                                0                              0
#> ASV899                                0                              0
#> ASV904                                0                              0
#> ASV906                                0                              0
#> ASV91                                 0                              0
#> ASV914                                0                              0
#> ASV926                                0                              0
#> ASV927                                0                              0
#> ASV930                                0                              0
#> ASV942                                0                              0
#> ASV946                                0                              0
#> ASV950                                0                              0
#> ASV953                                0                              0
#> ASV964                                0                              0
#> ASV969                                0                              0
#> ASV973                                0                              0
#> ASV987                                0                              0
#> ASV988                                0                              0
#>         ADABM30X.B_S14_MERGED.fastq.gz ADABM30X.H_S15_MERGED.fastq.gz
#> ASV1004                              0                              0
#> ASV1007                              0                              0
#> ASV101                               0                              0
#> ASV1017                              3                              0
#> ASV1022                              0                              0
#> ASV1030                              0                              0
#> ASV104                               3                              0
#> ASV105                               0                              0
#> ASV1059                              0                              0
#> ASV1069                              0                              0
#> ASV107                               0                              0
#> ASV1070                              0                              0
#> ASV1074                              0                              0
#> ASV1085                              0                              0
#> ASV1102                              0                              0
#> ASV1108                              3                              2
#> ASV111                              37                              0
#> ASV1118                              0                              0
#> ASV113                               0                             12
#> ASV1132                              0                              0
#> ASV1144                              0                              0
#> ASV1146                              1                              0
#> ASV1152                              0                              0
#> ASV1155                              0                              0
#> ASV1164                              0                              0
#> ASV1166                              0                              0
#> ASV1167                              0                              0
#> ASV1179                              0                              0
#> ASV1192                              0                              0
#> ASV12                                0                              0
#> ASV1200                              0                              0
#> ASV1212                              0                              0
#> ASV1218                              0                              0
#> ASV1223                              0                              0
#> ASV1230                              0                              8
#> ASV1233                              0                              0
#> ASV1234                              0                             26
#> ASV1236                              0                              0
#> ASV1239                              0                              0
#> ASV1242                              0                              0
#> ASV1253                              0                              0
#> ASV1270                              0                              0
#> ASV1278                              0                              0
#> ASV1279                              3                             14
#> ASV1300                              0                              0
#> ASV1302                              0                              0
#> ASV131                             235                           2751
#> ASV1311                              0                              0
#> ASV1314                              0                              0
#> ASV1321                              0                              0
#> ASV1332                              0                              1
#> ASV1340                              0                              0
#> ASV1359                              0                              0
#> ASV1365                              0                              0
#> ASV1369                              0                              0
#> ASV137                               0                              0
#> ASV1370                              0                              0
#> ASV1379                              0                              0
#> ASV139                               0                              0
#> ASV1396                              0                              0
#> ASV1409                              0                              0
#> ASV1419                              0                              0
#> ASV1420                              0                              0
#> ASV1421                              0                              0
#> ASV1422                              0                              0
#> ASV1428                              0                              0
#> ASV1434                              0                              0
#> ASV144                               0                           2926
#> ASV1458                              0                              0
#> ASV1459                              0                              0
#> ASV1460                              0                             15
#> ASV1461                              0                             72
#> ASV1462                              0                              0
#> ASV1468                              0                              0
#> ASV1483                              0                              1
#> ASV1485                              0                              0
#> ASV151                               0                              0
#> ASV1510                              0                              0
#> ASV1523                              0                              0
#> ASV1532                              0                              0
#> ASV1537                              0                              0
#> ASV1548                             15                              0
#> ASV1550                              0                              0
#> ASV1552                              0                             23
#> ASV1561                              0                              0
#> ASV1562                              0                              0
#> ASV1574                              0                              0
#> ASV1576                              0                              0
#> ASV1577                              0                              0
#> ASV1588                              0                              0
#> ASV159                               0                              0
#> ASV1609                              0                             20
#> ASV162                            2068                              0
#> ASV1624                              0                              0
#> ASV1628                              0                              0
#> ASV1630                              0                              0
#> ASV1632                              0                              0
#> ASV1635                              0                              0
#> ASV1642                              0                              0
#> ASV1644                              0                              0
#> ASV1650                              0                              6
#> ASV1653                              1                              0
#> ASV1674                              0                              0
#> ASV1683                              0                              0
#> ASV1690                              0                              2
#> ASV170                               0                              0
#> ASV171                               1                              0
#> ASV1712                              0                              1
#> ASV172                               0                              0
#> ASV1726                              0                              0
#> ASV175                             383                             15
#> ASV18                                0                              0
#> ASV19                                9                              8
#> ASV193                               0                              0
#> ASV198                               0                              0
#> ASV199                               0                              0
#> ASV2                                 3                              3
#> ASV201                               0                              0
#> ASV202                               0                              0
#> ASV210                               0                              0
#> ASV219                               1                              1
#> ASV221                               0                              0
#> ASV223                               0                            298
#> ASV24                                8                              0
#> ASV242                               0                              0
#> ASV244                               0                              0
#> ASV248                               0                              0
#> ASV25                                0                              2
#> ASV251                               0                              0
#> ASV253                            1594                              0
#> ASV255                               0                              0
#> ASV26                                0                              0
#> ASV264                               0                           1013
#> ASV266                               3                              0
#> ASV27                                0                              0
#> ASV28                                1                              0
#> ASV286                               0                              0
#> ASV29                                0                              1
#> ASV292                               0                             11
#> ASV297                               0                              0
#> ASV310                               0                              9
#> ASV314                               0                              0
#> ASV320                               0                              0
#> ASV322                              41                              0
#> ASV329                               0                              0
#> ASV334                               0                              0
#> ASV338                               0                              0
#> ASV339                               0                              0
#> ASV344                               0                              0
#> ASV357                               0                              0
#> ASV359                               0                              0
#> ASV365                               7                              4
#> ASV368                               0                              8
#> ASV377                               0                              0
#> ASV378                               1                              0
#> ASV38                               23                            282
#> ASV383                               0                              0
#> ASV404                               0                              0
#> ASV41                                0                              0
#> ASV42                                0                              0
#> ASV428                               0                              0
#> ASV43                                0                              0
#> ASV431                               0                              0
#> ASV443                               0                             27
#> ASV444                               2                              0
#> ASV45                                0                              0
#> ASV454                               0                              0
#> ASV46                             1238                              0
#> ASV464                               0                              0
#> ASV47                                0                              0
#> ASV477                               0                              0
#> ASV48                                0                             22
#> ASV483                              44                              0
#> ASV489                               0                              0
#> ASV49                                0                              0
#> ASV500                               0                              0
#> ASV504                               0                              0
#> ASV505                               0                              0
#> ASV507                               0                              0
#> ASV509                               0                              0
#> ASV515                               0                              0
#> ASV517                               0                              0
#> ASV52                                0                              0
#> ASV523                               0                              0
#> ASV524                               0                              0
#> ASV533                               0                              0
#> ASV543                               0                              0
#> ASV545                               0                              0
#> ASV546                               0                              3
#> ASV549                               0                              0
#> ASV55                                0                              0
#> ASV552                               0                              0
#> ASV559                               0                              0
#> ASV563                               0                              0
#> ASV569                               0                              0
#> ASV570                               0                              0
#> ASV577                               0                              0
#> ASV586                               0                              0
#> ASV59                                0                              0
#> ASV598                               0                              0
#> ASV60                                0                              0
#> ASV602                               0                              1
#> ASV618                               0                              0
#> ASV619                               0                              0
#> ASV631                               0                              0
#> ASV635                               3                              1
#> ASV637                               0                              0
#> ASV64                                0                              0
#> ASV643                               0                              3
#> ASV65                                7                             15
#> ASV662                               0                              0
#> ASV672                               0                              0
#> ASV673                               0                              0
#> ASV685                               0                              0
#> ASV69                                0                             14
#> ASV692                               0                              0
#> ASV694                              56                              5
#> ASV715                               0                              0
#> ASV717                               0                              0
#> ASV727                               1                              7
#> ASV731                               0                              0
#> ASV737                               0                              0
#> ASV741                               0                              9
#> ASV743                               0                              0
#> ASV746                               9                            157
#> ASV749                               0                              0
#> ASV753                               0                              0
#> ASV756                               0                              1
#> ASV758                               0                              0
#> ASV772                               0                              3
#> ASV78                                0                              0
#> ASV784                               0                              0
#> ASV787                               0                              0
#> ASV8                                 2                              2
#> ASV801                               1                              0
#> ASV81                                0                              0
#> ASV817                               0                              0
#> ASV823                               0                              0
#> ASV829                               0                              0
#> ASV834                               0                              0
#> ASV839                               0                              0
#> ASV854                               0                              0
#> ASV87                                0                              0
#> ASV89                                0                              0
#> ASV891                               0                              0
#> ASV892                               0                              0
#> ASV898                               0                              0
#> ASV899                               5                              0
#> ASV904                               0                              0
#> ASV906                               0                              0
#> ASV91                                1                              0
#> ASV914                               0                              0
#> ASV926                               0                              0
#> ASV927                               0                              0
#> ASV930                               0                              0
#> ASV942                               0                              0
#> ASV946                               0                              0
#> ASV950                               0                              0
#> ASV953                               0                              0
#> ASV964                               0                              0
#> ASV969                               0                              0
#> ASV973                               1                              0
#> ASV987                               0                              0
#> ASV988                               0                              0
#>         ADABM30X.M_S16_MERGED.fastq.gz AE30.ABM507_S17_MERGED.fastq.gz
#> ASV1004                              0                               0
#> ASV1007                              0                              19
#> ASV101                               0                               0
#> ASV1017                              0                               0
#> ASV1022                              0                               2
#> ASV1030                              0                               0
#> ASV104                               0                               0
#> ASV105                               0                               0
#> ASV1059                              0                              14
#> ASV1069                              0                               0
#> ASV107                               0                               0
#> ASV1070                              0                               0
#> ASV1074                              0                              12
#> ASV1085                              0                               0
#> ASV1102                              0                             144
#> ASV1108                              0                               0
#> ASV111                               0                               0
#> ASV1118                              0                               0
#> ASV113                               0                               0
#> ASV1132                              0                               0
#> ASV1144                              0                               0
#> ASV1146                              0                               0
#> ASV1152                              0                               0
#> ASV1155                              0                               0
#> ASV1164                              0                               0
#> ASV1166                              0                               0
#> ASV1167                              0                               0
#> ASV1179                              0                               0
#> ASV1192                              0                               6
#> ASV12                                0                               0
#> ASV1200                              0                               0
#> ASV1212                              0                               0
#> ASV1218                              0                               0
#> ASV1223                              0                               0
#> ASV1230                              0                               0
#> ASV1233                              0                               0
#> ASV1234                              0                               0
#> ASV1236                              0                               0
#> ASV1239                              0                               0
#> ASV1242                              0                               0
#> ASV1253                              0                               9
#> ASV1270                              0                               0
#> ASV1278                              0                               0
#> ASV1279                              0                               0
#> ASV1300                              0                               0
#> ASV1302                              0                               0
#> ASV131                               0                               0
#> ASV1311                              0                               0
#> ASV1314                              0                               0
#> ASV1321                              0                               0
#> ASV1332                              0                               0
#> ASV1340                              0                               0
#> ASV1359                              0                               0
#> ASV1365                              0                               0
#> ASV1369                              0                              18
#> ASV137                               0                               0
#> ASV1370                              0                              62
#> ASV1379                              0                               0
#> ASV139                               0                            2375
#> ASV1396                              0                               0
#> ASV1409                              0                               2
#> ASV1419                              0                               0
#> ASV1420                              0                               0
#> ASV1421                              0                               0
#> ASV1422                              0                               0
#> ASV1428                              0                               0
#> ASV1434                              0                               0
#> ASV144                               0                               0
#> ASV1458                              0                               0
#> ASV1459                              0                               3
#> ASV1460                              0                               0
#> ASV1461                              0                               0
#> ASV1462                              0                              13
#> ASV1468                              5                               0
#> ASV1483                              0                               0
#> ASV1485                              0                               0
#> ASV151                               0                               0
#> ASV1510                              0                               0
#> ASV1523                              0                               0
#> ASV1532                              0                               0
#> ASV1537                              0                               0
#> ASV1548                              0                               0
#> ASV1550                              0                               0
#> ASV1552                              0                               0
#> ASV1561                              0                               0
#> ASV1562                              0                               0
#> ASV1574                              0                               0
#> ASV1576                              0                               0
#> ASV1577                              0                               0
#> ASV1588                              0                               0
#> ASV159                               0                               0
#> ASV1609                              0                               0
#> ASV162                               0                               0
#> ASV1624                              0                               0
#> ASV1628                              0                               0
#> ASV1630                              0                               0
#> ASV1632                              0                               0
#> ASV1635                              0                               0
#> ASV1642                              0                               0
#> ASV1644                              0                               0
#> ASV1650                              0                               0
#> ASV1653                              0                               0
#> ASV1674                              0                               0
#> ASV1683                              0                               0
#> ASV1690                              0                               0
#> ASV170                               0                               0
#> ASV171                               0                               0
#> ASV1712                              0                               0
#> ASV172                               0                             265
#> ASV1726                              0                               0
#> ASV175                               2                               0
#> ASV18                                0                              44
#> ASV19                                0                               1
#> ASV193                               0                               0
#> ASV198                               0                               0
#> ASV199                               0                               0
#> ASV2                                 3                               7
#> ASV201                               0                               6
#> ASV202                               0                               0
#> ASV210                               0                               0
#> ASV219                               0                               0
#> ASV221                               0                               0
#> ASV223                               2                               0
#> ASV24                                1                               0
#> ASV242                               0                               0
#> ASV244                               0                             577
#> ASV248                               0                               0
#> ASV25                                0                               0
#> ASV251                               0                               0
#> ASV253                               0                               0
#> ASV255                               3                               0
#> ASV26                                0                               0
#> ASV264                               0                               0
#> ASV266                               1                               0
#> ASV27                                0                               0
#> ASV28                                0                               1
#> ASV286                               0                               0
#> ASV29                                0                               0
#> ASV292                               0                               0
#> ASV297                               0                               0
#> ASV310                               0                               0
#> ASV314                               0                               0
#> ASV320                               0                               0
#> ASV322                               0                               0
#> ASV329                               0                               0
#> ASV334                               0                               0
#> ASV338                             213                               0
#> ASV339                               0                               5
#> ASV344                               0                               0
#> ASV357                               0                               0
#> ASV359                               0                               0
#> ASV365                               0                               0
#> ASV368                               0                               0
#> ASV377                               0                               0
#> ASV378                               0                              34
#> ASV38                               39                               0
#> ASV383                               0                               0
#> ASV404                               0                               0
#> ASV41                                0                               0
#> ASV42                                0                               0
#> ASV428                               0                               0
#> ASV43                                0                              36
#> ASV431                               0                               0
#> ASV443                               0                               0
#> ASV444                               0                               0
#> ASV45                                0                               0
#> ASV454                               0                               0
#> ASV46                                0                               0
#> ASV464                               0                               0
#> ASV47                                0                             315
#> ASV477                               0                               0
#> ASV48                                0                               0
#> ASV483                             703                               0
#> ASV489                               0                               0
#> ASV49                                0                               0
#> ASV500                               0                               0
#> ASV504                               0                               0
#> ASV505                               0                               0
#> ASV507                               0                             167
#> ASV509                               0                               0
#> ASV515                               0                               0
#> ASV517                               0                               3
#> ASV52                                0                               0
#> ASV523                               0                               0
#> ASV524                               0                               0
#> ASV533                               0                               0
#> ASV543                               3                               0
#> ASV545                               0                               0
#> ASV546                               0                               0
#> ASV549                               0                               0
#> ASV55                                0                               0
#> ASV552                               0                               0
#> ASV559                               0                              65
#> ASV563                               0                               0
#> ASV569                               0                             357
#> ASV570                               0                               0
#> ASV577                               0                               0
#> ASV586                               0                               0
#> ASV59                                0                            1207
#> ASV598                               0                               0
#> ASV60                                0                               0
#> ASV602                               0                             327
#> ASV618                               0                               0
#> ASV619                               0                             419
#> ASV631                               0                               0
#> ASV635                               0                               0
#> ASV637                               0                               0
#> ASV64                                0                               0
#> ASV643                               0                               0
#> ASV65                                0                               0
#> ASV662                               0                               0
#> ASV672                               0                               0
#> ASV673                               0                               0
#> ASV685                               1                               0
#> ASV69                                0                               0
#> ASV692                               0                               0
#> ASV694                               2                               0
#> ASV715                               0                               3
#> ASV717                               0                               3
#> ASV727                               0                               0
#> ASV731                               1                               0
#> ASV737                               0                              75
#> ASV741                               0                             210
#> ASV743                               0                               0
#> ASV746                              49                               0
#> ASV749                               0                               0
#> ASV753                               0                               0
#> ASV756                               0                               0
#> ASV758                               0                               0
#> ASV772                               0                               0
#> ASV78                                0                               0
#> ASV784                               0                              12
#> ASV787                               0                               0
#> ASV8                                 0                              37
#> ASV801                               0                               0
#> ASV81                                0                               5
#> ASV817                               0                               0
#> ASV823                               2                               0
#> ASV829                               0                               0
#> ASV834                               0                               1
#> ASV839                               0                               0
#> ASV854                               0                               0
#> ASV87                                0                               0
#> ASV89                                0                               8
#> ASV891                               2                               0
#> ASV892                               0                               0
#> ASV898                               0                               0
#> ASV899                               2                               0
#> ASV904                               0                               0
#> ASV906                               0                               2
#> ASV91                                0                               0
#> ASV914                               0                               0
#> ASV926                               0                               0
#> ASV927                               0                               8
#> ASV930                               0                               0
#> ASV942                               0                               0
#> ASV946                               0                               0
#> ASV950                               0                             174
#> ASV953                               0                               0
#> ASV964                             133                               0
#> ASV969                               0                               0
#> ASV973                               5                               0
#> ASV987                               0                               0
#> ASV988                               0                               0
#> 
#> $res_lulu$curated_count
#> [1] 264
#> 
#> $res_lulu$curated_otus
#>   [1] "ASV1004" "ASV1007" "ASV101"  "ASV1017" "ASV1022" "ASV1030" "ASV104" 
#>   [8] "ASV105"  "ASV1059" "ASV1069" "ASV107"  "ASV1070" "ASV1074" "ASV1085"
#>  [15] "ASV1102" "ASV1108" "ASV111"  "ASV1118" "ASV113"  "ASV1132" "ASV1144"
#>  [22] "ASV1146" "ASV1152" "ASV1155" "ASV1164" "ASV1166" "ASV1167" "ASV1179"
#>  [29] "ASV1192" "ASV12"   "ASV1200" "ASV1212" "ASV1218" "ASV1223" "ASV1230"
#>  [36] "ASV1233" "ASV1234" "ASV1236" "ASV1239" "ASV1242" "ASV1253" "ASV1270"
#>  [43] "ASV1278" "ASV1279" "ASV1300" "ASV1302" "ASV131"  "ASV1311" "ASV1314"
#>  [50] "ASV1321" "ASV1332" "ASV1340" "ASV1359" "ASV1365" "ASV1369" "ASV137" 
#>  [57] "ASV1370" "ASV1379" "ASV139"  "ASV1396" "ASV1409" "ASV1419" "ASV1420"
#>  [64] "ASV1421" "ASV1422" "ASV1428" "ASV1434" "ASV144"  "ASV1458" "ASV1459"
#>  [71] "ASV1460" "ASV1461" "ASV1462" "ASV1468" "ASV1483" "ASV1485" "ASV151" 
#>  [78] "ASV1510" "ASV1523" "ASV1532" "ASV1537" "ASV1548" "ASV1550" "ASV1552"
#>  [85] "ASV1561" "ASV1562" "ASV1574" "ASV1576" "ASV1577" "ASV1588" "ASV159" 
#>  [92] "ASV1609" "ASV162"  "ASV1624" "ASV1628" "ASV1630" "ASV1632" "ASV1635"
#>  [99] "ASV1642" "ASV1644" "ASV1650" "ASV1653" "ASV1674" "ASV1683" "ASV1690"
#> [106] "ASV170"  "ASV171"  "ASV1712" "ASV172"  "ASV1726" "ASV175"  "ASV18"  
#> [113] "ASV19"   "ASV193"  "ASV198"  "ASV199"  "ASV2"    "ASV201"  "ASV202" 
#> [120] "ASV210"  "ASV219"  "ASV221"  "ASV223"  "ASV24"   "ASV242"  "ASV244" 
#> [127] "ASV248"  "ASV25"   "ASV251"  "ASV253"  "ASV255"  "ASV26"   "ASV264" 
#> [134] "ASV266"  "ASV27"   "ASV28"   "ASV286"  "ASV29"   "ASV292"  "ASV297" 
#> [141] "ASV310"  "ASV314"  "ASV320"  "ASV322"  "ASV329"  "ASV334"  "ASV338" 
#> [148] "ASV339"  "ASV344"  "ASV357"  "ASV359"  "ASV365"  "ASV368"  "ASV377" 
#> [155] "ASV378"  "ASV38"   "ASV383"  "ASV404"  "ASV41"   "ASV42"   "ASV428" 
#> [162] "ASV43"   "ASV431"  "ASV443"  "ASV444"  "ASV45"   "ASV454"  "ASV46"  
#> [169] "ASV464"  "ASV47"   "ASV477"  "ASV48"   "ASV483"  "ASV489"  "ASV49"  
#> [176] "ASV500"  "ASV504"  "ASV505"  "ASV507"  "ASV509"  "ASV515"  "ASV517" 
#> [183] "ASV52"   "ASV523"  "ASV524"  "ASV533"  "ASV543"  "ASV545"  "ASV546" 
#> [190] "ASV549"  "ASV55"   "ASV552"  "ASV559"  "ASV563"  "ASV569"  "ASV570" 
#> [197] "ASV577"  "ASV586"  "ASV59"   "ASV598"  "ASV60"   "ASV602"  "ASV618" 
#> [204] "ASV619"  "ASV631"  "ASV635"  "ASV637"  "ASV64"   "ASV643"  "ASV65"  
#> [211] "ASV662"  "ASV672"  "ASV673"  "ASV685"  "ASV69"   "ASV692"  "ASV694" 
#> [218] "ASV715"  "ASV717"  "ASV727"  "ASV731"  "ASV737"  "ASV741"  "ASV743" 
#> [225] "ASV746"  "ASV749"  "ASV753"  "ASV756"  "ASV758"  "ASV772"  "ASV78"  
#> [232] "ASV784"  "ASV787"  "ASV8"    "ASV801"  "ASV81"   "ASV817"  "ASV823" 
#> [239] "ASV829"  "ASV834"  "ASV839"  "ASV854"  "ASV87"   "ASV89"   "ASV891" 
#> [246] "ASV892"  "ASV898"  "ASV899"  "ASV904"  "ASV906"  "ASV91"   "ASV914" 
#> [253] "ASV926"  "ASV927"  "ASV930"  "ASV942"  "ASV946"  "ASV950"  "ASV953" 
#> [260] "ASV964"  "ASV969"  "ASV973"  "ASV987"  "ASV988" 
#> 
#> $res_lulu$discarded_count
#> [1] 147
#> 
#> $res_lulu$discarded_otus
#>   [1] "ASV666"  "ASV989"  "ASV94"   "ASV178"  "ASV209"  "ASV1058" "ASV261" 
#>   [8] "ASV816"  "ASV542"  "ASV85"   "ASV870"  "ASV333"  "ASV348"  "ASV98"  
#>  [15] "ASV493"  "ASV724"  "ASV566"  "ASV831"  "ASV1078" "ASV859"  "ASV1356"
#>  [22] "ASV313"  "ASV462"  "ASV594"  "ASV534"  "ASV915"  "ASV346"  "ASV1268"
#>  [29] "ASV1387" "ASV592"  "ASV766"  "ASV975"  "ASV1493" "ASV1131" "ASV1101"
#>  [36] "ASV474"  "ASV744"  "ASV1261" "ASV963"  "ASV880"  "ASV1052" "ASV1128"
#>  [43] "ASV358"  "ASV420"  "ASV693"  "ASV884"  "ASV1319" "ASV1557" "ASV1315"
#>  [50] "ASV576"  "ASV641"  "ASV1198" "ASV466"  "ASV1036" "ASV1247" "ASV1427"
#>  [57] "ASV1176" "ASV168"  "ASV762"  "ASV461"  "ASV833"  "ASV998"  "ASV254" 
#>  [64] "ASV1092" "ASV1117" "ASV832"  "ASV1182" "ASV624"  "ASV1320" "ASV999" 
#>  [71] "ASV1313" "ASV847"  "ASV722"  "ASV1526" "ASV367"  "ASV1386" "ASV1035"
#>  [78] "ASV1037" "ASV711"  "ASV33"   "ASV1410" "ASV853"  "ASV1067" "ASV440" 
#>  [85] "ASV580"  "ASV858"  "ASV526"  "ASV355"  "ASV1010" "ASV1542" "ASV1262"
#>  [92] "ASV962"  "ASV903"  "ASV984"  "ASV1011" "ASV1224" "ASV1276" "ASV996" 
#>  [99] "ASV1082" "ASV1109" "ASV1467" "ASV203"  "ASV491"  "ASV626"  "ASV650" 
#> [106] "ASV790"  "ASV1014" "ASV1625" "ASV249"  "ASV272"  "ASV488"  "ASV604" 
#> [113] "ASV616"  "ASV736"  "ASV748"  "ASV760"  "ASV828"  "ASV961"  "ASV979" 
#> [120] "ASV986"  "ASV1133" "ASV1603" "ASV67"   "ASV116"  "ASV405"  "ASV415" 
#> [127] "ASV796"  "ASV840"  "ASV860"  "ASV875"  "ASV940"  "ASV941"  "ASV959" 
#> [134] "ASV981"  "ASV1015" "ASV1020" "ASV1027" "ASV1031" "ASV1066" "ASV1107"
#> [141] "ASV1111" "ASV1263" "ASV1303" "ASV1341" "ASV1380" "ASV1388" "ASV1484"
#> 
#> $res_lulu$runtime
#> Time difference of 2.096411 secs
#> 
#> $res_lulu$minimum_match
#> [1] 84
#> 
#> $res_lulu$minimum_relative_cooccurence
#> [1] 0.95
#> 
#> $res_lulu$otu_map
#>         total spread parent_id curated rank
#> ASV2      662     19      ASV2  parent   34
#> ASV8     3318     15      ASV8  parent    9
#> ASV38    1124      7     ASV38  parent   27
#> ASV756     70      7    ASV756  parent  132
#> ASV12    4528      6     ASV12  parent    7
#> ASV19    3035      6     ASV19  parent   12
#> ASV175   1933      6    ASV175  parent   18
#> ASV18    1188      6     ASV18  parent   26
#> ASV694    185      6    ASV694  parent   77
#> ASV113    118      6    ASV113  parent  105
#> ASV378     46      6    ASV378  parent  156
#> ASV717     34      6    ASV717  parent  178
#> ASV41   10606      5     ASV41  parent    2
#> ASV65    8546      5     ASV65  parent    3
#> ASV89    1429      5     ASV89  parent   24
#> ASV618    369      5    ASV618  parent   52
#> ASV170    350      5    ASV170  parent   56
#> ASV746    224      5    ASV746  parent   70
#> ASV953    171      5    ASV953  parent   82
#> ASV28      82      5     ASV28  parent  122
#> ASV310     82      5    ASV310  parent  123
#> ASV666     47      5      ASV8  merged  154
#> ASV1022    43      5   ASV1022  parent  163
#> ASV1192    20      5   ASV1192  parent  195
#> ASV727     18      5    ASV727  parent  204
#> ASV989     16      5     ASV38  merged  208
#> ASV635     15      5    ASV635  parent  210
#> ASV1074    12      5   ASV1074  parent  230
#> ASV1108     8      5   ASV1108  parent  256
#> ASV1483     6      5   ASV1483  parent  266
#> ASV43   11075      4     ASV43  parent    1
#> ASV46    3058      4     ASV46  parent   10
#> ASV365    730      4    ASV365  parent   32
#> ASV94     678      4      ASV8  merged   33
#> ASV602    366      4    ASV602  parent   53
#> ASV178    308      4     ASV19  merged   61
#> ASV172    284      4    ASV172  parent   65
#> ASV209    264      4      ASV8  merged   67
#> ASV741    221      4    ASV741  parent   72
#> ASV1058   139      4     ASV89  merged   97
#> ASV261    102      4      ASV8  merged  109
#> ASV329     77      4    ASV329  parent  126
#> ASV906     33      4    ASV906  parent  180
#> ASV546     29      4    ASV546  parent  184
#> ASV1253    22      4   ASV1253  parent  192
#> ASV24      19      4     ASV24  parent  200
#> ASV1421    19      4   ASV1421  parent  201
#> ASV816     16      4     ASV38  merged  209
#> ASV662     13      4    ASV662  parent  220
#> ASV266      6      4    ASV266  parent  267
#> ASV162   3055      3    ASV162  parent   11
#> ASV139   2056      3    ASV139  parent   16
#> ASV219   1963      3    ASV219  parent   17
#> ASV201   1791      3    ASV201  parent   20
#> ASV251   1004      3    ASV251  parent   29
#> ASV244    621      3    ASV244  parent   37
#> ASV507    520      3    ASV507  parent   41
#> ASV483    510      3    ASV483  parent   42
#> ASV500    477      3    ASV500  parent   44
#> ASV577    459      3    ASV577  parent   45
#> ASV344    451      3    ASV344  parent   46
#> ASV542    428      3    ASV507  merged   48
#> ASV85     371      3    ASV483  merged   51
#> ASV870    224      3    ASV577  merged   71
#> ASV899    217      3    ASV899  parent   73
#> ASV52     165      3     ASV52  parent   86
#> ASV926    165      3    ASV926  parent   87
#> ASV159    160      3    ASV159  parent   89
#> ASV333    149      3     ASV19  merged   91
#> ASV631    112      3    ASV631  parent  107
#> ASV348    108      3      ASV8  merged  108
#> ASV817     93      3    ASV817  parent  111
#> ASV1007    91      3   ASV1007  parent  114
#> ASV543     86      3    ASV543  parent  118
#> ASV60      83      3     ASV60  parent  121
#> ASV559     69      3    ASV559  parent  135
#> ASV1370    68      3   ASV1370  parent  136
#> ASV1459    67      3   ASV1459  parent  138
#> ASV1359    56      3   ASV1359  parent  148
#> ASV292     55      3    ASV292  parent  149
#> ASV1059    54      3   ASV1059  parent  150
#> ASV443     49      3    ASV443  parent  152
#> ASV505     44      3    ASV505  parent  160
#> ASV509     40      3    ASV509  parent  164
#> ASV1624    40      3   ASV1624  parent  165
#> ASV98      39      3    ASV443  merged  166
#> ASV368     38      3    ASV368  parent  170
#> ASV493     36      3      ASV8  merged  173
#> ASV1279    36      3   ASV1279  parent  174
#> ASV904     35      3    ASV904  parent  177
#> ASV563     32      3    ASV563  parent  182
#> ASV1644    29      3   ASV1644  parent  185
#> ASV1369    25      3   ASV1369  parent  189
#> ASV1409    20      3   ASV1409  parent  196
#> ASV724     18      3    ASV509  merged  205
#> ASV255     17      3    ASV255  parent  207
#> ASV69      15      3     ASV69  parent  211
#> ASV151     15      3    ASV151  parent  212
#> ASV566     15      3    ASV443  merged  213
#> ASV831     12      3      ASV8  merged  231
#> ASV1078    12      3     ASV19  merged  232
#> ASV1230    12      3   ASV1230  parent  233
#> ASV1004     9      3   ASV1004  parent  250
#> ASV586      8      3    ASV586  parent  257
#> ASV859      5      3    ASV509  merged  286
#> ASV1278     5      3   ASV1278  parent  287
#> ASV1356     5      3   ASV1279  merged  288
#> ASV801      4      3    ASV801  parent  298
#> ASV64    8303      2     ASV64  parent    4
#> ASV144   3596      2    ASV144  parent    8
#> ASV131   2986      2    ASV131  parent   13
#> ASV202   2288      2    ASV202  parent   15
#> ASV107   1682      2    ASV107  parent   22
#> ASV322    998      2    ASV322  parent   30
#> ASV444    631      2    ASV444  parent   36
#> ASV313    574      2     ASV89  merged   39
#> ASV569    361      2    ASV569  parent   54
#> ASV672    342      2    ASV672  parent   58
#> ASV462    334      2     ASV89  merged   59
#> ASV223    300      2    ASV223  parent   62
#> ASV594    293      2     ASV19  merged   63
#> ASV534    277      2     ASV89  merged   66
#> ASV950    175      2    ASV950  parent   79
#> ASV428    172      2    ASV428  parent   81
#> ASV915    167      2     ASV43  merged   84
#> ASV346    165      2    ASV618  merged   88
#> ASV1146   133      2   ASV1146  parent   99
#> ASV942    130      2    ASV942  parent  101
#> ASV1164   127      2   ASV1164  parent  104
#> ASV1268    90      2    ASV139  merged  117
#> ASV1152    84      2   ASV1152  parent  119
#> ASV1387    82      2     ASV43  merged  124
#> ASV592     78      2    ASV618  merged  125
#> ASV737     77      2    ASV737  parent  127
#> ASV1458    76      2   ASV1458  parent  128
#> ASV772     74      2    ASV772  parent  129
#> ASV523     60      2    ASV523  parent  140
#> ASV766     60      2    ASV139  merged  141
#> ASV55      58      2     ASV55  parent  143
#> ASV1236    58      2   ASV1236  parent  144
#> ASV749     57      2    ASV749  parent  146
#> ASV854     48      2    ASV854  parent  153
#> ASV743     47      2    ASV743  parent  155
#> ASV975     46      2     ASV38  merged  157
#> ASV111     44      2    ASV111  parent  161
#> ASV1532    44      2   ASV1532  parent  162
#> ASV1493    38      2    ASV344  merged  171
#> ASV1131    33      2     ASV19  merged  181
#> ASV1101    32      2     ASV89  merged  183
#> ASV474     20      2      ASV8  merged  197
#> ASV744     20      2    ASV139  merged  198
#> ASV1396    19      2   ASV1396  parent  202
#> ASV758     15      2    ASV758  parent  214
#> ASV784     15      2    ASV784  parent  215
#> ASV297     14      2    ASV297  parent  218
#> ASV91      13      2     ASV91  parent  221
#> ASV199     13      2    ASV199  parent  222
#> ASV643     13      2    ASV643  parent  223
#> ASV787     13      2    ASV787  parent  224
#> ASV1261    13      2    ASV618  merged  225
#> ASV1434    13      2   ASV1434  parent  226
#> ASV357     12      2    ASV357  parent  234
#> ASV927     12      2    ASV927  parent  235
#> ASV248     10      2    ASV248  parent  242
#> ASV963     10      2    ASV172  merged  243
#> ASV1212    10      2   ASV1212  parent  244
#> ASV1321    10      2   ASV1321  parent  245
#> ASV880      9      2     ASV52  merged  251
#> ASV1052     9      2    ASV509  merged  252
#> ASV1155     8      2   ASV1155  parent  258
#> ASV685      7      2    ASV685  parent  263
#> ASV1128     7      2    ASV509  merged  264
#> ASV81       6      2     ASV81  parent  268
#> ASV339      6      2    ASV339  parent  269
#> ASV358      6      2     ASV69  merged  270
#> ASV420      6      2    ASV509  merged  271
#> ASV693      6      2   ASV1409  merged  272
#> ASV715      6      2    ASV715  parent  273
#> ASV884      6      2    ASV509  merged  274
#> ASV973      6      2    ASV973  parent  275
#> ASV1218     6      2   ASV1218  parent  276
#> ASV1319     6      2   ASV1074  merged  277
#> ASV1557     6      2     ASV89  merged  278
#> ASV1315     5      2    ASV784  merged  289
#> ASV576      4      2      ASV2  merged  299
#> ASV641      4      2    ASV505  merged  300
#> ASV1198     4      2    ASV509  merged  301
#> ASV29       3      2     ASV29  parent  312
#> ASV286      3      2    ASV286  parent  313
#> ASV464      3      2    ASV464  parent  314
#> ASV466      3      2     ASV38  merged  315
#> ASV1036     3      2    ASV113  merged  316
#> ASV1247     3      2      ASV2  merged  317
#> ASV1300     3      2   ASV1300  parent  318
#> ASV1427     3      2    ASV509  merged  319
#> ASV1690     3      2   ASV1690  parent  320
#> ASV834      2      2    ASV834  parent  337
#> ASV1176     2      2    ASV292  merged  338
#> ASV1223     2      2   ASV1223  parent  339
#> ASV1332     2      2   ASV1332  parent  340
#> ASV1365     2      2   ASV1365  parent  341
#> ASV78    6005      1     ASV78  parent    5
#> ASV101   4669      1    ASV101  parent    6
#> ASV168   2877      1     ASV64  merged   14
#> ASV221   1927      1    ASV221  parent   19
#> ASV242   1692      1    ASV242  parent   21
#> ASV253   1594      1    ASV253  parent   23
#> ASV59    1207      1     ASV59  parent   25
#> ASV264   1013      1    ASV264  parent   28
#> ASV45     836      1     ASV45  parent   31
#> ASV454    646      1    ASV454  parent   35
#> ASV489    583      1    ASV489  parent   38
#> ASV524    531      1    ASV524  parent   40
#> ASV552    489      1    ASV552  parent   43
#> ASV533    434      1    ASV533  parent   47
#> ASV619    419      1    ASV619  parent   49
#> ASV504    404      1    ASV504  parent   50
#> ASV26     351      1     ASV26  parent   55
#> ASV692    346      1    ASV692  parent   57
#> ASV47     315      1     ASV47  parent   60
#> ASV762    293      1    ASV489  merged   64
#> ASV461    258      1    ASV139  merged   68
#> ASV833    248      1    ASV221  merged   69
#> ASV338    213      1    ASV338  parent   74
#> ASV898    201      1    ASV898  parent   75
#> ASV946    197      1    ASV946  parent   76
#> ASV987    176      1    ASV987  parent   78
#> ASV998    173      1    ASV444  merged   80
#> ASV254    168      1     ASV38  merged   83
#> ASV637    166      1    ASV637  parent   85
#> ASV1070   153      1   ASV1070  parent   90
#> ASV1092   146      1    ASV251  merged   92
#> ASV1102   144      1   ASV1102  parent   93
#> ASV49     141      1     ASV49  parent   94
#> ASV1117   140      1     ASV43  merged   95
#> ASV1118   140      1   ASV1118  parent   96
#> ASV1132   136      1   ASV1132  parent   98
#> ASV964    133      1    ASV964  parent  100
#> ASV832    130      1    ASV504  merged  102
#> ASV1166   129      1   ASV1166  parent  103
#> ASV1182   114      1      ASV8  merged  106
#> ASV1314    97      1   ASV1314  parent  110
#> ASV969     92      1    ASV969  parent  112
#> ASV1200    92      1   ASV1200  parent  113
#> ASV624     91      1     ASV12  merged  115
#> ASV1320    91      1    ASV898  merged  116
#> ASV999     84      1    ASV139  merged  120
#> ASV1485    73      1   ASV1485  parent  130
#> ASV1461    72      1   ASV1461  parent  131
#> ASV1313    70      1     ASV78  merged  133
#> ASV1420    70      1   ASV1420  parent  134
#> ASV105     68      1    ASV105  parent  137
#> ASV431     65      1    ASV431  parent  139
#> ASV377     60      1    ASV377  parent  142
#> ASV1561    58      1   ASV1561  parent  145
#> ASV477     57      1    ASV477  parent  147
#> ASV1537    53      1   ASV1537  parent  151
#> ASV847     46      1    ASV107  merged  158
#> ASV722     45      1    ASV344  merged  159
#> ASV404     39      1    ASV404  parent  167
#> ASV1510    39      1   ASV1510  parent  168
#> ASV1526    39      1     ASV55  merged  169
#> ASV1379    38      1   ASV1379  parent  172
#> ASV367     36      1     ASV28  merged  175
#> ASV1642    36      1   ASV1642  parent  176
#> ASV1233    34      1   ASV1233  parent  179
#> ASV1386    29      1    ASV969  merged  186
#> ASV1234    26      1   ASV1234  parent  187
#> ASV1340    26      1   ASV1340  parent  188
#> ASV1035    25      1   ASV1007  merged  190
#> ASV1037    23      1    ASV692  merged  191
#> ASV48      22      1     ASV48  parent  193
#> ASV711     21      1     ASV45  merged  194
#> ASV1609    20      1   ASV1609  parent  199
#> ASV33      19      1     ASV78  merged  203
#> ASV1683    18      1   ASV1683  parent  206
#> ASV673     15      1    ASV673  parent  216
#> ASV1548    15      1   ASV1548  parent  217
#> ASV1460    14      1   ASV1460  parent  219
#> ASV545     13      1    ASV545  parent  227
#> ASV1410    13      1    ASV505  merged  228
#> ASV1462    13      1   ASV1462  parent  229
#> ASV853     12      1     ASV38  merged  236
#> ASV1067    12      1     ASV45  merged  237
#> ASV1311    12      1   ASV1311  parent  238
#> ASV1552    12      1   ASV1552  parent  239
#> ASV1588    12      1   ASV1588  parent  240
#> ASV359     11      1    ASV359  parent  241
#> ASV440     10      1    ASV545  merged  246
#> ASV580     10      1     ASV64  merged  247
#> ASV858     10      1     ASV19  merged  248
#> ASV1632    10      1   ASV1632  parent  249
#> ASV526      9      1   ASV1552  merged  253
#> ASV1179     9      1   ASV1179  parent  254
#> ASV1428     9      1   ASV1428  parent  255
#> ASV355      8      1    ASV151  merged  259
#> ASV1010     8      1     ASV19  merged  260
#> ASV1542     8      1     ASV12  merged  261
#> ASV1577     8      1   ASV1577  parent  262
#> ASV1262     7      1     ASV12  merged  265
#> ASV27       6      1     ASV27  parent  279
#> ASV137      6      1    ASV137  parent  280
#> ASV930      6      1    ASV930  parent  281
#> ASV962      6      1      ASV2  merged  282
#> ASV1242     6      1   ASV1242  parent  283
#> ASV1628     6      1   ASV1628  parent  284
#> ASV1650     6      1   ASV1650  parent  285
#> ASV903      5      1    ASV637  merged  290
#> ASV984      5      1    ASV113  merged  291
#> ASV1011     5      1   ASV1532  merged  292
#> ASV1224     5      1   ASV1359  merged  293
#> ASV1276     5      1     ASV52  merged  294
#> ASV1468     5      1   ASV1468  parent  295
#> ASV1574     5      1   ASV1574  parent  296
#> ASV1630     5      1   ASV1630  parent  297
#> ASV210      4      1    ASV210  parent  302
#> ASV383      4      1    ASV383  parent  303
#> ASV753      4      1    ASV753  parent  304
#> ASV996      4      1     ASV45  merged  305
#> ASV1069     4      1   ASV1069  parent  306
#> ASV1082     4      1      ASV2  merged  307
#> ASV1109     4      1   ASV1233  merged  308
#> ASV1144     4      1   ASV1144  parent  309
#> ASV1419     4      1   ASV1419  parent  310
#> ASV1467     4      1     ASV19  merged  311
#> ASV42       3      1     ASV42  parent  321
#> ASV104      3      1    ASV104  parent  322
#> ASV203      3      1     ASV64  merged  323
#> ASV491      3      1    ASV210  merged  324
#> ASV517      3      1    ASV517  parent  325
#> ASV598      3      1    ASV598  parent  326
#> ASV626      3      1    ASV357  merged  327
#> ASV650      3      1    ASV378  merged  328
#> ASV790      3      1    ASV505  merged  329
#> ASV839      3      1    ASV839  parent  330
#> ASV1014     3      1    ASV969  merged  331
#> ASV1017     3      1   ASV1017  parent  332
#> ASV1030     3      1   ASV1030  parent  333
#> ASV1085     3      1   ASV1085  parent  334
#> ASV1167     3      1   ASV1167  parent  335
#> ASV1625     3      1     ASV78  merged  336
#> ASV25       2      1     ASV25  parent  342
#> ASV249      2      1     ASV64  merged  343
#> ASV272      2      1    ASV839  merged  344
#> ASV320      2      1    ASV320  parent  345
#> ASV334      2      1    ASV334  parent  346
#> ASV488      2      1    ASV404  merged  347
#> ASV515      2      1    ASV515  parent  348
#> ASV549      2      1    ASV549  parent  349
#> ASV570      2      1    ASV570  parent  350
#> ASV604      2      1   ASV1311  merged  351
#> ASV616      2      1    ASV144  merged  352
#> ASV736      2      1     ASV19  merged  353
#> ASV748      2      1    ASV753  merged  354
#> ASV760      2      1     ASV45  merged  355
#> ASV823      2      1    ASV823  parent  356
#> ASV828      2      1   ASV1192  merged  357
#> ASV829      2      1    ASV829  parent  358
#> ASV891      2      1    ASV891  parent  359
#> ASV914      2      1    ASV914  parent  360
#> ASV961      2      1    ASV292  merged  361
#> ASV979      2      1     ASV19  merged  362
#> ASV986      2      1     ASV19  merged  363
#> ASV1133     2      1   ASV1552  merged  364
#> ASV1270     2      1   ASV1270  parent  365
#> ASV1422     2      1   ASV1422  parent  366
#> ASV1523     2      1   ASV1523  parent  367
#> ASV1576     2      1   ASV1576  parent  368
#> ASV1603     2      1    ASV113  merged  369
#> ASV1726     2      1   ASV1726  parent  370
#> ASV67       1      1     ASV12  merged  371
#> ASV87       1      1     ASV87  parent  372
#> ASV116      1      1     ASV38  merged  373
#> ASV171      1      1    ASV171  parent  374
#> ASV193      1      1    ASV193  parent  375
#> ASV198      1      1    ASV198  parent  376
#> ASV314      1      1    ASV314  parent  377
#> ASV405      1      1     ASV69  merged  378
#> ASV415      1      1    ASV113  merged  379
#> ASV731      1      1    ASV731  parent  380
#> ASV796      1      1   ASV1532  merged  381
#> ASV840      1      1    ASV404  merged  382
#> ASV860      1      1    ASV504  merged  383
#> ASV875      1      1    ASV787  merged  384
#> ASV892      1      1    ASV892  parent  385
#> ASV940      1      1     ASV45  merged  386
#> ASV941      1      1     ASV78  merged  387
#> ASV959      1      1     ASV81  merged  388
#> ASV981      1      1    ASV172  merged  389
#> ASV988      1      1    ASV988  parent  390
#> ASV1015     1      1    ASV637  merged  391
#> ASV1020     1      1    ASV113  merged  392
#> ASV1027     1      1      ASV2  merged  393
#> ASV1031     1      1   ASV1300  merged  394
#> ASV1066     1      1    ASV113  merged  395
#> ASV1107     1      1     ASV38  merged  396
#> ASV1111     1      1    ASV443  merged  397
#> ASV1239     1      1   ASV1239  parent  398
#> ASV1263     1      1   ASV1144  merged  399
#> ASV1302     1      1   ASV1302  parent  400
#> ASV1303     1      1    ASV746  merged  401
#> ASV1341     1      1   ASV1460  merged  402
#> ASV1380     1      1    ASV113  merged  403
#> ASV1388     1      1    ASV618  merged  404
#> ASV1484     1      1     ASV45  merged  405
#> ASV1550     1      1   ASV1550  parent  406
#> ASV1562     1      1   ASV1562  parent  407
#> ASV1635     1      1   ASV1635  parent  408
#> ASV1653     1      1   ASV1653  parent  409
#> ASV1674     1      1   ASV1674  parent  410
#> ASV1712     1      1   ASV1712  parent  411
#> 
#> $res_lulu$original_table
#>         A10.005.B_S188_MERGED.fastq.gz A10.005.H_S189_MERGED.fastq.gz
#> ASV2                                 1                             23
#> ASV8                               578                              9
#> ASV38                                0                              9
#> ASV756                              49                              0
#> ASV12                                0                              1
#> ASV19                                0                              0
#> ASV175                               0                              0
#> ASV18                               30                              0
#> ASV694                               0                              0
#> ASV113                               0                              4
#> ASV378                               0                              0
#> ASV717                               0                              0
#> ASV41                                6                              0
#> ASV65                                0                              0
#> ASV89                                0                              0
#> ASV618                               0                             11
#> ASV170                             104                              0
#> ASV746                               0                              0
#> ASV953                               0                              7
#> ASV28                                1                              0
#> ASV310                               0                              0
#> ASV666                              32                              0
#> ASV1022                              0                              0
#> ASV1192                              5                              0
#> ASV727                               0                              0
#> ASV989                               0                              0
#> ASV635                               0                              0
#> ASV1074                              0                              0
#> ASV1108                              0                              0
#> ASV1483                              0                              1
#> ASV43                                0                              0
#> ASV46                                0                              0
#> ASV365                               0                              0
#> ASV94                              153                              0
#> ASV602                               0                              0
#> ASV178                               0                              0
#> ASV172                               0                              0
#> ASV209                              52                              0
#> ASV741                               0                              0
#> ASV1058                              0                              0
#> ASV261                              25                              0
#> ASV329                               0                              0
#> ASV906                               0                              0
#> ASV546                               0                              0
#> ASV1253                              0                              0
#> ASV24                                0                              0
#> ASV1421                              0                              0
#> ASV816                               0                              0
#> ASV662                               0                              0
#> ASV266                               0                              0
#> ASV162                               0                              0
#> ASV139                               0                              0
#> ASV219                               0                              0
#> ASV201                               0                              0
#> ASV251                               3                              0
#> ASV244                               0                              0
#> ASV507                               0                              0
#> ASV483                               0                              0
#> ASV500                               0                             41
#> ASV577                               0                              0
#> ASV344                               0                              0
#> ASV542                               0                              0
#> ASV85                                0                              0
#> ASV870                               0                              0
#> ASV899                               0                              0
#> ASV52                                0                              0
#> ASV926                               0                              0
#> ASV159                               0                              0
#> ASV333                               0                              0
#> ASV631                               0                              0
#> ASV348                              24                              0
#> ASV817                               0                              0
#> ASV1007                              0                              0
#> ASV543                               0                              0
#> ASV60                                0                              0
#> ASV559                               0                              0
#> ASV1370                              0                              0
#> ASV1459                              0                              0
#> ASV1359                              0                              0
#> ASV292                               0                              0
#> ASV1059                              0                              0
#> ASV443                               0                              0
#> ASV505                               0                              0
#> ASV509                               0                              0
#> ASV1624                              0                              0
#> ASV98                                0                              0
#> ASV368                               0                              0
#> ASV493                               7                              0
#> ASV1279                              0                              0
#> ASV904                               0                              0
#> ASV563                              13                              0
#> ASV1644                             11                              0
#> ASV1369                              0                              0
#> ASV1409                              0                              0
#> ASV724                               0                              0
#> ASV255                               0                              0
#> ASV69                                0                              0
#> ASV151                               0                              0
#> ASV566                               0                              0
#> ASV831                               3                              0
#> ASV1078                              0                              0
#> ASV1230                              0                              3
#> ASV1004                              1                              0
#> ASV586                               0                              0
#> ASV859                               0                              0
#> ASV1278                              1                              0
#> ASV1356                              0                              0
#> ASV801                               0                              0
#> ASV64                                0                              0
#> ASV144                               0                              0
#> ASV131                               0                              0
#> ASV202                               0                              0
#> ASV107                               0                              0
#> ASV322                               0                              0
#> ASV444                               0                              0
#> ASV313                               0                              0
#> ASV569                               0                              0
#> ASV672                             319                              0
#> ASV462                               0                              0
#> ASV223                               0                              0
#> ASV594                               0                              0
#> ASV534                               0                              0
#> ASV950                               0                              0
#> ASV428                               0                              0
#> ASV915                               0                              0
#> ASV346                               0                             10
#> ASV1146                              0                              0
#> ASV942                               0                              0
#> ASV1164                              0                              0
#> ASV1268                              0                              0
#> ASV1152                              0                             83
#> ASV1387                              0                              0
#> ASV592                               0                              4
#> ASV737                               0                              0
#> ASV1458                              0                              0
#> ASV772                               0                              0
#> ASV523                               0                              0
#> ASV766                               0                              0
#> ASV55                                0                              0
#> ASV1236                             11                              0
#> ASV749                               0                              0
#> ASV854                               0                              0
#> ASV743                               0                              0
#> ASV975                               0                              0
#> ASV111                               0                              7
#> ASV1532                              0                              0
#> ASV1493                              0                              0
#> ASV1131                              0                              0
#> ASV1101                              0                              0
#> ASV474                               8                              0
#> ASV744                               0                              0
#> ASV1396                              0                              0
#> ASV758                               0                              0
#> ASV784                               0                              0
#> ASV297                               0                              0
#> ASV91                                0                              0
#> ASV199                               0                              0
#> ASV643                               0                              0
#> ASV787                               0                              0
#> ASV1261                              0                              1
#> ASV1434                              0                              0
#> ASV357                               0                              0
#> ASV927                               0                              0
#> ASV248                               0                              0
#> ASV963                               0                              0
#> ASV1212                              0                              0
#> ASV1321                              0                              0
#> ASV880                               0                              0
#> ASV1052                              0                              0
#> ASV1155                              0                              0
#> ASV685                               0                              0
#> ASV1128                              0                              0
#> ASV81                                0                              0
#> ASV339                               0                              0
#> ASV358                               0                              0
#> ASV420                               0                              0
#> ASV693                               0                              0
#> ASV715                               0                              0
#> ASV884                               0                              0
#> ASV973                               0                              0
#> ASV1218                              0                              0
#> ASV1319                              0                              0
#> ASV1557                              0                              0
#> ASV1315                              0                              0
#> ASV576                               0                              0
#> ASV641                               0                              0
#> ASV1198                              0                              0
#> ASV29                                0                              0
#> ASV286                               0                              0
#> ASV464                               0                              0
#> ASV466                               0                              0
#> ASV1036                              0                              0
#> ASV1247                              0                              0
#> ASV1300                              0                              0
#> ASV1427                              0                              0
#> ASV1690                              0                              0
#> ASV834                               0                              0
#> ASV1176                              0                              0
#> ASV1223                              0                              0
#> ASV1332                              0                              0
#> ASV1365                              0                              1
#> ASV78                                0                              0
#> ASV101                               0                              0
#> ASV168                               0                              0
#> ASV221                               0                              0
#> ASV242                               0                              0
#> ASV253                               0                              0
#> ASV59                                0                              0
#> ASV264                               0                              0
#> ASV45                                0                              0
#> ASV454                               0                              0
#> ASV489                               0                              0
#> ASV524                               0                              0
#> ASV552                               0                            489
#> ASV533                               0                              0
#> ASV619                               0                              0
#> ASV504                               0                              0
#> ASV26                              351                              0
#> ASV692                               0                              0
#> ASV47                                0                              0
#> ASV762                               0                              0
#> ASV461                               0                              0
#> ASV833                               0                              0
#> ASV338                               0                              0
#> ASV898                               0                              0
#> ASV946                             197                              0
#> ASV987                               0                              0
#> ASV998                               0                              0
#> ASV254                               0                              0
#> ASV637                               0                              0
#> ASV1070                              0                              0
#> ASV1092                              0                              0
#> ASV1102                              0                              0
#> ASV49                                0                              0
#> ASV1117                              0                              0
#> ASV1118                              0                              0
#> ASV1132                              0                              0
#> ASV964                               0                              0
#> ASV832                               0                              0
#> ASV1166                              0                              0
#> ASV1182                            114                              0
#> ASV1314                              0                              0
#> ASV969                               0                             92
#> ASV1200                              0                              0
#> ASV624                               0                              0
#> ASV1320                              0                              0
#> ASV999                               0                              0
#> ASV1485                              0                              0
#> ASV1461                              0                              0
#> ASV1313                              0                              0
#> ASV1420                              0                              0
#> ASV105                               0                              0
#> ASV431                               0                              0
#> ASV377                               0                              0
#> ASV1561                              0                              0
#> ASV477                               0                              0
#> ASV1537                              0                              0
#> ASV847                               0                              0
#> ASV722                               0                              0
#> ASV404                               0                              0
#> ASV1510                              0                              0
#> ASV1526                              0                              0
#> ASV1379                              0                              0
#> ASV367                               0                              0
#> ASV1642                              0                              0
#> ASV1233                              0                              0
#> ASV1386                              0                             29
#> ASV1234                              0                              0
#> ASV1340                              0                              0
#> ASV1035                              0                              0
#> ASV1037                              0                              0
#> ASV48                                0                              0
#> ASV711                               0                              0
#> ASV1609                              0                              0
#> ASV33                                0                              0
#> ASV1683                              0                              0
#> ASV673                               0                              0
#> ASV1548                              0                              0
#> ASV1460                              0                              0
#> ASV545                               0                              0
#> ASV1410                              0                              0
#> ASV1462                              0                              0
#> ASV853                               0                              0
#> ASV1067                              0                              0
#> ASV1311                              0                              0
#> ASV1552                              0                              0
#> ASV1588                              0                              0
#> ASV359                               0                              0
#> ASV440                               0                              0
#> ASV580                               0                              0
#> ASV858                               0                              0
#> ASV1632                              0                              0
#> ASV526                               0                              0
#> ASV1179                              0                              0
#> ASV1428                              0                              0
#> ASV355                               0                              0
#> ASV1010                              0                              0
#> ASV1542                              0                              0
#> ASV1577                              0                              0
#> ASV1262                              0                              0
#> ASV27                                0                              0
#> ASV137                               0                              0
#> ASV930                               0                              0
#> ASV962                               0                              0
#> ASV1242                              6                              0
#> ASV1628                              0                              0
#> ASV1650                              0                              0
#> ASV903                               0                              0
#> ASV984                               0                              0
#> ASV1011                              0                              0
#> ASV1224                              0                              0
#> ASV1276                              0                              0
#> ASV1468                              0                              0
#> ASV1574                              0                              0
#> ASV1630                              0                              0
#> ASV210                               0                              0
#> ASV383                               0                              0
#> ASV753                               0                              0
#> ASV996                               0                              0
#> ASV1069                              0                              0
#> ASV1082                              0                              0
#> ASV1109                              0                              0
#> ASV1144                              0                              0
#> ASV1419                              0                              0
#> ASV1467                              0                              0
#> ASV42                                0                              0
#> ASV104                               0                              0
#> ASV203                               0                              0
#> ASV491                               0                              0
#> ASV517                               0                              0
#> ASV598                               0                              0
#> ASV626                               0                              0
#> ASV650                               0                              0
#> ASV790                               0                              0
#> ASV839                               0                              0
#> ASV1014                              0                              3
#> ASV1017                              0                              0
#> ASV1030                              0                              0
#> ASV1085                              0                              0
#> ASV1167                              0                              0
#> ASV1625                              0                              0
#> ASV25                                0                              0
#> ASV249                               0                              0
#> ASV272                               0                              0
#> ASV320                               0                              0
#> ASV334                               0                              0
#> ASV488                               0                              0
#> ASV515                               0                              0
#> ASV549                               0                              0
#> ASV570                               0                              2
#> ASV604                               0                              0
#> ASV616                               0                              0
#> ASV736                               0                              0
#> ASV748                               0                              0
#> ASV760                               0                              0
#> ASV823                               0                              0
#> ASV828                               0                              0
#> ASV829                               2                              0
#> ASV891                               0                              0
#> ASV914                               0                              0
#> ASV961                               0                              0
#> ASV979                               0                              0
#> ASV986                               0                              0
#> ASV1133                              0                              0
#> ASV1270                              0                              0
#> ASV1422                              0                              0
#> ASV1523                              0                              0
#> ASV1576                              0                              0
#> ASV1603                              0                              0
#> ASV1726                              0                              0
#> ASV67                                0                              0
#> ASV87                                0                              0
#> ASV116                               0                              1
#> ASV171                               0                              0
#> ASV193                               0                              0
#> ASV198                               0                              0
#> ASV314                               0                              0
#> ASV405                               0                              0
#> ASV415                               0                              0
#> ASV731                               0                              0
#> ASV796                               0                              0
#> ASV840                               0                              0
#> ASV860                               0                              0
#> ASV875                               0                              0
#> ASV892                               0                              0
#> ASV940                               0                              0
#> ASV941                               0                              0
#> ASV959                               0                              0
#> ASV981                               0                              0
#> ASV988                               0                              1
#> ASV1015                              0                              0
#> ASV1020                              0                              0
#> ASV1027                              0                              0
#> ASV1031                              0                              0
#> ASV1066                              0                              1
#> ASV1107                              0                              0
#> ASV1111                              0                              0
#> ASV1239                              0                              0
#> ASV1263                              0                              0
#> ASV1302                              0                              0
#> ASV1303                              0                              0
#> ASV1341                              0                              0
#> ASV1380                              0                              0
#> ASV1388                              0                              0
#> ASV1484                              0                              0
#> ASV1550                              0                              0
#> ASV1562                              0                              0
#> ASV1635                              0                              0
#> ASV1653                              0                              0
#> ASV1674                              0                              0
#> ASV1712                              0                              0
#>         A10.005.M_S190_MERGED.fastq.gz A12.007_S191_MERGED.fastq.gz
#> ASV2                                 1                          281
#> ASV8                               399                          107
#> ASV38                                0                          337
#> ASV756                               0                            4
#> ASV12                                0                         3186
#> ASV19                                0                         2973
#> ASV175                               0                            0
#> ASV18                                3                           29
#> ASV694                               0                            0
#> ASV113                               0                            0
#> ASV378                               0                            0
#> ASV717                               0                            0
#> ASV41                                5                           22
#> ASV65                                0                            0
#> ASV89                                0                            0
#> ASV618                             131                           56
#> ASV170                              22                            4
#> ASV746                               0                            8
#> ASV953                              53                           26
#> ASV28                                0                            0
#> ASV310                               0                            0
#> ASV666                               3                            0
#> ASV1022                              0                            0
#> ASV1192                              0                            0
#> ASV727                               0                            0
#> ASV989                               0                            3
#> ASV635                               0                            0
#> ASV1074                              1                            2
#> ASV1108                              1                            1
#> ASV1483                              0                            1
#> ASV43                                0                            0
#> ASV46                                0                            0
#> ASV365                               0                            0
#> ASV94                               43                            8
#> ASV602                               0                            2
#> ASV178                               0                          292
#> ASV172                               0                            0
#> ASV209                              13                            2
#> ASV741                               0                            1
#> ASV1058                              0                            0
#> ASV261                               3                            3
#> ASV329                               0                            0
#> ASV906                               2                            5
#> ASV546                               0                            0
#> ASV1253                              0                            0
#> ASV24                                0                            9
#> ASV1421                              0                            0
#> ASV816                               0                            8
#> ASV662                               0                            0
#> ASV266                               0                            1
#> ASV162                               0                            0
#> ASV139                               0                            0
#> ASV219                               0                            0
#> ASV201                               0                            0
#> ASV251                             622                            0
#> ASV244                               0                            0
#> ASV507                               0                            0
#> ASV483                               0                           79
#> ASV500                             435                            0
#> ASV577                               0                          157
#> ASV344                               0                            0
#> ASV542                               0                            0
#> ASV85                                0                           55
#> ASV870                               0                           63
#> ASV899                               0                            0
#> ASV52                                0                            0
#> ASV926                               0                            0
#> ASV159                               0                            0
#> ASV333                               0                          147
#> ASV631                              63                            0
#> ASV348                               1                            0
#> ASV817                               0                            0
#> ASV1007                              0                           68
#> ASV543                               0                            0
#> ASV60                                0                            0
#> ASV559                               0                            0
#> ASV1370                              0                            0
#> ASV1459                              0                            0
#> ASV1359                              0                           46
#> ASV292                               0                            0
#> ASV1059                              0                            0
#> ASV443                               0                            0
#> ASV505                               0                            0
#> ASV509                               0                            0
#> ASV1624                              0                            2
#> ASV98                                0                            0
#> ASV368                               0                            0
#> ASV493                               1                            0
#> ASV1279                              0                            0
#> ASV904                               0                            1
#> ASV563                               0                           17
#> ASV1644                              0                            0
#> ASV1369                              0                            0
#> ASV1409                              0                            4
#> ASV724                               0                            0
#> ASV255                               0                            0
#> ASV69                                0                            0
#> ASV151                               0                            0
#> ASV566                               0                            0
#> ASV831                               0                            0
#> ASV1078                              0                           10
#> ASV1230                              0                            0
#> ASV1004                              0                            0
#> ASV586                               0                            0
#> ASV859                               0                            0
#> ASV1278                              0                            0
#> ASV1356                              0                            0
#> ASV801                               0                            1
#> ASV64                                0                            0
#> ASV144                               0                            0
#> ASV131                               0                            0
#> ASV202                               0                         2283
#> ASV107                               0                            0
#> ASV322                               0                            0
#> ASV444                               0                            0
#> ASV313                               0                            0
#> ASV569                               0                            0
#> ASV672                               0                           23
#> ASV462                               0                            0
#> ASV223                               0                            0
#> ASV594                               0                          292
#> ASV534                               0                            0
#> ASV950                               0                            0
#> ASV428                               0                            0
#> ASV915                               0                            0
#> ASV346                               0                            0
#> ASV1146                              0                            0
#> ASV942                             128                            0
#> ASV1164                              0                            0
#> ASV1268                              0                            0
#> ASV1152                              0                            1
#> ASV1387                              0                            0
#> ASV592                               0                            0
#> ASV737                               0                            0
#> ASV1458                             73                            3
#> ASV772                               0                           71
#> ASV523                               0                            0
#> ASV766                               0                            0
#> ASV55                                0                            0
#> ASV1236                              0                            0
#> ASV749                               0                            0
#> ASV854                               0                            0
#> ASV743                               0                            0
#> ASV975                               0                           44
#> ASV111                               0                            0
#> ASV1532                              0                            0
#> ASV1493                              0                            0
#> ASV1131                              0                           32
#> ASV1101                              0                            0
#> ASV474                               0                            0
#> ASV744                               0                            0
#> ASV1396                              0                            0
#> ASV758                               0                            0
#> ASV784                               0                            0
#> ASV297                               0                            0
#> ASV91                                0                            0
#> ASV199                               0                            0
#> ASV643                               0                           10
#> ASV787                               0                            0
#> ASV1261                              0                            0
#> ASV1434                              0                           12
#> ASV357                               0                            0
#> ASV927                               0                            0
#> ASV248                               0                            9
#> ASV963                               0                            0
#> ASV1212                              0                            0
#> ASV1321                              0                            0
#> ASV880                               0                            0
#> ASV1052                              0                            0
#> ASV1155                              0                            7
#> ASV685                               0                            0
#> ASV1128                              0                            0
#> ASV81                                0                            0
#> ASV339                               0                            0
#> ASV358                               0                            0
#> ASV420                               0                            0
#> ASV693                               0                            1
#> ASV715                               0                            0
#> ASV884                               0                            0
#> ASV973                               0                            0
#> ASV1218                              0                            0
#> ASV1319                              0                            1
#> ASV1557                              0                            0
#> ASV1315                              0                            0
#> ASV576                               0                            3
#> ASV641                               0                            0
#> ASV1198                              0                            0
#> ASV29                                0                            0
#> ASV286                               0                            0
#> ASV464                               0                            0
#> ASV466                               0                            0
#> ASV1036                              0                            0
#> ASV1247                              0                            1
#> ASV1300                              0                            0
#> ASV1427                              0                            0
#> ASV1690                              0                            1
#> ASV834                               0                            0
#> ASV1176                              0                            0
#> ASV1223                              0                            0
#> ASV1332                              0                            0
#> ASV1365                              0                            0
#> ASV78                                0                            0
#> ASV101                               0                            0
#> ASV168                               0                            0
#> ASV221                               0                            0
#> ASV242                               0                            0
#> ASV253                               0                            0
#> ASV59                                0                            0
#> ASV264                               0                            0
#> ASV45                                0                            0
#> ASV454                               0                            0
#> ASV489                               0                            0
#> ASV524                               0                            0
#> ASV552                               0                            0
#> ASV533                               0                            0
#> ASV619                               0                            0
#> ASV504                               0                            0
#> ASV26                                0                            0
#> ASV692                               0                            0
#> ASV47                                0                            0
#> ASV762                               0                            0
#> ASV461                               0                            0
#> ASV833                               0                            0
#> ASV338                               0                            0
#> ASV898                               0                            0
#> ASV946                               0                            0
#> ASV987                               0                            0
#> ASV998                               0                            0
#> ASV254                               0                            0
#> ASV637                               0                          166
#> ASV1070                              0                            0
#> ASV1092                            146                            0
#> ASV1102                              0                            0
#> ASV49                                0                            0
#> ASV1117                              0                            0
#> ASV1118                              0                            0
#> ASV1132                              0                            0
#> ASV964                               0                            0
#> ASV832                               0                            0
#> ASV1166                              0                            0
#> ASV1182                              0                            0
#> ASV1314                              0                            0
#> ASV969                               0                            0
#> ASV1200                              0                            0
#> ASV624                               0                           91
#> ASV1320                              0                            0
#> ASV999                               0                            0
#> ASV1485                              0                            0
#> ASV1461                              0                            0
#> ASV1313                              0                            0
#> ASV1420                              0                            0
#> ASV105                               0                            0
#> ASV431                               0                            0
#> ASV377                               0                            0
#> ASV1561                              0                            0
#> ASV477                               0                            0
#> ASV1537                              0                           53
#> ASV847                               0                            0
#> ASV722                               0                            0
#> ASV404                               0                            0
#> ASV1510                              0                            0
#> ASV1526                              0                            0
#> ASV1379                              0                            0
#> ASV367                               0                            0
#> ASV1642                              0                            0
#> ASV1233                              0                            0
#> ASV1386                              0                            0
#> ASV1234                              0                            0
#> ASV1340                              0                            0
#> ASV1035                              0                           25
#> ASV1037                              0                            0
#> ASV48                                0                            0
#> ASV711                               0                            0
#> ASV1609                              0                            0
#> ASV33                                0                            0
#> ASV1683                              0                            0
#> ASV673                               0                            0
#> ASV1548                              0                            0
#> ASV1460                              0                            0
#> ASV545                               0                            0
#> ASV1410                              0                            0
#> ASV1462                              0                            0
#> ASV853                               0                           12
#> ASV1067                              0                            0
#> ASV1311                              0                            0
#> ASV1552                              0                            0
#> ASV1588                              0                            0
#> ASV359                               0                            0
#> ASV440                               0                            0
#> ASV580                               0                            0
#> ASV858                               0                           10
#> ASV1632                              0                            0
#> ASV526                               0                            0
#> ASV1179                              0                            0
#> ASV1428                              0                            0
#> ASV355                               0                            0
#> ASV1010                              0                            8
#> ASV1542                              0                            8
#> ASV1577                              0                            0
#> ASV1262                              0                            7
#> ASV27                                0                            0
#> ASV137                               0                            0
#> ASV930                               0                            0
#> ASV962                               0                            6
#> ASV1242                              0                            0
#> ASV1628                              0                            0
#> ASV1650                              0                            0
#> ASV903                               0                            5
#> ASV984                               0                            0
#> ASV1011                              0                            0
#> ASV1224                              0                            0
#> ASV1276                              0                            0
#> ASV1468                              0                            0
#> ASV1574                              0                            5
#> ASV1630                              0                            0
#> ASV210                               0                            0
#> ASV383                               0                            0
#> ASV753                               0                            0
#> ASV996                               0                            0
#> ASV1069                              0                            0
#> ASV1082                              0                            4
#> ASV1109                              0                            0
#> ASV1144                              0                            0
#> ASV1419                              0                            0
#> ASV1467                              0                            4
#> ASV42                                0                            0
#> ASV104                               0                            0
#> ASV203                               0                            0
#> ASV491                               0                            0
#> ASV517                               0                            0
#> ASV598                               0                            0
#> ASV626                               0                            0
#> ASV650                               0                            0
#> ASV790                               0                            0
#> ASV839                               0                            0
#> ASV1014                              0                            0
#> ASV1017                              0                            0
#> ASV1030                              0                            0
#> ASV1085                              0                            3
#> ASV1167                              0                            0
#> ASV1625                              0                            0
#> ASV25                                0                            0
#> ASV249                               0                            0
#> ASV272                               0                            0
#> ASV320                               0                            2
#> ASV334                               0                            0
#> ASV488                               0                            0
#> ASV515                               0                            0
#> ASV549                               0                            0
#> ASV570                               0                            0
#> ASV604                               0                            0
#> ASV616                               0                            0
#> ASV736                               0                            2
#> ASV748                               0                            0
#> ASV760                               0                            0
#> ASV823                               0                            0
#> ASV828                               0                            0
#> ASV829                               0                            0
#> ASV891                               0                            0
#> ASV914                               0                            0
#> ASV961                               0                            0
#> ASV979                               0                            2
#> ASV986                               0                            2
#> ASV1133                              0                            0
#> ASV1270                              0                            0
#> ASV1422                              0                            0
#> ASV1523                              0                            0
#> ASV1576                              0                            0
#> ASV1603                              0                            0
#> ASV1726                              0                            0
#> ASV67                                0                            1
#> ASV87                                0                            0
#> ASV116                               0                            0
#> ASV171                               0                            0
#> ASV193                               0                            0
#> ASV198                               0                            0
#> ASV314                               0                            1
#> ASV405                               0                            0
#> ASV415                               0                            0
#> ASV731                               0                            0
#> ASV796                               0                            0
#> ASV840                               0                            0
#> ASV860                               0                            0
#> ASV875                               0                            0
#> ASV892                               0                            0
#> ASV940                               0                            0
#> ASV941                               0                            0
#> ASV959                               0                            0
#> ASV981                               0                            0
#> ASV988                               0                            0
#> ASV1015                              0                            1
#> ASV1020                              0                            0
#> ASV1027                              0                            1
#> ASV1031                              0                            0
#> ASV1066                              0                            0
#> ASV1107                              0                            0
#> ASV1111                              0                            0
#> ASV1239                              0                            0
#> ASV1263                              0                            0
#> ASV1302                              0                            0
#> ASV1303                              0                            0
#> ASV1341                              0                            0
#> ASV1380                              0                            0
#> ASV1388                              0                            0
#> ASV1484                              0                            0
#> ASV1550                              0                            0
#> ASV1562                              0                            0
#> ASV1635                              0                            0
#> ASV1653                              0                            0
#> ASV1674                              0                            0
#> ASV1712                              0                            0
#>         A12.007.B_S2_MERGED.fastq.gz A15.004_S3_MERGED.fastq.gz
#> ASV2                             260                          0
#> ASV8                             829                          2
#> ASV38                              0                          0
#> ASV756                             0                          0
#> ASV12                              0                          0
#> ASV19                              0                          0
#> ASV175                             0                          0
#> ASV18                           1040                          0
#> ASV694                             0                          0
#> ASV113                             0                          4
#> ASV378                             0                          0
#> ASV717                             0                          0
#> ASV41                              0                          0
#> ASV65                              0                          0
#> ASV89                              0                          0
#> ASV618                             0                          0
#> ASV170                             0                          0
#> ASV746                             0                          0
#> ASV953                             0                          0
#> ASV28                              0                          0
#> ASV310                             0                          0
#> ASV666                             4                          0
#> ASV1022                            0                          0
#> ASV1192                            0                          0
#> ASV727                             0                          0
#> ASV989                             0                          0
#> ASV635                             0                          0
#> ASV1074                            0                          0
#> ASV1108                            0                          0
#> ASV1483                            0                          0
#> ASV43                              0                          0
#> ASV46                              0                          0
#> ASV365                             0                          0
#> ASV94                              0                          0
#> ASV602                             0                          0
#> ASV178                             0                          0
#> ASV172                             0                          0
#> ASV209                             0                          0
#> ASV741                             0                          0
#> ASV1058                            0                          0
#> ASV261                             0                          0
#> ASV329                             0                          1
#> ASV906                             0                          0
#> ASV546                             0                          0
#> ASV1253                            0                          0
#> ASV24                              0                          1
#> ASV1421                            0                          0
#> ASV816                             0                          0
#> ASV662                             6                          1
#> ASV266                             0                          0
#> ASV162                             0                          0
#> ASV139                             0                          0
#> ASV219                             0                          0
#> ASV201                             0                          0
#> ASV251                             0                          0
#> ASV244                             0                          0
#> ASV507                             0                          0
#> ASV483                             0                          0
#> ASV500                             0                          0
#> ASV577                             0                          0
#> ASV344                             0                          0
#> ASV542                             0                          0
#> ASV85                              0                          0
#> ASV870                             0                          0
#> ASV899                             0                          0
#> ASV52                              0                          0
#> ASV926                             0                          0
#> ASV159                             0                          0
#> ASV333                             0                          0
#> ASV631                             0                          0
#> ASV348                             0                          0
#> ASV817                             0                          0
#> ASV1007                            0                          0
#> ASV543                             0                          0
#> ASV60                              0                         21
#> ASV559                             0                          0
#> ASV1370                            0                          0
#> ASV1459                            0                          0
#> ASV1359                            0                          0
#> ASV292                             0                          0
#> ASV1059                            0                          0
#> ASV443                             0                          0
#> ASV505                             0                         40
#> ASV509                             0                          7
#> ASV1624                           37                          0
#> ASV98                              0                          0
#> ASV368                             0                          0
#> ASV493                             0                          0
#> ASV1279                            0                          0
#> ASV904                             0                          0
#> ASV563                             0                          0
#> ASV1644                            1                          0
#> ASV1369                            0                          0
#> ASV1409                            0                          0
#> ASV724                             0                          2
#> ASV255                             0                          0
#> ASV69                              0                          2
#> ASV151                             0                          0
#> ASV566                             0                          0
#> ASV831                             8                          0
#> ASV1078                            0                          0
#> ASV1230                            0                          0
#> ASV1004                            7                          0
#> ASV586                             0                          1
#> ASV859                             0                          2
#> ASV1278                            0                          0
#> ASV1356                            0                          0
#> ASV801                             0                          0
#> ASV64                              0                          0
#> ASV144                             0                          0
#> ASV131                             0                          0
#> ASV202                             0                          0
#> ASV107                             0                          0
#> ASV322                             0                          0
#> ASV444                             0                          0
#> ASV313                             0                          0
#> ASV569                             0                          0
#> ASV672                             0                          0
#> ASV462                             0                          0
#> ASV223                             0                          0
#> ASV594                             0                          0
#> ASV534                             0                          0
#> ASV950                             0                          0
#> ASV428                             0                          0
#> ASV915                             0                          0
#> ASV346                             0                          0
#> ASV1146                            0                          0
#> ASV942                             0                          0
#> ASV1164                            0                          0
#> ASV1268                            0                          0
#> ASV1152                            0                          0
#> ASV1387                            0                          0
#> ASV592                             0                          0
#> ASV737                             0                          0
#> ASV1458                            0                          0
#> ASV772                             0                          0
#> ASV523                             0                          2
#> ASV766                             0                          0
#> ASV55                              0                          0
#> ASV1236                           47                          0
#> ASV749                             7                          0
#> ASV854                             0                          0
#> ASV743                             0                          0
#> ASV975                             0                          0
#> ASV111                             0                          0
#> ASV1532                           36                          0
#> ASV1493                            0                          0
#> ASV1131                            0                          0
#> ASV1101                            0                          0
#> ASV474                             0                          0
#> ASV744                             0                          0
#> ASV1396                            0                          0
#> ASV758                             0                          6
#> ASV784                             0                          0
#> ASV297                             0                          0
#> ASV91                              0                         12
#> ASV199                             0                          0
#> ASV643                             0                          0
#> ASV787                             0                          0
#> ASV1261                            0                          0
#> ASV1434                            0                          0
#> ASV357                             0                          0
#> ASV927                             0                          0
#> ASV248                             0                          0
#> ASV963                             0                          0
#> ASV1212                            0                          0
#> ASV1321                            0                          0
#> ASV880                             0                          0
#> ASV1052                            0                          6
#> ASV1155                            0                          0
#> ASV685                             0                          0
#> ASV1128                            0                          5
#> ASV81                              0                          0
#> ASV339                             0                          0
#> ASV358                             0                          0
#> ASV420                             0                          4
#> ASV693                             0                          0
#> ASV715                             0                          0
#> ASV884                             0                          3
#> ASV973                             0                          0
#> ASV1218                            0                          0
#> ASV1319                            0                          0
#> ASV1557                            0                          0
#> ASV1315                            0                          0
#> ASV576                             1                          0
#> ASV641                             0                          3
#> ASV1198                            0                          2
#> ASV29                              0                          0
#> ASV286                             0                          0
#> ASV464                             0                          0
#> ASV466                             0                          0
#> ASV1036                            0                          2
#> ASV1247                            2                          0
#> ASV1300                            1                          0
#> ASV1427                            0                          1
#> ASV1690                            0                          0
#> ASV834                             0                          0
#> ASV1176                            0                          0
#> ASV1223                            1                          0
#> ASV1332                            0                          0
#> ASV1365                            0                          0
#> ASV78                              0                          0
#> ASV101                             0                          0
#> ASV168                             0                          0
#> ASV221                             0                          0
#> ASV242                             0                          0
#> ASV253                             0                          0
#> ASV59                              0                          0
#> ASV264                             0                          0
#> ASV45                              0                          0
#> ASV454                             0                          0
#> ASV489                             0                          0
#> ASV524                             0                          0
#> ASV552                             0                          0
#> ASV533                             0                          0
#> ASV619                             0                          0
#> ASV504                             0                          0
#> ASV26                              0                          0
#> ASV692                             0                          0
#> ASV47                              0                          0
#> ASV762                             0                          0
#> ASV461                             0                          0
#> ASV833                             0                          0
#> ASV338                             0                          0
#> ASV898                             0                          0
#> ASV946                             0                          0
#> ASV987                             0                          0
#> ASV998                             0                          0
#> ASV254                             0                          0
#> ASV637                             0                          0
#> ASV1070                            0                          0
#> ASV1092                            0                          0
#> ASV1102                            0                          0
#> ASV49                              0                          0
#> ASV1117                            0                          0
#> ASV1118                            0                          0
#> ASV1132                            0                          0
#> ASV964                             0                          0
#> ASV832                             0                          0
#> ASV1166                          129                          0
#> ASV1182                            0                          0
#> ASV1314                            0                          0
#> ASV969                             0                          0
#> ASV1200                            0                          0
#> ASV624                             0                          0
#> ASV1320                            0                          0
#> ASV999                             0                          0
#> ASV1485                            0                          0
#> ASV1461                            0                          0
#> ASV1313                            0                          0
#> ASV1420                            0                          0
#> ASV105                             0                          0
#> ASV431                             0                          0
#> ASV377                             0                          0
#> ASV1561                            0                          0
#> ASV477                             0                          0
#> ASV1537                            0                          0
#> ASV847                             0                          0
#> ASV722                             0                          0
#> ASV404                             0                          0
#> ASV1510                            0                          0
#> ASV1526                            0                          0
#> ASV1379                            0                          0
#> ASV367                             0                          0
#> ASV1642                            0                          0
#> ASV1233                            0                          0
#> ASV1386                            0                          0
#> ASV1234                            0                          0
#> ASV1340                            0                          0
#> ASV1035                            0                          0
#> ASV1037                            0                          0
#> ASV48                              0                          0
#> ASV711                             0                          0
#> ASV1609                            0                          0
#> ASV33                              0                          0
#> ASV1683                            0                          0
#> ASV673                             0                          0
#> ASV1548                            0                          0
#> ASV1460                            0                          0
#> ASV545                             0                          0
#> ASV1410                            0                         13
#> ASV1462                            0                          0
#> ASV853                             0                          0
#> ASV1067                            0                          0
#> ASV1311                            0                         12
#> ASV1552                            0                          0
#> ASV1588                            0                          0
#> ASV359                             0                          0
#> ASV440                             0                          0
#> ASV580                             0                          0
#> ASV858                             0                          0
#> ASV1632                            0                          0
#> ASV526                             0                          0
#> ASV1179                            0                          0
#> ASV1428                            0                          0
#> ASV355                             0                          0
#> ASV1010                            0                          0
#> ASV1542                            0                          0
#> ASV1577                            0                          0
#> ASV1262                            0                          0
#> ASV27                              0                          0
#> ASV137                             0                          0
#> ASV930                             6                          0
#> ASV962                             0                          0
#> ASV1242                            0                          0
#> ASV1628                            0                          6
#> ASV1650                            0                          0
#> ASV903                             0                          0
#> ASV984                             0                          0
#> ASV1011                            5                          0
#> ASV1224                            0                          0
#> ASV1276                            0                          0
#> ASV1468                            0                          0
#> ASV1574                            0                          0
#> ASV1630                            0                          5
#> ASV210                             4                          0
#> ASV383                             0                          0
#> ASV753                             0                          4
#> ASV996                             0                          0
#> ASV1069                            0                          0
#> ASV1082                            0                          0
#> ASV1109                            0                          0
#> ASV1144                            4                          0
#> ASV1419                            0                          4
#> ASV1467                            0                          0
#> ASV42                              0                          3
#> ASV104                             0                          0
#> ASV203                             0                          0
#> ASV491                             3                          0
#> ASV517                             0                          0
#> ASV598                             0                          0
#> ASV626                             0                          0
#> ASV650                             0                          0
#> ASV790                             0                          3
#> ASV839                             0                          0
#> ASV1014                            0                          0
#> ASV1017                            0                          0
#> ASV1030                            0                          0
#> ASV1085                            0                          0
#> ASV1167                            0                          0
#> ASV1625                            0                          0
#> ASV25                              0                          0
#> ASV249                             0                          0
#> ASV272                             0                          0
#> ASV320                             0                          0
#> ASV334                             0                          0
#> ASV488                             0                          0
#> ASV515                             0                          0
#> ASV549                             0                          0
#> ASV570                             0                          0
#> ASV604                             0                          2
#> ASV616                             0                          0
#> ASV736                             0                          0
#> ASV748                             0                          2
#> ASV760                             0                          0
#> ASV823                             0                          0
#> ASV828                             0                          0
#> ASV829                             0                          0
#> ASV891                             0                          0
#> ASV914                             0                          2
#> ASV961                             0                          0
#> ASV979                             0                          0
#> ASV986                             0                          0
#> ASV1133                            0                          0
#> ASV1270                            0                          0
#> ASV1422                            0                          0
#> ASV1523                            0                          0
#> ASV1576                            0                          0
#> ASV1603                            0                          2
#> ASV1726                            0                          0
#> ASV67                              0                          0
#> ASV87                              0                          0
#> ASV116                             0                          0
#> ASV171                             0                          0
#> ASV193                             0                          0
#> ASV198                             0                          1
#> ASV314                             0                          0
#> ASV405                             0                          0
#> ASV415                             0                          0
#> ASV731                             0                          0
#> ASV796                             1                          0
#> ASV840                             0                          0
#> ASV860                             0                          0
#> ASV875                             0                          0
#> ASV892                             0                          1
#> ASV940                             0                          0
#> ASV941                             0                          0
#> ASV959                             0                          0
#> ASV981                             0                          0
#> ASV988                             0                          0
#> ASV1015                            0                          0
#> ASV1020                            0                          1
#> ASV1027                            0                          0
#> ASV1031                            0                          0
#> ASV1066                            0                          0
#> ASV1107                            0                          0
#> ASV1111                            0                          0
#> ASV1239                            0                          1
#> ASV1263                            1                          0
#> ASV1302                            0                          0
#> ASV1303                            0                          0
#> ASV1341                            0                          0
#> ASV1380                            0                          0
#> ASV1388                            0                          0
#> ASV1484                            0                          0
#> ASV1550                            0                          0
#> ASV1562                            0                          0
#> ASV1635                            0                          0
#> ASV1653                            0                          0
#> ASV1674                            0                          0
#> ASV1712                            0                          0
#>         A8.005_S4_MERGED.fastq.gz A9.012_S5_MERGED.fastq.gz
#> ASV2                            5                         1
#> ASV8                         1341                         0
#> ASV38                           0                         0
#> ASV756                          2                         0
#> ASV12                           1                         0
#> ASV19                           0                         0
#> ASV175                          0                         0
#> ASV18                          42                         0
#> ASV694                          0                         0
#> ASV113                          0                         0
#> ASV378                          0                         0
#> ASV717                          0                         0
#> ASV41                           0                         0
#> ASV65                           0                         0
#> ASV89                        1179                         0
#> ASV618                          0                         0
#> ASV170                        218                         0
#> ASV746                          0                         0
#> ASV953                          0                         0
#> ASV28                           0                         0
#> ASV310                          0                         0
#> ASV666                          7                         0
#> ASV1022                         0                         0
#> ASV1192                         0                         0
#> ASV727                          0                         0
#> ASV989                          0                         0
#> ASV635                          0                         0
#> ASV1074                         0                         0
#> ASV1108                         0                         0
#> ASV1483                         0                         0
#> ASV43                           0                         0
#> ASV46                           0                         0
#> ASV365                          0                         0
#> ASV94                         474                         0
#> ASV602                          0                         0
#> ASV178                          0                         0
#> ASV172                          0                         0
#> ASV209                        197                         0
#> ASV741                          0                         0
#> ASV1058                         0                         0
#> ASV261                         71                         0
#> ASV329                          5                         0
#> ASV906                          0                         0
#> ASV546                          0                         0
#> ASV1253                         0                         0
#> ASV24                           0                         0
#> ASV1421                         0                         0
#> ASV816                          0                         0
#> ASV662                          0                         4
#> ASV266                          0                         0
#> ASV162                          0                         0
#> ASV139                          0                         0
#> ASV219                          0                         0
#> ASV201                          0                         0
#> ASV251                          0                         0
#> ASV244                          0                         0
#> ASV507                          0                         0
#> ASV483                          0                         0
#> ASV500                          0                         1
#> ASV577                          0                         0
#> ASV344                          0                         0
#> ASV542                          0                         0
#> ASV85                           0                         0
#> ASV870                          0                         0
#> ASV899                          0                         0
#> ASV52                           0                         0
#> ASV926                          0                         0
#> ASV159                          0                         0
#> ASV333                          0                         0
#> ASV631                         28                        21
#> ASV348                         83                         0
#> ASV817                          0                         0
#> ASV1007                         0                         0
#> ASV543                          0                         0
#> ASV60                           0                         0
#> ASV559                          0                         0
#> ASV1370                         0                         0
#> ASV1459                         0                         0
#> ASV1359                         0                         0
#> ASV292                          0                         0
#> ASV1059                         0                         0
#> ASV443                          0                         0
#> ASV505                          0                         2
#> ASV509                          0                         8
#> ASV1624                         0                         0
#> ASV98                           0                         0
#> ASV368                          0                         0
#> ASV493                         28                         0
#> ASV1279                         0                         0
#> ASV904                          1                         0
#> ASV563                          2                         0
#> ASV1644                         0                         0
#> ASV1369                         0                         0
#> ASV1409                         0                         0
#> ASV724                          0                         4
#> ASV255                          0                         0
#> ASV69                           0                         0
#> ASV151                          0                         0
#> ASV566                          0                         0
#> ASV831                          1                         0
#> ASV1078                         0                         0
#> ASV1230                         1                         0
#> ASV1004                         0                         0
#> ASV586                          0                         1
#> ASV859                          0                         2
#> ASV1278                         0                         0
#> ASV1356                         0                         0
#> ASV801                          0                         0
#> ASV64                         116                         0
#> ASV144                          0                         0
#> ASV131                          0                         0
#> ASV202                          0                         0
#> ASV107                          0                         0
#> ASV322                          0                         0
#> ASV444                          0                         0
#> ASV313                        555                         0
#> ASV569                          0                         0
#> ASV672                          0                         0
#> ASV462                        322                         0
#> ASV223                          0                         0
#> ASV594                          0                         0
#> ASV534                        272                         0
#> ASV950                          0                         0
#> ASV428                          0                         0
#> ASV915                          0                         0
#> ASV346                          0                         0
#> ASV1146                         0                         0
#> ASV942                          0                         0
#> ASV1164                         0                         0
#> ASV1268                         0                         0
#> ASV1152                         0                         0
#> ASV1387                         0                         0
#> ASV592                          0                         0
#> ASV737                          0                         0
#> ASV1458                         0                         0
#> ASV772                          0                         0
#> ASV523                          0                         0
#> ASV766                          0                         0
#> ASV55                           0                         0
#> ASV1236                         0                         0
#> ASV749                          0                         0
#> ASV854                          0                         0
#> ASV743                          0                         0
#> ASV975                          0                         0
#> ASV111                          0                         0
#> ASV1532                         0                         0
#> ASV1493                         0                         0
#> ASV1131                         0                         0
#> ASV1101                        30                         0
#> ASV474                         12                         0
#> ASV744                          0                         0
#> ASV1396                         7                         0
#> ASV758                          0                         9
#> ASV784                          0                         0
#> ASV297                          0                         0
#> ASV91                           0                         0
#> ASV199                          0                         0
#> ASV643                          0                         0
#> ASV787                          0                         7
#> ASV1261                         0                         0
#> ASV1434                         0                         0
#> ASV357                          0                         0
#> ASV927                          0                         0
#> ASV248                          0                         0
#> ASV963                          0                         0
#> ASV1212                         0                         0
#> ASV1321                         0                         0
#> ASV880                          0                         0
#> ASV1052                         0                         3
#> ASV1155                         0                         0
#> ASV685                          0                         0
#> ASV1128                         0                         2
#> ASV81                           0                         0
#> ASV339                          0                         0
#> ASV358                          0                         0
#> ASV420                          0                         2
#> ASV693                          0                         0
#> ASV715                          3                         0
#> ASV884                          0                         3
#> ASV973                          0                         0
#> ASV1218                         0                         1
#> ASV1319                         0                         0
#> ASV1557                         0                         0
#> ASV1315                         0                         0
#> ASV576                          0                         0
#> ASV641                          0                         1
#> ASV1198                         0                         0
#> ASV29                           0                         0
#> ASV286                          0                         0
#> ASV464                          0                         0
#> ASV466                          0                         0
#> ASV1036                         0                         0
#> ASV1247                         0                         0
#> ASV1300                         0                         0
#> ASV1427                         0                         2
#> ASV1690                         0                         0
#> ASV834                          0                         0
#> ASV1176                         0                         0
#> ASV1223                         0                         0
#> ASV1332                         0                         0
#> ASV1365                         1                         0
#> ASV78                           0                         0
#> ASV101                          0                         0
#> ASV168                          0                         0
#> ASV221                          0                         0
#> ASV242                          0                      1692
#> ASV253                          0                         0
#> ASV59                           0                         0
#> ASV264                          0                         0
#> ASV45                           0                         0
#> ASV454                          0                         0
#> ASV489                          0                         0
#> ASV524                          0                         0
#> ASV552                          0                         0
#> ASV533                          0                         0
#> ASV619                          0                         0
#> ASV504                        404                         0
#> ASV26                           0                         0
#> ASV692                          0                         0
#> ASV47                           0                         0
#> ASV762                          0                         0
#> ASV461                          0                         0
#> ASV833                          0                         0
#> ASV338                          0                         0
#> ASV898                          0                       201
#> ASV946                          0                         0
#> ASV987                          0                         0
#> ASV998                          0                         0
#> ASV254                          0                         0
#> ASV637                          0                         0
#> ASV1070                       153                         0
#> ASV1092                         0                         0
#> ASV1102                         0                         0
#> ASV49                           0                         0
#> ASV1117                         0                         0
#> ASV1118                         0                         0
#> ASV1132                         0                         0
#> ASV964                          0                         0
#> ASV832                        130                         0
#> ASV1166                         0                         0
#> ASV1182                         0                         0
#> ASV1314                         0                         0
#> ASV969                          0                         0
#> ASV1200                         0                         0
#> ASV624                          0                         0
#> ASV1320                         0                        91
#> ASV999                          0                         0
#> ASV1485                         0                         0
#> ASV1461                         0                         0
#> ASV1313                         0                         0
#> ASV1420                         0                        70
#> ASV105                          0                         0
#> ASV431                         65                         0
#> ASV377                          0                         0
#> ASV1561                         0                         0
#> ASV477                         57                         0
#> ASV1537                         0                         0
#> ASV847                          0                         0
#> ASV722                          0                         0
#> ASV404                          0                         0
#> ASV1510                         0                        39
#> ASV1526                         0                         0
#> ASV1379                         0                         0
#> ASV367                          0                         0
#> ASV1642                         0                        36
#> ASV1233                         0                         0
#> ASV1386                         0                         0
#> ASV1234                         0                         0
#> ASV1340                         0                         0
#> ASV1035                         0                         0
#> ASV1037                         0                         0
#> ASV48                           0                         0
#> ASV711                          0                         0
#> ASV1609                         0                         0
#> ASV33                           0                         0
#> ASV1683                         0                         0
#> ASV673                          0                         0
#> ASV1548                         0                         0
#> ASV1460                         0                         0
#> ASV545                          0                         0
#> ASV1410                         0                         0
#> ASV1462                         0                         0
#> ASV853                          0                         0
#> ASV1067                         0                         0
#> ASV1311                         0                         0
#> ASV1552                         0                         0
#> ASV1588                         0                         0
#> ASV359                         11                         0
#> ASV440                          0                         0
#> ASV580                          0                         0
#> ASV858                          0                         0
#> ASV1632                         0                         0
#> ASV526                          0                         0
#> ASV1179                         0                         0
#> ASV1428                         0                         0
#> ASV355                          0                         0
#> ASV1010                         0                         0
#> ASV1542                         0                         0
#> ASV1577                         0                         0
#> ASV1262                         0                         0
#> ASV27                           0                         0
#> ASV137                          0                         0
#> ASV930                          0                         0
#> ASV962                          0                         0
#> ASV1242                         0                         0
#> ASV1628                         0                         0
#> ASV1650                         0                         0
#> ASV903                          0                         0
#> ASV984                          0                         0
#> ASV1011                         0                         0
#> ASV1224                         0                         0
#> ASV1276                         0                         0
#> ASV1468                         0                         0
#> ASV1574                         0                         0
#> ASV1630                         0                         0
#> ASV210                          0                         0
#> ASV383                          0                         0
#> ASV753                          0                         0
#> ASV996                          0                         0
#> ASV1069                         0                         0
#> ASV1082                         0                         0
#> ASV1109                         0                         0
#> ASV1144                         0                         0
#> ASV1419                         0                         0
#> ASV1467                         0                         0
#> ASV42                           0                         0
#> ASV104                          0                         0
#> ASV203                          0                         0
#> ASV491                          0                         0
#> ASV517                          0                         0
#> ASV598                          0                         0
#> ASV626                          0                         0
#> ASV650                          0                         0
#> ASV790                          0                         0
#> ASV839                          0                         0
#> ASV1014                         0                         0
#> ASV1017                         0                         0
#> ASV1030                         0                         3
#> ASV1085                         0                         0
#> ASV1167                         0                         0
#> ASV1625                         0                         0
#> ASV25                           0                         0
#> ASV249                          0                         0
#> ASV272                          0                         0
#> ASV320                          0                         0
#> ASV334                          0                         0
#> ASV488                          0                         0
#> ASV515                          0                         2
#> ASV549                          0                         0
#> ASV570                          0                         0
#> ASV604                          0                         0
#> ASV616                          0                         0
#> ASV736                          0                         0
#> ASV748                          0                         0
#> ASV760                          0                         0
#> ASV823                          0                         0
#> ASV828                          0                         0
#> ASV829                          0                         0
#> ASV891                          0                         0
#> ASV914                          0                         0
#> ASV961                          0                         0
#> ASV979                          0                         0
#> ASV986                          0                         0
#> ASV1133                         0                         0
#> ASV1270                         0                         0
#> ASV1422                         0                         0
#> ASV1523                         0                         0
#> ASV1576                         0                         2
#> ASV1603                         0                         0
#> ASV1726                         0                         0
#> ASV67                           0                         0
#> ASV87                           0                         0
#> ASV116                          0                         0
#> ASV171                          0                         0
#> ASV193                          0                         0
#> ASV198                          0                         0
#> ASV314                          0                         0
#> ASV405                          0                         0
#> ASV415                          0                         0
#> ASV731                          0                         0
#> ASV796                          0                         0
#> ASV840                          0                         0
#> ASV860                          1                         0
#> ASV875                          0                         1
#> ASV892                          0                         0
#> ASV940                          0                         0
#> ASV941                          0                         0
#> ASV959                          0                         0
#> ASV981                          0                         0
#> ASV988                          0                         0
#> ASV1015                         0                         0
#> ASV1020                         0                         0
#> ASV1027                         0                         0
#> ASV1031                         0                         0
#> ASV1066                         0                         0
#> ASV1107                         0                         0
#> ASV1111                         0                         0
#> ASV1239                         0                         0
#> ASV1263                         0                         0
#> ASV1302                         0                         0
#> ASV1303                         0                         0
#> ASV1341                         0                         0
#> ASV1380                         0                         0
#> ASV1388                         0                         0
#> ASV1484                         0                         0
#> ASV1550                         0                         0
#> ASV1562                         0                         0
#> ASV1635                         0                         0
#> ASV1653                         0                         0
#> ASV1674                         0                         0
#> ASV1712                         0                         0
#>         AB29.ABMX.H_S6_MERGED.fastq.gz AC27.013_S7_MERGED.fastq.gz
#> ASV2                                 5                           5
#> ASV8                                 1                           3
#> ASV38                                0                           0
#> ASV756                               0                           0
#> ASV12                                1                           1
#> ASV19                                0                           0
#> ASV175                            1295                           0
#> ASV18                                0                           0
#> ASV694                             101                           0
#> ASV113                               0                          14
#> ASV378                               0                           0
#> ASV717                               0                           0
#> ASV41                                0                           0
#> ASV65                             8505                           0
#> ASV89                                0                           0
#> ASV618                               0                           0
#> ASV170                               0                           0
#> ASV746                               0                           0
#> ASV953                               0                           0
#> ASV28                                0                           0
#> ASV310                              12                           0
#> ASV666                               0                           0
#> ASV1022                              0                           0
#> ASV1192                              0                           0
#> ASV727                               4                           0
#> ASV989                               0                           0
#> ASV635                               0                           0
#> ASV1074                              0                           0
#> ASV1108                              0                           0
#> ASV1483                              0                           0
#> ASV43                                0                           0
#> ASV46                             1230                           0
#> ASV365                             705                           0
#> ASV94                                0                           0
#> ASV602                               0                           0
#> ASV178                               0                           0
#> ASV172                               0                           0
#> ASV209                               0                           0
#> ASV741                               0                           0
#> ASV1058                              0                           0
#> ASV261                               0                           0
#> ASV329                               0                           4
#> ASV906                               0                           0
#> ASV546                               8                           0
#> ASV1253                              0                           0
#> ASV24                                0                           0
#> ASV1421                              3                           0
#> ASV816                               0                           0
#> ASV662                               0                           0
#> ASV266                               0                           0
#> ASV162                               2                           0
#> ASV139                               0                           0
#> ASV219                            1961                           0
#> ASV201                               0                           0
#> ASV251                               0                           0
#> ASV244                               0                           0
#> ASV507                               0                           0
#> ASV483                               0                           0
#> ASV500                               0                           0
#> ASV577                             292                           0
#> ASV344                               0                           0
#> ASV542                               0                           0
#> ASV85                                0                           0
#> ASV870                             152                           0
#> ASV899                             210                           0
#> ASV52                                0                           0
#> ASV926                               0                           0
#> ASV159                               0                           0
#> ASV333                               0                           0
#> ASV631                               0                           0
#> ASV348                               0                           0
#> ASV817                              13                           0
#> ASV1007                              0                           0
#> ASV543                               0                           0
#> ASV60                                0                           0
#> ASV559                               0                           0
#> ASV1370                              0                           0
#> ASV1459                              0                           0
#> ASV1359                              0                           0
#> ASV292                               0                           0
#> ASV1059                              0                           0
#> ASV443                               0                           0
#> ASV505                               0                           2
#> ASV509                               0                          25
#> ASV1624                              0                           0
#> ASV98                                0                           0
#> ASV368                               2                           0
#> ASV493                               0                           0
#> ASV1279                              0                           0
#> ASV904                               0                           0
#> ASV563                               0                           0
#> ASV1644                              0                           0
#> ASV1369                              0                           0
#> ASV1409                              0                           0
#> ASV724                               0                          12
#> ASV255                               0                           0
#> ASV69                                0                           0
#> ASV151                              13                           0
#> ASV566                               0                           0
#> ASV831                               0                           0
#> ASV1078                              0                           0
#> ASV1230                              0                           0
#> ASV1004                              0                           0
#> ASV586                               0                           0
#> ASV859                               0                           1
#> ASV1278                              0                           0
#> ASV1356                              0                           0
#> ASV801                               0                           0
#> ASV64                                0                        8187
#> ASV144                               0                           0
#> ASV131                               0                           0
#> ASV202                               0                           0
#> ASV107                               0                          50
#> ASV322                             957                           0
#> ASV444                               0                           0
#> ASV313                               0                           0
#> ASV569                               0                           0
#> ASV672                               0                           0
#> ASV462                               0                           0
#> ASV223                               0                           0
#> ASV594                               0                           0
#> ASV534                               0                           0
#> ASV950                               0                           0
#> ASV428                              26                           0
#> ASV915                               0                           0
#> ASV346                               0                           0
#> ASV1146                            132                           0
#> ASV942                               0                           0
#> ASV1164                              0                           0
#> ASV1268                              0                           0
#> ASV1152                              0                           0
#> ASV1387                              0                           0
#> ASV592                               0                           0
#> ASV737                               0                           0
#> ASV1458                              0                           0
#> ASV772                               0                           0
#> ASV523                               0                          58
#> ASV766                               0                           0
#> ASV55                                0                           0
#> ASV1236                              0                           0
#> ASV749                               0                          50
#> ASV854                               0                           0
#> ASV743                               0                           0
#> ASV975                               0                           0
#> ASV111                               0                           0
#> ASV1532                              0                           0
#> ASV1493                              0                           0
#> ASV1131                              0                           0
#> ASV1101                              0                           0
#> ASV474                               0                           0
#> ASV744                               0                           0
#> ASV1396                              0                           0
#> ASV758                               0                           0
#> ASV784                               0                           0
#> ASV297                               0                           1
#> ASV91                                0                           0
#> ASV199                               0                           0
#> ASV643                               0                           0
#> ASV787                               0                           0
#> ASV1261                              0                           0
#> ASV1434                              0                           0
#> ASV357                               0                           0
#> ASV927                               0                           0
#> ASV248                               0                           0
#> ASV963                               0                           0
#> ASV1212                              0                           0
#> ASV1321                              0                           0
#> ASV880                               0                           0
#> ASV1052                              0                           0
#> ASV1155                              0                           0
#> ASV685                               0                           0
#> ASV1128                              0                           0
#> ASV81                                0                           0
#> ASV339                               0                           0
#> ASV358                               0                           0
#> ASV420                               0                           0
#> ASV693                               0                           0
#> ASV715                               0                           0
#> ASV884                               0                           0
#> ASV973                               0                           0
#> ASV1218                              0                           0
#> ASV1319                              0                           0
#> ASV1557                              0                           0
#> ASV1315                              0                           0
#> ASV576                               0                           0
#> ASV641                               0                           0
#> ASV1198                              0                           2
#> ASV29                                0                           0
#> ASV286                               0                           0
#> ASV464                               0                           0
#> ASV466                               0                           0
#> ASV1036                              0                           0
#> ASV1247                              0                           0
#> ASV1300                              0                           0
#> ASV1427                              0                           0
#> ASV1690                              0                           0
#> ASV834                               0                           0
#> ASV1176                              0                           0
#> ASV1223                              0                           0
#> ASV1332                              0                           0
#> ASV1365                              0                           0
#> ASV78                                0                        6005
#> ASV101                               0                           0
#> ASV168                               0                        2877
#> ASV221                               0                           0
#> ASV242                               0                           0
#> ASV253                               0                           0
#> ASV59                                0                           0
#> ASV264                               0                           0
#> ASV45                              836                           0
#> ASV454                               0                           0
#> ASV489                               0                           0
#> ASV524                             531                           0
#> ASV552                               0                           0
#> ASV533                               0                           0
#> ASV619                               0                           0
#> ASV504                               0                           0
#> ASV26                                0                           0
#> ASV692                               0                           0
#> ASV47                                0                           0
#> ASV762                               0                           0
#> ASV461                               0                           0
#> ASV833                               0                           0
#> ASV338                               0                           0
#> ASV898                               0                           0
#> ASV946                               0                           0
#> ASV987                               0                         176
#> ASV998                               0                           0
#> ASV254                               0                           0
#> ASV637                               0                           0
#> ASV1070                              0                           0
#> ASV1092                              0                           0
#> ASV1102                              0                           0
#> ASV49                                0                           0
#> ASV1117                              0                           0
#> ASV1118                              0                           0
#> ASV1132                              0                         136
#> ASV964                               0                           0
#> ASV832                               0                           0
#> ASV1166                              0                           0
#> ASV1182                              0                           0
#> ASV1314                              0                           0
#> ASV969                               0                           0
#> ASV1200                              0                           0
#> ASV624                               0                           0
#> ASV1320                              0                           0
#> ASV999                               0                           0
#> ASV1485                              0                           0
#> ASV1461                              0                           0
#> ASV1313                              0                          70
#> ASV1420                              0                           0
#> ASV105                              68                           0
#> ASV431                               0                           0
#> ASV377                               0                           0
#> ASV1561                              0                           0
#> ASV477                               0                           0
#> ASV1537                              0                           0
#> ASV847                               0                          46
#> ASV722                               0                           0
#> ASV404                               0                           0
#> ASV1510                              0                           0
#> ASV1526                              0                           0
#> ASV1379                              0                           0
#> ASV367                               0                           0
#> ASV1642                              0                           0
#> ASV1233                              0                           0
#> ASV1386                              0                           0
#> ASV1234                              0                           0
#> ASV1340                             26                           0
#> ASV1035                              0                           0
#> ASV1037                              0                           0
#> ASV48                                0                           0
#> ASV711                              21                           0
#> ASV1609                              0                           0
#> ASV33                                0                          19
#> ASV1683                              0                           0
#> ASV673                               0                          15
#> ASV1548                              0                           0
#> ASV1460                              0                           0
#> ASV545                               0                           0
#> ASV1410                              0                           0
#> ASV1462                              0                           0
#> ASV853                               0                           0
#> ASV1067                             12                           0
#> ASV1311                              0                           0
#> ASV1552                              0                           0
#> ASV1588                              0                           0
#> ASV359                               0                           0
#> ASV440                               0                           0
#> ASV580                               0                          10
#> ASV858                               0                           0
#> ASV1632                              0                           0
#> ASV526                               0                           0
#> ASV1179                              9                           0
#> ASV1428                              9                           0
#> ASV355                               8                           0
#> ASV1010                              0                           0
#> ASV1542                              0                           0
#> ASV1577                              0                           8
#> ASV1262                              0                           0
#> ASV27                                0                           6
#> ASV137                               0                           0
#> ASV930                               0                           0
#> ASV962                               0                           0
#> ASV1242                              0                           0
#> ASV1628                              0                           0
#> ASV1650                              0                           0
#> ASV903                               0                           0
#> ASV984                               0                           0
#> ASV1011                              0                           0
#> ASV1224                              0                           0
#> ASV1276                              0                           0
#> ASV1468                              0                           0
#> ASV1574                              0                           0
#> ASV1630                              0                           0
#> ASV210                               0                           0
#> ASV383                               0                           0
#> ASV753                               0                           0
#> ASV996                               4                           0
#> ASV1069                              0                           0
#> ASV1082                              0                           0
#> ASV1109                              0                           0
#> ASV1144                              0                           0
#> ASV1419                              0                           0
#> ASV1467                              0                           0
#> ASV42                                0                           0
#> ASV104                               0                           0
#> ASV203                               0                           3
#> ASV491                               0                           0
#> ASV517                               0                           0
#> ASV598                               0                           0
#> ASV626                               0                           0
#> ASV650                               0                           0
#> ASV790                               0                           0
#> ASV839                               0                           0
#> ASV1014                              0                           0
#> ASV1017                              0                           0
#> ASV1030                              0                           0
#> ASV1085                              0                           0
#> ASV1167                              0                           0
#> ASV1625                              0                           3
#> ASV25                                0                           0
#> ASV249                               0                           2
#> ASV272                               0                           0
#> ASV320                               0                           0
#> ASV334                               0                           0
#> ASV488                               0                           0
#> ASV515                               0                           0
#> ASV549                               0                           0
#> ASV570                               0                           0
#> ASV604                               0                           0
#> ASV616                               0                           0
#> ASV736                               0                           0
#> ASV748                               0                           0
#> ASV760                               2                           0
#> ASV823                               0                           0
#> ASV828                               0                           0
#> ASV829                               0                           0
#> ASV891                               0                           0
#> ASV914                               0                           0
#> ASV961                               0                           0
#> ASV979                               0                           0
#> ASV986                               0                           0
#> ASV1133                              0                           0
#> ASV1270                              0                           0
#> ASV1422                              0                           0
#> ASV1523                              0                           0
#> ASV1576                              0                           0
#> ASV1603                              0                           0
#> ASV1726                              0                           0
#> ASV67                                0                           0
#> ASV87                                0                           0
#> ASV116                               0                           0
#> ASV171                               0                           0
#> ASV193                               0                           0
#> ASV198                               0                           0
#> ASV314                               0                           0
#> ASV405                               0                           0
#> ASV415                               0                           0
#> ASV731                               0                           0
#> ASV796                               0                           0
#> ASV840                               0                           0
#> ASV860                               0                           0
#> ASV875                               0                           0
#> ASV892                               0                           0
#> ASV940                               1                           0
#> ASV941                               0                           1
#> ASV959                               0                           0
#> ASV981                               0                           0
#> ASV988                               0                           0
#> ASV1015                              0                           0
#> ASV1020                              0                           0
#> ASV1027                              0                           0
#> ASV1031                              0                           0
#> ASV1066                              0                           0
#> ASV1107                              0                           0
#> ASV1111                              0                           0
#> ASV1239                              0                           0
#> ASV1263                              0                           0
#> ASV1302                              0                           0
#> ASV1303                              0                           0
#> ASV1341                              0                           0
#> ASV1380                              0                           0
#> ASV1388                              0                           0
#> ASV1484                              1                           0
#> ASV1550                              0                           0
#> ASV1562                              0                           0
#> ASV1635                              0                           0
#> ASV1653                              0                           0
#> ASV1674                              0                           1
#> ASV1712                              0                           0
#>         AC29033_S8_MERGED.fastq.gz AD26.005.B_S9_MERGED.fastq.gz
#> ASV2                            19                             8
#> ASV8                             1                             0
#> ASV38                            0                           437
#> ASV756                           0                             1
#> ASV12                            0                             0
#> ASV19                            0                             1
#> ASV175                           0                            73
#> ASV18                            0                             0
#> ASV694                           0                            16
#> ASV113                          83                             0
#> ASV378                           1                             0
#> ASV717                           3                             1
#> ASV41                            1                             0
#> ASV65                            1                             0
#> ASV89                          144                             0
#> ASV618                           0                             5
#> ASV170                           0                             0
#> ASV746                           0                             0
#> ASV953                           0                             8
#> ASV28                           78                             0
#> ASV310                          12                             2
#> ASV666                           0                             0
#> ASV1022                         14                             0
#> ASV1192                          3                             0
#> ASV727                           1                             0
#> ASV989                           0                             3
#> ASV635                           1                             9
#> ASV1074                          0                             0
#> ASV1108                          0                             0
#> ASV1483                          0                             1
#> ASV43                            0                         11042
#> ASV46                            0                             1
#> ASV365                           0                             0
#> ASV94                            0                             0
#> ASV602                           0                             0
#> ASV178                           0                             0
#> ASV172                           0                             0
#> ASV209                           0                             0
#> ASV741                           0                             0
#> ASV1058                         97                             0
#> ASV261                           0                             0
#> ASV329                          67                             0
#> ASV906                           0                             0
#> ASV546                          10                             0
#> ASV1253                          3                             0
#> ASV24                            0                             0
#> ASV1421                          1                             0
#> ASV816                           0                             6
#> ASV662                           2                             0
#> ASV266                           1                             0
#> ASV162                           0                             0
#> ASV139                           0                           141
#> ASV219                           0                             0
#> ASV201                           0                          1784
#> ASV251                           0                             0
#> ASV244                           0                            13
#> ASV507                           0                             5
#> ASV483                           0                             0
#> ASV500                           0                             0
#> ASV577                           0                             0
#> ASV344                           2                             0
#> ASV542                           0                             3
#> ASV85                            0                             0
#> ASV870                           0                             0
#> ASV899                           0                             0
#> ASV52                            0                            84
#> ASV926                           3                             0
#> ASV159                           0                            32
#> ASV333                           0                             0
#> ASV631                           0                             0
#> ASV348                           0                             0
#> ASV817                           0                             2
#> ASV1007                          0                             0
#> ASV543                           5                             0
#> ASV60                           53                             0
#> ASV559                           3                             0
#> ASV1370                          0                             0
#> ASV1459                          0                            48
#> ASV1359                          0                             1
#> ASV292                           1                             0
#> ASV1059                          0                             0
#> ASV443                           2                             0
#> ASV505                           0                             0
#> ASV509                           0                             0
#> ASV1624                          1                             0
#> ASV98                            1                             0
#> ASV368                           0                             0
#> ASV493                           0                             0
#> ASV1279                          0                            21
#> ASV904                           0                             0
#> ASV563                           0                             0
#> ASV1644                          0                             0
#> ASV1369                          6                             1
#> ASV1409                          0                             0
#> ASV724                           0                             0
#> ASV255                           0                             0
#> ASV69                            3                             0
#> ASV151                           0                             1
#> ASV566                           1                             0
#> ASV831                           0                             0
#> ASV1078                          0                             0
#> ASV1230                          0                             0
#> ASV1004                          1                             0
#> ASV586                           0                             0
#> ASV859                           0                             0
#> ASV1278                          0                             0
#> ASV1356                          0                             3
#> ASV801                           0                             0
#> ASV64                            0                             0
#> ASV144                           0                             0
#> ASV131                           0                             0
#> ASV202                           0                             5
#> ASV107                           0                             0
#> ASV322                           0                             0
#> ASV444                         629                             0
#> ASV313                          19                             0
#> ASV569                           0                             4
#> ASV672                           0                             0
#> ASV462                          12                             0
#> ASV223                           0                             0
#> ASV594                           0                             0
#> ASV534                           5                             0
#> ASV950                           0                             0
#> ASV428                           0                             0
#> ASV915                           0                           163
#> ASV346                           0                             0
#> ASV1146                          0                             0
#> ASV942                           0                             2
#> ASV1164                          0                            10
#> ASV1268                          0                             1
#> ASV1152                          0                             0
#> ASV1387                          0                            81
#> ASV592                           0                             0
#> ASV737                           2                             0
#> ASV1458                          0                             0
#> ASV772                           0                             0
#> ASV523                           0                             0
#> ASV766                           0                             0
#> ASV55                            0                            11
#> ASV1236                          0                             0
#> ASV749                           0                             0
#> ASV854                           0                             0
#> ASV743                           0                            25
#> ASV975                           0                             0
#> ASV111                           0                             0
#> ASV1532                          8                             0
#> ASV1493                          0                             0
#> ASV1131                          0                             0
#> ASV1101                          2                             0
#> ASV474                           0                             0
#> ASV744                           0                             0
#> ASV1396                          0                            12
#> ASV758                           0                             0
#> ASV784                           0                             0
#> ASV297                          13                             0
#> ASV91                            0                             0
#> ASV199                          12                             0
#> ASV643                           0                             0
#> ASV787                           0                             0
#> ASV1261                          0                             0
#> ASV1434                          0                             1
#> ASV357                           0                             0
#> ASV927                           4                             0
#> ASV248                           0                             0
#> ASV963                           0                             0
#> ASV1212                          0                             0
#> ASV1321                          0                             3
#> ASV880                           0                             2
#> ASV1052                          0                             0
#> ASV1155                          0                             0
#> ASV685                           0                             0
#> ASV1128                          0                             0
#> ASV81                            0                             0
#> ASV339                           0                             0
#> ASV358                           2                             0
#> ASV420                           0                             0
#> ASV693                           0                             0
#> ASV715                           0                             0
#> ASV884                           0                             0
#> ASV973                           0                             0
#> ASV1218                          0                             5
#> ASV1319                          0                             0
#> ASV1557                          2                             0
#> ASV1315                          0                             0
#> ASV576                           0                             0
#> ASV641                           0                             0
#> ASV1198                          0                             0
#> ASV29                            0                             0
#> ASV286                           2                             0
#> ASV464                           1                             0
#> ASV466                           0                             2
#> ASV1036                          1                             0
#> ASV1247                          0                             0
#> ASV1300                          2                             0
#> ASV1427                          0                             0
#> ASV1690                          0                             0
#> ASV834                           0                             0
#> ASV1176                          0                             0
#> ASV1223                          1                             0
#> ASV1332                          1                             0
#> ASV1365                          0                             0
#> ASV78                            0                             0
#> ASV101                           0                             0
#> ASV168                           0                             0
#> ASV221                        1927                             0
#> ASV242                           0                             0
#> ASV253                           0                             0
#> ASV59                            0                             0
#> ASV264                           0                             0
#> ASV45                            0                             0
#> ASV454                           0                             0
#> ASV489                           0                             0
#> ASV524                           0                             0
#> ASV552                           0                             0
#> ASV533                           0                           434
#> ASV619                           0                             0
#> ASV504                           0                             0
#> ASV26                            0                             0
#> ASV692                         346                             0
#> ASV47                            0                             0
#> ASV762                           0                             0
#> ASV461                           0                             0
#> ASV833                         248                             0
#> ASV338                           0                             0
#> ASV898                           0                             0
#> ASV946                           0                             0
#> ASV987                           0                             0
#> ASV998                         173                             0
#> ASV254                           0                           168
#> ASV637                           0                             0
#> ASV1070                          0                             0
#> ASV1092                          0                             0
#> ASV1102                          0                             0
#> ASV49                            0                             0
#> ASV1117                          0                           140
#> ASV1118                          0                             0
#> ASV1132                          0                             0
#> ASV964                           0                             0
#> ASV832                           0                             0
#> ASV1166                          0                             0
#> ASV1182                          0                             0
#> ASV1314                          0                             0
#> ASV969                           0                             0
#> ASV1200                          0                             0
#> ASV624                           0                             0
#> ASV1320                          0                             0
#> ASV999                           0                             0
#> ASV1485                          0                            73
#> ASV1461                          0                             0
#> ASV1313                          0                             0
#> ASV1420                          0                             0
#> ASV105                           0                             0
#> ASV431                           0                             0
#> ASV377                           0                             0
#> ASV1561                          0                             0
#> ASV477                           0                             0
#> ASV1537                          0                             0
#> ASV847                           0                             0
#> ASV722                           0                             0
#> ASV404                           0                             0
#> ASV1510                          0                             0
#> ASV1526                          0                             0
#> ASV1379                          0                            38
#> ASV367                          36                             0
#> ASV1642                          0                             0
#> ASV1233                         34                             0
#> ASV1386                          0                             0
#> ASV1234                          0                             0
#> ASV1340                          0                             0
#> ASV1035                          0                             0
#> ASV1037                         23                             0
#> ASV48                            0                             0
#> ASV711                           0                             0
#> ASV1609                          0                             0
#> ASV33                            0                             0
#> ASV1683                          0                             0
#> ASV673                           0                             0
#> ASV1548                          0                             0
#> ASV1460                          0                             0
#> ASV545                           0                             0
#> ASV1410                          0                             0
#> ASV1462                          0                             0
#> ASV853                           0                             0
#> ASV1067                          0                             0
#> ASV1311                          0                             0
#> ASV1552                          0                             0
#> ASV1588                          0                             0
#> ASV359                           0                             0
#> ASV440                           0                             0
#> ASV580                           0                             0
#> ASV858                           0                             0
#> ASV1632                          0                            10
#> ASV526                           0                             0
#> ASV1179                          0                             0
#> ASV1428                          0                             0
#> ASV355                           0                             0
#> ASV1010                          0                             0
#> ASV1542                          0                             0
#> ASV1577                          0                             0
#> ASV1262                          0                             0
#> ASV27                            0                             0
#> ASV137                           0                             0
#> ASV930                           0                             0
#> ASV962                           0                             0
#> ASV1242                          0                             0
#> ASV1628                          0                             0
#> ASV1650                          0                             0
#> ASV903                           0                             0
#> ASV984                           5                             0
#> ASV1011                          0                             0
#> ASV1224                          0                             0
#> ASV1276                          0                             0
#> ASV1468                          0                             0
#> ASV1574                          0                             0
#> ASV1630                          0                             0
#> ASV210                           0                             0
#> ASV383                           4                             0
#> ASV753                           0                             0
#> ASV996                           0                             0
#> ASV1069                          0                             0
#> ASV1082                          0                             0
#> ASV1109                          4                             0
#> ASV1144                          0                             0
#> ASV1419                          0                             0
#> ASV1467                          0                             0
#> ASV42                            0                             0
#> ASV104                           0                             0
#> ASV203                           0                             0
#> ASV491                           0                             0
#> ASV517                           0                             0
#> ASV598                           0                             0
#> ASV626                           0                             0
#> ASV650                           0                             0
#> ASV790                           0                             0
#> ASV839                           0                             0
#> ASV1014                          0                             0
#> ASV1017                          0                             0
#> ASV1030                          0                             0
#> ASV1085                          0                             0
#> ASV1167                          0                             3
#> ASV1625                          0                             0
#> ASV25                            0                             0
#> ASV249                           0                             0
#> ASV272                           0                             0
#> ASV320                           0                             0
#> ASV334                           0                             0
#> ASV488                           0                             0
#> ASV515                           0                             0
#> ASV549                           0                             2
#> ASV570                           0                             0
#> ASV604                           0                             0
#> ASV616                           0                             0
#> ASV736                           0                             0
#> ASV748                           0                             0
#> ASV760                           0                             0
#> ASV823                           0                             0
#> ASV828                           0                             0
#> ASV829                           0                             0
#> ASV891                           0                             0
#> ASV914                           0                             0
#> ASV961                           0                             0
#> ASV979                           0                             0
#> ASV986                           0                             0
#> ASV1133                          0                             0
#> ASV1270                          0                             0
#> ASV1422                          2                             0
#> ASV1523                          2                             0
#> ASV1576                          0                             0
#> ASV1603                          0                             0
#> ASV1726                          0                             2
#> ASV67                            0                             0
#> ASV87                            0                             1
#> ASV116                           0                             0
#> ASV171                           0                             0
#> ASV193                           1                             0
#> ASV198                           0                             0
#> ASV314                           0                             0
#> ASV405                           1                             0
#> ASV415                           1                             0
#> ASV731                           0                             0
#> ASV796                           0                             0
#> ASV840                           0                             0
#> ASV860                           0                             0
#> ASV875                           0                             0
#> ASV892                           0                             0
#> ASV940                           0                             0
#> ASV941                           0                             0
#> ASV959                           0                             0
#> ASV981                           0                             0
#> ASV988                           0                             0
#> ASV1015                          0                             0
#> ASV1020                          0                             0
#> ASV1027                          0                             0
#> ASV1031                          1                             0
#> ASV1066                          0                             0
#> ASV1107                          0                             1
#> ASV1111                          1                             0
#> ASV1239                          0                             0
#> ASV1263                          0                             0
#> ASV1302                          0                             1
#> ASV1303                          0                             0
#> ASV1341                          0                             0
#> ASV1380                          0                             0
#> ASV1388                          0                             0
#> ASV1484                          0                             0
#> ASV1550                          1                             0
#> ASV1562                          0                             0
#> ASV1635                          0                             0
#> ASV1653                          0                             0
#> ASV1674                          0                             0
#> ASV1712                          0                             0
#>         AD26.005.H_S10_MERGED.fastq.gz AD26.005.M_S11_MERGED.fastq.gz
#> ASV2                                 7                              5
#> ASV8                                 0                              0
#> ASV38                                0                             12
#> ASV756                               9                              4
#> ASV12                                0                           1338
#> ASV19                                0                              0
#> ASV175                               0                              0
#> ASV18                                0                              0
#> ASV694                               0                              0
#> ASV113                               0                              0
#> ASV378                               2                             10
#> ASV717                              20                              4
#> ASV41                                0                              0
#> ASV65                                0                              0
#> ASV89                               91                              0
#> ASV618                             166                              0
#> ASV170                               2                              0
#> ASV746                               0                              0
#> ASV953                              77                              0
#> ASV28                                0                              0
#> ASV310                               0                              0
#> ASV666                               0                              0
#> ASV1022                             23                              2
#> ASV1192                              3                              0
#> ASV727                               0                              0
#> ASV989                               0                              0
#> ASV635                               0                              0
#> ASV1074                              1                              0
#> ASV1108                              0                              0
#> ASV1483                              2                              0
#> ASV43                                0                              0
#> ASV46                                0                              0
#> ASV365                               0                              0
#> ASV94                                0                              0
#> ASV602                               0                             36
#> ASV178                               0                              0
#> ASV172                               3                             21
#> ASV209                               0                              0
#> ASV741                               0                              0
#> ASV1058                             40                              0
#> ASV261                               0                              0
#> ASV329                               0                              0
#> ASV906                              24                              0
#> ASV546                               0                              0
#> ASV1253                              2                              0
#> ASV24                                0                              0
#> ASV1421                              8                              7
#> ASV816                               0                              0
#> ASV662                               0                              0
#> ASV266                               0                              0
#> ASV162                               0                              0
#> ASV139                               0                              0
#> ASV219                               0                              0
#> ASV201                               0                              0
#> ASV251                               0                              0
#> ASV244                               0                              0
#> ASV507                               0                              0
#> ASV483                               0                              0
#> ASV500                               0                              0
#> ASV577                               0                              0
#> ASV344                             110                            339
#> ASV542                               0                              0
#> ASV85                                0                              0
#> ASV870                               0                              0
#> ASV899                               0                              0
#> ASV52                                0                              6
#> ASV926                             104                             58
#> ASV159                               0                             14
#> ASV333                               0                              0
#> ASV631                               0                              0
#> ASV348                               0                              0
#> ASV817                               0                             78
#> ASV1007                              0                              0
#> ASV543                               0                              0
#> ASV60                                9                              0
#> ASV559                               1                              0
#> ASV1370                              3                              3
#> ASV1459                             16                              0
#> ASV1359                              0                              0
#> ASV292                               0                             44
#> ASV1059                             38                              2
#> ASV443                               0                             34
#> ASV505                               0                              0
#> ASV509                               0                              0
#> ASV1624                              0                              0
#> ASV98                                0                             26
#> ASV368                               0                             28
#> ASV493                               0                              0
#> ASV1279                              0                              0
#> ASV904                               0                             33
#> ASV563                               0                              0
#> ASV1644                              0                             17
#> ASV1369                              0                              0
#> ASV1409                              0                              0
#> ASV724                               0                              0
#> ASV255                               0                              0
#> ASV69                                0                              0
#> ASV151                               0                              0
#> ASV566                               0                             12
#> ASV831                               0                              0
#> ASV1078                              0                              0
#> ASV1230                              0                              0
#> ASV1004                              0                              0
#> ASV586                               6                              0
#> ASV859                               0                              0
#> ASV1278                              0                              2
#> ASV1356                              0                              0
#> ASV801                               0                              2
#> ASV64                                0                              0
#> ASV144                               0                            672
#> ASV131                               0                              0
#> ASV202                               0                              0
#> ASV107                               0                           1632
#> ASV322                               0                              0
#> ASV444                               0                              0
#> ASV313                               0                              0
#> ASV569                               0                              0
#> ASV672                               0                              0
#> ASV462                               0                              0
#> ASV223                               0                              0
#> ASV594                               0                              0
#> ASV534                               0                              0
#> ASV950                               1                              0
#> ASV428                               0                              0
#> ASV915                               0                              0
#> ASV346                             155                              0
#> ASV1146                              0                              0
#> ASV942                               0                              0
#> ASV1164                            117                              0
#> ASV1268                              0                              0
#> ASV1152                              0                              0
#> ASV1387                              0                              0
#> ASV592                              74                              0
#> ASV737                               0                              0
#> ASV1458                              0                              0
#> ASV772                               0                              0
#> ASV523                               0                              0
#> ASV766                               0                              0
#> ASV55                                0                             47
#> ASV1236                              0                              0
#> ASV749                               0                              0
#> ASV854                               0                             13
#> ASV743                               0                             22
#> ASV975                               0                              0
#> ASV111                               0                              0
#> ASV1532                              0                              0
#> ASV1493                             37                              1
#> ASV1131                              0                              0
#> ASV1101                              0                              0
#> ASV474                               0                              0
#> ASV744                               0                              0
#> ASV1396                              0                              0
#> ASV758                               0                              0
#> ASV784                               0                              0
#> ASV297                               0                              0
#> ASV91                                0                              0
#> ASV199                               0                              0
#> ASV643                               0                              0
#> ASV787                               6                              0
#> ASV1261                             12                              0
#> ASV1434                              0                              0
#> ASV357                               0                              0
#> ASV927                               0                              0
#> ASV248                               0                              0
#> ASV963                               0                              0
#> ASV1212                              1                              0
#> ASV1321                              7                              0
#> ASV880                               0                              0
#> ASV1052                              0                              0
#> ASV1155                              0                              1
#> ASV685                               0                              0
#> ASV1128                              0                              0
#> ASV81                                0                              0
#> ASV339                               0                              0
#> ASV358                               0                              0
#> ASV420                               0                              0
#> ASV693                               0                              0
#> ASV715                               0                              0
#> ASV884                               0                              0
#> ASV973                               0                              0
#> ASV1218                              0                              0
#> ASV1319                              0                              0
#> ASV1557                              0                              0
#> ASV1315                              0                              0
#> ASV576                               0                              0
#> ASV641                               0                              0
#> ASV1198                              0                              0
#> ASV29                                0                              0
#> ASV286                               0                              0
#> ASV464                               0                              2
#> ASV466                               0                              0
#> ASV1036                              0                              0
#> ASV1247                              0                              0
#> ASV1300                              0                              0
#> ASV1427                              0                              0
#> ASV1690                              0                              0
#> ASV834                               0                              0
#> ASV1176                              0                              1
#> ASV1223                              0                              0
#> ASV1332                              0                              0
#> ASV1365                              0                              0
#> ASV78                                0                              0
#> ASV101                            4669                              0
#> ASV168                               0                              0
#> ASV221                               0                              0
#> ASV242                               0                              0
#> ASV253                               0                              0
#> ASV59                                0                              0
#> ASV264                               0                              0
#> ASV45                                0                              0
#> ASV454                               0                            646
#> ASV489                               0                              0
#> ASV524                               0                              0
#> ASV552                               0                              0
#> ASV533                               0                              0
#> ASV619                               0                              0
#> ASV504                               0                              0
#> ASV26                                0                              0
#> ASV692                               0                              0
#> ASV47                                0                              0
#> ASV762                               0                              0
#> ASV461                               0                              0
#> ASV833                               0                              0
#> ASV338                               0                              0
#> ASV898                               0                              0
#> ASV946                               0                              0
#> ASV987                               0                              0
#> ASV998                               0                              0
#> ASV254                               0                              0
#> ASV637                               0                              0
#> ASV1070                              0                              0
#> ASV1092                              0                              0
#> ASV1102                              0                              0
#> ASV49                                0                              0
#> ASV1117                              0                              0
#> ASV1118                              0                            140
#> ASV1132                              0                              0
#> ASV964                               0                              0
#> ASV832                               0                              0
#> ASV1166                              0                              0
#> ASV1182                              0                              0
#> ASV1314                              0                             97
#> ASV969                               0                              0
#> ASV1200                              0                              0
#> ASV624                               0                              0
#> ASV1320                              0                              0
#> ASV999                               0                              0
#> ASV1485                              0                              0
#> ASV1461                              0                              0
#> ASV1313                              0                              0
#> ASV1420                              0                              0
#> ASV105                               0                              0
#> ASV431                               0                              0
#> ASV377                               0                             60
#> ASV1561                              0                              0
#> ASV477                               0                              0
#> ASV1537                              0                              0
#> ASV847                               0                              0
#> ASV722                               0                             45
#> ASV404                               0                              0
#> ASV1510                              0                              0
#> ASV1526                              0                             39
#> ASV1379                              0                              0
#> ASV367                               0                              0
#> ASV1642                              0                              0
#> ASV1233                              0                              0
#> ASV1386                              0                              0
#> ASV1234                              0                              0
#> ASV1340                              0                              0
#> ASV1035                              0                              0
#> ASV1037                              0                              0
#> ASV48                                0                              0
#> ASV711                               0                              0
#> ASV1609                              0                              0
#> ASV33                                0                              0
#> ASV1683                              0                             18
#> ASV673                               0                              0
#> ASV1548                              0                              0
#> ASV1460                              0                              0
#> ASV545                              13                              0
#> ASV1410                              0                              0
#> ASV1462                              0                              0
#> ASV853                               0                              0
#> ASV1067                              0                              0
#> ASV1311                              0                              0
#> ASV1552                              0                              0
#> ASV1588                              0                             12
#> ASV359                               0                              0
#> ASV440                              10                              0
#> ASV580                               0                              0
#> ASV858                               0                              0
#> ASV1632                              0                              0
#> ASV526                               0                              0
#> ASV1179                              0                              0
#> ASV1428                              0                              0
#> ASV355                               0                              0
#> ASV1010                              0                              0
#> ASV1542                              0                              0
#> ASV1577                              0                              0
#> ASV1262                              0                              0
#> ASV27                                0                              0
#> ASV137                               6                              0
#> ASV930                               0                              0
#> ASV962                               0                              0
#> ASV1242                              0                              0
#> ASV1628                              0                              0
#> ASV1650                              0                              0
#> ASV903                               0                              0
#> ASV984                               0                              0
#> ASV1011                              0                              0
#> ASV1224                              0                              0
#> ASV1276                              0                              0
#> ASV1468                              0                              0
#> ASV1574                              0                              0
#> ASV1630                              0                              0
#> ASV210                               0                              0
#> ASV383                               0                              0
#> ASV753                               0                              0
#> ASV996                               0                              0
#> ASV1069                              0                              4
#> ASV1082                              0                              0
#> ASV1109                              0                              0
#> ASV1144                              0                              0
#> ASV1419                              0                              0
#> ASV1467                              0                              0
#> ASV42                                0                              0
#> ASV104                               0                              0
#> ASV203                               0                              0
#> ASV491                               0                              0
#> ASV517                               0                              0
#> ASV598                               3                              0
#> ASV626                               0                              0
#> ASV650                               0                              0
#> ASV790                               0                              0
#> ASV839                               0                              0
#> ASV1014                              0                              0
#> ASV1017                              0                              0
#> ASV1030                              0                              0
#> ASV1085                              0                              0
#> ASV1167                              0                              0
#> ASV1625                              0                              0
#> ASV25                                0                              0
#> ASV249                               0                              0
#> ASV272                               0                              0
#> ASV320                               0                              0
#> ASV334                               2                              0
#> ASV488                               0                              0
#> ASV515                               0                              0
#> ASV549                               0                              0
#> ASV570                               0                              0
#> ASV604                               0                              0
#> ASV616                               0                              0
#> ASV736                               0                              0
#> ASV748                               0                              0
#> ASV760                               0                              0
#> ASV823                               0                              0
#> ASV828                               0                              0
#> ASV829                               0                              0
#> ASV891                               0                              0
#> ASV914                               0                              0
#> ASV961                               0                              2
#> ASV979                               0                              0
#> ASV986                               0                              0
#> ASV1133                              0                              0
#> ASV1270                              0                              0
#> ASV1422                              0                              0
#> ASV1523                              0                              0
#> ASV1576                              0                              0
#> ASV1603                              0                              0
#> ASV1726                              0                              0
#> ASV67                                0                              0
#> ASV87                                0                              0
#> ASV116                               0                              0
#> ASV171                               0                              0
#> ASV193                               0                              0
#> ASV198                               0                              0
#> ASV314                               0                              0
#> ASV405                               0                              0
#> ASV415                               0                              0
#> ASV731                               0                              0
#> ASV796                               0                              0
#> ASV840                               0                              0
#> ASV860                               0                              0
#> ASV875                               0                              0
#> ASV892                               0                              0
#> ASV940                               0                              0
#> ASV941                               0                              0
#> ASV959                               0                              0
#> ASV981                               0                              0
#> ASV988                               0                              0
#> ASV1015                              0                              0
#> ASV1020                              0                              0
#> ASV1027                              0                              0
#> ASV1031                              0                              0
#> ASV1066                              0                              0
#> ASV1107                              0                              0
#> ASV1111                              0                              0
#> ASV1239                              0                              0
#> ASV1263                              0                              0
#> ASV1302                              0                              0
#> ASV1303                              0                              0
#> ASV1341                              0                              0
#> ASV1380                              0                              0
#> ASV1388                              1                              0
#> ASV1484                              0                              0
#> ASV1550                              0                              0
#> ASV1562                              0                              0
#> ASV1635                              1                              0
#> ASV1653                              0                              0
#> ASV1674                              0                              0
#> ASV1712                              0                              0
#>         AD30.ABMX.M_S12_MERGED.fastq.gz AD32.007.M_S13_MERGED.fastq.gz
#> ASV2                                 20                              5
#> ASV8                                  1                              7
#> ASV38                                 0                              0
#> ASV756                                0                              0
#> ASV12                                 0                              0
#> ASV19                                 0                             51
#> ASV175                              165                              0
#> ASV18                                 0                              0
#> ASV694                                5                              0
#> ASV113                                2                              0
#> ASV378                                0                              1
#> ASV717                                3                              0
#> ASV41                                 0                          10572
#> ASV65                                18                              0
#> ASV89                                 8                              0
#> ASV618                                0                              0
#> ASV170                                0                              0
#> ASV746                                0                              2
#> ASV953                                0                              0
#> ASV28                                 0                              1
#> ASV310                                0                             47
#> ASV666                                0                              0
#> ASV1022                               2                              0
#> ASV1192                               0                              5
#> ASV727                                5                              0
#> ASV989                                0                              0
#> ASV635                                1                              0
#> ASV1074                               0                              1
#> ASV1108                               0                              1
#> ASV1483                               0                              0
#> ASV43                                 1                              1
#> ASV46                               589                              0
#> ASV365                               14                              0
#> ASV94                                 0                              0
#> ASV602                                0                              0
#> ASV178                                0                             10
#> ASV172                                5                              0
#> ASV209                                0                              0
#> ASV741                                0                              1
#> ASV1058                               1                              0
#> ASV261                                0                              0
#> ASV329                                0                              0
#> ASV906                                0                              0
#> ASV546                                8                              0
#> ASV1253                               8                              0
#> ASV24                                 0                              0
#> ASV1421                               0                              0
#> ASV816                                0                              0
#> ASV662                                0                              0
#> ASV266                                0                              0
#> ASV162                              985                              0
#> ASV139                                0                             26
#> ASV219                                0                              0
#> ASV201                                1                              0
#> ASV251                                0                            379
#> ASV244                                0                             31
#> ASV507                                0                            406
#> ASV483                                0                              0
#> ASV500                                0                              0
#> ASV577                                0                             10
#> ASV344                                0                              0
#> ASV542                                0                            367
#> ASV85                                 0                              0
#> ASV870                                0                              9
#> ASV899                                0                              0
#> ASV52                                 0                             75
#> ASV926                                0                              0
#> ASV159                                0                            114
#> ASV333                                0                              1
#> ASV631                                0                              0
#> ASV348                                0                              0
#> ASV817                                0                              0
#> ASV1007                               0                              4
#> ASV543                                0                             78
#> ASV60                                 0                              0
#> ASV559                                0                              0
#> ASV1370                               0                              0
#> ASV1459                               0                              0
#> ASV1359                               0                              9
#> ASV292                                0                              0
#> ASV1059                               0                              0
#> ASV443                                0                              0
#> ASV505                                0                              0
#> ASV509                                0                              0
#> ASV1624                               0                              0
#> ASV98                                 0                              0
#> ASV368                                0                              0
#> ASV493                                0                              0
#> ASV1279                               0                              0
#> ASV904                                0                              0
#> ASV563                                0                              0
#> ASV1644                               0                              0
#> ASV1369                               0                              0
#> ASV1409                               0                             14
#> ASV724                                0                              0
#> ASV255                                2                             12
#> ASV69                                 0                              0
#> ASV151                                1                              0
#> ASV566                                0                              0
#> ASV831                                0                              0
#> ASV1078                               0                              1
#> ASV1230                               0                              0
#> ASV1004                               0                              0
#> ASV586                                0                              0
#> ASV859                                0                              0
#> ASV1278                               0                              2
#> ASV1356                               0                              0
#> ASV801                                0                              0
#> ASV64                                 0                              0
#> ASV144                                0                              0
#> ASV131                                0                              0
#> ASV202                                0                              0
#> ASV107                                0                              0
#> ASV322                                0                              0
#> ASV444                                0                              0
#> ASV313                                0                              0
#> ASV569                                0                              0
#> ASV672                                0                              0
#> ASV462                                0                              0
#> ASV223                                0                              0
#> ASV594                                0                              1
#> ASV534                                0                              0
#> ASV950                                0                              0
#> ASV428                                0                            146
#> ASV915                                0                              0
#> ASV346                                0                              0
#> ASV1146                               0                              0
#> ASV942                                0                              0
#> ASV1164                               0                              0
#> ASV1268                               0                              0
#> ASV1152                               0                              0
#> ASV1387                               0                              0
#> ASV592                                0                              0
#> ASV737                                0                              0
#> ASV1458                               0                              0
#> ASV772                                0                              0
#> ASV523                                0                              0
#> ASV766                                0                              6
#> ASV55                                 0                              0
#> ASV1236                               0                              0
#> ASV749                                0                              0
#> ASV854                                0                             35
#> ASV743                                0                              0
#> ASV975                                0                              0
#> ASV111                                0                              0
#> ASV1532                               0                              0
#> ASV1493                               0                              0
#> ASV1131                               0                              1
#> ASV1101                               0                              0
#> ASV474                                0                              0
#> ASV744                                0                             19
#> ASV1396                               0                              0
#> ASV758                                0                              0
#> ASV784                                0                              5
#> ASV297                                0                              0
#> ASV91                                 0                              0
#> ASV199                                0                              1
#> ASV643                                0                              0
#> ASV787                                0                              0
#> ASV1261                               0                              0
#> ASV1434                               0                              0
#> ASV357                                1                             11
#> ASV927                                0                              0
#> ASV248                                0                              1
#> ASV963                                1                              0
#> ASV1212                               0                              9
#> ASV1321                               0                              0
#> ASV880                                0                              7
#> ASV1052                               0                              0
#> ASV1155                               0                              0
#> ASV685                                6                              0
#> ASV1128                               0                              0
#> ASV81                                 2                              0
#> ASV339                                1                              0
#> ASV358                                0                              0
#> ASV420                                0                              0
#> ASV693                                0                              5
#> ASV715                                0                              0
#> ASV884                                0                              0
#> ASV973                                0                              0
#> ASV1218                               0                              0
#> ASV1319                               0                              0
#> ASV1557                               4                              0
#> ASV1315                               0                              3
#> ASV576                                0                              0
#> ASV641                                0                              0
#> ASV1198                               0                              0
#> ASV29                                 2                              0
#> ASV286                                1                              0
#> ASV464                                0                              0
#> ASV466                                0                              0
#> ASV1036                               0                              0
#> ASV1247                               0                              0
#> ASV1300                               0                              0
#> ASV1427                               0                              0
#> ASV1690                               0                              0
#> ASV834                                1                              0
#> ASV1176                               0                              0
#> ASV1223                               0                              0
#> ASV1332                               0                              0
#> ASV1365                               0                              0
#> ASV78                                 0                              0
#> ASV101                                0                              0
#> ASV168                                0                              0
#> ASV221                                0                              0
#> ASV242                                0                              0
#> ASV253                                0                              0
#> ASV59                                 0                              0
#> ASV264                                0                              0
#> ASV45                                 0                              0
#> ASV454                                0                              0
#> ASV489                                0                            583
#> ASV524                                0                              0
#> ASV552                                0                              0
#> ASV533                                0                              0
#> ASV619                                0                              0
#> ASV504                                0                              0
#> ASV26                                 0                              0
#> ASV692                                0                              0
#> ASV47                                 0                              0
#> ASV762                                0                            293
#> ASV461                                0                              0
#> ASV833                                0                              0
#> ASV338                                0                              0
#> ASV898                                0                              0
#> ASV946                                0                              0
#> ASV987                                0                              0
#> ASV998                                0                              0
#> ASV254                                0                              0
#> ASV637                                0                              0
#> ASV1070                               0                              0
#> ASV1092                               0                              0
#> ASV1102                               0                              0
#> ASV49                                 0                            141
#> ASV1117                               0                              0
#> ASV1118                               0                              0
#> ASV1132                               0                              0
#> ASV964                                0                              0
#> ASV832                                0                              0
#> ASV1166                               0                              0
#> ASV1182                               0                              0
#> ASV1314                               0                              0
#> ASV969                                0                              0
#> ASV1200                               0                             92
#> ASV624                                0                              0
#> ASV1320                               0                              0
#> ASV999                                0                              0
#> ASV1485                               0                              0
#> ASV1461                               0                              0
#> ASV1313                               0                              0
#> ASV1420                               0                              0
#> ASV105                                0                              0
#> ASV431                                0                              0
#> ASV377                                0                              0
#> ASV1561                              58                              0
#> ASV477                                0                              0
#> ASV1537                               0                              0
#> ASV847                                0                              0
#> ASV722                                0                              0
#> ASV404                                0                             39
#> ASV1510                               0                              0
#> ASV1526                               0                              0
#> ASV1379                               0                              0
#> ASV367                                0                              0
#> ASV1642                               0                              0
#> ASV1233                               0                              0
#> ASV1386                               0                              0
#> ASV1234                               0                              0
#> ASV1340                               0                              0
#> ASV1035                               0                              0
#> ASV1037                               0                              0
#> ASV48                                 0                              0
#> ASV711                                0                              0
#> ASV1609                               0                              0
#> ASV33                                 0                              0
#> ASV1683                               0                              0
#> ASV673                                0                              0
#> ASV1548                               0                              0
#> ASV1460                               0                              0
#> ASV545                                0                              0
#> ASV1410                               0                              0
#> ASV1462                               0                              0
#> ASV853                                0                              0
#> ASV1067                               0                              0
#> ASV1311                               0                              0
#> ASV1552                               0                              0
#> ASV1588                               0                              0
#> ASV359                                0                              0
#> ASV440                                0                              0
#> ASV580                                0                              0
#> ASV858                                0                              0
#> ASV1632                               0                              0
#> ASV526                                0                              0
#> ASV1179                               0                              0
#> ASV1428                               0                              0
#> ASV355                                0                              0
#> ASV1010                               0                              0
#> ASV1542                               0                              0
#> ASV1577                               0                              0
#> ASV1262                               0                              0
#> ASV27                                 0                              0
#> ASV137                                0                              0
#> ASV930                                0                              0
#> ASV962                                0                              0
#> ASV1242                               0                              0
#> ASV1628                               0                              0
#> ASV1650                               0                              0
#> ASV903                                0                              0
#> ASV984                                0                              0
#> ASV1011                               0                              0
#> ASV1224                               0                              5
#> ASV1276                               0                              5
#> ASV1468                               0                              0
#> ASV1574                               0                              0
#> ASV1630                               0                              0
#> ASV210                                0                              0
#> ASV383                                0                              0
#> ASV753                                0                              0
#> ASV996                                0                              0
#> ASV1069                               0                              0
#> ASV1082                               0                              0
#> ASV1109                               0                              0
#> ASV1144                               0                              0
#> ASV1419                               0                              0
#> ASV1467                               0                              0
#> ASV42                                 0                              0
#> ASV104                                0                              0
#> ASV203                                0                              0
#> ASV491                                0                              0
#> ASV517                                0                              0
#> ASV598                                0                              0
#> ASV626                                0                              3
#> ASV650                                0                              0
#> ASV790                                0                              0
#> ASV839                                3                              0
#> ASV1014                               0                              0
#> ASV1017                               0                              0
#> ASV1030                               0                              0
#> ASV1085                               0                              0
#> ASV1167                               0                              0
#> ASV1625                               0                              0
#> ASV25                                 0                              0
#> ASV249                                0                              0
#> ASV272                                2                              0
#> ASV320                                0                              0
#> ASV334                                0                              0
#> ASV488                                0                              2
#> ASV515                                0                              0
#> ASV549                                0                              0
#> ASV570                                0                              0
#> ASV604                                0                              0
#> ASV616                                0                              0
#> ASV736                                0                              0
#> ASV748                                0                              0
#> ASV760                                0                              0
#> ASV823                                0                              0
#> ASV828                                0                              0
#> ASV829                                0                              0
#> ASV891                                0                              0
#> ASV914                                0                              0
#> ASV961                                0                              0
#> ASV979                                0                              0
#> ASV986                                0                              0
#> ASV1133                               0                              0
#> ASV1270                               2                              0
#> ASV1422                               0                              0
#> ASV1523                               0                              0
#> ASV1576                               0                              0
#> ASV1603                               0                              0
#> ASV1726                               0                              0
#> ASV67                                 0                              0
#> ASV87                                 0                              0
#> ASV116                                0                              0
#> ASV171                                0                              0
#> ASV193                                0                              0
#> ASV198                                0                              0
#> ASV314                                0                              0
#> ASV405                                0                              0
#> ASV415                                0                              0
#> ASV731                                0                              0
#> ASV796                                0                              0
#> ASV840                                0                              1
#> ASV860                                0                              0
#> ASV875                                0                              0
#> ASV892                                0                              0
#> ASV940                                0                              0
#> ASV941                                0                              0
#> ASV959                                0                              0
#> ASV981                                0                              0
#> ASV988                                0                              0
#> ASV1015                               0                              0
#> ASV1020                               0                              0
#> ASV1027                               0                              0
#> ASV1031                               0                              0
#> ASV1066                               0                              0
#> ASV1107                               0                              0
#> ASV1111                               0                              0
#> ASV1239                               0                              0
#> ASV1263                               0                              0
#> ASV1302                               0                              0
#> ASV1303                               0                              0
#> ASV1341                               0                              0
#> ASV1380                               0                              0
#> ASV1388                               0                              0
#> ASV1484                               0                              0
#> ASV1550                               0                              0
#> ASV1562                               1                              0
#> ASV1635                               0                              0
#> ASV1653                               0                              0
#> ASV1674                               0                              0
#> ASV1712                               0                              0
#>         ADABM30X.B_S14_MERGED.fastq.gz ADABM30X.H_S15_MERGED.fastq.gz
#> ASV2                                 3                              3
#> ASV8                                 2                              2
#> ASV38                               19                            273
#> ASV756                               0                              1
#> ASV12                                0                              0
#> ASV19                                5                              4
#> ASV175                             383                             15
#> ASV18                                0                              0
#> ASV694                              56                              5
#> ASV113                               0                             11
#> ASV378                               1                              0
#> ASV717                               0                              0
#> ASV41                                0                              0
#> ASV65                                7                             15
#> ASV89                                0                              0
#> ASV618                               0                              0
#> ASV170                               0                              0
#> ASV746                               9                            156
#> ASV953                               0                              0
#> ASV28                                1                              0
#> ASV310                               0                              9
#> ASV666                               0                              0
#> ASV1022                              0                              0
#> ASV1192                              0                              0
#> ASV727                               1                              7
#> ASV989                               3                              6
#> ASV635                               3                              1
#> ASV1074                              0                              0
#> ASV1108                              3                              2
#> ASV1483                              0                              1
#> ASV43                                0                              0
#> ASV46                             1238                              0
#> ASV365                               7                              4
#> ASV94                                0                              0
#> ASV602                               0                              1
#> ASV178                               4                              2
#> ASV172                               0                              0
#> ASV209                               0                              0
#> ASV741                               0                              9
#> ASV1058                              0                              0
#> ASV261                               0                              0
#> ASV329                               0                              0
#> ASV906                               0                              0
#> ASV546                               0                              3
#> ASV1253                              0                              0
#> ASV24                                8                              0
#> ASV1421                              0                              0
#> ASV816                               1                              0
#> ASV662                               0                              0
#> ASV266                               3                              0
#> ASV162                            2068                              0
#> ASV139                               0                              0
#> ASV219                               1                              1
#> ASV201                               0                              0
#> ASV251                               0                              0
#> ASV244                               0                              0
#> ASV507                               0                              0
#> ASV483                              25                              0
#> ASV500                               0                              0
#> ASV577                               0                              0
#> ASV344                               0                              0
#> ASV542                               0                              0
#> ASV85                               19                              0
#> ASV870                               0                              0
#> ASV899                               5                              0
#> ASV52                                0                              0
#> ASV926                               0                              0
#> ASV159                               0                              0
#> ASV333                               0                              1
#> ASV631                               0                              0
#> ASV348                               0                              0
#> ASV817                               0                              0
#> ASV1007                              0                              0
#> ASV543                               0                              0
#> ASV60                                0                              0
#> ASV559                               0                              0
#> ASV1370                              0                              0
#> ASV1459                              0                              0
#> ASV1359                              0                              0
#> ASV292                               0                             10
#> ASV1059                              0                              0
#> ASV443                               0                             13
#> ASV505                               0                              0
#> ASV509                               0                              0
#> ASV1624                              0                              0
#> ASV98                                0                             12
#> ASV368                               0                              8
#> ASV493                               0                              0
#> ASV1279                              2                             13
#> ASV904                               0                              0
#> ASV563                               0                              0
#> ASV1644                              0                              0
#> ASV1369                              0                              0
#> ASV1409                              0                              0
#> ASV724                               0                              0
#> ASV255                               0                              0
#> ASV69                                0                             10
#> ASV151                               0                              0
#> ASV566                               0                              2
#> ASV831                               0                              0
#> ASV1078                              0                              1
#> ASV1230                              0                              8
#> ASV1004                              0                              0
#> ASV586                               0                              0
#> ASV859                               0                              0
#> ASV1278                              0                              0
#> ASV1356                              1                              1
#> ASV801                               1                              0
#> ASV64                                0                              0
#> ASV144                               0                           2924
#> ASV131                             235                           2751
#> ASV202                               0                              0
#> ASV107                               0                              0
#> ASV322                              41                              0
#> ASV444                               2                              0
#> ASV313                               0                              0
#> ASV569                               0                              0
#> ASV672                               0                              0
#> ASV462                               0                              0
#> ASV223                               0                            298
#> ASV594                               0                              0
#> ASV534                               0                              0
#> ASV950                               0                              0
#> ASV428                               0                              0
#> ASV915                               0                              0
#> ASV346                               0                              0
#> ASV1146                              1                              0
#> ASV942                               0                              0
#> ASV1164                              0                              0
#> ASV1268                              0                              0
#> ASV1152                              0                              0
#> ASV1387                              0                              0
#> ASV592                               0                              0
#> ASV737                               0                              0
#> ASV1458                              0                              0
#> ASV772                               0                              3
#> ASV523                               0                              0
#> ASV766                               0                              0
#> ASV55                                0                              0
#> ASV1236                              0                              0
#> ASV749                               0                              0
#> ASV854                               0                              0
#> ASV743                               0                              0
#> ASV975                               0                              2
#> ASV111                              37                              0
#> ASV1532                              0                              0
#> ASV1493                              0                              0
#> ASV1131                              0                              0
#> ASV1101                              0                              0
#> ASV474                               0                              0
#> ASV744                               0                              0
#> ASV1396                              0                              0
#> ASV758                               0                              0
#> ASV784                               0                              0
#> ASV297                               0                              0
#> ASV91                                1                              0
#> ASV199                               0                              0
#> ASV643                               0                              3
#> ASV787                               0                              0
#> ASV1261                              0                              0
#> ASV1434                              0                              0
#> ASV357                               0                              0
#> ASV927                               0                              0
#> ASV248                               0                              0
#> ASV963                               0                              0
#> ASV1212                              0                              0
#> ASV1321                              0                              0
#> ASV880                               0                              0
#> ASV1052                              0                              0
#> ASV1155                              0                              0
#> ASV685                               0                              0
#> ASV1128                              0                              0
#> ASV81                                0                              0
#> ASV339                               0                              0
#> ASV358                               0                              4
#> ASV420                               0                              0
#> ASV693                               0                              0
#> ASV715                               0                              0
#> ASV884                               0                              0
#> ASV973                               1                              0
#> ASV1218                              0                              0
#> ASV1319                              0                              0
#> ASV1557                              0                              0
#> ASV1315                              0                              0
#> ASV576                               0                              0
#> ASV641                               0                              0
#> ASV1198                              0                              0
#> ASV29                                0                              1
#> ASV286                               0                              0
#> ASV464                               0                              0
#> ASV466                               0                              1
#> ASV1036                              0                              0
#> ASV1247                              0                              0
#> ASV1300                              0                              0
#> ASV1427                              0                              0
#> ASV1690                              0                              2
#> ASV834                               0                              0
#> ASV1176                              0                              1
#> ASV1223                              0                              0
#> ASV1332                              0                              1
#> ASV1365                              0                              0
#> ASV78                                0                              0
#> ASV101                               0                              0
#> ASV168                               0                              0
#> ASV221                               0                              0
#> ASV242                               0                              0
#> ASV253                            1594                              0
#> ASV59                                0                              0
#> ASV264                               0                           1013
#> ASV45                                0                              0
#> ASV454                               0                              0
#> ASV489                               0                              0
#> ASV524                               0                              0
#> ASV552                               0                              0
#> ASV533                               0                              0
#> ASV619                               0                              0
#> ASV504                               0                              0
#> ASV26                                0                              0
#> ASV692                               0                              0
#> ASV47                                0                              0
#> ASV762                               0                              0
#> ASV461                               0                              0
#> ASV833                               0                              0
#> ASV338                               0                              0
#> ASV898                               0                              0
#> ASV946                               0                              0
#> ASV987                               0                              0
#> ASV998                               0                              0
#> ASV254                               0                              0
#> ASV637                               0                              0
#> ASV1070                              0                              0
#> ASV1092                              0                              0
#> ASV1102                              0                              0
#> ASV49                                0                              0
#> ASV1117                              0                              0
#> ASV1118                              0                              0
#> ASV1132                              0                              0
#> ASV964                               0                              0
#> ASV832                               0                              0
#> ASV1166                              0                              0
#> ASV1182                              0                              0
#> ASV1314                              0                              0
#> ASV969                               0                              0
#> ASV1200                              0                              0
#> ASV624                               0                              0
#> ASV1320                              0                              0
#> ASV999                               0                              0
#> ASV1485                              0                              0
#> ASV1461                              0                             72
#> ASV1313                              0                              0
#> ASV1420                              0                              0
#> ASV105                               0                              0
#> ASV431                               0                              0
#> ASV377                               0                              0
#> ASV1561                              0                              0
#> ASV477                               0                              0
#> ASV1537                              0                              0
#> ASV847                               0                              0
#> ASV722                               0                              0
#> ASV404                               0                              0
#> ASV1510                              0                              0
#> ASV1526                              0                              0
#> ASV1379                              0                              0
#> ASV367                               0                              0
#> ASV1642                              0                              0
#> ASV1233                              0                              0
#> ASV1386                              0                              0
#> ASV1234                              0                             26
#> ASV1340                              0                              0
#> ASV1035                              0                              0
#> ASV1037                              0                              0
#> ASV48                                0                             22
#> ASV711                               0                              0
#> ASV1609                              0                             20
#> ASV33                                0                              0
#> ASV1683                              0                              0
#> ASV673                               0                              0
#> ASV1548                             15                              0
#> ASV1460                              0                             14
#> ASV545                               0                              0
#> ASV1410                              0                              0
#> ASV1462                              0                              0
#> ASV853                               0                              0
#> ASV1067                              0                              0
#> ASV1311                              0                              0
#> ASV1552                              0                             12
#> ASV1588                              0                              0
#> ASV359                               0                              0
#> ASV440                               0                              0
#> ASV580                               0                              0
#> ASV858                               0                              0
#> ASV1632                              0                              0
#> ASV526                               0                              9
#> ASV1179                              0                              0
#> ASV1428                              0                              0
#> ASV355                               0                              0
#> ASV1010                              0                              0
#> ASV1542                              0                              0
#> ASV1577                              0                              0
#> ASV1262                              0                              0
#> ASV27                                0                              0
#> ASV137                               0                              0
#> ASV930                               0                              0
#> ASV962                               0                              0
#> ASV1242                              0                              0
#> ASV1628                              0                              0
#> ASV1650                              0                              6
#> ASV903                               0                              0
#> ASV984                               0                              0
#> ASV1011                              0                              0
#> ASV1224                              0                              0
#> ASV1276                              0                              0
#> ASV1468                              0                              0
#> ASV1574                              0                              0
#> ASV1630                              0                              0
#> ASV210                               0                              0
#> ASV383                               0                              0
#> ASV753                               0                              0
#> ASV996                               0                              0
#> ASV1069                              0                              0
#> ASV1082                              0                              0
#> ASV1109                              0                              0
#> ASV1144                              0                              0
#> ASV1419                              0                              0
#> ASV1467                              0                              0
#> ASV42                                0                              0
#> ASV104                               3                              0
#> ASV203                               0                              0
#> ASV491                               0                              0
#> ASV517                               0                              0
#> ASV598                               0                              0
#> ASV626                               0                              0
#> ASV650                               0                              0
#> ASV790                               0                              0
#> ASV839                               0                              0
#> ASV1014                              0                              0
#> ASV1017                              3                              0
#> ASV1030                              0                              0
#> ASV1085                              0                              0
#> ASV1167                              0                              0
#> ASV1625                              0                              0
#> ASV25                                0                              2
#> ASV249                               0                              0
#> ASV272                               0                              0
#> ASV320                               0                              0
#> ASV334                               0                              0
#> ASV488                               0                              0
#> ASV515                               0                              0
#> ASV549                               0                              0
#> ASV570                               0                              0
#> ASV604                               0                              0
#> ASV616                               0                              2
#> ASV736                               0                              0
#> ASV748                               0                              0
#> ASV760                               0                              0
#> ASV823                               0                              0
#> ASV828                               0                              0
#> ASV829                               0                              0
#> ASV891                               0                              0
#> ASV914                               0                              0
#> ASV961                               0                              0
#> ASV979                               0                              0
#> ASV986                               0                              0
#> ASV1133                              0                              2
#> ASV1270                              0                              0
#> ASV1422                              0                              0
#> ASV1523                              0                              0
#> ASV1576                              0                              0
#> ASV1603                              0                              0
#> ASV1726                              0                              0
#> ASV67                                0                              0
#> ASV87                                0                              0
#> ASV116                               0                              0
#> ASV171                               1                              0
#> ASV193                               0                              0
#> ASV198                               0                              0
#> ASV314                               0                              0
#> ASV405                               0                              0
#> ASV415                               0                              0
#> ASV731                               0                              0
#> ASV796                               0                              0
#> ASV840                               0                              0
#> ASV860                               0                              0
#> ASV875                               0                              0
#> ASV892                               0                              0
#> ASV940                               0                              0
#> ASV941                               0                              0
#> ASV959                               0                              0
#> ASV981                               0                              0
#> ASV988                               0                              0
#> ASV1015                              0                              0
#> ASV1020                              0                              0
#> ASV1027                              0                              0
#> ASV1031                              0                              0
#> ASV1066                              0                              0
#> ASV1107                              0                              0
#> ASV1111                              0                              0
#> ASV1239                              0                              0
#> ASV1263                              0                              0
#> ASV1302                              0                              0
#> ASV1303                              0                              1
#> ASV1341                              0                              1
#> ASV1380                              0                              1
#> ASV1388                              0                              0
#> ASV1484                              0                              0
#> ASV1550                              0                              0
#> ASV1562                              0                              0
#> ASV1635                              0                              0
#> ASV1653                              1                              0
#> ASV1674                              0                              0
#> ASV1712                              0                              1
#>         ADABM30X.M_S16_MERGED.fastq.gz AE30.ABM507_S17_MERGED.fastq.gz
#> ASV2                                 3                               7
#> ASV8                                 0                              36
#> ASV38                               37                               0
#> ASV756                               0                               0
#> ASV12                                0                               0
#> ASV19                                0                               1
#> ASV175                               2                               0
#> ASV18                                0                              44
#> ASV694                               2                               0
#> ASV113                               0                               0
#> ASV378                               0                              31
#> ASV717                               0                               3
#> ASV41                                0                               0
#> ASV65                                0                               0
#> ASV89                                0                               7
#> ASV618                               0                               0
#> ASV170                               0                               0
#> ASV746                              49                               0
#> ASV953                               0                               0
#> ASV28                                0                               1
#> ASV310                               0                               0
#> ASV666                               0                               1
#> ASV1022                              0                               2
#> ASV1192                              0                               4
#> ASV727                               0                               0
#> ASV989                               1                               0
#> ASV635                               0                               0
#> ASV1074                              0                               7
#> ASV1108                              0                               0
#> ASV1483                              0                               0
#> ASV43                                0                              31
#> ASV46                                0                               0
#> ASV365                               0                               0
#> ASV94                                0                               0
#> ASV602                               0                             327
#> ASV178                               0                               0
#> ASV172                               0                             255
#> ASV209                               0                               0
#> ASV741                               0                             210
#> ASV1058                              0                               1
#> ASV261                               0                               0
#> ASV329                               0                               0
#> ASV906                               0                               2
#> ASV546                               0                               0
#> ASV1253                              0                               9
#> ASV24                                1                               0
#> ASV1421                              0                               0
#> ASV816                               1                               0
#> ASV662                               0                               0
#> ASV266                               1                               0
#> ASV162                               0                               0
#> ASV139                               0                            1889
#> ASV219                               0                               0
#> ASV201                               0                               6
#> ASV251                               0                               0
#> ASV244                               0                             577
#> ASV507                               0                             109
#> ASV483                             406                               0
#> ASV500                               0                               0
#> ASV577                               0                               0
#> ASV344                               0                               0
#> ASV542                               0                              58
#> ASV85                              297                               0
#> ASV870                               0                               0
#> ASV899                               2                               0
#> ASV52                                0                               0
#> ASV926                               0                               0
#> ASV159                               0                               0
#> ASV333                               0                               0
#> ASV631                               0                               0
#> ASV348                               0                               0
#> ASV817                               0                               0
#> ASV1007                              0                              19
#> ASV543                               3                               0
#> ASV60                                0                               0
#> ASV559                               0                              65
#> ASV1370                              0                              62
#> ASV1459                              0                               3
#> ASV1359                              0                               0
#> ASV292                               0                               0
#> ASV1059                              0                              14
#> ASV443                               0                               0
#> ASV505                               0                               0
#> ASV509                               0                               0
#> ASV1624                              0                               0
#> ASV98                                0                               0
#> ASV368                               0                               0
#> ASV493                               0                               0
#> ASV1279                              0                               0
#> ASV904                               0                               0
#> ASV563                               0                               0
#> ASV1644                              0                               0
#> ASV1369                              0                              18
#> ASV1409                              0                               2
#> ASV724                               0                               0
#> ASV255                               3                               0
#> ASV69                                0                               0
#> ASV151                               0                               0
#> ASV566                               0                               0
#> ASV831                               0                               0
#> ASV1078                              0                               0
#> ASV1230                              0                               0
#> ASV1004                              0                               0
#> ASV586                               0                               0
#> ASV859                               0                               0
#> ASV1278                              0                               0
#> ASV1356                              0                               0
#> ASV801                               0                               0
#> ASV64                                0                               0
#> ASV144                               0                               0
#> ASV131                               0                               0
#> ASV202                               0                               0
#> ASV107                               0                               0
#> ASV322                               0                               0
#> ASV444                               0                               0
#> ASV313                               0                               0
#> ASV569                               0                             357
#> ASV672                               0                               0
#> ASV462                               0                               0
#> ASV223                               2                               0
#> ASV594                               0                               0
#> ASV534                               0                               0
#> ASV950                               0                             174
#> ASV428                               0                               0
#> ASV915                               0                               4
#> ASV346                               0                               0
#> ASV1146                              0                               0
#> ASV942                               0                               0
#> ASV1164                              0                               0
#> ASV1268                              0                              89
#> ASV1152                              0                               0
#> ASV1387                              0                               1
#> ASV592                               0                               0
#> ASV737                               0                              75
#> ASV1458                              0                               0
#> ASV772                               0                               0
#> ASV523                               0                               0
#> ASV766                               0                              54
#> ASV55                                0                               0
#> ASV1236                              0                               0
#> ASV749                               0                               0
#> ASV854                               0                               0
#> ASV743                               0                               0
#> ASV975                               0                               0
#> ASV111                               0                               0
#> ASV1532                              0                               0
#> ASV1493                              0                               0
#> ASV1131                              0                               0
#> ASV1101                              0                               0
#> ASV474                               0                               0
#> ASV744                               0                               1
#> ASV1396                              0                               0
#> ASV758                               0                               0
#> ASV784                               0                              10
#> ASV297                               0                               0
#> ASV91                                0                               0
#> ASV199                               0                               0
#> ASV643                               0                               0
#> ASV787                               0                               0
#> ASV1261                              0                               0
#> ASV1434                              0                               0
#> ASV357                               0                               0
#> ASV927                               0                               8
#> ASV248                               0                               0
#> ASV963                               0                               9
#> ASV1212                              0                               0
#> ASV1321                              0                               0
#> ASV880                               0                               0
#> ASV1052                              0                               0
#> ASV1155                              0                               0
#> ASV685                               1                               0
#> ASV1128                              0                               0
#> ASV81                                0                               4
#> ASV339                               0                               5
#> ASV358                               0                               0
#> ASV420                               0                               0
#> ASV693                               0                               0
#> ASV715                               0                               3
#> ASV884                               0                               0
#> ASV973                               5                               0
#> ASV1218                              0                               0
#> ASV1319                              0                               5
#> ASV1557                              0                               0
#> ASV1315                              0                               2
#> ASV576                               0                               0
#> ASV641                               0                               0
#> ASV1198                              0                               0
#> ASV29                                0                               0
#> ASV286                               0                               0
#> ASV464                               0                               0
#> ASV466                               0                               0
#> ASV1036                              0                               0
#> ASV1247                              0                               0
#> ASV1300                              0                               0
#> ASV1427                              0                               0
#> ASV1690                              0                               0
#> ASV834                               0                               1
#> ASV1176                              0                               0
#> ASV1223                              0                               0
#> ASV1332                              0                               0
#> ASV1365                              0                               0
#> ASV78                                0                               0
#> ASV101                               0                               0
#> ASV168                               0                               0
#> ASV221                               0                               0
#> ASV242                               0                               0
#> ASV253                               0                               0
#> ASV59                                0                            1207
#> ASV264                               0                               0
#> ASV45                                0                               0
#> ASV454                               0                               0
#> ASV489                               0                               0
#> ASV524                               0                               0
#> ASV552                               0                               0
#> ASV533                               0                               0
#> ASV619                               0                             419
#> ASV504                               0                               0
#> ASV26                                0                               0
#> ASV692                               0                               0
#> ASV47                                0                             315
#> ASV762                               0                               0
#> ASV461                               0                             258
#> ASV833                               0                               0
#> ASV338                             213                               0
#> ASV898                               0                               0
#> ASV946                               0                               0
#> ASV987                               0                               0
#> ASV998                               0                               0
#> ASV254                               0                               0
#> ASV637                               0                               0
#> ASV1070                              0                               0
#> ASV1092                              0                               0
#> ASV1102                              0                             144
#> ASV49                                0                               0
#> ASV1117                              0                               0
#> ASV1118                              0                               0
#> ASV1132                              0                               0
#> ASV964                             133                               0
#> ASV832                               0                               0
#> ASV1166                              0                               0
#> ASV1182                              0                               0
#> ASV1314                              0                               0
#> ASV969                               0                               0
#> ASV1200                              0                               0
#> ASV624                               0                               0
#> ASV1320                              0                               0
#> ASV999                               0                              84
#> ASV1485                              0                               0
#> ASV1461                              0                               0
#> ASV1313                              0                               0
#> ASV1420                              0                               0
#> ASV105                               0                               0
#> ASV431                               0                               0
#> ASV377                               0                               0
#> ASV1561                              0                               0
#> ASV477                               0                               0
#> ASV1537                              0                               0
#> ASV847                               0                               0
#> ASV722                               0                               0
#> ASV404                               0                               0
#> ASV1510                              0                               0
#> ASV1526                              0                               0
#> ASV1379                              0                               0
#> ASV367                               0                               0
#> ASV1642                              0                               0
#> ASV1233                              0                               0
#> ASV1386                              0                               0
#> ASV1234                              0                               0
#> ASV1340                              0                               0
#> ASV1035                              0                               0
#> ASV1037                              0                               0
#> ASV48                                0                               0
#> ASV711                               0                               0
#> ASV1609                              0                               0
#> ASV33                                0                               0
#> ASV1683                              0                               0
#> ASV673                               0                               0
#> ASV1548                              0                               0
#> ASV1460                              0                               0
#> ASV545                               0                               0
#> ASV1410                              0                               0
#> ASV1462                              0                              13
#> ASV853                               0                               0
#> ASV1067                              0                               0
#> ASV1311                              0                               0
#> ASV1552                              0                               0
#> ASV1588                              0                               0
#> ASV359                               0                               0
#> ASV440                               0                               0
#> ASV580                               0                               0
#> ASV858                               0                               0
#> ASV1632                              0                               0
#> ASV526                               0                               0
#> ASV1179                              0                               0
#> ASV1428                              0                               0
#> ASV355                               0                               0
#> ASV1010                              0                               0
#> ASV1542                              0                               0
#> ASV1577                              0                               0
#> ASV1262                              0                               0
#> ASV27                                0                               0
#> ASV137                               0                               0
#> ASV930                               0                               0
#> ASV962                               0                               0
#> ASV1242                              0                               0
#> ASV1628                              0                               0
#> ASV1650                              0                               0
#> ASV903                               0                               0
#> ASV984                               0                               0
#> ASV1011                              0                               0
#> ASV1224                              0                               0
#> ASV1276                              0                               0
#> ASV1468                              5                               0
#> ASV1574                              0                               0
#> ASV1630                              0                               0
#> ASV210                               0                               0
#> ASV383                               0                               0
#> ASV753                               0                               0
#> ASV996                               0                               0
#> ASV1069                              0                               0
#> ASV1082                              0                               0
#> ASV1109                              0                               0
#> ASV1144                              0                               0
#> ASV1419                              0                               0
#> ASV1467                              0                               0
#> ASV42                                0                               0
#> ASV104                               0                               0
#> ASV203                               0                               0
#> ASV491                               0                               0
#> ASV517                               0                               3
#> ASV598                               0                               0
#> ASV626                               0                               0
#> ASV650                               0                               3
#> ASV790                               0                               0
#> ASV839                               0                               0
#> ASV1014                              0                               0
#> ASV1017                              0                               0
#> ASV1030                              0                               0
#> ASV1085                              0                               0
#> ASV1167                              0                               0
#> ASV1625                              0                               0
#> ASV25                                0                               0
#> ASV249                               0                               0
#> ASV272                               0                               0
#> ASV320                               0                               0
#> ASV334                               0                               0
#> ASV488                               0                               0
#> ASV515                               0                               0
#> ASV549                               0                               0
#> ASV570                               0                               0
#> ASV604                               0                               0
#> ASV616                               0                               0
#> ASV736                               0                               0
#> ASV748                               0                               0
#> ASV760                               0                               0
#> ASV823                               2                               0
#> ASV828                               0                               2
#> ASV829                               0                               0
#> ASV891                               2                               0
#> ASV914                               0                               0
#> ASV961                               0                               0
#> ASV979                               0                               0
#> ASV986                               0                               0
#> ASV1133                              0                               0
#> ASV1270                              0                               0
#> ASV1422                              0                               0
#> ASV1523                              0                               0
#> ASV1576                              0                               0
#> ASV1603                              0                               0
#> ASV1726                              0                               0
#> ASV67                                0                               0
#> ASV87                                0                               0
#> ASV116                               0                               0
#> ASV171                               0                               0
#> ASV193                               0                               0
#> ASV198                               0                               0
#> ASV314                               0                               0
#> ASV405                               0                               0
#> ASV415                               0                               0
#> ASV731                               1                               0
#> ASV796                               0                               0
#> ASV840                               0                               0
#> ASV860                               0                               0
#> ASV875                               0                               0
#> ASV892                               0                               0
#> ASV940                               0                               0
#> ASV941                               0                               0
#> ASV959                               0                               1
#> ASV981                               0                               1
#> ASV988                               0                               0
#> ASV1015                              0                               0
#> ASV1020                              0                               0
#> ASV1027                              0                               0
#> ASV1031                              0                               0
#> ASV1066                              0                               0
#> ASV1107                              0                               0
#> ASV1111                              0                               0
#> ASV1239                              0                               0
#> ASV1263                              0                               0
#> ASV1302                              0                               0
#> ASV1303                              0                               0
#> ASV1341                              0                               0
#> ASV1380                              0                               0
#> ASV1388                              0                               0
#> ASV1484                              0                               0
#> ASV1550                              0                               0
#> ASV1562                              0                               0
#> ASV1635                              0                               0
#> ASV1653                              0                               0
#> ASV1674                              0                               0
#> ASV1712                              0                               0
#> 
#> 
#> $merged_ASV
#>         total spread parent_id curated rank
#> ASV666     47      5      ASV8  merged  154
#> ASV989     16      5     ASV38  merged  208
#> ASV94     678      4      ASV8  merged   33
#> ASV178    308      4     ASV19  merged   61
#> ASV209    264      4      ASV8  merged   67
#> ASV1058   139      4     ASV89  merged   97
#> ASV261    102      4      ASV8  merged  109
#> ASV816     16      4     ASV38  merged  209
#> ASV542    428      3    ASV507  merged   48
#> ASV85     371      3    ASV483  merged   51
#> ASV870    224      3    ASV577  merged   71
#> ASV333    149      3     ASV19  merged   91
#> ASV348    108      3      ASV8  merged  108
#> ASV98      39      3    ASV443  merged  166
#> ASV493     36      3      ASV8  merged  173
#> ASV724     18      3    ASV509  merged  205
#> ASV566     15      3    ASV443  merged  213
#> ASV831     12      3      ASV8  merged  231
#> ASV1078    12      3     ASV19  merged  232
#> ASV859      5      3    ASV509  merged  286
#> ASV1356     5      3   ASV1279  merged  288
#> ASV313    574      2     ASV89  merged   39
#> ASV462    334      2     ASV89  merged   59
#> ASV594    293      2     ASV19  merged   63
#> ASV534    277      2     ASV89  merged   66
#> ASV915    167      2     ASV43  merged   84
#> ASV346    165      2    ASV618  merged   88
#> ASV1268    90      2    ASV139  merged  117
#> ASV1387    82      2     ASV43  merged  124
#> ASV592     78      2    ASV618  merged  125
#> ASV766     60      2    ASV139  merged  141
#> ASV975     46      2     ASV38  merged  157
#> ASV1493    38      2    ASV344  merged  171
#> ASV1131    33      2     ASV19  merged  181
#> ASV1101    32      2     ASV89  merged  183
#> ASV474     20      2      ASV8  merged  197
#> ASV744     20      2    ASV139  merged  198
#> ASV1261    13      2    ASV618  merged  225
#> ASV963     10      2    ASV172  merged  243
#> ASV880      9      2     ASV52  merged  251
#> ASV1052     9      2    ASV509  merged  252
#> ASV1128     7      2    ASV509  merged  264
#> ASV358      6      2     ASV69  merged  270
#> ASV420      6      2    ASV509  merged  271
#> ASV693      6      2   ASV1409  merged  272
#> ASV884      6      2    ASV509  merged  274
#> ASV1319     6      2   ASV1074  merged  277
#> ASV1557     6      2     ASV89  merged  278
#> ASV1315     5      2    ASV784  merged  289
#> ASV576      4      2      ASV2  merged  299
#> ASV641      4      2    ASV505  merged  300
#> ASV1198     4      2    ASV509  merged  301
#> ASV466      3      2     ASV38  merged  315
#> ASV1036     3      2    ASV113  merged  316
#> ASV1247     3      2      ASV2  merged  317
#> ASV1427     3      2    ASV509  merged  319
#> ASV1176     2      2    ASV292  merged  338
#> ASV168   2877      1     ASV64  merged   14
#> ASV762    293      1    ASV489  merged   64
#> ASV461    258      1    ASV139  merged   68
#> ASV833    248      1    ASV221  merged   69
#> ASV998    173      1    ASV444  merged   80
#> ASV254    168      1     ASV38  merged   83
#> ASV1092   146      1    ASV251  merged   92
#> ASV1117   140      1     ASV43  merged   95
#> ASV832    130      1    ASV504  merged  102
#> ASV1182   114      1      ASV8  merged  106
#> ASV624     91      1     ASV12  merged  115
#> ASV1320    91      1    ASV898  merged  116
#> ASV999     84      1    ASV139  merged  120
#> ASV1313    70      1     ASV78  merged  133
#> ASV847     46      1    ASV107  merged  158
#> ASV722     45      1    ASV344  merged  159
#> ASV1526    39      1     ASV55  merged  169
#> ASV367     36      1     ASV28  merged  175
#> ASV1386    29      1    ASV969  merged  186
#> ASV1035    25      1   ASV1007  merged  190
#> ASV1037    23      1    ASV692  merged  191
#> ASV711     21      1     ASV45  merged  194
#> ASV33      19      1     ASV78  merged  203
#> ASV1410    13      1    ASV505  merged  228
#> ASV853     12      1     ASV38  merged  236
#> ASV1067    12      1     ASV45  merged  237
#> ASV440     10      1    ASV545  merged  246
#> ASV580     10      1     ASV64  merged  247
#> ASV858     10      1     ASV19  merged  248
#> ASV526      9      1   ASV1552  merged  253
#> ASV355      8      1    ASV151  merged  259
#> ASV1010     8      1     ASV19  merged  260
#> ASV1542     8      1     ASV12  merged  261
#> ASV1262     7      1     ASV12  merged  265
#> ASV962      6      1      ASV2  merged  282
#> ASV903      5      1    ASV637  merged  290
#> ASV984      5      1    ASV113  merged  291
#> ASV1011     5      1   ASV1532  merged  292
#> ASV1224     5      1   ASV1359  merged  293
#> ASV1276     5      1     ASV52  merged  294
#> ASV996      4      1     ASV45  merged  305
#> ASV1082     4      1      ASV2  merged  307
#> ASV1109     4      1   ASV1233  merged  308
#> ASV1467     4      1     ASV19  merged  311
#> ASV203      3      1     ASV64  merged  323
#> ASV491      3      1    ASV210  merged  324
#> ASV626      3      1    ASV357  merged  327
#> ASV650      3      1    ASV378  merged  328
#> ASV790      3      1    ASV505  merged  329
#> ASV1014     3      1    ASV969  merged  331
#> ASV1625     3      1     ASV78  merged  336
#> ASV249      2      1     ASV64  merged  343
#> ASV272      2      1    ASV839  merged  344
#> ASV488      2      1    ASV404  merged  347
#> ASV604      2      1   ASV1311  merged  351
#> ASV616      2      1    ASV144  merged  352
#> ASV736      2      1     ASV19  merged  353
#> ASV748      2      1    ASV753  merged  354
#> ASV760      2      1     ASV45  merged  355
#> ASV828      2      1   ASV1192  merged  357
#> ASV961      2      1    ASV292  merged  361
#> ASV979      2      1     ASV19  merged  362
#> ASV986      2      1     ASV19  merged  363
#> ASV1133     2      1   ASV1552  merged  364
#> ASV1603     2      1    ASV113  merged  369
#> ASV67       1      1     ASV12  merged  371
#> ASV116      1      1     ASV38  merged  373
#> ASV405      1      1     ASV69  merged  378
#> ASV415      1      1    ASV113  merged  379
#> ASV796      1      1   ASV1532  merged  381
#> ASV840      1      1    ASV404  merged  382
#> ASV860      1      1    ASV504  merged  383
#> ASV875      1      1    ASV787  merged  384
#> ASV940      1      1     ASV45  merged  386
#> ASV941      1      1     ASV78  merged  387
#> ASV959      1      1     ASV81  merged  388
#> ASV981      1      1    ASV172  merged  389
#> ASV1015     1      1    ASV637  merged  391
#> ASV1020     1      1    ASV113  merged  392
#> ASV1027     1      1      ASV2  merged  393
#> ASV1031     1      1   ASV1300  merged  394
#> ASV1066     1      1    ASV113  merged  395
#> ASV1107     1      1     ASV38  merged  396
#> ASV1111     1      1    ASV443  merged  397
#> ASV1263     1      1   ASV1144  merged  399
#> ASV1303     1      1    ASV746  merged  401
#> ASV1341     1      1   ASV1460  merged  402
#> ASV1380     1      1    ASV113  merged  403
#> ASV1388     1      1    ASV618  merged  404
#> ASV1484     1      1     ASV45  merged  405
#> 
# }
```
