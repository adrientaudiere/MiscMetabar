# Load sample data from file and rename samples using names of samples and an optional order

[![lifecycle-maturing](https://img.shields.io/badge/lifecycle-maturing-blue)](https://adrientaudiere.github.io/MiscMetabar/articles/Rules.html#lifecycle)

Useful for targets bioinformatic pipeline.

## Usage

``` r
sample_data_with_new_names(
  file_path,
  names_of_samples,
  samples_order = NULL,
  ...
)
```

## Arguments

- file_path:

  (required) a path to the sample_data file

- names_of_samples:

  (required) a vector of sample names

- samples_order:

  Optional numeric vector to sort sample names

- ...:

  Additional arguments passed on to
  [`utils::read.delim()`](https://rdrr.io/r/utils/read.table.html)
  function.

## Value

A data.frame from file_path and new names

## See also

[`rename_samples()`](https://adrientaudiere.github.io/MiscMetabar/reference/rename_samples.md)

## Author

Adrien Taudière

## Examples

``` r
sam_file <- system.file("extdata", "sam_data.csv", package = "MiscMetabar")
sample_data_with_new_names(sam_file, paste0("Samples_", seq(1, 185)))
#>                                             X      Sample_names   Tree_name
#> Samples_1      A10-005-B_S188_MERGED.fastq.gz    A10-005-B_S188     A10-005
#> Samples_2      A10-005-H_S189_MERGED.fastq.gz    A10-005-H_S189     A10-005
#> Samples_3      A10-005-M_S190_MERGED.fastq.gz    A10-005-M_S190     A10-005
#> Samples_4        A12-007_S191_MERGED.fastq.gz      A12-007_S191     A12-007
#> Samples_5        A12-007-B_S2_MERGED.fastq.gz      A12-007-B_S2     A12-007
#> Samples_6          A15-004_S3_MERGED.fastq.gz        A15-004_S3     A15-004
#> Samples_7           A8-005_S4_MERGED.fastq.gz         A8-005_S4      A8-005
#> Samples_8           A9-012_S5_MERGED.fastq.gz         A9-012_S5      A9-012
#> Samples_9      AB29-ABMX-H_S6_MERGED.fastq.gz    AB29-ABMX-H_S6  AB29-abm-X
#> Samples_10        AC27-013_S7_MERGED.fastq.gz       AC27-013_S7    AC27-013
#> Samples_11         AC29033_S8_MERGED.fastq.gz        AC29033_S8    AC29-033
#> Samples_12      AD26-005-B_S9_MERGED.fastq.gz     AD26-005-B_S9    AD26-005
#> Samples_13     AD26-005-H_S10_MERGED.fastq.gz    AD26-005-H_S10    AD26-005
#> Samples_14     AD26-005-M_S11_MERGED.fastq.gz    AD26-005-M_S11    AD26-005
#> Samples_15    AD30-ABMX-M_S12_MERGED.fastq.gz   AD30-ABMX-M_S12  AD30-abm-X
#> Samples_16     AD32-007-M_S13_MERGED.fastq.gz    AD32-007-M_S13    AD32-007
#> Samples_17     ADABM30X-B_S14_MERGED.fastq.gz    ADABM30X-B_S14  AD30-abm-X
#> Samples_18     ADABM30X-H_S15_MERGED.fastq.gz    ADABM30X-H_S15  AD30-abm-X
#> Samples_19     ADABM30X-M_S16_MERGED.fastq.gz    ADABM30X-M_S16  AD30-abm-X
#> Samples_20    AE30-ABM507_S17_MERGED.fastq.gz   AE30-ABM507_S17 AE30-abm507
#> Samples_21        B17-014_S18_MERGED.fastq.gz       B17-014_S18     BI7-014
#> Samples_22      B18-006-B_S19_MERGED.fastq.gz     B18-006-B_S19     B18-006
#> Samples_23    BA16-036bis_S20_MERGED.fastq.gz   BA16-036bis_S20    BA16-036
#> Samples_24     BA17-050-B_S21_MERGED.fastq.gz    BA17-050-B_S21    BA17-050
#> Samples_25     BB19-006-H_S22_MERGED.fastq.gz    BB19-006-H_S22    BB19-006
#> Samples_26      BB6-019-B_S23_MERGED.fastq.gz     BB6-019-B_S23     BB6-019
#> Samples_27      BB6-019-H_S24_MERGED.fastq.gz     BB6-019-H_S24     BB6-019
#> Samples_28      BB6-019-M_S25_MERGED.fastq.gz     BB6-019-M_S25     BB6-019
#> Samples_29       BD14-021_S26_MERGED.fastq.gz      BD14-021_S26    BD14-021
#> Samples_30      BE9-006-B_S27_MERGED.fastq.gz     BE9-006-B_S27     BE9-006
#> Samples_31      BE9-006-H_S28_MERGED.fastq.gz     BE9-006-H_S28     BE9-006
#> Samples_32      BE9-006-M_S29_MERGED.fastq.gz     BE9-006-M_S29     BE9-006
#> Samples_33      BG7-010-B_S30_MERGED.fastq.gz     BG7-010-B_S30     BG7-010
#> Samples_34      BG7-010-H_S31_MERGED.fastq.gz     BG7-010-H_S31     BG7-010
#> Samples_35      BG7-010-M_S32_MERGED.fastq.gz     BG7-010-M_S32     BG7-010
#> Samples_36        BH9-021_S33_MERGED.fastq.gz       BH9-021_S33     BH9-021
#> Samples_37     BJ17-007-M_S34_MERGED.fastq.gz    BJ17-007-M_S34    BJ17-007
#> Samples_38    BJ8-ABM-003_S35_MERGED.fastq.gz   BJ8-ABM-003_S35   BJ8abm003
#> Samples_39      BL7-006-B_S36_MERGED.fastq.gz     BL7-006-B_S36     BL7-006
#> Samples_40      BL7-006-H_S37_MERGED.fastq.gz     BL7-006-H_S37     BL7-006
#> Samples_41      BL7-006-M_S38_MERGED.fastq.gz     BL7-006-M_S38     BL7-006
#> Samples_42       BN11-041_S39_MERGED.fastq.gz      BN11-041_S39    BN11-041
#> Samples_43        BN9-002_S40_MERGED.fastq.gz       BN9-002_S40     BN9-002
#> Samples_44        BO8-002_S41_MERGED.fastq.gz       BO8-002_S41     BO8-002
#> Samples_45        BO8-005_S42_MERGED.fastq.gz       BO8-005_S42     BO8-005
#> Samples_46     BP11-001-B_S43_MERGED.fastq.gz    BP11-001-B_S43    BP11-001
#> Samples_47     BP11-001-H_S44_MERGED.fastq.gz    BP11-001-H_S44    BP11-001
#> Samples_48     BP11-001-M_S45_MERGED.fastq.gz    BP11-001-M_S45    BP11-001
#> Samples_49     BP12-025-B_S46_MERGED.fastq.gz    BP12-025-B_S46    BP12-025
#> Samples_50       BP14-006_S47_MERGED.fastq.gz      BP14-006_S47    BP14-006
#> Samples_51        BQ3-019_S48_MERGED.fastq.gz       BQ3-019_S48     BQ3-019
#> Samples_52      BQ4-018-B_S49_MERGED.fastq.gz     BQ4-018-B_S49     BQ4-018
#> Samples_53      BQ4-018-H_S50_MERGED.fastq.gz     BQ4-018-H_S50     BQ4-018
#> Samples_54      BQ4-018-M_S51_MERGED.fastq.gz     BQ4-018-M_S51     BQ4-018
#> Samples_55     BQ9ABM-002_S52_MERGED.fastq.gz    BQ9ABM-002_S52   BQ9abm002
#> Samples_56        BR8-005_S53_MERGED.fastq.gz       BR8-005_S53     BR8-005
#> Samples_57       BS14-006_S54_MERGED.fastq.gz      BS14-006_S54    BS14-006
#> Samples_58       BT-006-M_S55_MERGED.fastq.gz      BT-006-M_S55     BT7-006
#> Samples_59        BT7-006_S56_MERGED.fastq.gz       BT7-006_S56     BT7-006
#> Samples_60     BV11-002-B_S57_MERGED.fastq.gz    BV11-002-B_S57    BU11-002
#> Samples_61     BV11-002-H_S58_MERGED.fastq.gz    BV11-002-H_S58    BU11-002
#> Samples_62     BV11-002-M_S59_MERGED.fastq.gz    BV11-002-M_S59    BU11-002
#> Samples_63        BW8-003_S60_MERGED.fastq.gz       BW8-003_S60     BW8-003
#> Samples_64         C1-001_S61_MERGED.fastq.gz        C1-001_S61      C1-001
#> Samples_65      C21-NV1-B_S62_MERGED.fastq.gz     C21-NV1-B_S62    C21-nv-1
#> Samples_66      C21-NV1-H_S63_MERGED.fastq.gz     C21-NV1-H_S63    C21-nv-1
#> Samples_67      C21-NV1-M_S64_MERGED.fastq.gz     C21-NV1-M_S64    C21-nv-1
#> Samples_68         C9-005_S65_MERGED.fastq.gz        C9-005_S65      C9-005
#> Samples_69       CA12-024_S66_MERGED.fastq.gz      CA12-024_S66    CA12-024
#> Samples_70        CA9-027_S67_MERGED.fastq.gz       CA9-027_S67     CA9-027
#> Samples_71          CA9-X_S68_MERGED.fastq.gz         CA9-X_S68       CA9-X
#> Samples_72      CB8-019-B_S69_MERGED.fastq.gz     CB8-019-B_S69     CB8-019
#> Samples_73      CB8-019-H_S70_MERGED.fastq.gz     CB8-019-H_S70     CB8-019
#> Samples_74      CB8-019-M_S71_MERGED.fastq.gz     CB8-019-M_S71     CB8-019
#> Samples_75        CB9-013_S72_MERGED.fastq.gz       CB9-013_S72     CB9-013
#> Samples_76        CC3-044_S73_MERGED.fastq.gz       CC3-044_S73     CC3-044
#> Samples_77        CC8-003_S74_MERGED.fastq.gz       CC8-003_S74     CC8-003
#> Samples_78        D17-011_S77_MERGED.fastq.gz       D17-011_S77     D17-011
#> Samples_79      D18-003-B_S78_MERGED.fastq.gz     D18-003-B_S78     D18-003
#> Samples_80      D18-003-H_S79_MERGED.fastq.gz     D18-003-H_S79     D18-003
#> Samples_81      D18-003-M_S80_MERGED.fastq.gz     D18-003-M_S80     D18-003
#> Samples_82     D22-NVABM1_S81_MERGED.fastq.gz    D22-NVABM1_S81 D22-nvabm-1
#> Samples_83      D61-010-B_S82_MERGED.fastq.gz     D61-010-B_S82     DG1-010
#> Samples_84       D9-027-B_S83_MERGED.fastq.gz      D9-027-B_S83      D9-027
#> Samples_85       D9-027-H_S84_MERGED.fastq.gz      D9-027-H_S84      D9-027
#> Samples_86       D9-027-M_S85_MERGED.fastq.gz      D9-027-M_S85      D9-027
#> Samples_87    DBM-ABM-001_S86_MERGED.fastq.gz   DBM-ABM-001_S86   DB1abm001
#> Samples_88      DJ2-008-B_S87_MERGED.fastq.gz     DJ2-008-B_S87     DJ2-008
#> Samples_89      DJ2-008-H_S88_MERGED.fastq.gz     DJ2-008-H_S88     DJ2-008
#> Samples_90      DJ2-008-M_S89_MERGED.fastq.gz     DJ2-008-M_S89     DJ2-008
#> Samples_91     DP4-ABM001_S90_MERGED.fastq.gz    DP4-ABM001_S90   DP4abm001
#> Samples_92   DS1-ABM002-B_S91_MERGED.fastq.gz  DS1-ABM002-B_S91   DS1abm002
#> Samples_93   DS1-ABM002-H_S92_MERGED.fastq.gz  DS1-ABM002-H_S92   DS1abm002
#> Samples_94   DS1-ABM002-M_S93_MERGED.fastq.gz  DS1-ABM002-M_S93   DS1abm002
#> Samples_95      DU3-045-B_S94_MERGED.fastq.gz     DU3-045-B_S94     DU3-045
#> Samples_96        DW4-007_S95_MERGED.fastq.gz       DW4-007_S95     DW4-007
#> Samples_97      DY5-004-B_S96_MERGED.fastq.gz     DY5-004-B_S96     DY5-004
#> Samples_98      DY5-004-H_S97_MERGED.fastq.gz     DY5-004-H_S97     DY5-004
#> Samples_99      DY5-004-M_S98_MERGED.fastq.gz     DY5-004-M_S98     DY5-004
#> Samples_100   DZ6-ABM-001_S99_MERGED.fastq.gz   DZ6-ABM-001_S99   DZ6abm001
#> Samples_101     E9-009-B_S100_MERGED.fastq.gz     E9-009-B_S100      E9-009
#> Samples_102     E9-009-H_S101_MERGED.fastq.gz     E9-009-H_S101      E9-009
#> Samples_103     E9-009-M_S102_MERGED.fastq.gz     E9-009-M_S102      E9-009
#> Samples_104  EA5-ABM-001_S103_MERGED.fastq.gz  EA5-ABM-001_S103    EAabm001
#> Samples_105    EC2-013-B_S104_MERGED.fastq.gz    EC2-013-B_S104     EC2-013
#> Samples_106   F6-ABM-001_S105_MERGED.fastq.gz   F6-ABM-001_S105    F6abm001
#> Samples_107     F7-015-M_S106_MERGED.fastq.gz     F7-015-M_S106      F7-015
#> Samples_108    FOMES19-H_S108_MERGED.fastq.gz    FOMES19-H_S108     FOMES19
#> Samples_109    FOMES19-M_S109_MERGED.fastq.gz    FOMES19-M_S109     FOMES19
#> Samples_110    H10-018-M_S110_MERGED.fastq.gz    H10-018-M_S110     H10-018
#> Samples_111 H24-NVABM1-H_S111_MERGED.fastq.gz H24-NVABM1-H_S111 H24-nvabm-1
#> Samples_112    J18-004-B_S114_MERGED.fastq.gz    J18-004-B_S114     J18-004
#> Samples_113    J18-004-H_S115_MERGED.fastq.gz    J18-004-H_S115     J18-004
#> Samples_114    J18-004-M_S116_MERGED.fastq.gz    J18-004-M_S116     J18-004
#> Samples_115    K18-002-H_S117_MERGED.fastq.gz    K18-002-H_S117     K18-002
#> Samples_116   K26-NVABM1_S118_MERGED.fastq.gz   K26-NVABM1_S118 K26-nvamb-1
#> Samples_117       L19X-B_S119_MERGED.fastq.gz       L19X-B_S119       L19-X
#> Samples_118       L19X-H_S120_MERGED.fastq.gz       L19X-H_S120       L19-X
#> Samples_119       L19X-M_S121_MERGED.fastq.gz       L19X-M_S121       L19-X
#> Samples_120    L23-002-B_S122_MERGED.fastq.gz    L23-002-B_S122     L23-002
#> Samples_121    L23-002-H_S123_MERGED.fastq.gz    L23-002-H_S123     L23-002
#> Samples_122    L23-002-M_S124_MERGED.fastq.gz    L23-002-M_S124     L23-002
#> Samples_123      M22-001_S125_MERGED.fastq.gz      M22-001_S125     M22-001
#> Samples_124       N19X-B_S126_MERGED.fastq.gz       N19X-B_S126       N19-X
#> Samples_125       N19X-H_S127_MERGED.fastq.gz       N19X-H_S127       N19-X
#> Samples_126       N19X-M_S128_MERGED.fastq.gz       N19X-M_S128       N19-X
#> Samples_127    N22-001-B_S129_MERGED.fastq.gz    N22-001-B_S129     M22-001
#> Samples_128    N23-002-B_S130_MERGED.fastq.gz    N23-002-B_S130     N23-002
#> Samples_129    N23-002-H_S131_MERGED.fastq.gz    N23-002-H_S131     N23-002
#> Samples_130    N23-002-M_S132_MERGED.fastq.gz    N23-002-M_S132     N23-002
#> Samples_131     N25-ABMX_S133_MERGED.fastq.gz     N25-ABMX_S133   N25-abm-X
#> Samples_132   NVABM-0058_S134_MERGED.fastq.gz   NVABM-0058_S134   nvabm0058
#> Samples_133 NVABM-0163-H_S135_MERGED.fastq.gz NVABM-0163-H_S135   nvabm0163
#> Samples_134   NVABM-0397_S138_MERGED.fastq.gz   NVABM-0397_S138   nvabm0397
#> Samples_135    NVABM0216_S136_MERGED.fastq.gz    NVABM0216_S136   nvabm0216
#> Samples_136  NVABM0244-M_S137_MERGED.fastq.gz  NVABM0244-M_S137   nvabm0244
#> Samples_137      O20-X-B_S139_MERGED.fastq.gz      O20-X-B_S139       O20-X
#> Samples_138      O20-X-H_S140_MERGED.fastq.gz      O20-X-H_S140       O20-X
#> Samples_139      O20-X-M_S141_MERGED.fastq.gz      O20-X-M_S141       O20-X
#> Samples_140    O21-007-B_S142_MERGED.fastq.gz    O21-007-B_S142     O21-007
#> Samples_141    O21-007-H_S143_MERGED.fastq.gz    O21-007-H_S143     O21-007
#> Samples_142    O21-007-M_S144_MERGED.fastq.gz    O21-007-M_S144     O21-007
#> Samples_143    O24-003-B_S145_MERGED.fastq.gz    O24-003-B_S145     O24-003
#> Samples_144    O24-003-H_S146_MERGED.fastq.gz    O24-003-H_S146     O24-003
#> Samples_145    O24-003-M_S147_MERGED.fastq.gz    O24-003-M_S147     O24-003
#> Samples_146    O26-004-B_S148_MERGED.fastq.gz    O26-004-B_S148     O26-004
#> Samples_147    O26-004-H_S149_MERGED.fastq.gz    O26-004-H_S149     O26-004
#> Samples_148    O26-004-M_S150_MERGED.fastq.gz    O26-004-M_S150     O26-004
#> Samples_149      O27-012_S151_MERGED.fastq.gz      O27-012_S151     O27-012
#> Samples_150     O9-005-B_S152_MERGED.fastq.gz     O9-005-B_S152      C9-005
#> Samples_151    P19-023-M_S153_MERGED.fastq.gz    P19-023-M_S153     P19-023
#> Samples_152    P27-015-M_S154_MERGED.fastq.gz    P27-015-M_S154     P27-015
#> Samples_153   P27-ABM001_S155_MERGED.fastq.gz   P27-ABM001_S155  P27-abm001
#> Samples_154 Q27-ABM003-B_S156_MERGED.fastq.gz Q27-ABM003-B_S156  Q27-abm003
#> Samples_155     R25-ABMX_S157_MERGED.fastq.gz     R25-ABMX_S157   R25-abm-X
#> Samples_156    R28-008-B_S158_MERGED.fastq.gz    R28-008-B_S158     R28-008
#> Samples_157    R28-008-H_S159_MERGED.fastq.gz    R28-008-H_S159     R28-008
#> Samples_158    R28-008-M_S160_MERGED.fastq.gz    R28-008-M_S160     R28-008
#> Samples_159      T28-011_S161_MERGED.fastq.gz      T28-011_S161     T28-011
#> Samples_160 T28-ABM602-B_S162_MERGED.fastq.gz T28-ABM602-B_S162  T28-abm602
#> Samples_161   U27-ABM002_S163_MERGED.fastq.gz   U27-ABM002_S163  U27-abm002
#> Samples_162     W25-ABMX_S164_MERGED.fastq.gz     W25-ABMX_S164   W25-abm-X
#> Samples_163    W26-001-B_S165_MERGED.fastq.gz    W26-001-B_S165     W26-001
#> Samples_164    W26-001-H_S166_MERGED.fastq.gz    W26-001-H_S166     W26-001
#> Samples_165    W26-001-M_S167_MERGED.fastq.gz    W26-001-M_S167     W26-001
#> Samples_166      W30-006_S168_MERGED.fastq.gz      W30-006_S168     W30-006
#> Samples_167     W9-025-M_S169_MERGED.fastq.gz     W9-025-M_S169     W29-025
#> Samples_168    X24-009-B_S170_MERGED.fastq.gz    X24-009-B_S170     X24-009
#> Samples_169    X24-009-H_S171_MERGED.fastq.gz    X24-009-H_S171     X24-009
#> Samples_170    X24-009-M_S172_MERGED.fastq.gz    X24-009-M_S172     X24-009
#> Samples_171      X24-010_S173_MERGED.fastq.gz      X24-010_S173     X24-010
#> Samples_172    X29-004-B_S174_MERGED.fastq.gz    X29-004-B_S174     X29-004
#> Samples_173    X29-004-H_S175_MERGED.fastq.gz    X29-004-H_S175     X29-004
#> Samples_174    X29-004-M_S176_MERGED.fastq.gz    X29-004-M_S176     X29-004
#> Samples_175 Y21-ABM484-H_S177_MERGED.fastq.gz Y21-ABM484-H_S177  Y31-abm484
#> Samples_176    Y28-002-B_S178_MERGED.fastq.gz    Y28-002-B_S178     Y28-002
#> Samples_177    Y28-002-H_S179_MERGED.fastq.gz    Y28-002-H_S179     Y28-002
#> Samples_178    Y28-002-M_S180_MERGED.fastq.gz    Y28-002-M_S180     Y28-002
#> Samples_179    Y29-007-B_S181_MERGED.fastq.gz    Y29-007-B_S181     Y29-007
#> Samples_180    Y29-007-H_S182_MERGED.fastq.gz    Y29-007-H_S182     Y29-007
#> Samples_181    Y29-007-M_S183_MERGED.fastq.gz    Y29-007-M_S183     Y29-007
#> Samples_182 Y31-ABM484-B_S184_MERGED.fastq.gz Y31-ABM484-B_S184  Y31-abm484
#> Samples_183    Z29-001-H_S185_MERGED.fastq.gz    Z29-001-H_S185     Z29-001
#> Samples_184      Z30-002_S186_MERGED.fastq.gz      Z30-002_S186     Z30-002
#> Samples_185 Z30-ABM560-M_S187_MERGED.fastq.gz Z30-ABM560-M_S187  Z30-abm560
#>             Sample_id Height Diameter Time
#> Samples_1         188    Low       52   15
#> Samples_2         189   High       52   15
#> Samples_3         190 Middle       52   15
#> Samples_4         191   <NA>     28,4    0
#> Samples_5           2    Low     28,4    0
#> Samples_6           3   <NA>     30,7    0
#> Samples_7           4   <NA>     32,8    0
#> Samples_8           5   <NA>     33,3    0
#> Samples_9           6   High       99    5
#> Samples_10          7   <NA>       32   15
#> Samples_11          8   <NA>     55,4   15
#> Samples_12          9    Low    115,5   15
#> Samples_13         10   High    115,5   15
#> Samples_14         11 Middle    115,5   15
#> Samples_15         12 Middle        -    5
#> Samples_16         13 Middle       52   NA
#> Samples_17         14    Low        -    5
#> Samples_18         15   High        -    5
#> Samples_19         16 Middle        -    5
#> Samples_20         17   <NA>       10    5
#> Samples_21         18   <NA>     21,5   NA
#> Samples_22         19    Low     32,5   10
#> Samples_23         20   <NA>     40,3   15
#> Samples_24         21    Low       30   10
#> Samples_25         22   High     37,5   15
#> Samples_26         23    Low     33,8    5
#> Samples_27         24   High     33,8    5
#> Samples_28         25 Middle     33,8    5
#> Samples_29         26   <NA>     32,4    5
#> Samples_30         27    Low     23,1    0
#> Samples_31         28   High     23,1    0
#> Samples_32         29 Middle     23,1    0
#> Samples_33         30    Low       32    0
#> Samples_34         31   High       32    0
#> Samples_35         32 Middle       32    0
#> Samples_36         33   <NA>     44,3    0
#> Samples_37         34 Middle       32    5
#> Samples_38         35   <NA>       40   15
#> Samples_39         36    Low     35,6    0
#> Samples_40         37   High     35,6    0
#> Samples_41         38 Middle     35,6    0
#> Samples_42         39   <NA>     38,3   NA
#> Samples_43         40   <NA>     32,5    0
#> Samples_44         41   <NA>     26,5   10
#> Samples_45         42   <NA>       23    0
#> Samples_46         43    Low     47,5   15
#> Samples_47         44   High     47,5   15
#> Samples_48         45 Middle     47,5   15
#> Samples_49         46    Low     64,5   10
#> Samples_50         47   <NA>     25,3   10
#> Samples_51         48   <NA>       80   15
#> Samples_52         49    Low       34    0
#> Samples_53         50   High       34    0
#> Samples_54         51 Middle       34    0
#> Samples_55         52   <NA>       40   15
#> Samples_56         53   <NA>     36,8    0
#> Samples_57         54   <NA>       35   15
#> Samples_58         55 Middle       41    0
#> Samples_59         56   <NA>       41    0
#> Samples_60         57    Low       33    5
#> Samples_61         58   High       33    5
#> Samples_62         59 Middle       33    5
#> Samples_63         60   <NA>     38,2    0
#> Samples_64         61   <NA>       43   10
#> Samples_65         62    Low       30    0
#> Samples_66         63   High       30    0
#> Samples_67         64 Middle       30    0
#> Samples_68         65   <NA>       39    5
#> Samples_69         66   <NA>     18,3   NA
#> Samples_70         67   <NA>     23,4    0
#> Samples_71         68   <NA>        -   15
#> Samples_72         69    Low     33,3    0
#> Samples_73         70   High     33,3    0
#> Samples_74         71 Middle     33,3    0
#> Samples_75         72   <NA>     32,3    5
#> Samples_76         73   <NA>     35,7    5
#> Samples_77         74   <NA>     36,5    0
#> Samples_78         77   <NA>     59,8    5
#> Samples_79         78    Low     45,5    5
#> Samples_80         79   High     45,5    5
#> Samples_81         80 Middle     45,5    5
#> Samples_82         81   <NA>       30    5
#> Samples_83         82    Low     36,3   15
#> Samples_84         83    Low     18,9   NA
#> Samples_85         84   High     18,9   NA
#> Samples_86         85 Middle     18,9   NA
#> Samples_87         86   <NA>       40   15
#> Samples_88         87    Low     76,4    0
#> Samples_89         88   High     76,4    0
#> Samples_90         89 Middle     76,4    0
#> Samples_91         90   <NA>       60   15
#> Samples_92         91    Low       25   15
#> Samples_93         92   High       25   15
#> Samples_94         93 Middle       25   15
#> Samples_95         94    Low       38   15
#> Samples_96         95   <NA>       34   10
#> Samples_97         96    Low     60,4    0
#> Samples_98         97   High     60,4    0
#> Samples_99         98 Middle     60,4    0
#> Samples_100        99   <NA>       35   15
#> Samples_101       100    Low     29,5    0
#> Samples_102       101   High     29,5    0
#> Samples_103       102 Middle     29,5    0
#> Samples_104       103   <NA>       25   15
#> Samples_105       104    Low       52   15
#> Samples_106       105   <NA>       30   15
#> Samples_107       106 Middle       50   10
#> Samples_108       108   High     <NA>   NA
#> Samples_109       109 Middle     <NA>   NA
#> Samples_110       110 Middle     88,5   15
#> Samples_111       111   High       30    5
#> Samples_112       114    Low     27,6    0
#> Samples_113       115   High     27,6    0
#> Samples_114       116 Middle     27,6    0
#> Samples_115       117   High     33,5   10
#> Samples_116       118   <NA>       30    5
#> Samples_117       119    Low       30   NA
#> Samples_118       120   High       30   NA
#> Samples_119       121 Middle       30   NA
#> Samples_120       122    Low     40,5   10
#> Samples_121       123   High     40,5   10
#> Samples_122       124 Middle     40,5   10
#> Samples_123       125   <NA>       34    0
#> Samples_124       126    Low     33,5   NA
#> Samples_125       127   High     33,5   NA
#> Samples_126       128 Middle     33,5   NA
#> Samples_127       129    Low       34    0
#> Samples_128       130    Low     21,5   15
#> Samples_129       131   High     21,5   15
#> Samples_130       132 Middle     21,5   15
#> Samples_131       133   <NA>       20    5
#> Samples_132       134   <NA>       30   10
#> Samples_133       135   High       40   10
#> Samples_134       138   <NA>       35   10
#> Samples_135       136   <NA>       30   15
#> Samples_136       137 Middle       50   10
#> Samples_137       139    Low     <NA>   NA
#> Samples_138       140   High     <NA>   NA
#> Samples_139       141 Middle     <NA>   NA
#> Samples_140       142    Low       25    0
#> Samples_141       143   High       25    0
#> Samples_142       144 Middle       25    0
#> Samples_143       145    Low       75    0
#> Samples_144       146   High       75    0
#> Samples_145       147 Middle       75    0
#> Samples_146       148    Low     30,4    0
#> Samples_147       149   High     30,4    0
#> Samples_148       150 Middle     30,4    0
#> Samples_149       151   <NA>     69,5   15
#> Samples_150       152    Low       39    5
#> Samples_151       153 Middle     18,4   NA
#> Samples_152       154 Middle       96    5
#> Samples_153       155   <NA>       25    5
#> Samples_154       156    Low       10    5
#> Samples_155       157   <NA>        -    5
#> Samples_156       158    Low     37,2    0
#> Samples_157       159   High     37,2    0
#> Samples_158       160 Middle     37,2    0
#> Samples_159       161   <NA>     66,8   10
#> Samples_160       162    Low       64    5
#> Samples_161       163   <NA>       11    5
#> Samples_162       164   <NA>       40    5
#> Samples_163       165    Low     71,7    0
#> Samples_164       166   High     71,7    0
#> Samples_165       167 Middle     71,7    0
#> Samples_166       168   <NA>     59,6   10
#> Samples_167       169 Middle     71,5   NA
#> Samples_168       170    Low     38,8    0
#> Samples_169       171   High     38,8    0
#> Samples_170       172 Middle     38,8    0
#> Samples_171       173   <NA>     37,8    5
#> Samples_172       174    Low     52,7   10
#> Samples_173       175   High     52,7   10
#> Samples_174       176 Middle     52,7   10
#> Samples_175       177   High       20    5
#> Samples_176       178    Low     63,6   NA
#> Samples_177       179   High     63,6   NA
#> Samples_178       180 Middle     63,6   NA
#> Samples_179       181    Low     68,6    0
#> Samples_180       182   High     68,6    0
#> Samples_181       183 Middle     68,6    0
#> Samples_182       184    Low       20    5
#> Samples_183       185   High     35,2    5
#> Samples_184       186   <NA>     43,4   15
#> Samples_185       187 Middle       11    5
```
