# Blast all sequence of `refseq` slot of a [`phyloseq-class`](https://rdrr.io/pkg/phyloseq/man/phyloseq-class.html) object against a custom database.

[![lifecycle-experimental](https://img.shields.io/badge/lifecycle-experimental-orange)](https://adrientaudiere.github.io/MiscMetabar/articles/Rules.html#lifecycle)

Use the blast software.

## Usage

``` r
blast_pq(
  physeq,
  fasta_for_db = NULL,
  database = NULL,
  blastpath = NULL,
  id_cut = 90,
  bit_score_cut = 50,
  min_cover_cut = 50,
  e_value_cut = 1e-30,
  unique_per_seq = FALSE,
  score_filter = TRUE,
  nproc = 1,
  args_makedb = NULL,
  args_blastn = NULL,
  keep_temporary_files = FALSE
)
```

## Arguments

- physeq:

  (required) a
  [`phyloseq-class`](https://rdrr.io/pkg/phyloseq/man/phyloseq-class.html)
  object obtained using the `phyloseq` package.

- fasta_for_db:

  Either a DNAStringSet object or a path to a fasta file to make the
  blast database.

- database:

  path to a blast database

- blastpath:

  path to blast program

- id_cut:

  (default: 90) cut of in identity percent to keep result

- bit_score_cut:

  (default: 50) cut of in bit score to keep result The higher the
  bit-score, the better the sequence similarity. The bit-score is the
  requires size of a sequence database in which the current match could
  be found just by chance. The bit-score is a log2 scaled and normalized
  raw-score. Each increase by one doubles the required database size
  (2bit-score).

- min_cover_cut:

  (default: 50) cut of in query cover (%) to keep result

- e_value_cut:

  (default: 1e-30) cut of in e-value (%) to keep result The BLAST
  E-value is the number of expected hits of similar quality (score) that
  could be found just by chance.

- unique_per_seq:

  (logical, default FALSE) if TRUE only return the better match (higher
  **bit score**) for each sequence

- score_filter:

  (logical, default TRUE) does results are filter by score? If FALSE,
  `id_cut`,`bit_score_cut`, `e_value_cut` and `min_cover_cut` are
  ignored

- nproc:

  (default: 1) Set to number of cpus/processors to use for blast (args
  -num_threads for blastn command)

- args_makedb:

  Additional arguments passed on to makeblastdb command

- args_blastn:

  Additional arguments passed on to blastn command

- keep_temporary_files:

  (logical, default: FALSE) Do we keep temporary files

  - db.fasta (refseq transformed into a database)

  - dbase list of files (output of blastn)

  - blast_result.txt the summary result of blastn using
    `-outfmt "6 qseqid qlen sseqid slen length pident evalue bitscore qcovs"`

  - temp.fasta if fasta_for_db was a DNAStringSet object.

## Value

a blast table

## See also

[`blast_to_phyloseq()`](https://adrientaudiere.github.io/MiscMetabar/dev/reference/blast_to_phyloseq.md)
to use `refseq` slot as a database

## Examples

``` r
blast_pq(data_fungi_mini,
  fasta_for_db = system.file("extdata", "mini_UNITE_fungi.fasta.gz",
    package = "MiscMetabar"
  )
)
#> Build the database from fasta_for_db
#> Blast refseq from physeq object against the database
#>       Query name Query seq. length
#> 5501       ASV41               364
#> 5502       ASV41               364
#> 15001      ASV75               337
#> 5503       ASV41               364
#> 15002      ASV75               337
#> 15003      ASV75               337
#> 5504       ASV41               364
#> 5505       ASV41               364
#> 5506       ASV41               364
#> 5507       ASV41               364
#> 5508       ASV41               364
#> 5509       ASV41               364
#> 5510       ASV41               364
#> 5511       ASV41               364
#> 5512       ASV41               364
#> 15004      ASV75               337
#> 5513       ASV41               364
#> 5514       ASV41               364
#> 5515       ASV41               364
#> 15005      ASV75               337
#> 5516       ASV41               364
#> 5517       ASV41               364
#> 15006      ASV75               337
#> 15007      ASV75               337
#> 15008      ASV75               337
#> 15009      ASV75               337
#> 15011      ASV75               337
#> 15014      ASV75               337
#> 14002      ASV71               307
#> 14003      ASV71               307
#> 20001     ASV101               387
#>                                                                                                                                                                                                                                             Taxa name
#> 5501                                                                     Mycenaceae_sp_OP470041_SH0013594.10FU_reps_singleton_;tax=k:Fungi,p:Basidiomycota,c:Agaricomycetes,o:Agaricales,f:Mycenaceae,g:Mycenaceae_gen_Incertae_sedis,s:Mycenaceae_sp
#> 5502                                                                                                 Mycena_sp_UDB05068035_SH0013561.10FU_reps_singleton_;tax=k:Fungi,p:Basidiomycota,c:Agaricomycetes,o:Agaricales,f:Mycenaceae,g:Mycena,s:Mycena_sp
#> 15001                                                                           Peniophora_reidii_MH861616_SH0003674.10FU_reps_singleton_;tax=k:Fungi,p:Basidiomycota,c:Agaricomycetes,o:Russulales,f:Peniophoraceae,g:Peniophora,s:Peniophora_reidii
#> 5503                                                                                                 Mycena_sp_UDB05033057_SH0013579.10FU_reps_singleton_;tax=k:Fungi,p:Basidiomycota,c:Agaricomycetes,o:Agaricales,f:Mycenaceae,g:Mycena,s:Mycena_sp
#> 15002                                                                                   Peniophora_sp_OR149912_SH0003663.10FU_reps_singleton_;tax=k:Fungi,p:Basidiomycota,c:Agaricomycetes,o:Russulales,f:Peniophoraceae,g:Peniophora,s:Peniophora_sp
#> 15003                                                                                   Peniophora_sp_MK588791_SH0003657.10FU_reps_singleton_;tax=k:Fungi,p:Basidiomycota,c:Agaricomycetes,o:Russulales,f:Peniophoraceae,g:Peniophora,s:Peniophora_sp
#> 5504                                                                                      Mycena_venustula_OK239713_SH0008654.10FU_reps_singleton_;tax=k:Fungi,p:Basidiomycota,c:Agaricomycetes,o:Agaricales,f:Mycenaceae,g:Mycena,s:Mycena_venustula
#> 5505                                                                                                 Mycena_sp_UDB05006031_SH0013576.10FU_reps_singleton_;tax=k:Fungi,p:Basidiomycota,c:Agaricomycetes,o:Agaricales,f:Mycenaceae,g:Mycena,s:Mycena_sp
#> 5506                                                                  Mycenaceae_sp_UDB05167898_SH0013566.10FU_reps_singleton_;tax=k:Fungi,p:Basidiomycota,c:Agaricomycetes,o:Agaricales,f:Mycenaceae,g:Mycenaceae_gen_Incertae_sedis,s:Mycenaceae_sp
#> 5507                                                                                                    Mycena_sp_OR987288_SH0013580.10FU_reps_singleton_;tax=k:Fungi,p:Basidiomycota,c:Agaricomycetes,o:Agaricales,f:Mycenaceae,g:Mycena,s:Mycena_sp
#> 5508                                             Tricholomataceae_sp_LC574741_SH0013558.10FU_reps_singleton_;tax=k:Fungi,p:Basidiomycota,c:Agaricomycetes,o:Agaricales,f:Tricholomataceae,g:Tricholomataceae_gen_Incertae_sedis,s:Tricholomataceae_sp
#> 5509                                                                                Mycena_acanthophila_OQ633198_SH0008640.10FU_reps_singleton_;tax=k:Fungi,p:Basidiomycota,c:Agaricomycetes,o:Agaricales,f:Mycenaceae,g:Mycena,s:Mycena_acanthophila
#> 5510                                                                            Mycena_quercus-ilicis_OQ633195_SH0008657.10FU_reps_singleton_;tax=k:Fungi,p:Basidiomycota,c:Agaricomycetes,o:Agaricales,f:Mycenaceae,g:Mycena,s:Mycena_quercus-ilicis
#> 5511                                                                                                    Mycena_sp_OR883535_SH0013563.10FU_reps_singleton_;tax=k:Fungi,p:Basidiomycota,c:Agaricomycetes,o:Agaricales,f:Mycenaceae,g:Mycena,s:Mycena_sp
#> 5512                                                                                                   Mycena_sp_UDB023583_SH0008647.10FU_reps_singleton_;tax=k:Fungi,p:Basidiomycota,c:Agaricomycetes,o:Agaricales,f:Mycenaceae,g:Mycena,s:Mycena_sp
#> 15004                                                                               Peniophora_pini_KT210102_SH0003677.10FU_reps_singleton_;tax=k:Fungi,p:Basidiomycota,c:Agaricomycetes,o:Russulales,f:Peniophoraceae,g:Peniophora,s:Peniophora_pini
#> 5513                                                                                                 Mycena_sp_UDB05068051_SH0013586.10FU_reps_singleton_;tax=k:Fungi,p:Basidiomycota,c:Agaricomycetes,o:Agaricales,f:Mycenaceae,g:Mycena,s:Mycena_sp
#> 5514                                                                                                 Mycena_sp_UDB05018702_SH0013573.10FU_reps_singleton_;tax=k:Fungi,p:Basidiomycota,c:Agaricomycetes,o:Agaricales,f:Mycenaceae,g:Mycena,s:Mycena_sp
#> 5515                                                                                                 Mycena_sp_UDB05154274_SH0008652.10FU_reps_singleton_;tax=k:Fungi,p:Basidiomycota,c:Agaricomycetes,o:Agaricales,f:Mycenaceae,g:Mycena,s:Mycena_sp
#> 15005                                                                                   Peniophora_sp_KJ832048_SH0003664.10FU_reps_singleton_;tax=k:Fungi,p:Basidiomycota,c:Agaricomycetes,o:Russulales,f:Peniophoraceae,g:Peniophora,s:Peniophora_sp
#> 5516                                                                                                 Mycena_sp_UDB05027734_SH0013564.10FU_reps_singleton_;tax=k:Fungi,p:Basidiomycota,c:Agaricomycetes,o:Agaricales,f:Mycenaceae,g:Mycena,s:Mycena_sp
#> 5517                                                                                                 Mycena_sp_UDB05115971_SH0013554.10FU_reps_singleton_;tax=k:Fungi,p:Basidiomycota,c:Agaricomycetes,o:Agaricales,f:Mycenaceae,g:Mycena,s:Mycena_sp
#> 15006                                  Fungi_sp_KU515816_SH0003692.10FU_reps_singleton_;tax=k:Fungi,p:Fungi_phy_Incertae_sedis,c:Fungi_cls_Incertae_sedis,o:Fungi_ord_Incertae_sedis,f:Fungi_fam_Incertae_sedis,g:Fungi_gen_Incertae_sedis,s:Fungi_sp
#> 15007                                              Subulicystidium_longisporum_OK310732_SH0003688.10FU_reps_singleton_;tax=k:Fungi,p:Basidiomycota,c:Agaricomycetes,o:Trechisporales,f:Hydnodontaceae,g:Subulicystidium,s:Subulicystidium_longisporum
#> 15008                                  Fungi_sp_GU721554_SH0003687.10FU_reps_singleton_;tax=k:Fungi,p:Fungi_phy_Incertae_sedis,c:Fungi_cls_Incertae_sedis,o:Fungi_ord_Incertae_sedis,f:Fungi_fam_Incertae_sedis,g:Fungi_gen_Incertae_sedis,s:Fungi_sp
#> 15009                                  Fungi_sp_GQ999289_SH0003662.10FU_reps_singleton_;tax=k:Fungi,p:Fungi_phy_Incertae_sedis,c:Fungi_cls_Incertae_sedis,o:Fungi_ord_Incertae_sedis,f:Fungi_fam_Incertae_sedis,g:Fungi_gen_Incertae_sedis,s:Fungi_sp
#> 15011                                                                                   Peniophora_sp_AB672622_SH0003656.10FU_reps_singleton_;tax=k:Fungi,p:Basidiomycota,c:Agaricomycetes,o:Russulales,f:Peniophoraceae,g:Peniophora,s:Peniophora_sp
#> 15014                                  Fungi_sp_JX135952_SH0003681.10FU_reps_singleton_;tax=k:Fungi,p:Fungi_phy_Incertae_sedis,c:Fungi_cls_Incertae_sedis,o:Fungi_ord_Incertae_sedis,f:Fungi_fam_Incertae_sedis,g:Fungi_gen_Incertae_sedis,s:Fungi_sp
#> 14002                                        Tremellales_sp_UDB06404467_SH0007343.10FU_reps_singleton_;tax=k:Fungi,p:Basidiomycota,c:Tremellomycetes,o:Tremellales,f:Tremellales_fam_Incertae_sedis,g:Tremellales_gen_Incertae_sedis,s:Tremellales_sp
#> 14003 Tremellomycetes_sp_UDB06410903_SH0007341.10FU_reps_singleton_;tax=k:Fungi,p:Basidiomycota,c:Tremellomycetes,o:Tremellomycetes_ord_Incertae_sedis,f:Tremellomycetes_fam_Incertae_sedis,g:Tremellomycetes_gen_Incertae_sedis,s:Tremellomycetes_sp
#> 20001                                            Corticiales_sp_KY462668_SH0001937.10FU_reps_singleton_;tax=k:Fungi,p:Basidiomycota,c:Agaricomycetes,o:Corticiales,f:Corticiales_fam_Incertae_sedis,g:Corticiales_gen_Incertae_sedis,s:Corticiales_sp
#>       Taxa seq. length Alignment length % id. match   e-value bit score
#> 5501               625              326      94.479 3.12e-142       499
#> 5502               619              302      96.026 1.88e-139       490
#> 15001              551              297      95.286 2.25e-133       470
#> 5503               651              304      93.750 2.46e-128       453
#> 15002              578              297      94.276 2.27e-128       453
#> 15003              571              304      93.092 6.35e-124       438
#> 5504               628              307      92.508 8.92e-123       435
#> 5505               657              304      92.434 1.15e-121       431
#> 5506               611              302      93.046 1.15e-121       431
#> 5507               577              305      92.459 1.93e-119       424
#> 5508               527              303      92.079 1.93e-119       424
#> 5509               635              313      91.693 1.93e-119       424
#> 5510               645              312      91.667 6.94e-119       422
#> 5511               650              298      92.282 2.50e-118       420
#> 5512               628              304      91.776 2.50e-118       420
#> 15004              563              305      91.803 2.98e-117       416
#> 5513               633              305      91.803 1.16e-116       414
#> 5514               597              303      92.079 1.16e-116       414
#> 5515               642              308      91.234 4.18e-116       412
#> 15005              547              302      90.728 1.08e-111       398
#> 5516               603              304      90.789 1.51e-110       394
#> 5517               599              306      90.523 5.44e-110       392
#> 15006              551              300      90.667 6.49e-109       388
#> 15007              530              250      94.800 6.49e-109       388
#> 15008              544              308      90.260 6.49e-109       388
#> 15009              526              303      90.099 2.33e-108       387
#> 15011              530              286      90.559 1.82e-104       374
#> 15014              346              210      94.286  8.63e-88       318
#> 14002              408              186      93.548  2.23e-73       270
#> 14003              406              182      92.857  1.34e-70       261
#> 20001              603              201      90.050  2.86e-68       254
#>       Query cover
#> 5501           89
#> 5502           83
#> 15001          88
#> 5503           83
#> 15002          88
#> 15003          89
#> 5504           83
#> 5505           83
#> 5506           83
#> 5507           84
#> 5508           82
#> 5509           83
#> 5510           83
#> 5511           82
#> 5512           83
#> 15004          89
#> 5513           84
#> 5514           82
#> 5515           83
#> 15005          89
#> 5516           82
#> 5517           82
#> 15006          86
#> 15007          74
#> 15008          89
#> 15009          89
#> 15011          83
#> 15014          62
#> 14002          59
#> 14003          59
#> 20001          51
```
