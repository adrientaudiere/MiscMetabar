# Summarize a tax_table (taxonomic slot of phyloseq object) using gtsummary

[![lifecycle-maturing](https://img.shields.io/badge/lifecycle-maturing-blue)](https://adrientaudiere.github.io/MiscMetabar/articles/Rules.html#lifecycle)

Mainly a wrapper for the
[`gtsummary::tbl_summary()`](https://www.danieldsjoberg.com/gtsummary/reference/tbl_summary.html)
function in the case of `physeq` object.

## Usage

``` r
tbl_sum_taxtable(physeq, taxonomic_ranks = NULL, ...)
```

## Arguments

- physeq:

  (required) a
  [`phyloseq-class`](https://rdrr.io/pkg/phyloseq/man/phyloseq-class.html)
  object obtained using the `phyloseq` package.

- taxonomic_ranks:

  A list of taxonomic ranks we want to summarized.

- ...:

  Additional arguments passed on to
  [`gtsummary::tbl_summary()`](https://www.danieldsjoberg.com/gtsummary/reference/tbl_summary.html)

## Value

A table of class c('tbl_summary', 'gtsummary')

## Author

Adrien Taudière

## Examples

``` r
tbl_sum_taxtable(data_fungi_mini)


  

Characteristic
```

**N = 45**¹

Domain

  

    Fungi

45 (100%)

Phylum

  

    Basidiomycota

45 (100%)

Class

  

    Agaricomycetes

41 (93%)

    Atractiellomycetes

1 (2.3%)

    Tremellomycetes

2 (4.5%)

    Unknown

1

Order

  

    Agaricales

7 (17%)

    Atractiellales

1 (2.4%)

    Auriculariales

6 (14%)

    Cantharellales

1 (2.4%)

    Corticiales

1 (2.4%)

    Hymenochaetales

5 (12%)

    Polyporales

11 (26%)

    Russulales

9 (21%)

    Tremellales

1 (2.4%)

    Unknown

3

Family

  

    Aporpiaceae

1 (2.6%)

    Atractiellales_fam_Incertae_sedis

1 (2.6%)

    Auriculariaceae

2 (5.1%)

    Cantharellales_fam_Incertae_sedis

1 (2.6%)

    Corticiaceae

1 (2.6%)

    Entolomataceae

1 (2.6%)

    Exidiaceae

3 (7.7%)

    Hericiaceae

1 (2.6%)

    Hymenochaetales_fam_Incertae_sedis

1 (2.6%)

    Hyphodermataceae

2 (5.1%)

    Lyophyllaceae

4 (10%)

    Peniophoraceae

1 (2.6%)

    Phanerochaetaceae

1 (2.6%)

    Polyporaceae

5 (13%)

    Pterulaceae

1 (2.6%)

    Russulales_fam_Incertae_sedis

1 (2.6%)

    Schizoporaceae

4 (10%)

    Steccherinaceae

1 (2.6%)

    Stereaceae

6 (15%)

    Tricholomataceae

1 (2.6%)

    Unknown

6

Genus

  

    Antrodiella

1 (2.7%)

    Auricularia

2 (5.4%)

    Basidiodendron

1 (2.7%)

    Elmerina

1 (2.7%)

    Entocybe

1 (2.7%)

    Exidia

2 (5.4%)

    Fomes

4 (11%)

    Gloeohypochnicium

1 (2.7%)

    Helicogloea

1 (2.7%)

    Hericium

1 (2.7%)

    Hyphoderma

2 (5.4%)

    Marchandiomyces

1 (2.7%)

    Mycena

1 (2.7%)

    Ossicaulis

4 (11%)

    Peniophora

1 (2.7%)

    Peniophorella

1 (2.7%)

    Phanerochaete

1 (2.7%)

    Radulomyces

1 (2.7%)

    Sistotrema

1 (2.7%)

    Stereum

4 (11%)

    Trametes

1 (2.7%)

    Xylodon

4 (11%)

    Unknown

8

Species

  

    analogum

1 (2.9%)

    brasiliensis

1 (2.9%)

    buckii

1 (2.9%)

    caryae

1 (2.9%)

    coralloides

1 (2.9%)

    eyrei

1 (2.9%)

    flaviporus

1 (2.9%)

    fomentarius

4 (11%)

    glandulosa

2 (5.7%)

    hirsutum

1 (2.9%)

    lachnopus

4 (11%)

    livescens

1 (2.9%)

    mesenterica

1 (2.9%)

    molaris

1 (2.9%)

    oblongisporum

1 (2.9%)

    ostrea

3 (8.6%)

    pellucida

1 (2.9%)

    pubera

1 (2.9%)

    raduloides

3 (8.6%)

    renati

1 (2.9%)

    roseocremeum

1 (2.9%)

    setigerum

1 (2.9%)

    versicolor

1 (2.9%)

    versiformis

1 (2.9%)

    Unknown

10

Trophic.Mode

  

    -

5 (11%)

    Pathotroph

1 (2.2%)

    Pathotroph-Saprotroph

2 (4.4%)

    Pathotroph-Saprotroph-Symbiotroph

1 (2.2%)

    Saprotroph

33 (73%)

    Saprotroph-Symbiotroph

3 (6.7%)

Guild

  

    -

5 (11%)

    Ectomycorrhizal-Wood Saprotroph

1 (2.2%)

    Endophyte-Undefined Saprotroph

2 (4.4%)

    Fungal Parasite-Undefined Saprotroph

1 (2.2%)

    Leaf Saprotroph-Plant Pathogen-Undefined Saprotroph-Wood Saprotroph

1 (2.2%)

    Lichen Parasite

1 (2.2%)

    Plant Pathogen-Wood Saprotroph

1 (2.2%)

    Undefined Saprotroph

20 (44%)

    Wood Saprotroph

11 (24%)

    Wood Saprotroph-Undefined Saprotroph

2 (4.4%)

Trait

  

    -

13 (29%)

    Brown Rot

13 (29%)

    Brown Rot-White Rot

13 (29%)

    White Rot

13 (29%)

Confidence.Ranking

  

    -

5 (11%)

    Highly Probable

3 (6.7%)

    Possible

2 (4.4%)

    Probable

35 (78%)

Genus_species

  

    Antrodiella_brasiliensis

1 (2.2%)

    Auricularia_NA

1 (2.2%)

    Auricularia_mesenterica

1 (2.2%)

    Basidiodendron_eyrei

1 (2.2%)

    Elmerina_caryae

1 (2.2%)

    Entocybe_NA

1 (2.2%)

    Exidia_glandulosa

2 (4.4%)

    Fomes_fomentarius

4 (8.9%)

    Gloeohypochnicium_analogum

1 (2.2%)

    Helicogloea_pellucida

1 (2.2%)

    Hericium_coralloides

1 (2.2%)

    Hyphoderma_roseocremeum

1 (2.2%)

    Hyphoderma_setigerum

1 (2.2%)

    Marchandiomyces_buckii

1 (2.2%)

    Mycena_renati

1 (2.2%)

    NA_NA

8 (18%)

    Ossicaulis_lachnopus

4 (8.9%)

    Peniophora_versiformis

1 (2.2%)

    Peniophorella_pubera

1 (2.2%)

    Phanerochaete_livescens

1 (2.2%)

    Radulomyces_molaris

1 (2.2%)

    Sistotrema_oblongisporum

1 (2.2%)

    Stereum_hirsutum

1 (2.2%)

    Stereum_ostrea

3 (6.7%)

    Trametes_versicolor

1 (2.2%)

    Xylodon_flaviporus

1 (2.2%)

    Xylodon_raduloides

3 (6.7%)

¹ n (%)

data_fungi_mini \|\>
[filt_taxa_pq](https://adrientaudiere.github.io/MiscMetabar/reference/filt_taxa_pq.md)(min_occurence
= 2) \|\> tbl_sum_taxtable(taxonomic_rank =
[c](https://rdrr.io/r/base/c.html)("Species", "Genus")) \#\> Cleaning
suppress 0 taxa ( ) and 0 sample(s) ( ). \#\> Number of non-matching ASV
0 \#\> Number of matching ASV 45 \#\> Number of filtered-out ASV 2 \#\>
Number of kept ASV 43 \#\> Number of kept samples 137

[TABLE]
