# Assign taxonomy using blastn algorithm and the blast software

[![lifecycle-experimental](https://img.shields.io/badge/lifecycle-experimental-orange)](https://adrientaudiere.github.io/MiscMetabar/articles/Rules.html#lifecycle)

Use the blast software.

## Usage

``` r
assign_blastn(
  physeq = NULL,
  ref_fasta = NULL,
  database = NULL,
  blastpath = NULL,
  seq2search = NULL,
  behavior = c("return_matrix", "return_taxtab", "add_to_phyloseq"),
  method_algo = c("vote", "top-hit"),
  suffix = "_blastn",
  min_id = 95,
  min_bit_score = 50,
  min_cover = 95,
  min_e_value = 1e-30,
  nb_voting = NULL,
  column_names = c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species"),
  vote_algorithm = c("consensus", "rel_majority", "abs_majority", "unanimity"),
  strict = FALSE,
  nb_agree_threshold = 1,
  preference_index = NULL,
  collapse_string = "/",
  replace_collapsed_rank_by_NA = TRUE,
  simplify_taxo = TRUE,
  keep_blast_metrics = FALSE,
  blast_table = NULL,
  ...
)
```

## Arguments

- physeq:

  (optional) A phyloseq object with a `@refseq` slot. If `seq2search` is
  provided, `physeq` can be `NULL`.

- ref_fasta:

  Either a DNAStringSet object or a path to a fasta file to make the
  blast database. It must be in sintax format. See
  [`assign_sintax()`](https://adrientaudiere.github.io/MiscMetabar/dev/reference/assign_sintax.md).

- database:

  path to a blast database. Only used if ref_fasta is not set.

- blastpath:

  path to blast program.

- seq2search:

  A DNAStringSet object of sequences to search for, or a matrix whose
  `colnames` are the sequences to search for (e.g. a dada2 sequence
  table where columns are ASV sequences). In the latter case the
  `colnames` are converted to a DNAStringSet and used as the taxa names.
  Replace the physeq object.

- behavior:

  Either "return_matrix" (default), "return_taxtab", or
  "add_to_phyloseq":

  - "return_matrix" return a list of two matrix with taxonomic value in
    the first element of the list and bootstrap value in the second one.

  - "return_taxtab" return a character matrix of taxonomic values (rows
    are taxa, columns are the `column_names`) with rownames set to the
    taxa names. This is convenient as a direct input to
    [`phyloseq::tax_table()`](https://rdrr.io/pkg/phyloseq/man/tax_table-methods.html).

  - "add_to_phyloseq" return a phyloseq object with amended slot
    `@taxtable`. Only available if using physeq input and not seq2search
    input.

- method_algo:

  (One of "vote" or "top-hit"). If top-hit, only the better match is
  used to assign taxonomy. If vote, the algorithm takes all (or
  `nb_voting` if `nb_voting` is not null) select assignation and resolve
  the conflict using the function
  [`resolve_vector_ranks()`](https://adrientaudiere.github.io/MiscMetabar/dev/reference/resolve_vector_ranks.md).

- suffix:

  (character) The suffix to name the new columns. If set to "" (the
  default), the taxa_ranks algorithm is used without suffix.

- min_id:

  (default: 95) the identity percent to take into account a references
  taxa

- min_bit_score:

  (default: 50) the minimum bit score to take into account a references
  taxa

- min_cover:

  (default: 95) cut of in query cover (%) to keep result

- min_e_value:

  (default: 1e-30) cut of in e-value (%) to keep result The BLAST
  E-value is the number of expected hits of similar quality (score) that
  could be found just by chance.

- nb_voting:

  (Int, default NULL). The number of taxa to keep before apply a vote to
  resolve conflict. If NULL all taxa passing the filters (min_id,
  min_bit_score, min_cover and min_e_value) are selected.

- column_names:

  A vector of names for taxonomic ranks. Must correspond to names in the
  ref_fasta files.

- vote_algorithm:

  the method to vote among "consensus", "rel_majority", "abs_majority"
  and "unanimity". See
  [`resolve_vector_ranks()`](https://adrientaudiere.github.io/MiscMetabar/dev/reference/resolve_vector_ranks.md)
  for more details.

- strict:

  (Logical, default FALSE). See
  [`resolve_vector_ranks()`](https://adrientaudiere.github.io/MiscMetabar/dev/reference/resolve_vector_ranks.md)
  for more details.

- nb_agree_threshold:

  See
  [`resolve_vector_ranks()`](https://adrientaudiere.github.io/MiscMetabar/dev/reference/resolve_vector_ranks.md)
  for more details.

- preference_index:

  See
  [`resolve_vector_ranks()`](https://adrientaudiere.github.io/MiscMetabar/dev/reference/resolve_vector_ranks.md)
  for more details.

- collapse_string:

  See
  [`resolve_vector_ranks()`](https://adrientaudiere.github.io/MiscMetabar/dev/reference/resolve_vector_ranks.md)
  for more details.

- replace_collapsed_rank_by_NA:

  (Logical, default TRUE) See
  [`resolve_vector_ranks()`](https://adrientaudiere.github.io/MiscMetabar/dev/reference/resolve_vector_ranks.md)
  for more details.

- simplify_taxo:

  (logical default TRUE). Do we apply the function
  [`simplify_taxo()`](https://adrientaudiere.github.io/MiscMetabar/dev/reference/simplify_taxo.md)
  to the phyloseq object?

- keep_blast_metrics:

  (Logical, default FALSE). If TRUE, the blast metrics ("Query seq.
  length", "Taxa seq. length", "Alignment length", "% id. match",
  "e-value", "bit score" and "Query cover") are stored in the tax_table.

- blast_table:

  (data.frame, default NULL). A raw blast table as returned by
  [`blast_pq()`](https://adrientaudiere.github.io/MiscMetabar/dev/reference/blast_pq.md)
  with `unique_per_seq = FALSE` and `score_filter = FALSE`. If set,
  blastn is not run again: the score filters (`min_id`, `min_bit_score`,
  `min_cover` and `min_e_value`) and the vote (or top hit) are applied
  to this table, so one search can serve several `vote_algorithm` or
  filter values. `ref_fasta` and `database` are then ignored.

- ...:

  Additional arguments passed on to
  [`blast_pq()`](https://adrientaudiere.github.io/MiscMetabar/dev/reference/blast_pq.md)

## Value

- If behavior == "return_matrix" :

  - If method_algo = "top-hit" a matrix of taxonomic assignation

  - If method_algo = "vote", a list of two matrix, the first is the raw
    taxonomic assignation (before vote). The second one is the taxonomic
    assignation in which conflicts are resolved using vote.

- If behavior == "return_taxtab", a character matrix of taxonomic values
  ready for
  [`phyloseq::tax_table()`](https://rdrr.io/pkg/phyloseq/man/tax_table-methods.html)

- If behavior == "add_to_phyloseq", return a new phyloseq object

## Author

Adrien Taudière

## Examples

``` r
ref_fasta <- Biostrings::readDNAStringSet(system.file("extdata",
  "mini_UNITE_fungi.fasta.gz",
  package = "MiscMetabar", mustWork = TRUE
))

mat <- assign_blastn(data_fungi_mini,
  ref_fasta = ref_fasta,
  method_algo = "top-hit", min_id = 70, min_e_value = 1e-3, min_cover = 50,
  min_bit_score = 20
)
#> Build the database from fasta_for_db
#> Blast refseq from physeq object against the database
head(mat)
#>       taxa_names Kingdom_blastn   Phylum_blastn     Class_blastn Order_blastn
#> 5501       ASV41        k:Fungi p:Basidiomycota c:Agaricomycetes o:Agaricales
#> 15001      ASV75        k:Fungi p:Basidiomycota c:Agaricomycetes o:Russulales
#> 2501       ASV26        k:Fungi p:Basidiomycota c:Agaricomycetes o:Russulales
#> 2001       ASV25        k:Fungi p:Basidiomycota c:Agaricomycetes o:Agaricales
#> 13501      ASV68        k:Fungi p:Basidiomycota c:Agaricomycetes o:Agaricales
#> 6001       ASV42        k:Fungi p:Basidiomycota c:Agaricomycetes o:Agaricales
#>            Family_blastn                    Genus_blastn      Species_blastn
#> 5501        f:Mycenaceae g:Mycenaceae_gen_Incertae_sedis     s:Mycenaceae_sp
#> 15001   f:Peniophoraceae                    g:Peniophora s:Peniophora_reidii
#> 2501        f:Stereaceae                     g:Xylobolus      s:Xylobolus_sp
#> 2001  f:Tricholomataceae                    g:Tricholoma     s:Tricholoma_sp
#> 13501 f:Tricholomataceae                    g:Tricholoma     s:Tricholoma_sp
#> 6001  f:Tricholomataceae                    g:Tricholoma     s:Tricholoma_sp
#>                                             Taxa_name_db_blastn
#> 5501      Mycenaceae_sp_OP470041_SH0013594.10FU_reps_singleton_
#> 15001 Peniophora_reidii_MH861616_SH0003674.10FU_reps_singleton_
#> 2501       Xylobolus_sp_MH121216_SH0004343.10FU_reps_singleton_
#> 2001     Tricholoma_sp_UDB007415_SH0014274.10FU_reps_singleton_
#> 13501    Tricholoma_sp_UDB007415_SH0014274.10FU_reps_singleton_
#> 6001     Tricholoma_sp_UDB007415_SH0014274.10FU_reps_singleton_

assign_blastn(data_fungi_mini,
  ref_fasta = ref_fasta, method_algo = "vote",
  vote_algorithm = "rel_majority", min_id = 90, min_cover = 50,
  behavior = "add_to_phyloseq"
)@tax_table
#> Build the database from fasta_for_db
#> Blast refseq from physeq object against the database
#> Taxonomy Table:     [45 taxa by 20 taxonomic ranks]:
#>        Domain  Phylum          Class                Order            
#> ASV7   "Fungi" "Basidiomycota" "Agaricomycetes"     "Russulales"     
#> ASV8   "Fungi" "Basidiomycota" "Agaricomycetes"     "Russulales"     
#> ASV12  "Fungi" "Basidiomycota" "Agaricomycetes"     "Hymenochaetales"
#> ASV18  "Fungi" "Basidiomycota" "Agaricomycetes"     "Russulales"     
#> ASV25  "Fungi" "Basidiomycota" "Agaricomycetes"     "Agaricales"     
#> ASV26  "Fungi" "Basidiomycota" "Agaricomycetes"     "Russulales"     
#> ASV27  "Fungi" "Basidiomycota" "Agaricomycetes"     "Polyporales"    
#> ASV29  "Fungi" "Basidiomycota" "Agaricomycetes"     "Auriculariales" 
#> ASV32  "Fungi" "Basidiomycota" "Agaricomycetes"     "Cantharellales" 
#> ASV34  "Fungi" "Basidiomycota" "Agaricomycetes"     "Agaricales"     
#> ASV35  "Fungi" "Basidiomycota" "Agaricomycetes"     "Polyporales"    
#> ASV41  "Fungi" "Basidiomycota" "Agaricomycetes"     "Agaricales"     
#> ASV42  "Fungi" "Basidiomycota" "Agaricomycetes"     "Agaricales"     
#> ASV46  "Fungi" "Basidiomycota" "Atractiellomycetes" "Atractiellales" 
#> ASV47  "Fungi" "Basidiomycota" "Agaricomycetes"     "Agaricales"     
#> ASV48  "Fungi" "Basidiomycota" "Agaricomycetes"     "Auriculariales" 
#> ASV49  "Fungi" "Basidiomycota" "Agaricomycetes"     "Polyporales"    
#> ASV50  "Fungi" "Basidiomycota" "Agaricomycetes"     "Russulales"     
#> ASV53  "Fungi" "Basidiomycota" "Agaricomycetes"     "Polyporales"    
#> ASV54  "Fungi" "Basidiomycota" "Agaricomycetes"     "Auriculariales" 
#> ASV58  "Fungi" "Basidiomycota" "Agaricomycetes"     "Polyporales"    
#> ASV59  "Fungi" "Basidiomycota" "Agaricomycetes"     "Polyporales"    
#> ASV61  "Fungi" "Basidiomycota" "Agaricomycetes"     "Polyporales"    
#> ASV62  "Fungi" "Basidiomycota" "Agaricomycetes"     NA               
#> ASV63  "Fungi" "Basidiomycota" "Agaricomycetes"     "Polyporales"    
#> ASV64  "Fungi" "Basidiomycota" "Agaricomycetes"     "Polyporales"    
#> ASV67  "Fungi" "Basidiomycota" "Agaricomycetes"     "Hymenochaetales"
#> ASV68  "Fungi" "Basidiomycota" "Agaricomycetes"     "Agaricales"     
#> ASV71  "Fungi" "Basidiomycota" "Tremellomycetes"    NA               
#> ASV72  "Fungi" "Basidiomycota" "Agaricomycetes"     "Polyporales"    
#> ASV75  "Fungi" "Basidiomycota" "Agaricomycetes"     "Russulales"     
#> ASV77  "Fungi" "Basidiomycota" "Agaricomycetes"     "Agaricales"     
#> ASV82  "Fungi" "Basidiomycota" "Agaricomycetes"     "Auriculariales" 
#> ASV83  "Fungi" "Basidiomycota" "Tremellomycetes"    "Tremellales"    
#> ASV85  "Fungi" "Basidiomycota" "Agaricomycetes"     "Hymenochaetales"
#> ASV91  "Fungi" "Basidiomycota" "Agaricomycetes"     "Auriculariales" 
#> ASV93  "Fungi" "Basidiomycota" "Agaricomycetes"     "Russulales"     
#> ASV94  "Fungi" "Basidiomycota" "Agaricomycetes"     "Russulales"     
#> ASV99  "Fungi" "Basidiomycota" "Agaricomycetes"     "Polyporales"    
#> ASV100 "Fungi" "Basidiomycota" NA                   NA               
#> ASV101 "Fungi" "Basidiomycota" "Agaricomycetes"     "Corticiales"    
#> ASV104 "Fungi" "Basidiomycota" "Agaricomycetes"     "Russulales"     
#> ASV105 "Fungi" "Basidiomycota" "Agaricomycetes"     "Hymenochaetales"
#> ASV107 "Fungi" "Basidiomycota" "Agaricomycetes"     "Hymenochaetales"
#> ASV108 "Fungi" "Basidiomycota" "Agaricomycetes"     "Auriculariales" 
#>        Family                               Genus               Species        
#> ASV7   "Stereaceae"                         NA                  NA             
#> ASV8   "Stereaceae"                         "Stereum"           "ostrea"       
#> ASV12  "Schizoporaceae"                     "Xylodon"           "raduloides"   
#> ASV18  "Stereaceae"                         "Stereum"           "ostrea"       
#> ASV25  "Lyophyllaceae"                      "Ossicaulis"        "lachnopus"    
#> ASV26  "Stereaceae"                         "Stereum"           "hirsutum"     
#> ASV27  "Steccherinaceae"                    "Antrodiella"       "brasiliensis" 
#> ASV29  "Exidiaceae"                         "Basidiodendron"    "eyrei"        
#> ASV32  "Cantharellales_fam_Incertae_sedis"  "Sistotrema"        "oblongisporum"
#> ASV34  "Entolomataceae"                     "Entocybe"          NA             
#> ASV35  "Polyporaceae"                       "Fomes"             "fomentarius"  
#> ASV41  "Tricholomataceae"                   "Mycena"            "renati"       
#> ASV42  "Lyophyllaceae"                      "Ossicaulis"        "lachnopus"    
#> ASV46  "Atractiellales_fam_Incertae_sedis"  "Helicogloea"       "pellucida"    
#> ASV47  "Pterulaceae"                        "Radulomyces"       "molaris"      
#> ASV48  "Aporpiaceae"                        "Elmerina"          "caryae"       
#> ASV49  "Phanerochaetaceae"                  "Phanerochaete"     "livescens"    
#> ASV50  "Russulales_fam_Incertae_sedis"      "Gloeohypochnicium" "analogum"     
#> ASV53  "Polyporaceae"                       "Fomes"             "fomentarius"  
#> ASV54  "Auriculariaceae"                    "Auricularia"       NA             
#> ASV58  "Polyporaceae"                       "Fomes"             "fomentarius"  
#> ASV59  "Hyphodermataceae"                   "Hyphoderma"        "roseocremeum" 
#> ASV61  "Hyphodermataceae"                   "Hyphoderma"        "setigerum"    
#> ASV62  NA                                   NA                  NA             
#> ASV63  NA                                   NA                  NA             
#> ASV64  "Polyporaceae"                       "Trametes"          "versicolor"   
#> ASV67  "Schizoporaceae"                     "Xylodon"           "raduloides"   
#> ASV68  "Lyophyllaceae"                      "Ossicaulis"        "lachnopus"    
#> ASV71  NA                                   NA                  NA             
#> ASV72  NA                                   NA                  NA             
#> ASV75  "Peniophoraceae"                     "Peniophora"        "versiformis"  
#> ASV77  "Lyophyllaceae"                      "Ossicaulis"        "lachnopus"    
#> ASV82  "Exidiaceae"                         "Exidia"            "glandulosa"   
#> ASV83  NA                                   NA                  NA             
#> ASV85  "Hymenochaetales_fam_Incertae_sedis" "Peniophorella"     "pubera"       
#> ASV91  "Auriculariaceae"                    "Auricularia"       "mesenterica"  
#> ASV93  "Stereaceae"                         NA                  NA             
#> ASV94  "Stereaceae"                         "Stereum"           "ostrea"       
#> ASV99  "Polyporaceae"                       "Fomes"             "fomentarius"  
#> ASV100 NA                                   NA                  NA             
#> ASV101 "Corticiaceae"                       "Marchandiomyces"   "buckii"       
#> ASV104 "Hericiaceae"                        "Hericium"          "coralloides"  
#> ASV105 "Schizoporaceae"                     "Xylodon"           "flaviporus"   
#> ASV107 "Schizoporaceae"                     "Xylodon"           "raduloides"   
#> ASV108 "Exidiaceae"                         "Exidia"            "glandulosa"   
#>        Trophic.Mode                       
#> ASV7   "Saprotroph"                       
#> ASV8   "Saprotroph"                       
#> ASV12  "Saprotroph"                       
#> ASV18  "Saprotroph"                       
#> ASV25  "Saprotroph"                       
#> ASV26  "Saprotroph"                       
#> ASV27  "Saprotroph"                       
#> ASV29  "Saprotroph"                       
#> ASV32  "Saprotroph-Symbiotroph"           
#> ASV34  "Saprotroph"                       
#> ASV35  "Saprotroph"                       
#> ASV41  "Pathotroph-Saprotroph"            
#> ASV42  "Saprotroph"                       
#> ASV46  "Saprotroph"                       
#> ASV47  "Saprotroph"                       
#> ASV48  "Saprotroph"                       
#> ASV49  "Saprotroph"                       
#> ASV50  "Saprotroph"                       
#> ASV53  "Saprotroph"                       
#> ASV54  "Saprotroph"                       
#> ASV58  "Saprotroph"                       
#> ASV59  "Saprotroph"                       
#> ASV61  "Saprotroph"                       
#> ASV62  "-"                                
#> ASV63  "-"                                
#> ASV64  "Saprotroph"                       
#> ASV67  "Saprotroph"                       
#> ASV68  "Saprotroph"                       
#> ASV71  "-"                                
#> ASV72  "-"                                
#> ASV75  "Pathotroph-Saprotroph"            
#> ASV77  "Saprotroph"                       
#> ASV82  "Saprotroph-Symbiotroph"           
#> ASV83  "Pathotroph-Saprotroph-Symbiotroph"
#> ASV85  "Saprotroph"                       
#> ASV91  "Saprotroph"                       
#> ASV93  "Saprotroph"                       
#> ASV94  "Saprotroph"                       
#> ASV99  "Saprotroph"                       
#> ASV100 "-"                                
#> ASV101 "Pathotroph"                       
#> ASV104 "Saprotroph"                       
#> ASV105 "Saprotroph"                       
#> ASV107 "Saprotroph"                       
#> ASV108 "Saprotroph-Symbiotroph"           
#>        Guild                                                            
#> ASV7   "WoodSaprotroph-UndefinedSaprotroph"                             
#> ASV8   "UndefinedSaprotroph"                                            
#> ASV12  "UndefinedSaprotroph"                                            
#> ASV18  "UndefinedSaprotroph"                                            
#> ASV25  "WoodSaprotroph"                                                 
#> ASV26  "UndefinedSaprotroph"                                            
#> ASV27  "WoodSaprotroph"                                                 
#> ASV29  "UndefinedSaprotroph"                                            
#> ASV32  "Ectomycorrhizal-WoodSaprotroph"                                 
#> ASV34  "UndefinedSaprotroph"                                            
#> ASV35  "WoodSaprotroph"                                                 
#> ASV41  "LeafSaprotroph-PlantPathogen-UndefinedSaprotroph-WoodSaprotroph"
#> ASV42  "WoodSaprotroph"                                                 
#> ASV46  "UndefinedSaprotroph"                                            
#> ASV47  "UndefinedSaprotroph"                                            
#> ASV48  "UndefinedSaprotroph"                                            
#> ASV49  "WoodSaprotroph"                                                 
#> ASV50  "UndefinedSaprotroph"                                            
#> ASV53  "WoodSaprotroph"                                                 
#> ASV54  "UndefinedSaprotroph"                                            
#> ASV58  "WoodSaprotroph"                                                 
#> ASV59  "UndefinedSaprotroph"                                            
#> ASV61  "UndefinedSaprotroph"                                            
#> ASV62  "-"                                                              
#> ASV63  "-"                                                              
#> ASV64  "WoodSaprotroph"                                                 
#> ASV67  "UndefinedSaprotroph"                                            
#> ASV68  "WoodSaprotroph"                                                 
#> ASV71  "-"                                                              
#> ASV72  "-"                                                              
#> ASV75  "PlantPathogen-WoodSaprotroph"                                   
#> ASV77  "WoodSaprotroph"                                                 
#> ASV82  "Endophyte-UndefinedSaprotroph"                                  
#> ASV83  "FungalParasite-UndefinedSaprotroph"                             
#> ASV85  "UndefinedSaprotroph"                                            
#> ASV91  "UndefinedSaprotroph"                                            
#> ASV93  "WoodSaprotroph-UndefinedSaprotroph"                             
#> ASV94  "UndefinedSaprotroph"                                            
#> ASV99  "WoodSaprotroph"                                                 
#> ASV100 "-"                                                              
#> ASV101 "LichenParasite"                                                 
#> ASV104 "UndefinedSaprotroph"                                            
#> ASV105 "UndefinedSaprotroph"                                            
#> ASV107 "UndefinedSaprotroph"                                            
#> ASV108 "Endophyte-UndefinedSaprotroph"                                  
#>        Trait               Confidence.Ranking Genus_species               
#> ASV7   "NULL"              "Probable"         "_"                         
#> ASV8   "WhiteRot"          "Probable"         "Stereum_ostrea"            
#> ASV12  "WhiteRot"          "Probable"         "Xylodon_raduloides"        
#> ASV18  "WhiteRot"          "Probable"         "Stereum_ostrea"            
#> ASV25  "BrownRot"          "Probable"         "Ossicaulis_lachnopus"      
#> ASV26  "WhiteRot"          "Probable"         "Stereum_hirsutum"          
#> ASV27  "WhiteRot"          "HighlyProbable"   "Antrodiella_brasiliensis"  
#> ASV29  "NULL"              "Probable"         "Basidiodendron_eyrei"      
#> ASV32  "WhiteRot"          "Possible"         "Sistotrema_oblongisporum"  
#> ASV34  "NULL"              "Probable"         "Entocybe_"                 
#> ASV35  "BrownRot-WhiteRot" "Probable"         "Fomes_fomentarius"         
#> ASV41  "NULL"              "Probable"         "Mycena_renati"             
#> ASV42  "BrownRot"          "Probable"         "Ossicaulis_lachnopus"      
#> ASV46  "NULL"              "Probable"         "Helicogloea_pellucida"     
#> ASV47  "WhiteRot"          "Probable"         "Radulomyces_molaris"       
#> ASV48  "NULL"              "Probable"         "Elmerina_caryae"           
#> ASV49  "WhiteRot"          "HighlyProbable"   "Phanerochaete_livescens"   
#> ASV50  "WhiteRot"          "Probable"         "Gloeohypochnicium_analogum"
#> ASV53  "BrownRot-WhiteRot" "Probable"         "Fomes_fomentarius"         
#> ASV54  "NULL"              "Probable"         "Auricularia_"              
#> ASV58  "BrownRot-WhiteRot" "Probable"         "Fomes_fomentarius"         
#> ASV59  "WhiteRot"          "Probable"         "Hyphoderma_roseocremeum"   
#> ASV61  "WhiteRot"          "Probable"         "Hyphoderma_setigerum"      
#> ASV62  "-"                 "-"                "_"                         
#> ASV63  "-"                 "-"                "_"                         
#> ASV64  "WhiteRot"          "HighlyProbable"   "Trametes_versicolor"       
#> ASV67  "WhiteRot"          "Probable"         "Xylodon_raduloides"        
#> ASV68  "BrownRot"          "Probable"         "Ossicaulis_lachnopus"      
#> ASV71  "-"                 "-"                "_"                         
#> ASV72  "-"                 "-"                "_"                         
#> ASV75  "WhiteRot"          "Probable"         "Peniophora_versiformis"    
#> ASV77  "BrownRot"          "Probable"         "Ossicaulis_lachnopus"      
#> ASV82  "NULL"              "Probable"         "Exidia_glandulosa"         
#> ASV83  "NULL"              "Possible"         "_"                         
#> ASV85  "WhiteRot"          "Probable"         "Peniophorella_pubera"      
#> ASV91  "NULL"              "Probable"         "Auricularia_mesenterica"   
#> ASV93  "NULL"              "Probable"         "_"                         
#> ASV94  "WhiteRot"          "Probable"         "Stereum_ostrea"            
#> ASV99  "BrownRot-WhiteRot" "Probable"         "Fomes_fomentarius"         
#> ASV100 "-"                 "-"                "_"                         
#> ASV101 "NULL"              "Probable"         "Marchandiomyces_buckii"    
#> ASV104 "WhiteRot"          "Probable"         "Hericium_coralloides"      
#> ASV105 "WhiteRot"          "Probable"         "Xylodon_flaviporus"        
#> ASV107 "WhiteRot"          "Probable"         "Xylodon_raduloides"        
#> ASV108 "NULL"              "Probable"         "Exidia_glandulosa"         
#>        Kingdom_blastn Phylum_blastn   Class_blastn      Order_blastn 
#> ASV7   NA             NA              NA                NA           
#> ASV8   NA             NA              NA                NA           
#> ASV12  NA             NA              NA                NA           
#> ASV18  NA             NA              NA                NA           
#> ASV25  NA             NA              NA                NA           
#> ASV26  NA             NA              NA                NA           
#> ASV27  NA             NA              NA                NA           
#> ASV29  NA             NA              NA                NA           
#> ASV32  NA             NA              NA                NA           
#> ASV34  NA             NA              NA                NA           
#> ASV35  NA             NA              NA                NA           
#> ASV41  "Fungi"        "Basidiomycota" "Agaricomycetes"  "Agaricales" 
#> ASV42  NA             NA              NA                NA           
#> ASV46  NA             NA              NA                NA           
#> ASV47  NA             NA              NA                NA           
#> ASV48  NA             NA              NA                NA           
#> ASV49  NA             NA              NA                NA           
#> ASV50  NA             NA              NA                NA           
#> ASV53  NA             NA              NA                NA           
#> ASV54  NA             NA              NA                NA           
#> ASV58  NA             NA              NA                NA           
#> ASV59  NA             NA              NA                NA           
#> ASV61  NA             NA              NA                NA           
#> ASV62  NA             NA              NA                NA           
#> ASV63  NA             NA              NA                NA           
#> ASV64  NA             NA              NA                NA           
#> ASV67  NA             NA              NA                NA           
#> ASV68  NA             NA              NA                NA           
#> ASV71  "Fungi"        "Basidiomycota" "Tremellomycetes" NA           
#> ASV72  NA             NA              NA                NA           
#> ASV75  "Fungi"        "Basidiomycota" "Agaricomycetes"  "Russulales" 
#> ASV77  NA             NA              NA                NA           
#> ASV82  NA             NA              NA                NA           
#> ASV83  NA             NA              NA                NA           
#> ASV85  NA             NA              NA                NA           
#> ASV91  NA             NA              NA                NA           
#> ASV93  NA             NA              NA                NA           
#> ASV94  NA             NA              NA                NA           
#> ASV99  NA             NA              NA                NA           
#> ASV100 NA             NA              NA                NA           
#> ASV101 "Fungi"        "Basidiomycota" "Agaricomycetes"  "Corticiales"
#> ASV104 NA             NA              NA                NA           
#> ASV105 NA             NA              NA                NA           
#> ASV107 NA             NA              NA                NA           
#> ASV108 NA             NA              NA                NA           
#>        Family_blastn                    Genus_blastn                    
#> ASV7   NA                               NA                              
#> ASV8   NA                               NA                              
#> ASV12  NA                               NA                              
#> ASV18  NA                               NA                              
#> ASV25  NA                               NA                              
#> ASV26  NA                               NA                              
#> ASV27  NA                               NA                              
#> ASV29  NA                               NA                              
#> ASV32  NA                               NA                              
#> ASV34  NA                               NA                              
#> ASV35  NA                               NA                              
#> ASV41  "Mycenaceae"                     "Mycena"                        
#> ASV42  NA                               NA                              
#> ASV46  NA                               NA                              
#> ASV47  NA                               NA                              
#> ASV48  NA                               NA                              
#> ASV49  NA                               NA                              
#> ASV50  NA                               NA                              
#> ASV53  NA                               NA                              
#> ASV54  NA                               NA                              
#> ASV58  NA                               NA                              
#> ASV59  NA                               NA                              
#> ASV61  NA                               NA                              
#> ASV62  NA                               NA                              
#> ASV63  NA                               NA                              
#> ASV64  NA                               NA                              
#> ASV67  NA                               NA                              
#> ASV68  NA                               NA                              
#> ASV71  NA                               NA                              
#> ASV72  NA                               NA                              
#> ASV75  "Peniophoraceae"                 "Peniophora"                    
#> ASV77  NA                               NA                              
#> ASV82  NA                               NA                              
#> ASV83  NA                               NA                              
#> ASV85  NA                               NA                              
#> ASV91  NA                               NA                              
#> ASV93  NA                               NA                              
#> ASV94  NA                               NA                              
#> ASV99  NA                               NA                              
#> ASV100 NA                               NA                              
#> ASV101 "Corticiales_fam_Incertae_sedis" "Corticiales_gen_Incertae_sedis"
#> ASV104 NA                               NA                              
#> ASV105 NA                               NA                              
#> ASV107 NA                               NA                              
#> ASV108 NA                               NA                              
#>        Species_blastn  
#> ASV7   NA              
#> ASV8   NA              
#> ASV12  NA              
#> ASV18  NA              
#> ASV25  NA              
#> ASV26  NA              
#> ASV27  NA              
#> ASV29  NA              
#> ASV32  NA              
#> ASV34  NA              
#> ASV35  NA              
#> ASV41  "Mycena_sp"     
#> ASV42  NA              
#> ASV46  NA              
#> ASV47  NA              
#> ASV48  NA              
#> ASV49  NA              
#> ASV50  NA              
#> ASV53  NA              
#> ASV54  NA              
#> ASV58  NA              
#> ASV59  NA              
#> ASV61  NA              
#> ASV62  NA              
#> ASV63  NA              
#> ASV64  NA              
#> ASV67  NA              
#> ASV68  NA              
#> ASV71  NA              
#> ASV72  NA              
#> ASV75  NA              
#> ASV77  NA              
#> ASV82  NA              
#> ASV83  NA              
#> ASV85  NA              
#> ASV91  NA              
#> ASV93  NA              
#> ASV94  NA              
#> ASV99  NA              
#> ASV100 NA              
#> ASV101 "Corticiales_sp"
#> ASV104 NA              
#> ASV105 NA              
#> ASV107 NA              
#> ASV108 NA              
#>        Taxa_name_db_blastn                                     
#> ASV7   NA                                                      
#> ASV8   NA                                                      
#> ASV12  NA                                                      
#> ASV18  NA                                                      
#> ASV25  NA                                                      
#> ASV26  NA                                                      
#> ASV27  NA                                                      
#> ASV29  NA                                                      
#> ASV32  NA                                                      
#> ASV34  NA                                                      
#> ASV35  NA                                                      
#> ASV41  NA                                                      
#> ASV42  NA                                                      
#> ASV46  NA                                                      
#> ASV47  NA                                                      
#> ASV48  NA                                                      
#> ASV49  NA                                                      
#> ASV50  NA                                                      
#> ASV53  NA                                                      
#> ASV54  NA                                                      
#> ASV58  NA                                                      
#> ASV59  NA                                                      
#> ASV61  NA                                                      
#> ASV62  NA                                                      
#> ASV63  NA                                                      
#> ASV64  NA                                                      
#> ASV67  NA                                                      
#> ASV68  NA                                                      
#> ASV71  NA                                                      
#> ASV72  NA                                                      
#> ASV75  NA                                                      
#> ASV77  NA                                                      
#> ASV82  NA                                                      
#> ASV83  NA                                                      
#> ASV85  NA                                                      
#> ASV91  NA                                                      
#> ASV93  NA                                                      
#> ASV94  NA                                                      
#> ASV99  NA                                                      
#> ASV100 NA                                                      
#> ASV101 "Corticiales_sp_KY462668_SH0001937.10FU_reps_singleton_"
#> ASV104 NA                                                      
#> ASV105 NA                                                      
#> ASV107 NA                                                      
#> ASV108 NA                                                      

assign_blastn(data_fungi_mini,
  ref_fasta = ref_fasta, method_algo = "vote",
  vote_algorithm = "consensus", replace_collapsed_rank_by_NA = FALSE,
  min_id = 90, min_cover = 50, behavior = "add_to_phyloseq"
)@tax_table
#> Build the database from fasta_for_db
#> Blast refseq from physeq object against the database
#> Taxonomy Table:     [45 taxa by 20 taxonomic ranks]:
#>        Domain  Phylum          Class                Order            
#> ASV7   "Fungi" "Basidiomycota" "Agaricomycetes"     "Russulales"     
#> ASV8   "Fungi" "Basidiomycota" "Agaricomycetes"     "Russulales"     
#> ASV12  "Fungi" "Basidiomycota" "Agaricomycetes"     "Hymenochaetales"
#> ASV18  "Fungi" "Basidiomycota" "Agaricomycetes"     "Russulales"     
#> ASV25  "Fungi" "Basidiomycota" "Agaricomycetes"     "Agaricales"     
#> ASV26  "Fungi" "Basidiomycota" "Agaricomycetes"     "Russulales"     
#> ASV27  "Fungi" "Basidiomycota" "Agaricomycetes"     "Polyporales"    
#> ASV29  "Fungi" "Basidiomycota" "Agaricomycetes"     "Auriculariales" 
#> ASV32  "Fungi" "Basidiomycota" "Agaricomycetes"     "Cantharellales" 
#> ASV34  "Fungi" "Basidiomycota" "Agaricomycetes"     "Agaricales"     
#> ASV35  "Fungi" "Basidiomycota" "Agaricomycetes"     "Polyporales"    
#> ASV41  "Fungi" "Basidiomycota" "Agaricomycetes"     "Agaricales"     
#> ASV42  "Fungi" "Basidiomycota" "Agaricomycetes"     "Agaricales"     
#> ASV46  "Fungi" "Basidiomycota" "Atractiellomycetes" "Atractiellales" 
#> ASV47  "Fungi" "Basidiomycota" "Agaricomycetes"     "Agaricales"     
#> ASV48  "Fungi" "Basidiomycota" "Agaricomycetes"     "Auriculariales" 
#> ASV49  "Fungi" "Basidiomycota" "Agaricomycetes"     "Polyporales"    
#> ASV50  "Fungi" "Basidiomycota" "Agaricomycetes"     "Russulales"     
#> ASV53  "Fungi" "Basidiomycota" "Agaricomycetes"     "Polyporales"    
#> ASV54  "Fungi" "Basidiomycota" "Agaricomycetes"     "Auriculariales" 
#> ASV58  "Fungi" "Basidiomycota" "Agaricomycetes"     "Polyporales"    
#> ASV59  "Fungi" "Basidiomycota" "Agaricomycetes"     "Polyporales"    
#> ASV61  "Fungi" "Basidiomycota" "Agaricomycetes"     "Polyporales"    
#> ASV62  "Fungi" "Basidiomycota" "Agaricomycetes"     NA               
#> ASV63  "Fungi" "Basidiomycota" "Agaricomycetes"     "Polyporales"    
#> ASV64  "Fungi" "Basidiomycota" "Agaricomycetes"     "Polyporales"    
#> ASV67  "Fungi" "Basidiomycota" "Agaricomycetes"     "Hymenochaetales"
#> ASV68  "Fungi" "Basidiomycota" "Agaricomycetes"     "Agaricales"     
#> ASV71  "Fungi" "Basidiomycota" "Tremellomycetes"    NA               
#> ASV72  "Fungi" "Basidiomycota" "Agaricomycetes"     "Polyporales"    
#> ASV75  "Fungi" "Basidiomycota" "Agaricomycetes"     "Russulales"     
#> ASV77  "Fungi" "Basidiomycota" "Agaricomycetes"     "Agaricales"     
#> ASV82  "Fungi" "Basidiomycota" "Agaricomycetes"     "Auriculariales" 
#> ASV83  "Fungi" "Basidiomycota" "Tremellomycetes"    "Tremellales"    
#> ASV85  "Fungi" "Basidiomycota" "Agaricomycetes"     "Hymenochaetales"
#> ASV91  "Fungi" "Basidiomycota" "Agaricomycetes"     "Auriculariales" 
#> ASV93  "Fungi" "Basidiomycota" "Agaricomycetes"     "Russulales"     
#> ASV94  "Fungi" "Basidiomycota" "Agaricomycetes"     "Russulales"     
#> ASV99  "Fungi" "Basidiomycota" "Agaricomycetes"     "Polyporales"    
#> ASV100 "Fungi" "Basidiomycota" NA                   NA               
#> ASV101 "Fungi" "Basidiomycota" "Agaricomycetes"     "Corticiales"    
#> ASV104 "Fungi" "Basidiomycota" "Agaricomycetes"     "Russulales"     
#> ASV105 "Fungi" "Basidiomycota" "Agaricomycetes"     "Hymenochaetales"
#> ASV107 "Fungi" "Basidiomycota" "Agaricomycetes"     "Hymenochaetales"
#> ASV108 "Fungi" "Basidiomycota" "Agaricomycetes"     "Auriculariales" 
#>        Family                               Genus               Species        
#> ASV7   "Stereaceae"                         NA                  NA             
#> ASV8   "Stereaceae"                         "Stereum"           "ostrea"       
#> ASV12  "Schizoporaceae"                     "Xylodon"           "raduloides"   
#> ASV18  "Stereaceae"                         "Stereum"           "ostrea"       
#> ASV25  "Lyophyllaceae"                      "Ossicaulis"        "lachnopus"    
#> ASV26  "Stereaceae"                         "Stereum"           "hirsutum"     
#> ASV27  "Steccherinaceae"                    "Antrodiella"       "brasiliensis" 
#> ASV29  "Exidiaceae"                         "Basidiodendron"    "eyrei"        
#> ASV32  "Cantharellales_fam_Incertae_sedis"  "Sistotrema"        "oblongisporum"
#> ASV34  "Entolomataceae"                     "Entocybe"          NA             
#> ASV35  "Polyporaceae"                       "Fomes"             "fomentarius"  
#> ASV41  "Tricholomataceae"                   "Mycena"            "renati"       
#> ASV42  "Lyophyllaceae"                      "Ossicaulis"        "lachnopus"    
#> ASV46  "Atractiellales_fam_Incertae_sedis"  "Helicogloea"       "pellucida"    
#> ASV47  "Pterulaceae"                        "Radulomyces"       "molaris"      
#> ASV48  "Aporpiaceae"                        "Elmerina"          "caryae"       
#> ASV49  "Phanerochaetaceae"                  "Phanerochaete"     "livescens"    
#> ASV50  "Russulales_fam_Incertae_sedis"      "Gloeohypochnicium" "analogum"     
#> ASV53  "Polyporaceae"                       "Fomes"             "fomentarius"  
#> ASV54  "Auriculariaceae"                    "Auricularia"       NA             
#> ASV58  "Polyporaceae"                       "Fomes"             "fomentarius"  
#> ASV59  "Hyphodermataceae"                   "Hyphoderma"        "roseocremeum" 
#> ASV61  "Hyphodermataceae"                   "Hyphoderma"        "setigerum"    
#> ASV62  NA                                   NA                  NA             
#> ASV63  NA                                   NA                  NA             
#> ASV64  "Polyporaceae"                       "Trametes"          "versicolor"   
#> ASV67  "Schizoporaceae"                     "Xylodon"           "raduloides"   
#> ASV68  "Lyophyllaceae"                      "Ossicaulis"        "lachnopus"    
#> ASV71  NA                                   NA                  NA             
#> ASV72  NA                                   NA                  NA             
#> ASV75  "Peniophoraceae"                     "Peniophora"        "versiformis"  
#> ASV77  "Lyophyllaceae"                      "Ossicaulis"        "lachnopus"    
#> ASV82  "Exidiaceae"                         "Exidia"            "glandulosa"   
#> ASV83  NA                                   NA                  NA             
#> ASV85  "Hymenochaetales_fam_Incertae_sedis" "Peniophorella"     "pubera"       
#> ASV91  "Auriculariaceae"                    "Auricularia"       "mesenterica"  
#> ASV93  "Stereaceae"                         NA                  NA             
#> ASV94  "Stereaceae"                         "Stereum"           "ostrea"       
#> ASV99  "Polyporaceae"                       "Fomes"             "fomentarius"  
#> ASV100 NA                                   NA                  NA             
#> ASV101 "Corticiaceae"                       "Marchandiomyces"   "buckii"       
#> ASV104 "Hericiaceae"                        "Hericium"          "coralloides"  
#> ASV105 "Schizoporaceae"                     "Xylodon"           "flaviporus"   
#> ASV107 "Schizoporaceae"                     "Xylodon"           "raduloides"   
#> ASV108 "Exidiaceae"                         "Exidia"            "glandulosa"   
#>        Trophic.Mode                       
#> ASV7   "Saprotroph"                       
#> ASV8   "Saprotroph"                       
#> ASV12  "Saprotroph"                       
#> ASV18  "Saprotroph"                       
#> ASV25  "Saprotroph"                       
#> ASV26  "Saprotroph"                       
#> ASV27  "Saprotroph"                       
#> ASV29  "Saprotroph"                       
#> ASV32  "Saprotroph-Symbiotroph"           
#> ASV34  "Saprotroph"                       
#> ASV35  "Saprotroph"                       
#> ASV41  "Pathotroph-Saprotroph"            
#> ASV42  "Saprotroph"                       
#> ASV46  "Saprotroph"                       
#> ASV47  "Saprotroph"                       
#> ASV48  "Saprotroph"                       
#> ASV49  "Saprotroph"                       
#> ASV50  "Saprotroph"                       
#> ASV53  "Saprotroph"                       
#> ASV54  "Saprotroph"                       
#> ASV58  "Saprotroph"                       
#> ASV59  "Saprotroph"                       
#> ASV61  "Saprotroph"                       
#> ASV62  "-"                                
#> ASV63  "-"                                
#> ASV64  "Saprotroph"                       
#> ASV67  "Saprotroph"                       
#> ASV68  "Saprotroph"                       
#> ASV71  "-"                                
#> ASV72  "-"                                
#> ASV75  "Pathotroph-Saprotroph"            
#> ASV77  "Saprotroph"                       
#> ASV82  "Saprotroph-Symbiotroph"           
#> ASV83  "Pathotroph-Saprotroph-Symbiotroph"
#> ASV85  "Saprotroph"                       
#> ASV91  "Saprotroph"                       
#> ASV93  "Saprotroph"                       
#> ASV94  "Saprotroph"                       
#> ASV99  "Saprotroph"                       
#> ASV100 "-"                                
#> ASV101 "Pathotroph"                       
#> ASV104 "Saprotroph"                       
#> ASV105 "Saprotroph"                       
#> ASV107 "Saprotroph"                       
#> ASV108 "Saprotroph-Symbiotroph"           
#>        Guild                                                            
#> ASV7   "WoodSaprotroph-UndefinedSaprotroph"                             
#> ASV8   "UndefinedSaprotroph"                                            
#> ASV12  "UndefinedSaprotroph"                                            
#> ASV18  "UndefinedSaprotroph"                                            
#> ASV25  "WoodSaprotroph"                                                 
#> ASV26  "UndefinedSaprotroph"                                            
#> ASV27  "WoodSaprotroph"                                                 
#> ASV29  "UndefinedSaprotroph"                                            
#> ASV32  "Ectomycorrhizal-WoodSaprotroph"                                 
#> ASV34  "UndefinedSaprotroph"                                            
#> ASV35  "WoodSaprotroph"                                                 
#> ASV41  "LeafSaprotroph-PlantPathogen-UndefinedSaprotroph-WoodSaprotroph"
#> ASV42  "WoodSaprotroph"                                                 
#> ASV46  "UndefinedSaprotroph"                                            
#> ASV47  "UndefinedSaprotroph"                                            
#> ASV48  "UndefinedSaprotroph"                                            
#> ASV49  "WoodSaprotroph"                                                 
#> ASV50  "UndefinedSaprotroph"                                            
#> ASV53  "WoodSaprotroph"                                                 
#> ASV54  "UndefinedSaprotroph"                                            
#> ASV58  "WoodSaprotroph"                                                 
#> ASV59  "UndefinedSaprotroph"                                            
#> ASV61  "UndefinedSaprotroph"                                            
#> ASV62  "-"                                                              
#> ASV63  "-"                                                              
#> ASV64  "WoodSaprotroph"                                                 
#> ASV67  "UndefinedSaprotroph"                                            
#> ASV68  "WoodSaprotroph"                                                 
#> ASV71  "-"                                                              
#> ASV72  "-"                                                              
#> ASV75  "PlantPathogen-WoodSaprotroph"                                   
#> ASV77  "WoodSaprotroph"                                                 
#> ASV82  "Endophyte-UndefinedSaprotroph"                                  
#> ASV83  "FungalParasite-UndefinedSaprotroph"                             
#> ASV85  "UndefinedSaprotroph"                                            
#> ASV91  "UndefinedSaprotroph"                                            
#> ASV93  "WoodSaprotroph-UndefinedSaprotroph"                             
#> ASV94  "UndefinedSaprotroph"                                            
#> ASV99  "WoodSaprotroph"                                                 
#> ASV100 "-"                                                              
#> ASV101 "LichenParasite"                                                 
#> ASV104 "UndefinedSaprotroph"                                            
#> ASV105 "UndefinedSaprotroph"                                            
#> ASV107 "UndefinedSaprotroph"                                            
#> ASV108 "Endophyte-UndefinedSaprotroph"                                  
#>        Trait               Confidence.Ranking Genus_species               
#> ASV7   "NULL"              "Probable"         "_"                         
#> ASV8   "WhiteRot"          "Probable"         "Stereum_ostrea"            
#> ASV12  "WhiteRot"          "Probable"         "Xylodon_raduloides"        
#> ASV18  "WhiteRot"          "Probable"         "Stereum_ostrea"            
#> ASV25  "BrownRot"          "Probable"         "Ossicaulis_lachnopus"      
#> ASV26  "WhiteRot"          "Probable"         "Stereum_hirsutum"          
#> ASV27  "WhiteRot"          "HighlyProbable"   "Antrodiella_brasiliensis"  
#> ASV29  "NULL"              "Probable"         "Basidiodendron_eyrei"      
#> ASV32  "WhiteRot"          "Possible"         "Sistotrema_oblongisporum"  
#> ASV34  "NULL"              "Probable"         "Entocybe_"                 
#> ASV35  "BrownRot-WhiteRot" "Probable"         "Fomes_fomentarius"         
#> ASV41  "NULL"              "Probable"         "Mycena_renati"             
#> ASV42  "BrownRot"          "Probable"         "Ossicaulis_lachnopus"      
#> ASV46  "NULL"              "Probable"         "Helicogloea_pellucida"     
#> ASV47  "WhiteRot"          "Probable"         "Radulomyces_molaris"       
#> ASV48  "NULL"              "Probable"         "Elmerina_caryae"           
#> ASV49  "WhiteRot"          "HighlyProbable"   "Phanerochaete_livescens"   
#> ASV50  "WhiteRot"          "Probable"         "Gloeohypochnicium_analogum"
#> ASV53  "BrownRot-WhiteRot" "Probable"         "Fomes_fomentarius"         
#> ASV54  "NULL"              "Probable"         "Auricularia_"              
#> ASV58  "BrownRot-WhiteRot" "Probable"         "Fomes_fomentarius"         
#> ASV59  "WhiteRot"          "Probable"         "Hyphoderma_roseocremeum"   
#> ASV61  "WhiteRot"          "Probable"         "Hyphoderma_setigerum"      
#> ASV62  "-"                 "-"                "_"                         
#> ASV63  "-"                 "-"                "_"                         
#> ASV64  "WhiteRot"          "HighlyProbable"   "Trametes_versicolor"       
#> ASV67  "WhiteRot"          "Probable"         "Xylodon_raduloides"        
#> ASV68  "BrownRot"          "Probable"         "Ossicaulis_lachnopus"      
#> ASV71  "-"                 "-"                "_"                         
#> ASV72  "-"                 "-"                "_"                         
#> ASV75  "WhiteRot"          "Probable"         "Peniophora_versiformis"    
#> ASV77  "BrownRot"          "Probable"         "Ossicaulis_lachnopus"      
#> ASV82  "NULL"              "Probable"         "Exidia_glandulosa"         
#> ASV83  "NULL"              "Possible"         "_"                         
#> ASV85  "WhiteRot"          "Probable"         "Peniophorella_pubera"      
#> ASV91  "NULL"              "Probable"         "Auricularia_mesenterica"   
#> ASV93  "NULL"              "Probable"         "_"                         
#> ASV94  "WhiteRot"          "Probable"         "Stereum_ostrea"            
#> ASV99  "BrownRot-WhiteRot" "Probable"         "Fomes_fomentarius"         
#> ASV100 "-"                 "-"                "_"                         
#> ASV101 "NULL"              "Probable"         "Marchandiomyces_buckii"    
#> ASV104 "WhiteRot"          "Probable"         "Hericium_coralloides"      
#> ASV105 "WhiteRot"          "Probable"         "Xylodon_flaviporus"        
#> ASV107 "WhiteRot"          "Probable"         "Xylodon_raduloides"        
#> ASV108 "NULL"              "Probable"         "Exidia_glandulosa"         
#>        Kingdom_blastn Phylum_blastn              Class_blastn              
#> ASV7   NA             NA                         NA                        
#> ASV8   NA             NA                         NA                        
#> ASV12  NA             NA                         NA                        
#> ASV18  NA             NA                         NA                        
#> ASV25  NA             NA                         NA                        
#> ASV26  NA             NA                         NA                        
#> ASV27  NA             NA                         NA                        
#> ASV29  NA             NA                         NA                        
#> ASV32  NA             NA                         NA                        
#> ASV34  NA             NA                         NA                        
#> ASV35  NA             NA                         NA                        
#> ASV41  "Fungi"        "Basidiomycota"            "Agaricomycetes"          
#> ASV42  NA             NA                         NA                        
#> ASV46  NA             NA                         NA                        
#> ASV47  NA             NA                         NA                        
#> ASV48  NA             NA                         NA                        
#> ASV49  NA             NA                         NA                        
#> ASV50  NA             NA                         NA                        
#> ASV53  NA             NA                         NA                        
#> ASV54  NA             NA                         NA                        
#> ASV58  NA             NA                         NA                        
#> ASV59  NA             NA                         NA                        
#> ASV61  NA             NA                         NA                        
#> ASV62  NA             NA                         NA                        
#> ASV63  NA             NA                         NA                        
#> ASV64  NA             NA                         NA                        
#> ASV67  NA             NA                         NA                        
#> ASV68  NA             NA                         NA                        
#> ASV71  "Fungi"        "Basidiomycota"            "Tremellomycetes"         
#> ASV72  NA             NA                         NA                        
#> ASV75  "Fungi"        "Fungi_phy_Incertae_sedis" "Fungi_cls_Incertae_sedis"
#> ASV77  NA             NA                         NA                        
#> ASV82  NA             NA                         NA                        
#> ASV83  NA             NA                         NA                        
#> ASV85  NA             NA                         NA                        
#> ASV91  NA             NA                         NA                        
#> ASV93  NA             NA                         NA                        
#> ASV94  NA             NA                         NA                        
#> ASV99  NA             NA                         NA                        
#> ASV100 NA             NA                         NA                        
#> ASV101 "Fungi"        "Basidiomycota"            "Agaricomycetes"          
#> ASV104 NA             NA                         NA                        
#> ASV105 NA             NA                         NA                        
#> ASV107 NA             NA                         NA                        
#> ASV108 NA             NA                         NA                        
#>        Order_blastn                        
#> ASV7   NA                                  
#> ASV8   NA                                  
#> ASV12  NA                                  
#> ASV18  NA                                  
#> ASV25  NA                                  
#> ASV26  NA                                  
#> ASV27  NA                                  
#> ASV29  NA                                  
#> ASV32  NA                                  
#> ASV34  NA                                  
#> ASV35  NA                                  
#> ASV41  "Agaricales"                        
#> ASV42  NA                                  
#> ASV46  NA                                  
#> ASV47  NA                                  
#> ASV48  NA                                  
#> ASV49  NA                                  
#> ASV50  NA                                  
#> ASV53  NA                                  
#> ASV54  NA                                  
#> ASV58  NA                                  
#> ASV59  NA                                  
#> ASV61  NA                                  
#> ASV62  NA                                  
#> ASV63  NA                                  
#> ASV64  NA                                  
#> ASV67  NA                                  
#> ASV68  NA                                  
#> ASV71  "Tremellomycetes_ord_Incertae_sedis"
#> ASV72  NA                                  
#> ASV75  "Trechisporales"                    
#> ASV77  NA                                  
#> ASV82  NA                                  
#> ASV83  NA                                  
#> ASV85  NA                                  
#> ASV91  NA                                  
#> ASV93  NA                                  
#> ASV94  NA                                  
#> ASV99  NA                                  
#> ASV100 NA                                  
#> ASV101 "Corticiales"                       
#> ASV104 NA                                  
#> ASV105 NA                                  
#> ASV107 NA                                  
#> ASV108 NA                                  
#>        Family_blastn                       
#> ASV7   NA                                  
#> ASV8   NA                                  
#> ASV12  NA                                  
#> ASV18  NA                                  
#> ASV25  NA                                  
#> ASV26  NA                                  
#> ASV27  NA                                  
#> ASV29  NA                                  
#> ASV32  NA                                  
#> ASV34  NA                                  
#> ASV35  NA                                  
#> ASV41  "Tricholomataceae"                  
#> ASV42  NA                                  
#> ASV46  NA                                  
#> ASV47  NA                                  
#> ASV48  NA                                  
#> ASV49  NA                                  
#> ASV50  NA                                  
#> ASV53  NA                                  
#> ASV54  NA                                  
#> ASV58  NA                                  
#> ASV59  NA                                  
#> ASV61  NA                                  
#> ASV62  NA                                  
#> ASV63  NA                                  
#> ASV64  NA                                  
#> ASV67  NA                                  
#> ASV68  NA                                  
#> ASV71  "Tremellomycetes_fam_Incertae_sedis"
#> ASV72  NA                                  
#> ASV75  "Hydnodontaceae"                    
#> ASV77  NA                                  
#> ASV82  NA                                  
#> ASV83  NA                                  
#> ASV85  NA                                  
#> ASV91  NA                                  
#> ASV93  NA                                  
#> ASV94  NA                                  
#> ASV99  NA                                  
#> ASV100 NA                                  
#> ASV101 "Corticiales_fam_Incertae_sedis"    
#> ASV104 NA                                  
#> ASV105 NA                                  
#> ASV107 NA                                  
#> ASV108 NA                                  
#>        Genus_blastn                          Species_blastn               
#> ASV7   NA                                    NA                           
#> ASV8   NA                                    NA                           
#> ASV12  NA                                    NA                           
#> ASV18  NA                                    NA                           
#> ASV25  NA                                    NA                           
#> ASV26  NA                                    NA                           
#> ASV27  NA                                    NA                           
#> ASV29  NA                                    NA                           
#> ASV32  NA                                    NA                           
#> ASV34  NA                                    NA                           
#> ASV35  NA                                    NA                           
#> ASV41  "Tricholomataceae_gen_Incertae_sedis" "Mycena_quercus-ilicis"      
#> ASV42  NA                                    NA                           
#> ASV46  NA                                    NA                           
#> ASV47  NA                                    NA                           
#> ASV48  NA                                    NA                           
#> ASV49  NA                                    NA                           
#> ASV50  NA                                    NA                           
#> ASV53  NA                                    NA                           
#> ASV54  NA                                    NA                           
#> ASV58  NA                                    NA                           
#> ASV59  NA                                    NA                           
#> ASV61  NA                                    NA                           
#> ASV62  NA                                    NA                           
#> ASV63  NA                                    NA                           
#> ASV64  NA                                    NA                           
#> ASV67  NA                                    NA                           
#> ASV68  NA                                    NA                           
#> ASV71  "Tremellomycetes_gen_Incertae_sedis"  "Tremellomycetes_sp"         
#> ASV72  NA                                    NA                           
#> ASV75  "Subulicystidium"                     "Subulicystidium_longisporum"
#> ASV77  NA                                    NA                           
#> ASV82  NA                                    NA                           
#> ASV83  NA                                    NA                           
#> ASV85  NA                                    NA                           
#> ASV91  NA                                    NA                           
#> ASV93  NA                                    NA                           
#> ASV94  NA                                    NA                           
#> ASV99  NA                                    NA                           
#> ASV100 NA                                    NA                           
#> ASV101 "Corticiales_gen_Incertae_sedis"      "Corticiales_sp"             
#> ASV104 NA                                    NA                           
#> ASV105 NA                                    NA                           
#> ASV107 NA                                    NA                           
#> ASV108 NA                                    NA                           
#>        Taxa_name_db_blastn                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                           
#> ASV7   NA                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                            
#> ASV8   NA                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                            
#> ASV12  NA                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                            
#> ASV18  NA                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                            
#> ASV25  NA                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                            
#> ASV26  NA                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                            
#> ASV27  NA                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                            
#> ASV29  NA                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                            
#> ASV32  NA                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                            
#> ASV34  NA                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                            
#> ASV35  NA                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                            
#> ASV41  "Mycenaceae_sp_OP470041_SH0013594.10FU_reps_singleton_/Mycena_sp_UDB05068035_SH0013561.10FU_reps_singleton_/Mycena_sp_UDB05033057_SH0013579.10FU_reps_singleton_/Mycena_venustula_OK239713_SH0008654.10FU_reps_singleton_/Mycena_sp_UDB05006031_SH0013576.10FU_reps_singleton_/Mycenaceae_sp_UDB05167898_SH0013566.10FU_reps_singleton_/Mycena_sp_OR987288_SH0013580.10FU_reps_singleton_/Tricholomataceae_sp_LC574741_SH0013558.10FU_reps_singleton_/Mycena_acanthophila_OQ633198_SH0008640.10FU_reps_singleton_/Mycena_quercus-ilicis_OQ633195_SH0008657.10FU_reps_singleton_/Mycena_sp_OR883535_SH0013563.10FU_reps_singleton_/Mycena_sp_UDB023583_SH0008647.10FU_reps_singleton_/Mycena_sp_UDB05068051_SH0013586.10FU_reps_singleton_/Mycena_sp_UDB05018702_SH0013573.10FU_reps_singleton_/Mycena_sp_UDB05154274_SH0008652.10FU_reps_singleton_/Mycena_sp_UDB05027734_SH0013564.10FU_reps_singleton_/Mycena_sp_UDB05115971_SH0013554.10FU_reps_singleton_"
#> ASV42  NA                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                            
#> ASV46  NA                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                            
#> ASV47  NA                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                            
#> ASV48  NA                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                            
#> ASV49  NA                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                            
#> ASV50  NA                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                            
#> ASV53  NA                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                            
#> ASV54  NA                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                            
#> ASV58  NA                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                            
#> ASV59  NA                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                            
#> ASV61  NA                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                            
#> ASV62  NA                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                            
#> ASV63  NA                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                            
#> ASV64  NA                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                            
#> ASV67  NA                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                            
#> ASV68  NA                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                            
#> ASV71  "Tremellales_sp_UDB06404467_SH0007343.10FU_reps_singleton_/Tremellomycetes_sp_UDB06410903_SH0007341.10FU_reps_singleton_"                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                     
#> ASV72  NA                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                            
#> ASV75  "Peniophora_reidii_MH861616_SH0003674.10FU_reps_singleton_/Peniophora_sp_OR149912_SH0003663.10FU_reps_singleton_/Peniophora_sp_MK588791_SH0003657.10FU_reps_singleton_/Peniophora_pini_KT210102_SH0003677.10FU_reps_singleton_/Peniophora_sp_KJ832048_SH0003664.10FU_reps_singleton_/Fungi_sp_KU515816_SH0003692.10FU_reps_singleton_/Subulicystidium_longisporum_OK310732_SH0003688.10FU_reps_singleton_/Fungi_sp_GU721554_SH0003687.10FU_reps_singleton_/Fungi_sp_GQ999289_SH0003662.10FU_reps_singleton_/Peniophora_sp_AB672622_SH0003656.10FU_reps_singleton_/Fungi_sp_JX135952_SH0003681.10FU_reps_singleton_"                                                                                                                                                                                                                                                                                                                                           
#> ASV77  NA                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                            
#> ASV82  NA                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                            
#> ASV83  NA                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                            
#> ASV85  NA                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                            
#> ASV91  NA                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                            
#> ASV93  NA                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                            
#> ASV94  NA                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                            
#> ASV99  NA                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                            
#> ASV100 NA                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                            
#> ASV101 "Corticiales_sp_KY462668_SH0001937.10FU_reps_singleton_"                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                      
#> ASV104 NA                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                            
#> ASV105 NA                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                            
#> ASV107 NA                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                            
#> ASV108 NA                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                            

## return_taxtab with seq2search (a DNAStringSet or a matrix)
seqs_to_assign <- refseq(data_fungi_mini)
tax_mat <- assign_blastn(seq2search = seqs_to_assign,
  ref_fasta = ref_fasta, method_algo = "top-hit", min_id = 70,
  min_e_value = 1e-3, min_cover = 50, min_bit_score = 20,
  behavior = "return_taxtab"
)
#> Build the database from fasta_for_db
#> Blast refseq from physeq object against the database
head(tax_mat)
#>       Kingdom   Phylum            Class              Order         
#> ASV41 "k:Fungi" "p:Basidiomycota" "c:Agaricomycetes" "o:Agaricales"
#> ASV75 "k:Fungi" "p:Basidiomycota" "c:Agaricomycetes" "o:Russulales"
#> ASV26 "k:Fungi" "p:Basidiomycota" "c:Agaricomycetes" "o:Russulales"
#> ASV25 "k:Fungi" "p:Basidiomycota" "c:Agaricomycetes" "o:Agaricales"
#> ASV68 "k:Fungi" "p:Basidiomycota" "c:Agaricomycetes" "o:Agaricales"
#> ASV42 "k:Fungi" "p:Basidiomycota" "c:Agaricomycetes" "o:Agaricales"
#>       Family               Genus                            
#> ASV41 "f:Mycenaceae"       "g:Mycenaceae_gen_Incertae_sedis"
#> ASV75 "f:Peniophoraceae"   "g:Peniophora"                   
#> ASV26 "f:Stereaceae"       "g:Xylobolus"                    
#> ASV25 "f:Tricholomataceae" "g:Tricholoma"                   
#> ASV68 "f:Tricholomataceae" "g:Tricholoma"                   
#> ASV42 "f:Tricholomataceae" "g:Tricholoma"                   
#>       Species              
#> ASV41 "s:Mycenaceae_sp"    
#> ASV75 "s:Peniophora_reidii"
#> ASV26 "s:Xylobolus_sp"     
#> ASV25 "s:Tricholoma_sp"    
#> ASV68 "s:Tricholoma_sp"    
#> ASV42 "s:Tricholoma_sp"    
```
