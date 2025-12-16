# Simplify taxonomy by removing some unused characters such as "k\_\_"

[![lifecycle-maturing](https://img.shields.io/badge/lifecycle-maturing-blue)](https://adrientaudiere.github.io/MiscMetabar/articles/Rules.html#lifecycle)

Internally used in
[`clean_pq()`](https://adrientaudiere.github.io/MiscMetabar/reference/clean_pq.md)

## Usage

``` r
simplify_taxo(
  physeq,
  pattern_to_remove = c(".__", ".*:"),
  remove_space = TRUE,
  remove_NA = FALSE
)
```

## Arguments

- physeq:

  (required) a
  [`phyloseq-class`](https://rdrr.io/pkg/phyloseq/man/phyloseq-class.html)
  object obtained using the `phyloseq` package.

- pattern_to_remove:

  (a vector of character) the pattern to remove using
  [`base::gsub()`](https://rdrr.io/r/base/grep.html) function.

- remove_space:

  (logical; default TRUE): do we remove space?

- remove_NA:

  (logical; default FALSE): do we remove NA (in majuscule)?

## Value

A
[`phyloseq-class`](https://rdrr.io/pkg/phyloseq/man/phyloseq-class.html)
object with simplified taxonomy

## Author

Adrien Taudière

## Examples

``` r
d_fm <- data_fungi_mini
d_fm@tax_table[, "Species"] <- paste0(rep(
  c("s__", "s:"),
  ntaxa(d_fm) / 2
), d_fm@tax_table[, "Species"])

# First column is the new vector of Species,
# second column is the column before simplification
cbind(
  simplify_taxo(d_fm)@tax_table[, "Species"],
  d_fm@tax_table[, "Species"]
)
#>        Species         Species           
#> ASV7   "NA"            "s__NA"           
#> ASV8   "ostrea"        "s:ostrea"        
#> ASV12  "raduloides"    "s__raduloides"   
#> ASV18  "ostrea"        "s:ostrea"        
#> ASV25  "lachnopus"     "s__lachnopus"    
#> ASV26  "hirsutum"      "s:hirsutum"      
#> ASV27  "brasiliensis"  "s__brasiliensis" 
#> ASV29  "eyrei"         "s:eyrei"         
#> ASV32  "oblongisporum" "s__oblongisporum"
#> ASV34  "NA"            "s:NA"            
#> ASV35  "fomentarius"   "s__fomentarius"  
#> ASV41  "renati"        "s:renati"        
#> ASV42  "lachnopus"     "s__lachnopus"    
#> ASV46  "pellucida"     "s:pellucida"     
#> ASV47  "molaris"       "s__molaris"      
#> ASV48  "caryae"        "s:caryae"        
#> ASV49  "livescens"     "s__livescens"    
#> ASV50  "analogum"      "s:analogum"      
#> ASV53  "fomentarius"   "s__fomentarius"  
#> ASV54  "NA"            "s:NA"            
#> ASV58  "fomentarius"   "s__fomentarius"  
#> ASV59  "roseocremeum"  "s:roseocremeum"  
#> ASV61  "setigerum"     "s__setigerum"    
#> ASV62  "NA"            "s:NA"            
#> ASV63  "NA"            "s__NA"           
#> ASV64  "versicolor"    "s:versicolor"    
#> ASV67  "raduloides"    "s__raduloides"   
#> ASV68  "lachnopus"     "s:lachnopus"     
#> ASV71  "NA"            "s__NA"           
#> ASV72  "NA"            "s:NA"            
#> ASV75  "versiformis"   "s__versiformis"  
#> ASV77  "lachnopus"     "s:lachnopus"     
#> ASV82  "glandulosa"    "s__glandulosa"   
#> ASV83  "NA"            "s:NA"            
#> ASV85  "pubera"        "s__pubera"       
#> ASV91  "mesenterica"   "s:mesenterica"   
#> ASV93  "NA"            "s__NA"           
#> ASV94  "ostrea"        "s:ostrea"        
#> ASV99  "fomentarius"   "s__fomentarius"  
#> ASV100 "NA"            "s:NA"            
#> ASV101 "buckii"        "s__buckii"       
#> ASV104 "coralloides"   "s:coralloides"   
#> ASV105 "flaviporus"    "s__flaviporus"   
#> ASV107 "raduloides"    "s:raduloides"    
#> ASV108 "glandulosa"    "s__glandulosa"   
cbind(
  simplify_taxo(d_fm, remove_NA = TRUE)@tax_table[, "Species"],
  d_fm@tax_table[, "Species"]
)
#>        Species         Species           
#> ASV7   ""              "s__NA"           
#> ASV8   "ostrea"        "s:ostrea"        
#> ASV12  "raduloides"    "s__raduloides"   
#> ASV18  "ostrea"        "s:ostrea"        
#> ASV25  "lachnopus"     "s__lachnopus"    
#> ASV26  "hirsutum"      "s:hirsutum"      
#> ASV27  "brasiliensis"  "s__brasiliensis" 
#> ASV29  "eyrei"         "s:eyrei"         
#> ASV32  "oblongisporum" "s__oblongisporum"
#> ASV34  ""              "s:NA"            
#> ASV35  "fomentarius"   "s__fomentarius"  
#> ASV41  "renati"        "s:renati"        
#> ASV42  "lachnopus"     "s__lachnopus"    
#> ASV46  "pellucida"     "s:pellucida"     
#> ASV47  "molaris"       "s__molaris"      
#> ASV48  "caryae"        "s:caryae"        
#> ASV49  "livescens"     "s__livescens"    
#> ASV50  "analogum"      "s:analogum"      
#> ASV53  "fomentarius"   "s__fomentarius"  
#> ASV54  ""              "s:NA"            
#> ASV58  "fomentarius"   "s__fomentarius"  
#> ASV59  "roseocremeum"  "s:roseocremeum"  
#> ASV61  "setigerum"     "s__setigerum"    
#> ASV62  ""              "s:NA"            
#> ASV63  ""              "s__NA"           
#> ASV64  "versicolor"    "s:versicolor"    
#> ASV67  "raduloides"    "s__raduloides"   
#> ASV68  "lachnopus"     "s:lachnopus"     
#> ASV71  ""              "s__NA"           
#> ASV72  ""              "s:NA"            
#> ASV75  "versiformis"   "s__versiformis"  
#> ASV77  "lachnopus"     "s:lachnopus"     
#> ASV82  "glandulosa"    "s__glandulosa"   
#> ASV83  ""              "s:NA"            
#> ASV85  "pubera"        "s__pubera"       
#> ASV91  "mesenterica"   "s:mesenterica"   
#> ASV93  ""              "s__NA"           
#> ASV94  "ostrea"        "s:ostrea"        
#> ASV99  "fomentarius"   "s__fomentarius"  
#> ASV100 ""              "s:NA"            
#> ASV101 "buckii"        "s__buckii"       
#> ASV104 "coralloides"   "s:coralloides"   
#> ASV105 "flaviporus"    "s__flaviporus"   
#> ASV107 "raduloides"    "s:raduloides"    
#> ASV108 "glandulosa"    "s__glandulosa"   
```
