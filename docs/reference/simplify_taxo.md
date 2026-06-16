# Simplify taxonomy by removing some unused characters such as "k\_\_"

[![lifecycle-maturing](https://img.shields.io/badge/lifecycle-maturing-blue)](https://adrientaudiere.github.io/MiscMetabar/articles/Rules.html#lifecycle)

Internally used in
[`clean_pq()`](https://adrientaudiere.github.io/MiscMetabar/reference/clean_pq.md)

## Usage

``` r
simplify_taxo(
  physeq,
  pattern_to_remove = c(".__", ".*:"),
  ranks_for_pattern_to_remove = phyloseq::rank_names(physeq),
  ranks_to_remove_space = phyloseq::rank_names(physeq)[!grepl("Species",
    phyloseq::rank_names(physeq))],
  ranks_to_remove_NA = phyloseq::rank_names(physeq),
  pattern_to_NA = NULL,
  ranks_for_pattern_to_NA = phyloseq::rank_names(physeq)
)
```

## Arguments

- physeq:

  (required) a
  [`phyloseq-class`](https://rdrr.io/pkg/phyloseq/man/phyloseq-class.html)
  object obtained using the `phyloseq` package.

- pattern_to_remove:

  (a vector of character) regex patterns passed to
  [`base::gsub()`](https://rdrr.io/r/base/grep.html): the matched
  *substring* is deleted from the cell value and the rest of the string
  is kept (e.g. `".__"` turns `"k__Fungi"` into `"Fungi"`).

- ranks_for_pattern_to_remove:

  (character vector or NULL; default all ranks) column names in
  `tax_table` to which `pattern_to_remove` is applied. Pass `NULL` to
  skip this operation on all columns.

- ranks_to_remove_space:

  (character vector or NULL; default all ranks whose name does not
  contain `"Species"`) column names from which ASCII spaces and
  non-breaking spaces (U+00A0) are stripped. Pass `NULL` to skip space
  removal entirely.

- ranks_to_remove_NA:

  (character vector or NULL; default all ranks) column names from which
  the literal string `"NA"` (case-sensitive) is removed. Pass `NULL` to
  skip this operation. **Breaking change from v0.16:** the old
  `remove_NA = FALSE` default is now `ranks_to_remove_NA` defaulting to
  all ranks; pass `NULL` to reproduce the old behaviour.

- pattern_to_NA:

  (character; default NULL): a regex; if an entire cell value matches,
  the *whole cell* is replaced with `NA` (nothing from the original
  value is kept). Designed for PR2-style placeholder unknowns such as
  `Embryophyceae_X`, `Embryophyceae_XX`, `Embryophyceae_XXX`,
  `Embryophyceae_XXX_sp.`, or `Mortierella_sp.`. Use `"_X+$|_sp\\.$"` to
  cover all such patterns: `_X+$` catches rank-filler X's; `_sp\\.$`
  catches any genus-only species placeholder.

- ranks_for_pattern_to_NA:

  (character vector or NULL; default all ranks) column names to which
  `pattern_to_NA` is applied. Pass `NULL` to skip this operation on all
  columns.

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
# Apply pattern_to_remove only to Genus and Species columns
cbind(
  simplify_taxo(d_fm,
    ranks_for_pattern_to_remove = c("Genus", "Species")
  )@tax_table[, "Species"],
  d_fm@tax_table[, "Species"]
)
#> Error in simplify_taxo(d_fm, ranks_for_pattern_to_remove = c("Genus",     "Species")): unused argument (ranks_for_pattern_to_remove = c("Genus", "Species"))
if (FALSE) { # \dontrun{
# Replace PR2 placeholder unknowns (_X, _XX, _XXX, _XXX_sp., Genus_sp.) with NA
simplify_taxo(pq_pr2, pattern_to_NA = "_X+$|_sp\\.$")
# Apply pattern_to_NA only to the Species column
simplify_taxo(pq_pr2,
  pattern_to_NA = "_X+$|_sp\\.$",
  ranks_for_pattern_to_NA = "Species"
)
} # }
```
