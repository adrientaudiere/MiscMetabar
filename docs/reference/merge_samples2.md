# Merge samples by a sample variable or factor

[![lifecycle-stable](https://img.shields.io/badge/lifecycle-stable-green)](https://adrientaudiere.github.io/MiscMetabar/articles/Rules.html#lifecycle)

Firstly release in the [speedyseq](https://github.com/mikemc/speedyseq/)
R package by Michael R. McLaren.

This function provides an alternative to
[`phyloseq::merge_samples()`](https://rdrr.io/pkg/phyloseq/man/merge_samples-methods.html)
that better handles sample variables of different types, especially
categorical sample variables. It combines the samples in `x` defined by
the sample variable or factor `group` by summing the abundances in
`otu_table(x)` and combines sample variables by the summary functions in
`funs`. The default summary function,
[`unique_or_na()`](https://adrientaudiere.github.io/MiscMetabar/reference/unique_or_na.md),
collapses the values within a group to a single unique value if it
exists and otherwise returns NA. The new (merged) samples are named by
the values in `group`.

## Usage

``` r
merge_samples2(
  x,
  group,
  fun_otu = sum,
  funs = list(),
  reorder = FALSE,
  default_fun = unique_or_na
)

# S4 method for class 'phyloseq'
merge_samples2(
  x,
  group,
  fun_otu = sum,
  funs = list(),
  reorder = FALSE,
  default_fun = unique_or_na
)

# S4 method for class 'otu_table'
merge_samples2(
  x,
  group,
  fun_otu = sum,
  reorder = FALSE,
  default_fun = unique_or_na
)

# S4 method for class 'sample_data'
merge_samples2(
  x,
  group,
  funs = list(),
  reorder = FALSE,
  default_fun = unique_or_na
)
```

## Arguments

- x:

  A `phyloseq`, `otu_table`, or `sample_data` object

- group:

  A sample variable or a vector of length `nsamples(x)` defining the
  sample grouping. A vector must be supplied if x is an otu_table

- fun_otu:

  Function for combining abundances in the otu_table; default is `sum`.
  Can be a formula to be converted to a function by
  [`purrr::as_mapper()`](https://purrr.tidyverse.org/reference/as_mapper.html)

- funs:

  Named list of merge functions for sample variables; default is
  `unique_or_na`

- reorder:

  Logical specifying whether to reorder the new (merged) samples by name

- default_fun:

  Default functions if funs is not set. Per default the function
  unique_or_na is used. See
  [`diff_fct_diff_class()`](https://adrientaudiere.github.io/MiscMetabar/reference/diff_fct_diff_class.md)
  for a useful alternative.

## Value

A new phyloseq-class, otu_table or sam_data object depending on the
class of the x param

## Author

Michael R. McLaren (orcid:
[0000-0003-1575-473X](https://orcid.org/0000-0003-1575-473X)) modified
by Adrien Taudiere

## Examples

``` r
data(enterotype)

# Merge samples with the same project and clinical status
ps <- enterotype
sample_data(ps) <- sample_data(ps) %>%
  transform(Project.ClinicalStatus = Project:ClinicalStatus)
sample_data(ps) %>% head()
#>           Enterotype Sample_ID SeqTech  SampleID     Project Nationality Gender
#> AM.AD.1         <NA>   AM.AD.1  Sanger   AM.AD.1      gill06    american      F
#> AM.AD.2         <NA>   AM.AD.2  Sanger   AM.AD.2      gill06    american      M
#> AM.F10.T1       <NA> AM.F10.T1  Sanger AM.F10.T1 turnbaugh09    american      F
#> AM.F10.T2          3 AM.F10.T2  Sanger AM.F10.T2 turnbaugh09    american      F
#> DA.AD.1            2   DA.AD.1  Sanger   DA.AD.1     MetaHIT      danish      F
#> DA.AD.1T        <NA>  DA.AD.1T  Sanger      <NA>        <NA>        <NA>   <NA>
#>           Age ClinicalStatus Project.ClinicalStatus
#> AM.AD.1    28        healthy         gill06:healthy
#> AM.AD.2    37        healthy         gill06:healthy
#> AM.F10.T1  NA          obese      turnbaugh09:obese
#> AM.F10.T2  NA          obese      turnbaugh09:obese
#> DA.AD.1    59        healthy        MetaHIT:healthy
#> DA.AD.1T   NA           <NA>                   <NA>
ps0 <- merge_samples2(ps, "Project.ClinicalStatus",
  fun_otu = mean,
  funs = list(Age = mean)
)
#> Warning: `group` has missing values; corresponding samples will be dropped
sample_data(ps0) %>% head()
#>                     Enterotype Sample_ID SeqTech SampleID     Project
#> gill06:healthy            <NA>      <NA>  Sanger     <NA>      gill06
#> kurokawa07:healthy        <NA>      <NA>  Sanger     <NA>  kurokawa07
#> kurokawa07 :healthy       <NA>      <NA>  Sanger     <NA> kurokawa07 
#> MetaHIT:CD                   1   ES.AD.1  Sanger  ES.AD.1     MetaHIT
#> MetaHIT:healthy           <NA>      <NA>  Sanger     <NA>     MetaHIT
#> MetaHIT:obese             <NA>      <NA>  Sanger     <NA>     MetaHIT
#>                     Nationality Gender      Age ClinicalStatus
#> gill06:healthy         american   <NA> 32.50000        healthy
#> kurokawa07:healthy     japanese      M 15.12500        healthy
#> kurokawa07 :healthy    japanese   <NA> 19.17364        healthy
#> MetaHIT:CD              spanish      F 25.00000             CD
#> MetaHIT:healthy            <NA>   <NA> 50.00000        healthy
#> MetaHIT:obese            danish   <NA> 54.00000          obese
#>                     Project.ClinicalStatus
#> gill06:healthy              gill06:healthy
#> kurokawa07:healthy      kurokawa07:healthy
#> kurokawa07 :healthy    kurokawa07 :healthy
#> MetaHIT:CD                      MetaHIT:CD
#> MetaHIT:healthy            MetaHIT:healthy
#> MetaHIT:obese                MetaHIT:obese
```
