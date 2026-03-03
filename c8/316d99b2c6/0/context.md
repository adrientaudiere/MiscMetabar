# Session Context

## User Prompts

### Prompt 1

entire status

### Prompt 2

[Request interrupted by user]

### Prompt 3

in tax_bar_pq with option fact = "a_factor_for_exemple_TIME" label_taxa=TRUE and add_ribbon = TRUE, if there is taxa that are present in other level(s) than the last level, we can't see the value. Please add a way to label taxa that are not present in other level than the last one.

### Prompt 4

I want the code you use to test for this new feature

### Prompt 5

[Request interrupted by user for tool use]

### Prompt 6

I want an example with taxa exclusive to the first bar get left-side labels (extra layers)

### Prompt 7

Ok still improving the tax_bar_pq. If their is taxonomic level(s) that are absent of both the right and the left side of the plot, add a warning that conseil to don't use label_taxa = TRUE

### Prompt 8

Style, document, and run R CMD check on the current R package, then propose fixes for any issues found.

## Steps

Run the following three commands **in sequence**, capturing all output:

1. Format code:
```r
air format .
```

2. Regenerate documentation:
```r
Rscript -e "devtools::document()"
```

3. Run package check:
```r
Rscript -e "rcmdcheck::rcmdcheck(args = c('--no-manual', '--as-cran'))"
```

## After each step

- Show a summary of what changed or was reported.
- For every ERROR, WARN...

### Prompt 9

For the note 2 how to avoid the creation of Rplots.pdf when running exemple

### Prompt 10

Yes please add

### Prompt 11

No, I did

### Prompt 12

Style, document, and run R CMD check on the current R package, then propose fixes for any issues found.

## Steps

Run the following three commands **in sequence**, capturing all output:

1. Format code:
```r
air format .
```

2. Regenerate documentation:
```r
Rscript -e "devtools::document()"
```

3. Run package check:
```r
Rscript -e "rcmdcheck::rcmdcheck(args = c('--no-manual', '--as-cran'))"
```

## After each step

- Show a summary of what changed or was reported.
- For every ERROR, WARN...

