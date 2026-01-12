# Multiple plot function

[![lifecycle-stable](https://img.shields.io/badge/lifecycle-stable-green)](https://adrientaudiere.github.io/MiscMetabar/articles/Rules.html#lifecycle)

ggplot objects can be passed in ..., or to plotlist (as a list of ggplot
objects)

If the layout is something like matrix(c(1,2,3,3), nrow=2, byrow=TRUE),
then plot 1 will go in the upper left, 2 will go in the upper right, and
3 will go all the way across the bottom.

## Usage

``` r
multiplot(..., plotlist = NULL, cols = 1, layout = NULL)
```

## Arguments

- ...:

  list of ggplot objects

- plotlist:

  list of ggplot objects

- cols:

  number of columns

- layout:

  A matrix specifying the layout. If present, 'cols' is ignored.

## Value

Nothing. Print the list of ggplot objects
