# Translates a factor into colors.

Translates a factor into colors.

## Usage

``` r
fac2col(x, col.pal = funky_color, na.col = "grey", seed = NULL)
```

## Arguments

- x:

  a numeric vector (for num2col) or a vector converted to a factor (for
  fac2col).

- col.pal:

  (default funky_color) a function generating colors according to a
  given palette.

- na.col:

  (default grey) the color to be used for missing values (NAs)

- seed:

  (default NULL) a seed for R's random number generated, used to fix the
  random permutation of colors in the palette used; if NULL, no
  randomization is used and the colors are taken from the palette
  according to the ordering of the levels

## Value

a color vector

## See also

The R package RColorBrewer, proposing a nice selection of color
palettes. The viridis package, with many excellent palettes

## Author

Thibaut Jombart in `adegenet` package
