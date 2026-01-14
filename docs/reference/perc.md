# Convert a value (or a fraction x/y) in percentage

[![lifecycle-maturing](https://img.shields.io/badge/lifecycle-maturing-blue)](https://adrientaudiere.github.io/MiscMetabar/articles/Rules.html#lifecycle)

Mostly for internal use.

## Usage

``` r
perc(x, y = NULL, accuracy = 0, add_symbol = FALSE)
```

## Arguments

- x:

  (required) value

- y:

  if y is set, compute the division of x by y

- accuracy:

  number of digits (number of digits after zero)

- add_symbol:

  if set to TRUE add the % symbol to the value

## Value

The percentage value (number or character if add_symbol is set to TRUE)

## Author

Adrien Taudière
