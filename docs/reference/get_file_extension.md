# Get the extension of a file

[![lifecycle-maturing](https://img.shields.io/badge/lifecycle-maturing-blue)](https://adrientaudiere.github.io/MiscMetabar/articles/Rules.html#lifecycle)

Internally used in
[`count_seq()`](https://adrientaudiere.github.io/MiscMetabar/reference/count_seq.md).
Warning: don't work when there is '.' in the name of the file before the
extension

## Usage

``` r
get_file_extension(file_path)
```

## Arguments

- file_path:

  (required) path to a file

## Value

The extension of a file.

## Author

Adrien Taudière
