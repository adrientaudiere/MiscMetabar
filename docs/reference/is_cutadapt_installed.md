# Test if cutadapt is installed.

[![lifecycle-maturing](https://img.shields.io/badge/lifecycle-maturing-blue)](https://adrientaudiere.github.io/MiscMetabar/articles/Rules.html#lifecycle)

Useful for testthat and examples compilation for R CMD CHECK and test
coverage

## Usage

``` r
is_cutadapt_installed(
  args_before_cutadapt =
    "source ~/miniconda3/etc/profile.d/conda.sh && conda activate cutadaptenv && "
)
```

## Arguments

- args_before_cutadapt:

  : (String) A one line bash command to run before to run cutadapt. For
  examples, "source ~/miniconda3/etc/profile.d/conda.sh && conda
  activate cutadaptenv &&" allow to bypass the conda init which asks to
  restart the shell

## Value

A logical that say if cutadapt is install in

## Author

Adrien Taudière

## Examples

``` r
MiscMetabar::is_cutadapt_installed()
#> Warning: running command 'bash /tmp/Rtmpr7txlM/script_cutadapt.sh 2>&1' had status 1
#> [1] TRUE
```
