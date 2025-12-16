# Install a package if not present

**\[experimental\]**

## Usage

``` r
install_pkg_needed(
  pkg,
  use_pak = TRUE,
  bioconductor_pkg = FALSE,
  github_pkg = FALSE,
  verbose = FALSE
)
```

## Arguments

- pkg:

  The name of the package

- use_pak:

  (logical, default TRUE) Use of
  [`pak::pkg_install()`](https://pak.r-lib.org/reference/pkg_install.html).
  If FALSE use the base `install.package()` function or the function
  [`BiocManager::install()`](https://bioconductor.github.io/BiocManager/reference/install.html)
  if bioconductor_pkg is true or the function

- bioconductor_pkg:

  (logical, default FALSE). If use_pak is TRUE, do nothing, else use
  [`BiocManager::install()`](https://bioconductor.github.io/BiocManager/reference/install.html)
  to install the package.

- github_pkg:

  (logical, default FALSE). If use_pak is TRUE, do nothing, else use
  [`devtools::install_github`](https://remotes.r-lib.org/reference/install_github.html)
  to install the package.

- verbose:

  (logical, default FALSE) Does the function print message?

## Value

Nothing

## Author

Adrien Taudière
