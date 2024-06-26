% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/targets_misc.R
\name{sample_data_with_new_names}
\alias{sample_data_with_new_names}
\title{Load sample data from file and rename samples using names of samples and an
optional order}
\usage{
sample_data_with_new_names(
  file_path,
  names_of_samples,
  samples_order = NULL,
  ...
)
}
\arguments{
\item{file_path}{(required) a path to the sample_data file}

\item{names_of_samples}{(required) a vector of sample names}

\item{samples_order}{Optional numeric vector to sort sample names}

\item{...}{Other arguments passed on to \code{\link[utils:read.table]{utils::read.delim()}} function.}
}
\value{
A data.frame from file_path and new names
}
\description{
\if{html}{\out{
<a href="https://adrientaudiere.github.io/MiscMetabar/articles/Rules.html#lifecycle">
<img src="https://img.shields.io/badge/lifecycle-maturing-blue" alt="lifecycle-maturing"></a>
}}


Useful for targets bioinformatic pipeline.
}
\examples{
sam_file <- system.file("extdata", "sam_data.csv", package = "MiscMetabar")
sample_data_with_new_names(sam_file, paste0("Samples_", seq(1, 185)))

}
\seealso{
\code{\link[=rename_samples]{rename_samples()}}
}
\author{
Adrien Taudière
}
