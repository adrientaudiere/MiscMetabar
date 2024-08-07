% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/plot_functions.R
\name{diff_fct_diff_class}
\alias{diff_fct_diff_class}
\title{Compute different functions for different class of vector.}
\usage{
diff_fct_diff_class(
  x,
  numeric_fonction = mean,
  logical_method = "TRUE_if_one",
  character_method = "unique_or_na",
  ...
)
}
\arguments{
\item{x}{: a vector}

\item{numeric_fonction}{: a function for numeric vector. For ex. \code{sum} or \code{mean}}

\item{logical_method}{: A method for logical vector. One of :
\itemize{
\item TRUE_if_one (default)
\item NA_if_not_all_TRUE
\item FALSE_if_not_all_TRUE
}}

\item{character_method}{: A method for character vector (and factor). One of :
\itemize{
\item unique_or_na (default)
\item more_frequent
\item more_frequent_without_equality
}}

\item{...}{Other arguments passed on to the numeric function (ex. na.rm=TRUE)}
}
\value{
a single value
}
\description{
\if{html}{\out{
<a href="https://adrientaudiere.github.io/MiscMetabar/articles/Rules.html#lifecycle">
<img src="https://img.shields.io/badge/lifecycle-experimental-orange" alt="lifecycle-experimental"></a>
}}


Mainly an internal function useful in "sapply(..., tapply)" methods
}
\examples{

diff_fct_diff_class(
  data_fungi@sam_data$Sample_id,
  numeric_fonction = sum,
  na.rm = TRUE
)
diff_fct_diff_class(
  data_fungi@sam_data$Time,
  numeric_fonction = mean,
  na.rm = TRUE
)
diff_fct_diff_class(
  data_fungi@sam_data$Height == "Low",
  logical_method = "TRUE_if_one"
)
diff_fct_diff_class(
  data_fungi@sam_data$Height == "Low",
  logical_method = "NA_if_not_all_TRUE"
)
diff_fct_diff_class(
  data_fungi@sam_data$Height == "Low",
  logical_method = "FALSE_if_not_all_TRUE"
)
diff_fct_diff_class(
  data_fungi@sam_data$Height,
  character_method = "unique_or_na"
)
diff_fct_diff_class(
  c("IE", "IE"),
  character_method = "unique_or_na"
)
diff_fct_diff_class(
  c("IE", "IE", "TE", "TE"),
  character_method = "more_frequent"
)
diff_fct_diff_class(
  c("IE", "IE", "TE", "TE"),
  character_method = "more_frequent_without_equality"
)
}
\author{
Adrien Taudière
}
