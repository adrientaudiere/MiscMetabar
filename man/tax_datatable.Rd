% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/plot_functions.R
\name{tax_datatable}
\alias{tax_datatable}
\title{Make a datatable with the taxonomy of a \code{\link{phyloseq-class}} object}
\usage{
tax_datatable(physeq, abundance = TRUE, ...)
}
\arguments{
\item{physeq}{(required): A \code{\link{phyloseq-class}} object}

\item{abundance}{(Default: TRUE): Does the number of sequences is print}

\item{...}{Other argument for the datatable function}
}
\value{
A datatable
}
\description{
Make a datatable with the taxonomy of a \code{\link{phyloseq-class}} object
}
\examples{
data("GlobalPatterns")
tax_datatable(subset_taxa(GlobalPatterns,
              rowSums(GlobalPatterns@otu_table)>10000))
}
\author{
Adrien Taudière
}
