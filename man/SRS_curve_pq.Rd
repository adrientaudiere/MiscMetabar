% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/plot_functions.R
\name{SRS_curve_pq}
\alias{SRS_curve_pq}
\title{Scaling with ranked subsampling (SRS) curve of phyloseq object}
\usage{
SRS_curve_pq(physeq, clean_pq = FALSE, ...)
}
\arguments{
\item{physeq}{(required): a \code{\link{phyloseq-class}} object obtained
using the \code{phyloseq} package.}

\item{clean_pq}{(logical): Does the phyloseq
object is cleaned using the \code{\link[=clean_pq]{clean_pq()}} function?}

\item{...}{Other arguments passed on to \code{SRS::SRScurve()}}
}
\value{
A plot
}
\description{
\verb{r lifecycle::badge("experimental")}
}