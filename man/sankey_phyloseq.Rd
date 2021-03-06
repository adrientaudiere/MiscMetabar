% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/plot_functions.R
\name{sankey_phyloseq}
\alias{sankey_phyloseq}
\title{Sankey plot of \code{\link{phyloseq-class}} object}
\usage{
sankey_phyloseq(
  physeq = NULL,
  fact = NULL,
  taxa = c(1:4),
  nb_seq = FALSE,
  min_prop_tax = 0,
  tax2remove = NULL,
  units = NULL,
  symbol2sub = c("\\\\.", "-"),
  ...
)
}
\arguments{
\item{physeq}{(required): a \code{\link{phyloseq-class}} object.}

\item{fact}{(Optional): Name of the factor to cluster samples by modalities.
Need to be in \code{physeq@sam_data}.}

\item{taxa}{(Default: c(1:4)): a vector of taxonomic rank to plot}

\item{nb_seq}{(Default: FALSE): Represent the number of sequences or the
number of OTUs (nb_seq = FALSE). Note that ploting the number of sequences
is slower.}

\item{min_prop_tax}{(Default: 0): The minimum proportion for taxon to be
ploted. EXPERIMENTAL. For the moment each links below the min.prop.
tax is discard from the sankey network resulting in sometimes weird plot.}

\item{tax2remove}{: a vector of taxonomic groups to remove from the analysis
(e.g. \code{c('Incertae sedis', 'unidentified')})}

\item{units}{: character string describing physical units (if any) for Value}

\item{symbol2sub}{(default = c('\\.', '-')): vector of symbol to delete in
the taxonomy}

\item{...}{Additional arguments passed on to
\code{\link[networkD3]{sankeyNetwork}}}
}
\value{
A \code{\link[networkD3]{sankeyNetwork}} plot representing the
 taxonomic distribution of OTUs or sequences. If \code{fact} is set,
 represent the distribution of the last taxonomic level in the modalities
 of \code{fact}
}
\description{
Sankey plot of \code{\link{phyloseq-class}} object
}
\examples{
data("GlobalPatterns")
GP <- subset_taxa(GlobalPatterns, GlobalPatterns@tax_table[,1] == "Archaea")
sankey_phyloseq(GP, fact = 'SampleType')
sankey_phyloseq(GP, taxa = c(1:4), min_prop_tax = 0.01)
sankey_phyloseq(GP, taxa = c(1:4), min_prop_tax = 0.01, nb_seq = TRUE)
}
\seealso{
\code{\link[networkD3]{sankeyNetwork}}
}
\author{
Adrien Taudière
}
