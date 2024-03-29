% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/colocalisation.R
\name{coloc.RedRibbonColoc}
\alias{coloc.RedRibbonColoc}
\title{Compute a colocalisation}
\usage{
\method{coloc}{RedRibbonColoc}(self, n.reduce = min)
}
\arguments{
\item{self}{a colocRedRibbon object}

\item{n.reduce}{function to reduce the number of sample columns for a and b into a number (Default: min)
This is needed as some eQTL/GWAS toolchains output an effective number of samples by SNP.}
}
\value{
RedRibbonColoc object
}
\description{
Compute a colocalisation
}
