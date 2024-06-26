% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/colocalisation.R
\name{ggRedRibbonColoc}
\alias{ggRedRibbonColoc}
\title{Plot a colocalisation with ggplot}
\usage{
ggRedRibbonColoc(
  self,
  plot.order = 1:4,
  show.title = TRUE,
  labels = NULL,
  tss = NULL,
  shortid = NULL,
  title = NULL,
  highlight = c(),
  .log10 = FALSE
)
}
\arguments{
\item{self}{a colocRedRibbon object}

\item{plot.order}{a vector specifying the plot order (default =1:4,  1 = RedRibbon plot, 2 =  manhantan plot for \code{a}, 3 = plot for \code{a}and 'b\verb{, 4 = manhantan plot for }b`).}

\item{show.title}{shows the title (default = TRUE)}

\item{labels}{axis labels}

\item{tss}{transcription start site}

\item{shortid}{name of the gene}

\item{title}{the title of the plot}

\item{highlight}{a list of SNPs to highlight}

\item{.log10}{output log10 pval (default = FALSE)}
}
\value{
ggplot object
}
\description{
Plot a colocalisation with ggplot
}
