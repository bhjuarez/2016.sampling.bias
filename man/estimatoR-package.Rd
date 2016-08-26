\name{estimatoR-package}
\alias{estimatoR-package}
\alias{estimatoR}
\docType{package}
\title{
Estimates Number of Samples Needed to Change Significance for Phylogenetic and Non-Phylogenetic Statistical Tests
}
\description{
Tests for effect of known missing samples on significance of regressions. It does so by iteratively simulating additional data (and species on the phylogeny, if included), until significance of the regression changes.
}
\details{
\tabular{ll}{
Package: \tab estimatoR\cr
Type: \tab Package\cr
Version: \tab 1.0\cr
Date: \tab 2016-08-25\cr
License: \tab GPL-3\cr
}
Tests for effect of known missing samples on significance of regressions. The function 'estimate.lm' will do so without incorporating phylogeny and just simulate new data, whereas 'estimate.pgls' will include the phylogeny and also randomply add new tips for the new data points.
}
\author{
Bryan H. Juarez, Orlando Schwery, Giulio Valentino Dalla Riva

Maintainer: Bryan H. Juarez <bhjuarez@iastate.edu>
}
\keyword{ package }
\seealso{
\code{\link{estimate.lm}}, \code{\link{estimate.pgls}}
}
