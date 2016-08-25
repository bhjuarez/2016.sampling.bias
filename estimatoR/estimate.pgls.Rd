\name{estimate.pgls}
\alias{estimate.pgls}
\title{
Estimates Number of Samples Needed to Change Significance for Phylogenetic Regression
}
\description{
Tests for effect of known missing samples on significance of phylogenetic regression. It does so by iteratively simulating additional data and species on the phylogeny, until the significance of the regression changes.
}
\usage{
estimate.pgls(dat, phy, est)
}
\arguments{
  \item{dat}{
Data frame containing, data in first column and group assignment (using integers) in second column. Rownames corresponding to tip labels of tree.
}
  \item{phy}{
An object of type 'phy', a rooted phylogeny. Can be ultrametric or not. Tip labels corresponding to rownames of data table.
}
  \item{est}{
Number of missing taxa of the focal group (if known). This value will act as a maximum limit for iterations adding new species.
}
}
\value{
Matrix with number of iterations (which corresponds to missing taxa to be sampled) in first column and new p-values for each iteration in second column.
}
\author{
Bryan H. Juarez, Orlando Schwery, Giulio Valentino Dalla Riva
}
\seealso{
\code{\link{estimate.lm}}, \code{\link{estimatoR-package}}
}
\examples{
  fun <- dat[, 1] ~ as.vector(dat[, 2])  # define function
  dat <- cbind(rnorm(50, 50, 105), rep(1, 50))  # simulate data group 1
  dat <- rbind(dat, cbind(rnorm(50, 110, 90), rep(2, 50)))  # simulate data group 2
  colnames(dat) <- c("length", "group")  # give colnames
  rownames(dat) <- paste("t", seq(1:100), sep = "")  # name species to match tree tip labels later
  dat <- as.data.frame(dat)  # convert to data frame
  est <- 10  # set number of estimates (best set to number of known unsampled taxa)
  tree <- rtree(100, rooted = T, tip.label = paste("t", seq(1:100), sep = ""))  # simulate tree with matching tiplabels
  estimate.pgls(dat, tree, est) # test running
}
