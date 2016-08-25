\name{estimate.lm}
\alias{estimate.lm}
\title{
Estimates Number of Samples Needed to Change Significance for Regression
}
\description{
Tests for effect of known missing samples on significance of regression, without incorporating phylogeny. It does so by iteratively simulating additional data until the significance of the regression changes.
}
\usage{
estimate.lm(mod.obj, est)
}
\arguments{
  \item{mod.obj}{
Regression object made from data frame containing, data in first column and group assignment (using integers) in second column.
}
  \item{est}{
Number of missing taxa of the focal group (if known). This value will act as a maximum limit for iterations adding new species.
}
}
\value{
Number of iterations (which corresponds to missing taxa to be sampled).
}
\author{
Bryan H. Juarez, Orlando Schwery, Giulio Valentino Dalla Riva
}
\seealso{
\code{\link{estimate.pgls}}, \code{\link{estimatoR-package}}
}
\examples{
  dat <- cbind(rnorm(50, 50, 125), rep(1, 50))  #Generate data 1
  dat <- rbind(dat, cbind(rnorm(50, 110, 90), rep(2, 50)))  #Generate data 2 and bind
  colnames(dat) <- c("length", "group")  #Assing column names
  dat <- as.data.frame(dat)
  test <- lm(dat[, 1] ~ dat [, 2])  #Run initial regression
  estimate.lm(test, 200)  # test running
}
