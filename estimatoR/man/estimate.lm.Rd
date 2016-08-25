\name{estimate.lm}
\alias{estimate.lm}
\title{
%%  ~~function to do ... ~~
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
estimate.lm(mod.obj, est)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{mod.obj}{
%%     ~~Describe \code{mod.obj} here~~
}
  \item{est}{
%%     ~~Describe \code{est} here~~
}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
%%  ~Describe the value returned
%%  If it is a LIST, use
%%  \item{comp1 }{Description of 'comp1'}
%%  \item{comp2 }{Description of 'comp2'}
%% ...
}
\references{
%% ~put references to the literature/web site here ~
}
\author{
%%  ~~who you are~~
}
\note{
%%  ~~further notes~~
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
}
\examples{
##---- Should be DIRECTLY executable !! ----
##-- ==>  Define data, use random,
##--	or do  help(data=index)  for the standard data sets.

## The function is currently defined as
function (mod.obj, est)
{
    if (class(mod.obj) != "lm" & class(mod.obj) != "procD.lm") {
        stop("Please provide a lm or procD.lm model object")
    }
    if (class(mod.obj) == "lm") {
        p <- broom::glance(mod.obj)$p.value
        if (p > 0.05) {
            stop("Regression is not significant")
        }
        dat <- mod.obj$model
        mean <- test[[1]][[1]]
        sd <- sqrt(var(dat[, 1][which(dat[, 2] == 1, )]))
        n <- 0
        stats <- matrix(ncol = 1, nrow = 0)
        while (n < est & p < 0.05) {
            n <- n + 1
            sim <- rnorm(1, mean, sd)
            sub.dat <- data.frame(sim, 1)
            colnames(sub.dat) <- colnames(dat)
            dat <- rbind(dat, sub.dat)
            sub.p <- summary(lm(dat[, 1] ~ dat[, 2]))$coeff[2,
                4]
            stats <- rbind(stats, sub.p)
        }
        return(dim(stats)[1])
    }
    else {
        reg <- procD.lm(fun)
        p <- reg$aov.table[1, 7]
        if (p > 0.05) {
            stop("Regression is not significant")
        }
        mean <- reg$coeff[1, 1]
        sd <- sqrt(var(dat[, 1][which(dat[, 2] == 1, )]))
        n <- 0
        i <- 0
        stats <- matrix(ncol = 1, nrow = 0)
        while (n <= est & p < 0.05) {
            n <- n + 1
            i <- i + 1
            sim <- rnorm(1, mean, sd)
            sub.dat <- data.frame(sim, 1)
            colnames(sub.dat) <- colnames(dat)
            dat <- rbind(dat, sub.dat)
            sub.reg <- procD.lm(length ~ group, data = dat)
            p <- sub.reg$aov.table[1, 7]
            stats <- rbind(stats, p)
        }
        return(list(stats, dim(stats)[1]))
    }
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }% use one of  RShowDoc("KEYWORDS")
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
