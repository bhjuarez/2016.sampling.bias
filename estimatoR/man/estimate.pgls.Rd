\name{estimate.pgls}
\alias{estimate.pgls}
\title{
%%  ~~function to do ... ~~
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
estimate.pgls(dat, phy, est)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{dat}{
%%     ~~Describe \code{dat} here~~
}
  \item{phy}{
%%     ~~Describe \code{phy} here~~
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
function (dat, phy, est)
{
    reg <- summary(lm(length ~ as.factor(group), data = dat))
    mean <- reg$coeff[1, 1]
    sd <- sqrt(var(dat$length[which(dat$group == 1, )]))
    p <- reg$coeff[2, 4]
    pre.tree <- tree
    n <- 0
    stats <- matrix(ncol = 2, nrow = 0)
    while (n < est & p < 0.05) {
        n <- n + 1
        tree.labs <- paste("t", length(pre.tree$tip.label) +
            n, sep = "")
        node <- round(runif(1, 1, Nnode(pre.tree) + length(tree$tip.label)))
        position <- runif(1) * pre.tree$edge.length[which(pre.tree$edge[,
            2] == node)]
        if (is.ultrametric(pre.tree) == TRUE) {
            new.tree <- bind.tip(pre.tree, tree.labs, where = node,
                position = position)
        }
        else {
            edge.length <- runif(1, min = (min(pre.tree$edge.length)),
                max = (max(pre.tree$edge.length)))
            new.tree <- bind.tip(pre.tree, tree.labs, edge.length,
                where = node, position = position)
        }
        pre.tree <- new.tree
        sim <- rnorm(1, mean, sd)
        sub.dat <- data.frame(sim, 1)
        colnames(sub.dat) <- c("length", "group")
        dat <- rbind(dat, sub.dat)
        sub.reg <- summary(lm(length ~ as.factor(group), data = dat))
        p <- sub.reg$coeff[2, 4]
        stats <- rbind(stats, c(n, p))
    }
    print(stats)
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }% use one of  RShowDoc("KEYWORDS")
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
