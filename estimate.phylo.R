# load required packages
require(ape)
require(geiger)
require(geomorph)
require(phytools)

fun <- dat[, 1] ~ as.vector(dat[, 2])  # define function
dat <- cbind(rnorm(50, 50, 105), rep(1, 50))  # simulate data group 1
dat <- rbind(dat, cbind(rnorm(50, 110, 90), rep(2, 50)))  # simulate data group 2
colnames(dat) <- c("length", "group")  # give colnames
rownames(dat) <- paste("t", seq(1:100), sep = "")  # name species to match tree tip labels later
dat <- as.data.frame(dat)  # convert to data frame
est <- 10  # set number of estimates (best set to number of known unsampled taxa)

tree <- rtree(100, rooted = T, tip.label = paste("t", seq(1:100), sep = ""))  # simulate tree with matching tiplabels
ultratree <- chronoMPL(tree, se = TRUE, test = TRUE)  # create ultrametric version


estimate.pgls <- function(dat, phy, est){  # main function checking number of taxa to be added
  ###Add data to new tips
  pre.tree <- tree #define temporary tree to preserve original
  for( i in 1:est) {  # loop adding random tips to tree
    tree.labs <- paste("t", length(pre.tree$tip.label)+i, sep = "")  # generate label that corresponds to format of existing ones (continuing the numeration)
    node <- round(runif(1, 1, Nnode(pre.tree)+length(tree$tip.label)))  # pick node at random
    position <- runif(1)*pre.tree$edge.length[which(pre.tree$edge[, 2]==node)]  # pick location on node at random
    if (is.ultrametric(pre.tree) == TRUE) {  # go without BL if tree is ultrametric
      new.tree <- bind.tip(pre.tree, tree.labs, where = node, position = position )  # add new tip to tree at selected node and position
    } else {  # go with assigned BL if tree is not ultrametric
      edge.length <- runif(1,  min=(min(pre.tree$edge.length)), max=(max(pre.tree$edge.length)))  # create BL for new node, picked from uniform distribution with existing BLs as bounds
      new.tree <- bind.tip(pre.tree, tree.labs, edge.length, where = node, position = position )  # add new tip to tree at selected node and position
    }
    pre.tree <- new.tree  # overwrite temporary tree with the one which has added node (so it carries over to next iteration)
  }

  reg <- procD.pgls(dat[, 1] ~ dat[, 2], phy = tree)
  mean <- reg$coeff[1, 1]
  sd <- sqrt(var(dat$length[which(dat$group == 1, )]))
  p <- reg$coeff[2, 4]
  n <- 0  ______CHECKIFNEEDSCHANGE
  stats <- matrix(ncol = 2, nrow = 0)
  while(n < est & p < 0.05){
    n <- n + 1
    sim <- rnorm(1, mean, sd)
    sub.dat <- data.frame(sim, 1)
    colnames(sub.dat) <- c("length", "group")
    dat <- rbind(dat, sub.dat)
    sub.reg <- summary(lm(length~as.factor(group), data = dat))
    p <- sub.reg$coeff[2, 4]
    stats <- rbind(stats, c(n, p))
  }
  print(stats)
}

###Simulations
BM.stats <- fitContinuous(tree, dat[, 1])$opt
fastBM(tree, BM.stats$z0, sig2 = BM.stats$sigsq)
reg <- summary(lm(length~as.factor(group), data = as.data.frame(dat)))

#treedat <- treedata(tree, as.matrix(dat))
#reg <- procD.pgls(treedat$data[, 1] ~ treedat$data[, 2], phy = treedat$phy)
