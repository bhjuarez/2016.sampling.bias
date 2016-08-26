estimate.pgls <- function(dat, phy, est) {  # main function checking number of taxa to be added

  #run regression on initial tree & data set
  dat <- geomorph.data.frame(dat)
  reg <- procD.pgls(dat[[1]][, 1] ~ dat[[1]][, 2], phy)  # run regression
  #reg <- summary(lm(length~as.factor(group), data = dat))
  mean <- reg$coeff[1, 1]  # save initial mean
  sd <- sqrt(var(dat$length[which(dat$group == 1, )]))  # save initial SD
  p <- reg$coeff[2, 4]  # save initial p-value

  pre.tree <- tree #define temporary tree to preserve original
  n <- 0   #set initial iteration to 0
  stats <- matrix(ncol = 2, nrow = 0)  # initiate stats table
#  for( i in 1:est) {  # loop adding random tips to tree
  while(n < est & p < 0.05) {
    n <- n + 1
    tree.labs <- paste("t", length(pre.tree$tip.label) + n, sep = "")  # generate label that corresponds to format of existing ones (continuing the numeration)
    node <- round(runif(1, 1, Nnode(pre.tree)+length(tree$tip.label)))  # pick node at random
    position <- runif(1)*pre.tree$edge.length[which(pre.tree$edge[, 2]==node)]  # pick location on node at random
    if (is.ultrametric(pre.tree) == TRUE) {  # go without BL if tree is ultrametric
      new.tree <- bind.tip(pre.tree, tree.labs, where = node, position = position )  # add new tip to tree at selected node and position
    } else {  # go with assigned BL if tree is not ultrametric
      edge.length <- runif(1,  min=(min(pre.tree$edge.length)), max=(max(pre.tree$edge.length)))  # create BL for new node, picked from uniform distribution with existing BLs as bounds
      new.tree <- bind.tip(pre.tree, tree.labs, edge.length, where = node, position = position )  # add new tip to tree at selected node and position
    }
    pre.tree <- new.tree  # overwrite temporary tree with the one which has added node (so it carries over to next iteration)

    sim <- rnorm(1, mean, sd)
    sub.dat <- data.frame(sim, 1)
    colnames(sub.dat) <- c("length", "group")
    dat <- rbind(dat, sub.dat)
    sub.reg <- summary(lm(length~as.factor(group), data = dat))
    #sub.reg <- procD.pgls(dat[, 1] ~ dat[, 2], phy = new.tree)  # run regression
    p <- sub.reg$coeff[2, 4]
    stats <- rbind(stats, c(n, p))
  }
  print(stats)
}