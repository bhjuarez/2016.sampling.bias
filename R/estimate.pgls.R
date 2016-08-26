estimate.pgls <- function(f1, gdf, est){  # main function checking number of taxa to be added
  
  #run regression on initial tree & data set
  num.phy <- which(lapply(gdf, class)=='phylo')
  phy <- gdf[[-num.phy]]

  reg <- procD.pgls(f1, phy = phy, data = gdf)  # run regression
  mean <- reg$pgls.coefficients[[1]]  # save initial mean
  sd <- sqrt(var(reg$Y[which(reg$X[, 2] == 0)]))  # save initial SD
  p <- reg$aov.table[nrow(reg$aov.table) - 2, 7]  # save initial p-value

  if(p > 0.05){
      stop("Regression is not significant")
  }
    
  pre.tree <- phy #define temporary tree to preserve original
  n <- 0   #set initial iteration to 0
  stats <- matrix(ncol = 2, nrow = 0)  # initiate stats table

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
    dat <- gdf[-num.phy]
    colnames(sub.dat) <- names(dat)
    dat <- rbind(dat, sub.dat)
    rownames(dat) <- c(rownames(reg$Y), tree.labs)
    
    sub.reg <- procD.pgls(f1, phy = pre.tree, data = test)  # run regression
    p <- sub.reg$coeff[2, 4]
    stats <- rbind(stats, c(n, p))
  }
  print(stats)
}