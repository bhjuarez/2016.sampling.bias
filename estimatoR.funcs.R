require(ape)
require(geiger)
require(geomorph)
require(phytools)

fun <- dat[, 1] ~ as.vector(dat[, 2])
dat <- cbind(rnorm(50, 50, 125), rep(1, 50))
dat <- rbind(dat, cbind(rnorm(50, 110, 90), rep(2, 50)))
colnames(dat) <- c("length", "group")
dat <- as.data.frame(dat)
est <- 200

estimate.lm.norm(fun, dat, est)

estimate.lm <- function(fun, dat, est){
	reg <- procD.lm(fun)
	p <- reg$aov.table[1, 7]

	if(p > 0.05){
		stop("Regression is not significant")
	}
	
	mean <- reg$coeff[1, 1]
	sd <- sqrt(var(dat[, 1][which(dat[, 2] == 1, )]))
	n <- 0
	i <- 0
	stats <- matrix(ncol = 2, nrow = 0)
	
	while(n < est & p < 0.05){
  		n <- n + 1
  		i <- i +1
  		sim <- rnorm(1, mean, sd)
  		sub.dat <- data.frame(sim, 1)
  		colnames(sub.dat) <- colnames(dat)
  		dat <- rbind(dat, sub.dat)
  		sub.reg <- procD.lm(length ~ group, data = dat)
  		p <- sub.reg$aov.table[1, 7]
  		stats <- rbind(stats, c(n, p))	
	}
	return(stats)
}

estimate.lm.norm <- function(fun, dat, est){
  reg <- summary(lm(fun))
  p <- reg$coeff[2, 4]
  
  if(p > 0.05){
    stop("Regression is not significant")
  }
  
  mean <- reg$coeff[1, 1]
  sd <- sqrt(var(dat[, 1][which(dat[, 2] == 1, )]))
  n <- 0
  i <- 0
  stats <- matrix(ncol = 2, nrow = 0)
  
  while(n < est & p < 0.05){
    n <- n + 1
    i <- i +1
    sim <- rnorm(1, mean, sd)
    sub.dat <- data.frame(sim, 1)
    colnames(sub.dat) <- colnames(dat)
    dat <- rbind(dat, sub.dat)
    sub.reg <- summary(lm(dat[, 1] ~ as.factor(dat[, 2])))
    p <- sub.reg$coeff[2, 4]
    stats <- rbind(stats, c(n, p))	
  }
  return(stats)
}

estimate.lm.norm(fun, dat, est)






fun <- dat[, 1] ~ as.vector(dat[, 2])
dat <- cbind(rnorm(50, 50, 105), rep(1, 50))
dat <- rbind(dat, cbind(rnorm(50, 110, 90), rep(2, 50)))
colnames(dat) <- c("length", "group")
rownames(dat) <- paste("t", seq(1:100), sep = "")
dat <- as.data.frame(dat)
est <- 10

tree <- rtree(100, rooted = T, tip.label = paste("t", seq(1:100), sep = ""))
estimate.pgls <- function(dat, phy, est){
	###Add data to new tips
	pre.tree <- tree
	n <- 100
	for( i in 1:n){
  	tree.labs <- paste("t", length(pre.tree$tip.label)+i, sep = "")
  	node <- round(runif(1, 1, Nnode(pre.tree)+length(tree$tip.label)))
  	position <- runif(1)*pre.tree$edge.length[which(pre.tree$edge[, 2]==node)]
  	new.tree <- bind.tip(pre.tree, tree.labs, where = node, position = position )
  	pre.tree <- new.tree
	}
	
		reg <- procD.pgls(dat[, 1] ~ dat[, 2], phy = tree)
	mean <- reg$coeff[1, 1]
	sd <- sqrt(var(dat$length[which(dat$group == 1, )]))
	p <- reg$coeff[2, 4]
	n <- 0
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









