# Main script

## load required packages
require(ape)
require(broom)
require(geiger)
require(geomorph)
require(phytools)

## initialize objects and run analysis

### for standard lm
source("./R/estimate.lm.R")

dat <- cbind(rnorm(50, 50, 125), rep(1, 50))  #Generate data 1
dat <- rbind(dat, cbind(rnorm(50, 110, 90), rep(2, 50)))  #Generate data 2 and bind
colnames(dat) <- c("length", "group")  #Assing column names
dat <- as.data.frame(dat)
test <- lm(dat[, 1] ~ dat [, 2])  #Run initial regression

#' test running
estimate.lm(test, 200)  #Run function

### for phylogenetic lm
source("./R/estimate.pgls.R")

f1 <- length ~ as.factor(group)  # define function
dat <- cbind(rnorm(50, 50, 105), rep(1, 50))  # simulate data group 1
dat <- rbind(dat, cbind(rnorm(50, 110, 90), rep(2, 50)))  # simulate data group 2
colnames(dat) <- c("length", "group")  # give colnames
rownames(dat) <- paste("t", seq(1:100), sep = "")  # name species to match tree tip labels later
tree <- ape::rtree(100, rooted = T, tip.label = paste("t", seq(1:100), sep = ""))  # simulate tree with matching tiplabels
gdf <- geomorph.data.frame(length = dat[, 1], group = dat[, 2], phy = tree)
est <- 10  # set number of estimates (best set to number of known unsampled taxa)

#' test running
estimate.pgls(f1, gdf, est)



###Simulations
data_vec <- dat[,1]
names(data_vec) <- rownames(dat)

BM.stats <- geiger::fitContinuous(tree, data_vec)$opt
brownian_simul <- phytools::fastBM(tree, BM.stats$z0, sig2 = BM.stats$sigsq)
reg <- summary(lm(length~as.factor(group), data = as.data.frame(dat)))

#treedat <- treedata(tree, as.matrix(dat))
#reg <- procD.pgls(treedat$data[, 1] ~ treedat$data[, 2], phy = treedat$phy)
