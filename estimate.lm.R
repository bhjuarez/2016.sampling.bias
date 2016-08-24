require(ape)
require(geiger)
require(geomorph)
require(phytools)

fun <- dat[, 1] ~ as.vector(dat[, 2])
dat <- cbind(rnorm(50, 50, 125), rep(1, 50))
dat <- rbind(dat, cbind(rnorm(50, 110, 90), rep(2, 50)))
colnames(dat) <- c("length", "group")
dat <- as.data.frame(dat)
estimate.lm(fun, dat, 200)

estimate.lm <- function(fun, dat, est, method){
  method <- "lm"
  if (method == "procD"){
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
  else {
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

}















