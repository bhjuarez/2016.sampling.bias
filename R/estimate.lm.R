estimate.lm <- function(mod.obj, est){
  if (class(mod.obj) != "lm" & class(mod.obj) != "procD.lm"){  #Check if user provided correct object
    stop("Please provide a lm or procD.lm model object")
  }
  
  if (class(mod.obj) == "lm"){  #Implement lm
    p <- broom::glance(mod.obj)$p.value
    if(p > 0.05){
      stop("Regression is not significant")
    }
    
    dat <- mod.obj$model
    mean <- test[[1]][[1]]
    sd <- sqrt(var(dat[, 1][which(dat[, 2] == 1, )]))  #Make multivariate ind and dep
    n <- 0
    stats <- matrix(ncol = 1, nrow = 0)
    
    while(n < est & p < 0.05){
      n <- n + 1
      sim <- rnorm(1, mean, sd)
      sub.dat <- data.frame(sim, 1)
      colnames(sub.dat) <- colnames(dat)
      dat <- rbind(dat, sub.dat)
      sub.p <- summary(lm(dat[, 1] ~ dat[, 2]))$coeff[2, 4] #make multivariate
      stats <- rbind(stats, sub.p)	
    }
    return(dim(stats)[1])
  }
  
  else {  #Implement procD.lm
    reg <- procD.lm(fun)
    p <- reg$aov.table[1, 7]
    
    if(p > 0.05){
      stop("Regression is not significant")
    }
    
    mean <- reg$coeff[1, 1]
    sd <- sqrt(var(dat[, 1][which(dat[, 2] == 1, )]))
    n <- 0
    i <- 0
    stats <- matrix(ncol = 1, nrow = 0)
    
    while(n <= est & p < 0.05){
      n <- n + 1
      i <- i +1
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