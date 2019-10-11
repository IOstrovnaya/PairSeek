PairSeekLogRatio<-function(y.anal,  data.anal,LassoIterations= 50) #sean's method for ratios
{
  
  
  ####################################
  #################################### first get distribution of lambda 
  
  
  
  
  
  pairsmatrixnames <-  NULL 
  for(i in 1:(dim(data.anal)[2]-1)){
    for(j in (i+1):dim(data.anal)[2]){
      pairsmatrixnames <- c(pairsmatrixnames, paste(colnames(data.anal)[i], ".X.",colnames(data.anal)[j], sep=""))
    }
  }
  
  data.pairs <- log(RatioX(data.anal, choose(dim(data.anal)[2],2))+1)
  colnames(data.pairs)  <- pairsmatrixnames
  
  
  grpA <- sample((1:dim(data.anal)[1]), floor(dim(data.anal)[1]/2)) ### Grp A: fitting lasso on these
  grpB <- (1:dim(data.anal)[1])[!((1:dim(data.anal)[1]) %in% grpA)] ### Grp B: no lasso, just used to find c
  
  datagrpA <- data.anal[grpA,]
  datagrpB <- data.anal[grpB,]
  
  ygrpB <-  1*y.anal[grpB]
  ygrpA <-  1*y.anal[grpA]
  mod1 <- glmnet(data.pairs[grpA,], ygrpA,family="binomial",  alpha = 1)
  lambdaforRato <- mod1$lambda
  
  ####################################
  #################################### first get distribution of lambda 
  
  retmatrixRatio <- NULL 
  
  for(simsim in 1:LassoIterations){
    
    grpA <- sample((1:dim(data.pairs)[1]), floor(dim(data.pairs)[1]/2)) ### Grp A: fitting lasso on these
    grpB <- (1:dim(data.pairs)[1])[!((1:dim(data.pairs)[1]) %in% grpA)] ### Grp B: no lasso, just used to find c
    
    datagrpA <- data.pairs[grpA,]
    datagrpB <- data.pairs[grpB,]
    
    ygrpB <-  1*y.anal[grpB]
    ygrpA <-  1*y.anal[grpA]
    
    
    
    
    ######### looking at just the Rank
    
    
    mod1 <- glmnet(datagrpA, ygrpA,family="binomial",  alpha = 1)
    
    selected <- predict(mod1, type = "nonzero", s= sample(lambdaforRato,1))
    selected <- selected[[length(selected)]]
    
    ret <- logical(ncol(datagrpA))
    ret[selected] <- TRUE
    names(ret) <- colnames(datagrpA)
    
    retmatrixRatio  <- rbind(retmatrixRatio, ret)
    
    
    
  }
  
  return(apply(retmatrixRatio != 0, 2, mean))
}
