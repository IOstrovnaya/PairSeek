
PairSeekBinary<-function(y.anal,  data.anal,LassoIterations= 50) #seans method for dichotomized
{
  ####################################
  #################################### first get distribution of lambda 
  
  
  grpA <- sample((1:dim(data.anal)[1]), floor(dim(data.anal)[1]/2)) ### Grp A: fitting lasso on these
  grpB <- (1:dim(data.anal)[1])[!((1:dim(data.anal)[1]) %in% grpA)] ### Grp B: no lasso, just used to find c
  
  datagrpA <- data.anal[grpA,]
  datagrpB <- data.anal[grpB,]
  
  ygrpB <-  1*y.anal[grpB]
  ygrpA <-  1*y.anal[grpA]
  
  
  pairsmatrixnames <-  NULL 
  for(i in 1:(dim(data.anal)[2]-1)){
    for(j in (i+1):dim(data.anal)[2]){
      pairsmatrixnames <- c(pairsmatrixnames, paste(colnames(data.anal)[i], ".X.",colnames(data.anal)[j], sep=""))
    }
  }
  
  xtest <- getpairsOddXOut( datagrpB, datagrpA, ygrpB, choose(dim(datagrpB)[2],2))
  colnames(xtest)  <- pairsmatrixnames
  
  mod1 <- glmnet(xtest, ygrpA,family="binomial",  alpha = 1)
  lambdaforC2 <- mod1$lambda
  
  
  
  retmatrix <-  NULL 
  
  for(simsim in 1:LassoIterations){
    
    grpA <- sample((1:dim(data.anal)[1]), floor(dim(data.anal)[1]/2)) ### Grp A: fitting lasso on these
    grpB <- (1:dim(data.anal)[1])[!((1:dim(data.anal)[1]) %in% grpA)] ### Grp B: no lasso, just used to find c
    
    datagrpA <- data.anal[grpA,]
    datagrpB <- data.anal[grpB,]
    
    ygrpB <-  1*y.anal[grpB]
    ygrpA <-  1*y.anal[grpA]
    
    xtest <- getpairsOddXOut( datagrpB, datagrpA, ygrpB, choose(dim(datagrpB)[2],2))
    colnames(xtest)  <- pairsmatrixnames
    
    mod1 <- glmnet(xtest, ygrpA,family="binomial",  alpha = 1)
    
    selected <- predict(mod1, type = "nonzero", s= sample(lambdaforC2,1))
    selected <- selected[[length(selected)]]
    
    ret <- logical(ncol(xtest))
    ret[selected] <- TRUE
    names(ret) <- colnames(xtest)
    
    retmatrix  <- rbind(retmatrix, ret)
    
    
    
    
  }
  
  
  return( apply(retmatrix != 0, 2, mean))
  
  
}
