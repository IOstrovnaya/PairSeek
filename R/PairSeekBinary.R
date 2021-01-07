PairSeekBinary<-function(y,  data,LassoIterations= 50)
{library(Rcpp)
  library(glmnet)
  ####################################
  #################################### first get distribution of lambda
  sourceCpp("src/newpairingOddXoutEntropy.cpp")
  grpA <- sample((1:dim(data)[1]), floor(dim(data)[1]/2)) ### Grp A: fitting lasso on these
  grpB <- (1:dim(data)[1])[!((1:dim(data)[1]) %in% grpA)] ### Grp B: no lasso, just used to find c
  datagrpA <- data[grpA,]
  datagrpB <- data[grpB,]
  ygrpB <-  1*y[grpB]
  ygrpA <-  1*y[grpA]
  pairsmatrixnames <-  NULL
  for(i in 1:(dim(data)[2]-1)){
    for(j in (i+1):dim(data)[2]){
      pairsmatrixnames <- c(pairsmatrixnames, paste(colnames(data)[i], ".X.",colnames(data)[j], sep=""))
    }
  }
  xtest <-getpairsOddXOut(datagrpB, datagrpA, ygrpB, choose(dim(datagrpB)[2],2))
  colnames(xtest)  <- pairsmatrixnames
  mod1 <- glmnet(xtest, ygrpA,family="binomial",  alpha = 1)
  lambdaforC2 <- mod1$lambda
  retmatrix <-  NULL
  cat("Percent computed: ")
  for(simsim in 1:LassoIterations){
    if (simsim %in% round(c(1:10)*LassoIterations/10)) cat(paste(round(100*simsim/LassoIterations),".. ",sep=""))
    grpA <- sample((1:dim(data)[1]), floor(dim(data)[1]/2)) ### Grp A: fitting lasso on these
    grpB <- (1:dim(data)[1])[!((1:dim(data)[1]) %in% grpA)] ### Grp B: no lasso, just used to find c
    datagrpA <- data[grpA,]
    datagrpB <- data[grpB,]
    ygrpB <-  1*y[grpB]
    ygrpA <-  1*y[grpA]
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
