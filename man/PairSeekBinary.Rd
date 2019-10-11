\name{PairSeekBinary}
\alias{PairSeekBinary}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
PairSeekBinary(y.anal, data.anal, LassoIterations = 50)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{y.anal}{
%%     ~~Describe \code{y.anal} here~~
}
  \item{data.anal}{
%%     ~~Describe \code{data.anal} here~~
}
  \item{LassoIterations}{
%%     ~~Describe \code{LassoIterations} here~~
}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
%%  ~Describe the value returned
%%  If it is a LIST, use
%%  \item{comp1 }{Description of 'comp1'}
%%  \item{comp2 }{Description of 'comp2'}
%% ...
}
\references{
%% ~put references to the literature/web site here ~
}
\author{
%%  ~~who you are~~
}
\note{
%%  ~~further notes~~
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
}
\examples{
##---- Should be DIRECTLY executable !! ----
##-- ==>  Define data, use random,
##--	or do  help(data=index)  for the standard data sets.

## The function is currently defined as
function (y.anal, data.anal, LassoIterations = 50) 
{
    grpA <- sample((1:dim(data.anal)[1]), floor(dim(data.anal)[1]/2))
    grpB <- (1:dim(data.anal)[1])[!((1:dim(data.anal)[1]) \%in\% 
        grpA)]
    datagrpA <- data.anal[grpA, ]
    datagrpB <- data.anal[grpB, ]
    ygrpB <- 1 * y.anal[grpB]
    ygrpA <- 1 * y.anal[grpA]
    pairsmatrixnames <- NULL
    for (i in 1:(dim(data.anal)[2] - 1)) {
        for (j in (i + 1):dim(data.anal)[2]) {
            pairsmatrixnames <- c(pairsmatrixnames, paste(colnames(data.anal)[i], 
                ".X.", colnames(data.anal)[j], sep = ""))
        }
    }
    xtest <- getpairsOddXOut(datagrpB, datagrpA, ygrpB, choose(dim(datagrpB)[2], 
        2))
    colnames(xtest) <- pairsmatrixnames
    mod1 <- glmnet(xtest, ygrpA, family = "binomial", alpha = 1)
    lambdaforC2 <- mod1$lambda
    retmatrix <- NULL
    for (simsim in 1:LassoIterations) {
        grpA <- sample((1:dim(data.anal)[1]), floor(dim(data.anal)[1]/2))
        grpB <- (1:dim(data.anal)[1])[!((1:dim(data.anal)[1]) \%in\% 
            grpA)]
        datagrpA <- data.anal[grpA, ]
        datagrpB <- data.anal[grpB, ]
        ygrpB <- 1 * y.anal[grpB]
        ygrpA <- 1 * y.anal[grpA]
        xtest <- getpairsOddXOut(datagrpB, datagrpA, ygrpB, choose(dim(datagrpB)[2], 
            2))
        colnames(xtest) <- pairsmatrixnames
        mod1 <- glmnet(xtest, ygrpA, family = "binomial", alpha = 1)
        selected <- predict(mod1, type = "nonzero", s = sample(lambdaforC2, 
            1))
        selected <- selected[[length(selected)]]
        ret <- logical(ncol(xtest))
        ret[selected] <- TRUE
        names(ret) <- colnames(xtest)
        retmatrix <- rbind(retmatrix, ret)
    }
    return(apply(retmatrix != 0, 2, mean))
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }% use one of  RShowDoc("KEYWORDS")
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
