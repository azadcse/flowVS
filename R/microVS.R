## ========================================================================================
## Estimate optimum parameters for variance stabilization in microarray data.
## Input:
#     data: The microarray data in a Matrix. 
#     cfLow, cfHigh: lowest and highest possible values for cofactor (log scale)
#     frac: fraction of differentially expressed genes used in variance stabilization (<0 & >=1)
# output: optimum cofactor
## ========================================================================================

microVS = function(data, cfLow=0, cfHigh=10, frac=1)
{
  ## some error checking
  if(!is(data, "matrix"))
    stop(" The microarray data must be a Matrix.")
  if(frac>1 || frac <= 0)
    stop(" 0< frac<=1 ")
  if(cfLow>=cfHigh)
  {
    print("Warning: cfLow>=cfHigh, using default values")
    cfLow=0
    cfHigh=10
  }
     
  cat("====================================================================\n")
  cat("Finding optimum cofactor for asinh transformation\n")
  cat("====================================================================\n")
  cat(sprintf("%15s     %15s \n", "cofactor(log scale)", "Bartlett\'s stat"))
  cat("====================================================================\n")
  
  cofactors = seq(cfLow,cfHigh,1)
  bartlett = NULL
  for(cf in cofactors)
  {
    data.t = asinh(data/exp(cf))
    
    #find the non-differentially exprssed genes
    if(frac<1)
    {
      diff = quantile(abs(data.t[,1] - data.t[,2]), frac)
      keep.idx = which(abs(data.t[,1] - data.t[,2]) <= diff)
      data.t = data.t[keep.idx,]
    }
    bt=bartlettTestMicro(data.t)
    bartlett = c(bartlett, bt)
    cat(sprintf("%10d %25.2f \n", cf, bt))
  }
  
  minIdx = which.min(bartlett)
  cat("\n Optimum cofactor :", sprintf("exp(%d)",cofactors[minIdx]), "\n")
  cat("====================================================================\n\n")
  
  plot(cofactors, bartlett, type='o', pch=16, 
       xlab="Cofactors (log scale)", ylab="Bartlett's statistics",
       main = paste("Optimum cofactor: ", 
                    sprintf("exp(%d)",cofactors[minIdx]), sep=""))
  points(cofactors[minIdx], bartlett[minIdx], pch=16, col='red')
  
  data.t = asinh(data/exp(cofactors[minIdx]))
  return(data.t)
}




##============================================================
## Internal function
## Compute Bartlett's statistics from a transformed microarray data
##============================================================
bartlettTestMicro = function(dataMatrix) 
{
  means = rowMeans(dataMatrix)
  k = nrow(dataMatrix)
  n = ncol(dataMatrix)  # all genes have equal number of samples 
  N = n*k
  vars = apply(dataMatrix,1,var)
  var.pooled = sum(vars) / (k-1) 
  numerator = (N-k) * log(var.pooled) - (n-1) * sum(log(vars))
  denom = 1 + ( ( k/(n-1) - 1/(N-k) ) / (3* (k-1) ) )
  bt = numerator/denom 
  return(bt)
}



## ==========================================================================
## meanSdPlot method for matrix modified from vsn package 
## - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
plotMeanSd=function(x, ranks=TRUE, xlab = ifelse(ranks, "Rank of means (ascending order)", "mean"),
                     ylab = "Standard deviation", pch  = ".", plot = TRUE, ...) 
  {            
              if(!is(x, "matrix"))
                stop("\'x\' must be a mtrix\n")
              stopifnot(is.logical(ranks), length(ranks)==1, !is.na(ranks))
              
              n = nrow(x)
              if(n==0L) {
                warning("In 'meanSdPlot': matrix has 0 rows, there is nothing to be done.")
                return()
              }
              
              px   = rowMeans(x, na.rm=TRUE)
              
              sqr     = function(x)  x*x
              rs       = rowSums(!is.na(x))
              rs[rs<1]  = NA
              py = sqrt(rowSums(sqr(x-px),na.rm=TRUE)/(rs-1))
              #py   = sqrt(rowV(x, mean=px, na.rm=TRUE))
              rpx  = rank(px, na.last=FALSE, ties.method = "random")
              
              ## run median with centers at dm,2*dm,3*dm,... and width 2*dm
              dm        = 0.05
              midpoints = seq(dm, 1-dm, by=dm)
              within    = function(x, x1, x2) { x>=x1 & x<=x2 }
              mediwind  = function(mp) median(py[within(rpx/n, mp-dm, mp+dm)], na.rm=TRUE)
              rq.sds    = sapply(midpoints, mediwind)
              
              res = if(ranks) {
                list(rank=midpoints*n, sd=rq.sds, px=rpx, py=py)
              } else {
                list(quantile=quantile(px, probs=midpoints, na.rm=TRUE), sd=rq.sds, px=px, py=py)
              }
              
              if(plot) {
                plot(res$px, res$py, pch=pch, xlab=xlab, cex=2, ylab=ylab, ...)
                #smoothScatter(res$px, res$py, xlab=xlab, ylab=ylab, ...)
                lines(res[[1L]], res$sd, col="red", pch=19, lwd=5, cex=1)
              }
              
              return(invisible(res))
}
