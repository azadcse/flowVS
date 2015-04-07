## ========================================================================================
## Estimate optimum parameters for within-population variance stabilization.
## Input:
##    fs: A flowSet containing a collection of flow cyometry samples. 
##    channels: A character vector identifying the channels/dimensions to be transformed.
## Output: A numeric vector representing the optimum cofactors for 
##  the requested channels. The optimum cofactor for the input channels[i]
##  is stored in the ith entry of the returned vector
## ========================================================================================

estParamFlowVS = function(fs, channels)
{
  ## some error checking
  checkClass(fs, "flowSet")
  checkClass(channels, "character")
  nmatch = which(channels %in% colnames(fs))
  if(length(nmatch) != length(channels))
    stop(" At least one channel name is not present in the flowSet.")
  
  cofactors = NULL
  ## estimate optimum cofactors
  for(col in channels)
  {
    cat("====================================================================\n")
    cat("Channel ", col, ' : ')
    cat("Finding optimum cofactor for asinh transformation\n")
    cat("====================================================================\n")
    fs1D = fs[,col]
    cf = optimStat(fs1D)
    cofactors = c(cofactors, cf)
    cat("\n Optimum cofactor for ", col, " : ", cf, "\n")
    cat("====================================================================\n\n")
  }
  return (cofactors)
}


## ========================================================================================
## Transform a flowSet by asinh transformation.
## Input:
##    fs: A flowSet containing a collection of flow cyometry samples. 
##    channels: A character vector identifying the channels/dimensions to be transformed.
##    cofactor: cofactors used in the asinh transformation 
## Output: returns a new flowSet with the transformed channels
## ========================================================================================

transFlowVS = function(fs, channels, cofactors)
{
  ## error checking
  checkClass(fs, "flowSet")
  checkClass(channels, "character")
  nmatch = which(channels %in% colnames(fs))
  if(length(nmatch) != length(channels))
    stop(" At least one channel name is not present in the flowSet.")
  
  
  cat("Transforming flowSet by asinh function with supplied cofactors.\n")
  fs1 = fs[1:length(fs)] # forcing deep copy, memory can be a problem
  for(i in 1:length(fs))
  {
    mat = as.matrix(exprs(fs[[i]]))
    for(j in 1:length(channels))
    {
      mat[,channels[j]] = asinh(mat[,channels[j]]/cofactors[j])
    }

    fs1[[i]] = flowFrame(mat)
  }
  
  return(fs1)
}


    

## ========================================================================================
## Internal function
## Compute Bartlett's statistics for a set of 1D clusters (density peaks)
##  of a single channel across all samples of a flowset
## Input:
##    cofactor: the cofactor used with asinh transformation for this channel
##    fs1D: A 1-D flowSet containing a single channel of a collection of FC samples. 
##    MAX_BT: The maximum value for bartlett's stat when it can not be computed 
## Output: A numeric vector representing the optimum cofactors for 
## ========================================================================================

flowVS1D=function(cofactor, fs1D, bwFac = 2, plot=FALSE, MAX_BT=10^9, save.plot=FALSE, populationQuant=.01, borderQuant=0.001, 
                  annotation="",save.folder="peaks", ...)
{
    peakStats = matrix(nrow=0, ncol=5)
    colnames(peakStats) = c("n", "mean", "median", "variance", "sample")
    stain = colnames(fs1D)[1]
    for(i in 1:length(fs1D))  # this loop can be parallelized
    {
      y = exprs(fs1D[[i]])[,1]
      y = asinh(y/cofactor)
      peakStat = densityPeaks(y, stain, bwFac = bwFac, plot=plot, save.plot=save.plot, populationQuant=populationQuant, borderQuant=borderQuant, 
                              annotation=paste(i),save.folder="", ...)
      peakStat = cbind(peakStat, rep(i, nrow(peakStat))) 
      peakStats = rbind(peakStats, peakStat)        
    }
    # remove outlying samples
    # if most of the samples contains k peaks except a few
    # then it is possible that peaks are not identified properly for the few samples.
    # so, we can remove those samples from consideration 
    # we remove 10% outlier samples from consideration

    sampleIdx = peakStats[,5]
    nsamples = length(unique(sampleIdx))
    npeaks = table(sampleIdx) # number of peaks per sample
    
    sampleOrder = order(abs(npeaks - median(npeaks)), decreasing=TRUE)
    rm = round(nsamples/10) 
    if(rm > 0 && nsamples>(rm+1))
    {
      rmIdx = sampleOrder[1:rm]
      rmSamples = as.integer(names(npeaks[rmIdx]))
      
      selcted = !(sampleIdx %in% rmSamples)
      peakStats = peakStats[selcted,]
      
    }

    ## remove outlying peaks  
    ## remove 10% peaks with the highest deviation from median variance
    ## so that 90% peaks can be stabilized without outliers' influence
    pkVar = peakStats[,4]
    npeaks = length(pkVar)
    varOrder = order(abs(pkVar - median(pkVar)), decreasing=TRUE)
    rm = round(npeaks/10) 
    if(rm > 0 && npeaks>(rm+1))
    {
      rmIdx = varOrder[1:rm]
      peakStats = peakStats[-rmIdx,]
    }
    
    if(nrow(peakStats)<=1)
      bt = MAX_BT
    else
      bt = bartlettTest(peakStats)
    return (bt)
    #return (peakStats)
}




peakStats1D=function(fs1D, bwFac = 2, plot=FALSE, save.plot=FALSE, populationQuant=.01, borderQuant=0.001, 
                  annotation="",save.folder="peaks", ...)
{
  peakStats = matrix(nrow=0, ncol=5)
  colnames(peakStats) = c("n", "mean", "median", "variance", "sample")
  stain = colnames(fs1D)[1]
  for(i in 1:length(fs1D))  # this loop can be parallelized
  {
    y = exprs(fs1D[[i]])[,1]
    y = (y-min(y))/(max(y) - min(y))
    peakStat = densityPeaks(y, stain, bwFac = bwFac, plot=plot, save.plot=save.plot, populationQuant=populationQuant, borderQuant=borderQuant, 
                            annotation=paste(i),save.folder="", ...)
    peakStat = cbind(peakStat, rep(i, nrow(peakStat))) 
    peakStats = rbind(peakStats, peakStat)        
  }
  
  peakStats1 = peakStats
  # remove outlying samples
  # if most of the samples contains k peaks except a few
  # then it is possible that peaks are not identified properly for those few samples.
  # so, we can remove those samples from consideration 
  # we remove 10% outlier samples from consideration
  
  sampleIdx = peakStats[,5]
  nsamples = length(unique(sampleIdx))
  npeaks = table(sampleIdx) # number of peaks per sample
  
  sampleOrder = order(abs(npeaks - median(npeaks)), decreasing=TRUE)
  rm = round(nsamples/10) 
  if(rm > 0 && nsamples>(rm+1))
  {
    rmIdx = sampleOrder[1:rm]
    rmSamples = as.integer(names(npeaks[rmIdx]))
    
    selcted = !(sampleIdx %in% rmSamples)
    peakStats = peakStats[selcted,]
    
  }
  
  ## remove outlying peaks  
  ## remove 10% peaks with the highest deviation from median variance
  ## so that 90% peaks can be stabilized without outliers' influence
  pkVar = peakStats[,4]
  npeaks = length(pkVar)
  varOrder = order(abs(pkVar - median(pkVar)), decreasing=TRUE)
  rm = round(npeaks/10) 
  if(rm > 0 && npeaks>(rm+1))
  {
    rmIdx = varOrder[1:rm]
    peakStats = peakStats[-rmIdx,]
  }
  

  bt = bartlettTest(peakStats)
  return (list(peakStats=peakStats1,bartlettStat=bt))
}


##============================================================
## Internal function
## Compute Bartlett's statistics from the statistics of 
## a set of 1D clusters (density peaks)
## Arguments:
##      mean, variance, number of cells in each of k peaks from all samples
## can be passed as a matrix/data frame 
##============================================================
bartlettTest = function(peakStats) 
{
  
  peakStats = peakStats[peakStats[,1]>1,,drop=FALSE] # remove peaks with single point
  size = peakStats[,1]
  
  variance = peakStats[,4]
  sum1 = sum((size-1)*log(variance))
  sum2 = sum((size-1)*variance)
  N = sum(size)
  sum3 = sum(1/(size-1))
  k = nrow(peakStats)

  bt = ( (N-k)*log(sum2/(N-k)) - sum1 ) / (1+( sum3 - 1/(N-k))/(3*(k-1)) )
  return(bt)
}



##============================================================
# internal function
# Identifying the optimum cofactor for a selected channel
# Input: fs1D, a flowSet with a single stain
#        cfLow, cfHigh: lowest and highest possible values for cofactor
#                       (log scale)
# output: optimum cofactor 
##============================================================
optimStat = function(fs1D, cfLow=-1, cfHigh=10, MAX_BT=10^9)
{
  if(cfLow>=cfHigh)
  {
    print("Warning: cfLow>=cfHigh, using default values")
    cfLow=-1
    cfHigh=10
  }
  cf = cfLow:cfHigh
  ncf =  length(cf)
  cfopt = rep(0,ncf-1)
  btopt = rep(0,ncf-1)
  
  cat(sprintf("%18s       %10s    %15s  %8s \n", "cf range", "opt cf", "Bartlett\'s stat", "time"))
  cat("====================================================================\n")
  for(i in 1:(ncf-1))
  {
    ptm <- proc.time()
    tol = (exp(cf[i+1]) - exp(cf[i]))/10
    opt = suppressWarnings(optimize(f = flowVS1D, interval = c(exp(cf[i]),exp(cf[i+1])), fs1D, tol=tol, plot=FALSE, MAX_BT=MAX_BT))
    btopt[i] = opt$objective
    cfopt[i] = opt$minimum
    #cat('cf range= [',format(round(exp(cf[i]), 2), nsmall = 2), ',',
    #    format(round(exp(cf[i+1]) , 2), nsmall = 2),']: ', 
    #    'opt cf= ', format(round(opt$minimum, 2), nsmall = 2), sep="")
   
    cat(sprintf("[%9.2f, %9.2f ] %10.2f ", exp(cf[i]), exp(cf[i+1]), opt$minimum))
    if(opt$objective==MAX_BT)
    {
      cat(sprintf("%10s ", " MAX (10^9)"))
    }
    else
      cat(sprintf("%15.2f ", opt$objective))
    cat(sprintf("%13.2f \n", (proc.time() - ptm)[1]))
  }

  minIdx = which.min(btopt)
  #now perform a local search around cfopt[minIdx] and plot
  del = cfopt[minIdx]/10
  btLocal = rep(0,11)
  btLocal[6] = btopt[minIdx]
  cfLocal = c(5:1,cfopt[minIdx],1:5)
  cfLocal[1:5] = cfopt[minIdx] - 5:1 * del
  cfLocal[7:11] = cfopt[minIdx] + 1:5 * del
  for(i in c(1:5,7:11))
  {
    btLocal[i] = flowVS1D(cfLocal[i], fs1D)
  }
  
  minIdx = which.min(btLocal)
  plot(cfLocal, btLocal, type='o', pch=16, 
       xlab="Cofactors", ylab="Bartlett's statistics",
       main = paste("Optimum cofactor for ", colnames(fs1D), " : ", 
                    format(round(cfLocal[minIdx], 2), nsmall = 2), sep=""))
  points(cfLocal[minIdx],btLocal[minIdx], pch=16, col='red')
  return (cfLocal[minIdx])
}



