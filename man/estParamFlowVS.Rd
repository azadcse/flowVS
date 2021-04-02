\name{estParamFlowVS}
\alias{estParamFlowVS}
\title{
  Estimate optimum parameters for per-channel within-population variance stabilization.
}
\description{
This function estimates the variance stabilizing cofactors, one for each channel for the entire dataset. When a fluorescence channel \code{z} is transformed by asinh transformation with the optimum cofactor for \code{z}, the within-population variances of populations from all samples in the channel \code{z}  are approximately stabilized. 
}
\usage{
  estParamFlowVS(fs, channels)
}

\arguments{
    \item{fs }{A flowSet containing a collection of flow cyometry samples.}
    \item{channels }{A character vector identifying the channels/dimensions to be transformed. If any entry in this vector is not present in the flowSet, the function returns with an error.}
}
\details{
Let \code{z} be a fluorescence channel (column of a flowFrame). 
We consider transforming \code{z} by asinh transformation such that after transformation we obtain the tranformed channel asinh(\code{z/c}), where \code{c} is a normalizing cofactor.

The \code{estParamFlowVS} function estimates cofactors, one for each channel for the entire dataset such that the within-population variance is stabilized in each fluorescence channel. When a fluorescence channel \code{z} is transformed by asinh transformation with the optimum cofactor for \code{z}, the within-population variances of populations from all samples in the channel \code{z}  are approximately stabilized. 
}
\value{
  \code{estParamFlowVS} returns a numeric vector representing the optimum cofactors for the requested channels. The optimum cofactor for the input \code{channels[i]} is stored in the \code{i}th entry of the returned vector.
}
\references{
  Ariful Azad, Bartek Rajwa, and Alex Pothen (2015), "flowVS: Channel-Specific Variance Stabilization in Flow Cytometry", BMC Bioinformatics, vol 17, pp 1-14, 2016.
}
\author{Ariful Azad}


\seealso{
  \code{\link{transFlowVS}}
}
\examples{
data(HD)
## identify optimum cofactor for CD3 and CD4 channels (from five samples)
\donttest{
cofactors = estParamFlowVS(HD[1:5],channels=c('CD3','CD4'))}
# See detail examples in the documentation of the transFlowVS function.
  
}


\keyword{ transformation }
\keyword{ variance stabilization }