\name{transFlowVS}
\alias{transFlowVS}
\title{
  Transform a flowSet by asinh transformation.
}
\description{
This function transforms a flowSet by asinh transformation with the cofactors passed on to the function.
The optimum cofactors that stabilize within-population variances in different fluorescence channels are estimated beforehand and passed to this function for data transformation. 
}
\usage{
  transFlowVS(fs, channels, cofactors)
}

\arguments{
    \item{fs }{A flowSet containing a collection of flow cyometry samples.}
    \item{channels }{A character vector identifying the channels/dimensions to be transformed. If any entry in this vector is not present in the flowSet, the function returns with an error.}
    \item{cofactors }{A numeric vector. \code{cofactors[i]} is used with asinh function to transform the column with name specified by \code{channels[i]}.}
}
\details{
This function transforms a flowSet by asinh transformation with selected cofactors.
The column with name \code{channels[i]} of every flowFrame of the input flowSet is transformed by asinh transformation with \code{cofactors[i]}. 
For example, let \code{z_ij} be the \code{i}th column of \code{j}th flowFrame in the input flowSet \code{fs}.
Then after transformation \code{z_ij} would be converted to asinh(\code{z_ij/cofactors[i]}).

For variance stabiilzation, the optimum cofactors that stabilize within-population variances in different fluorescence channels are estimated beforehand and passed to this function for data transformation. 
Variance stabilizing cofactors can be estimatd by the \code{\link{estParamFlowVS}} function.
}
\value{
  \code{transFlowVS} returns a new flowSet with the transformed channels.
 % The \code{i}th column of every flowFrame of the input flowSet \code{fs} is transformed by asinh transformation with \code{cofactors[i]}. 
}
\references{
  Ariful Azad, Bartek Rajwa, and Alex Pothen (2015), "flowVS: Channel-Specific Variance Stabilization in Flow Cytometry", manuscript submitted for publication.
}
\author{Ariful Azad}


\seealso{
  \code{\link{estParamFlowVS}}
}
\examples{
## ------------------------------------------------
## Example 1: Healthy data from flowVS package
## ------------------------------------------------
data(HD)
\donttest{
## identify optimum cofactor for CD3 and CD4 channels (from five samples)
cofactors = estParamFlowVS(HD[1:5],channels=c('CD3','CD4'))
## transform CD3 and CD4 channels in all samples
HD.VS = transFlowVS(HD, c('CD3','CD4'), cofactors)
densityplot(~CD3+CD4, HD.VS, main="Transfromed CD3 and CD4 channels in HD data") }

## ------------------------------------------------
## Example 2: ITN data from flowStats package
## ------------------------------------------------
\donttest{
library(flowStats)
data(ITN)
# identify lymphocytes
ITN.lymphs = fsApply(ITN,lymphs, list("FS"=c(200, 600),"SS"=c(0, 400)), "FSC", "SSC",FALSE)
## identify optimum cofactor for CD3 and CD4 channels (from five samples)
cofactors = estParamFlowVS(ITN.lymphs[1:5],channels=c('CD3', 'CD4'))  
## transform CD3 and CD4 channels in all samples
ITN.VS = transFlowVS(ITN.lymphs, c('CD3','CD4'), cofactors)
densityplot(~CD3+CD4, ITN.VS, main="Transfromed CD3 and CD4 channels in ITN data")}  
}

\keyword{ transformation }
\keyword{ variance stabilization }