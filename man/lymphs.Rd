\name{lymphs}
\alias{lymphs}
\title{
  Identify lymphocyte cells from a flow cytometry sample.
}
\description{
Identify and retain lymphocytes from a flow cytometry sample based on the forward and side scatters. 
}
\usage{
  lymphs(ff, lymph.boundary, fsc, ssc, plot=FALSE)
}

\arguments{
    \item{ff }{A flowFrame containing a flow cyometry sample.}
    \item{lymph.boundary }{A list denoting an approximate rectangular boundary for lymphocytes. The first element of the list represents the lower and upper limit of forward scatter (FSC), and the 2nd element represents the lower and upper limit of side scatter (SSC). Example: list("FSC"=c(180000, 500000),"SSC"=c(0, 180000)).}
    \item{fsc}{name (or numeric index) of the forward scatter channel.}
    \item{ssc}{name (or numeric index) of the side scatter channel.}
    \item{plot}{true/false. If true then plots the rectangular and elliptical gates for the lymphocytes.}
}
\details{
At first a rectangular gate is created based on the \code{lymph.boundary}.
Then the \code{norm2Filter} function is used to identify lymphocytes.
}
\value{
  \code{lymphs} returns a new flowFrame containing the lymphocytes.
}
\references{
  Ariful Azad, Bartek Rajwa, and Alex Pothen (2015), "flowVS: Channel-Specific Variance Stabilization in Flow Cytometry", manuscript submitted for publication.
}
\author{Ariful Azad}


\seealso{
  \code{\link{estParamFlowVS}}
}
\examples{
library(flowStats)
data(ITN)
# identify lymphocytes
ITN.lymphs = lymphs(ITN[[1]], list("FS"=c(200, 600),"SS"=c(0, 400)), "FSC", "SSC",TRUE)
}

\keyword{ transformation }
\keyword{ variance stabilization }