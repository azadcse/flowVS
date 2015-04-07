\name{microVS}
\alias{microVS}
\title{
  Variance stabilization for microarray data.
}
\description{
Variance-stabilizing inverse hyperbolic sine (\code{asinh}) transformation for microarray data.
}
\usage{
  microVS(data, cfLow=0, cfHigh=10, frac=1)
}

\arguments{
    \item{data }{The microarray data in a Matrix.}
    \item{cfLow }{lowest possible value of cofactor (log scale).}
    \item{cfHigh }{highest possible value of cofactor (log scale).}
     \item{frac }{fraction of differentially expressed genes used in variance stabilization (0< frac <=1).}
}
\details{
This function transforms a microarray data matrix \code{z} by \code{asinh(z/c)} transformation where \code{c} is a normalizing cofactor.
The cofactor is searched in the range \code{[cfLow, cfHigh]} and an optimum cofactor is obtained for which the transformed data is variance stabilized.
The optimum cofactor is obtained by minimizing Bartlett's test statistics for homogeneity of variance.
If the parameter \code{frac} is less then one, a fraction of differentially expressed genes are used in estimating the cofactor.
}
\value{
  \code{microVS} returns a matrix of the variance-stabilizing microarray data.
}
\references{
  Ariful Azad, Bartek Rajwa, and Alex Pothen (2015), "flowVS: Channel-Specific Variance Stabilization in Flow Cytometry", manuscript submitted for publication.
}
\author{Ariful Azad}


\examples{

# stabilize variance of the Kidney microarray data from the vsn package
library(vsn)
data(kidney)
kidney.t = microVS(exprs(kidney))
plotMeanSd(kidney.t)

}

\keyword{ transformation }
\keyword{ variance stabilization }