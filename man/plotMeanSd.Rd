\name{plotMeanSd}
\docType{methods}
\alias{plotMeanSd}
\title{Plot row standard deviations versus row means (modified from vsn package)}

\description{Plot row standard deviations versus row means of a data matrix.}
\usage{
plotMeanSd(x, 
           ranks = TRUE,
           xlab  = ifelse(ranks, "Rank of means (ascending order)", "mean"),
           ylab  = "Standard deviation",
           pch   = ".",
           plot  = TRUE,
           \dots)} 
\arguments{
  \item{x}{An object of class
    \code{\link{matrix}}}
\item{ranks}{Logical, indicating whether the x-axis (means) should be plotted
  on the original scale (\code{FALSE}) or on the rank scale (\code{TRUE}). The latter
  distributes the data more evenly along the x-axis and allows a
  better visual assessment of the standard deviation as a function of 
  the mean.}
\item{xlab}{Character, label for the x-axis.}
\item{ylab}{Character, label for the y-axis.}
\item{pch}{Plot symbol.}
\item{plot}{Logical.  If \code{TRUE} (default), a plot is produced.
  Calling the function with \code{plot=FALSE} can be useful if only
  its return value is of interest.}
\item{\dots}{Further arguments that get passed to plot.default.}
}
\details{Standard deviation and mean are calculated row-wise from the
  expression matrix (in) \code{x}. The scatterplot of these versus each other
  allows to visually verify whether there is a dependence of the standard
  deviation (or variance) on the mean.
  The red dots depict the running median estimator (window-width 10\%).
  If there is no variance-mean dependence, then the line formed by the
  red dots should be approximately horizontal.
}
\value{
  A named list with four components: its elements \code{px} and
  \code{py} are the x- and y-coordinates of the individual data points
  in the plot; its first and second element are the x-coordinates and values of
  the running median estimator (the red dots in the plot).
  Depending on the value of \code{plot}, the method can also have a side effect, which is to create a plot on the
  active graphics device.
}

\examples{
  library(vsn)
  data(kidney)
  kidney.t = microVS(exprs(kidney))
  plotMeanSd(kidney.t)
}
\keyword{hplot}
\keyword{methods}
