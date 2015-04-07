\name{HD}
\alias{HD}
\docType{data}
\title{ Sample flow cytometry data from healthy individuals}
\description{
  A \code{\link[flowCore:flowSet-class]{flowSet}} cotaining 12 flow cytometry samples from three healthy individuals "A", "C", and "D".
From each individual, the samples were drawn on two different days and two technical replicates were created from each sample (i.e., 3 x 2 x 2 = 12 samples).
Each HD sample was stained using labeled antibodies against CD45, CD3, CD4, CD8, and CD19 protein markers.
Here, an HD sample "C_4_2" means that  it is collected on day 4 from individual "C" and it is the second replicate on that day.
We have identified lymphocytes in each sample of the HD dataset and apply the subsequent analysis on lymphocytes.
}
\usage{data(HD)}
\value{
  A \code{flowSet} containing 12
  \code{\link[flowCore:flowFrame-class]{flowFrames}}. There are 3
  subject groups with 4 samples each (2 days and 2 technical replicates per day). 
}

\references{
  Ariful Azad, Arif Khan, Bartek Rajwa, Saumyadipta Pyne, and Alex Pothen, 
  "Classifying immunophenotypes with templates from flow cytometry", In Proceedings of the International Conference on Bioinformatics, Computational Biology and Biomedical Informatics (ACM BCB), 2013.

}

\keyword{datasets}
