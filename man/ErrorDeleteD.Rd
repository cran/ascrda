\name{ErrorDeleteD}
\alias{ErrorDeleteD}

\title{error rate using delete-d method}

\description{
Uses the functions ErrorRow and DeleteD to compute the
error rate in one hold-out sample.
Usually about 1000 iterations of this function are
used for a complete hold-out cross-validation estimate.
}

\usage{ErrorDeleteD(Data, cl, FUN = FitRda, d = 0.33, ...)}

\arguments{
  \item{Data}{Combination of training and test data}
  \item{cl}{Combination of training and test class}
  \item{FUN}{Methods to be implemented, default is rda}
  \item{d}{Proportion to hold out for delete-d method}
  \item{\dots}{Options passed to FUN}
}

\value{Output from FUN}


\author{M. S. Islam and A. I. McLeod, email: mislam22@uwo.ca}

\seealso{
\code{\link{DeleteD}},
\code{\link{ErrorRow}}
}

\examples{
require("rda")
data(brain)
#time: about 51 seconds
startTime <- proc.time()[3]
ErrorDeleteD(brain.x, brain.y, FUN=ascrda, d=0.33)
endTime <- proc.time()[3]
TotalTime <- endTime - startTime
TotalTime
}


\keyword{ classif }
