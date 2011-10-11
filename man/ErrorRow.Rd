\name{ErrorRow}
\alias{ErrorRow}

\title{training/test errors for a row}
\description{
Finding training and test test error rate where the test data is
r-th row. Number of variables selected beforehand.
}

\usage{
ErrorRow(r, Data, cl, FUN, ...)
}

\arguments{
  \item{r}{The row numbers which correspond to test data}
  \item{Data}{Combination of training and test data}
  \item{cl}{Combination of training and test class}
  \item{FUN}{Methods to be implemented}
  \item{\dots}{Options passed to FUN}
}

\details{Used with delete-d cross-validation}

\value{If successful, the output from FUN is returned. If not successful, -1 
is returned.}

\author{M. S. Islam and A. I. McLeod, email: mislam22@uwo.ca}

\seealso{
\code{\link{DeleteD}},
\code{\link{ErrorDeleteD}}
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

