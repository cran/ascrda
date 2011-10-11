\name{DeleteD}
\alias{DeleteD}
\title{Balanced hold-back sample indicies}
\description{
Finds the indicies for the holdout sample using the delete-d method.
At least one observation is kept in both the training/test samples.
}
\usage{
DeleteD(y, d = 0.25)
}

\arguments{
  \item{y}{output, classes must be 1, 2, ...}
  \item{d}{fraction of data in hold-out or test sample}
}

\details{Used with delete-d cross-validation}

\value{vector of indicies for hold-out sample}

\author{M. S. Islam and A. I. McLeod, email: mislam22@uwo.ca}

\seealso{
\code{\link{ErrorRow}},
\code{\link{ErrorDeleteD}}
}

\examples{
DeleteD(rep(1:5, 10), d=0.5)
}

\keyword{ classif }
