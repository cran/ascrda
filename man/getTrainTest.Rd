\name{getTrainTest}
\alias{getTrainTest}
\title{
Form training and test samples using delete-d method stratified sampling.
}
\description{
Supports delete-d cross validation.
For specified d, the training and test samples are constructed using
a stratified sample to preserve the proportions in each class.
}
\usage{
getTrainTest(X, y, d = 0.25)
}
\arguments{
  \item{X}{
inputs in matrix form, n-by-p matrix
}
  \item{y}{
output, p-vector
}
  \item{d}{
proportion for hold-out sample,  0<d<1.
}
}

\details{
Kim (2009) compares various CV methods and recommends the delete-d also known as the hold-out method
for better accuracy than bootstrap or k-fold methods.
}

\value{
\item{X}{training inputs}
\item{y}{training output}
\item{Xt}{test inputs}
\item{yt}{test output}
}

\references{
Kim, J. H. (2009).
Estimating Classification Error Rate: Repeated Cross-validation, Repeated
Hold-out and Bootstrap.
Computational Statistics and Data Analysis, 53, 3735--3745.
}



\author{M. S. Islam and A. I. McLeod, email: aimcleod@uwo.ca}

\examples{
X <- matrix(rnorm(90), nrow=30)
y <- rep(1:3, 10)
M <- getTrainTest(X, y, d=0.25)
}
\keyword{ classif }
