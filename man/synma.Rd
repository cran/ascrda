\name{synma}
\alias{synma}
\title{
Generate training and test synthetic microarray datasets with Q=2 classes.
}
\description{
Training and test expression matrices of size G-by-n where G is the number
of genes and n is the number of samples are generated as well as
vectors of length n for the two output classes.
}
\usage{
synma(n = c(100, 100), nt = c(500, 500), rho = c(0.9, 0.9), B = c(20, 20), m = c(100, 100), fE = 0.02, E = 0.5, Si = 1)
}

\arguments{
  \item{n}{
vector length 2, c(number healthy, number diseases), training sample
}
  \item{nt}{
vector length 2, c(number healthy, number diseases), test sample
}
  \item{rho}{
vector length 2, c(correlation healthy, correlation diseases)
}
  \item{B}{
vector length 2, c(block size healthy, block size disease)
}
  \item{m}{
vector length 2, c(# blocks healthy, # blocks disease)
}
  \item{fE}{
fraction of genes expressed in 'diseased' patients
}
  \item{E}{
Expression level of diseased genes
}
  \item{Si}{
Variance inflation factor applied to the 'diseased' group
}

}
\details{
The number of genes is determined as G=B[1]*m[1] = B[2]*m[2].
If these are not equal an informative error message is given.
}

\value{
A list with four elements:

\item{X}{expression matrix, training, (B*m)-by-(n[1]+n[2]), where B=B[1] and m=m[1]}
\item{y}{response, training, n[1] healthy followed by n[2] diseased}
\item{Xt}{expression matrix, test, (B*m)-by-(nt[1]+nt[2]), where B=B[2] and m=m[2], 
rows are augmented if B[2]*m[2] not equal to G}
\item{yt}{response, test, nt[1] healthy followed by nt[2] diseased}
}

\references{ 
Guo Y, Hastie T, Tibshirani R (2007). 
Regularized linear discriminant analysis and its application in microarrays. 
Biostatistics 8, 86-100. 

Islam, M. S. (2008). 
Periodicity, Change Detection and Prediction in Microarrays. Ph.D. thesis, 
The University of Western Ontario, Canada.
}

\author{M. S. Islam and A. I. McLeod, email: aimcleod@uwo.ca}

\examples{
#A micro-sized microarray!
ans<-synma(n=c(3,3),nt=c(2,2),B=c(3,3),m=c(3,3),fE=0.5)
ans
}

\keyword{classif}

