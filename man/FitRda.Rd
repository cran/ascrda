\name{FitRda}
\alias{FitRda}
\title{ Regularized Discriminant Analysis }
\description{
Fits RDA, computes mis-classification rates on training and test data.
}
\usage{
FitRda(X, y, Xt, yt, alpha = seq(0, 0.99, 0.11), delta = seq(0, 3, 0.2), ...)
}

\arguments{
  \item{X}{training, expression matrix }
  \item{y}{response vector taking values 1, 2, ...}
  \item{Xt}{test, expression matrix }
  \item{yt}{ response vector for test data taking values 1, 2, ... }
  \item{alpha}{ tuning parameter }
  \item{delta}{turning parameter }
  \item{...}{ optional arguments, passed to rda }
}

\value{
Error rates for training and test data if rda.cv() works. Otherwise
no result returned. You need to use: out<-FitRda(...) and exists(out)
to determine if an error occurred. 
}


\references{ 
Guo Y, Hastie T, Tibshirani R (2007). 
Regularized linear discriminant analysis and its application in microarrays. 
Biostatistics 8, 86-100. 
}

\author{M. S. Islam and A. I. McLeod, email: aimcleod@uwo.ca}

\examples{
#Example 1
#Tiny test example - 1.4 seconds
M<-synma(n=c(10,10), nt=c(10,10), B=c(5,5), m=c(5,5))
startTime <- proc.time()[3]
ans<-FitRda(M$X, M$y, M$Xt, M$yt)
endTime <- proc.time()[3]
TotalTime <- endTime-startTime
out <- c(ans[2], TotalTime)
names(out)<-c("MisclassificationRate", "TotalTime")
out
#
\dontrun{
#Example 2.  large example. About 12.6 sec
#training sample size = 200
#test sample size = 1000
#number of genes = 1000
set.seed(132561)#for exact reproducibility
BScenario <- list(c(50, 10),c(10, 50), c(10, 10),  c(50, 50))
mScenario <- list(c(20,100),c(100, 20),c(100, 100),c(20, 20))
M <- synma(n=c(100,100), nt=c(500,500), rho = c(0.9, 0.9),
  B=BScenario[[1]], m=mScenario[[1]], fE=0.05)
StartTime <- proc.time()[3]
ANS <- FitRda(M$X, M$y, M$Xt, M$yt)
EndTime <- proc.time()[3]
TotalTime <-EndTime-StartTime
ans <- c(ANS[2], TotalTime)
names(ans) <- c("MisclassificationRate", "TotalTime")
ans
##> ans
##MisclassificationRate             TotalTime
##                 0.30                 12.61
}

\dontrun{
#Example 3. With 3 classes.  About 76.5 sec
#Spira Dataset
y <- as.numeric(Spira$Patients$STATUS)
X <- Spira$X
set.seed(132561)#for exact reproducibility
M <-  getTrainTest(X, y, d=0.25)
StartTime <- proc.time()[3]
ANS <- FitRda(M$X, M$y, M$Xt, M$yt)
EndTime <- proc.time()[3]
TotalTime <- EndTime-StartTime
out <- c(ANS[2], TotalTime)
names(out) <- c("MisclassificationRate", "TotalTime")
out
#MisclassificationRate             TotalTime 
#            0.5555556            76.5500000 
}

}

\keyword{ classif }
