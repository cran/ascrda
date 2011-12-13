\name{ascrda}
\alias{ascrda}

\title{Test error rate for ASCRDA}

\description{The misclassifcation rate on the test data is computed for the following methods:
ASCRDAk - SCRDA with nearest neightbor autocovariates, k=kHat;
ASCRDA1 - SCRDA with nearest neightbor autocovariate, k=1;
SCRDA - regular SCRDA. shrunken centroid regularized discriminant analysis, k=0;
}

\usage{
ascrda(X, Y, Xt, Yt, alpha = seq(0, 0.99, 0.11), delta = seq(0, 3, 0.2), 
	SCRDAmethod=c("ASCRDAk","ASCRDA1","SCRDA"), ...)
}

\arguments{
  \item{X}{training, data matrix, n-by-G}
  \item{Y}{training output}
  \item{Xt}{test data matrix}
  \item{Yt}{test output}
  \item{alpha}{ alpha in rda function}
  \item{delta}{ delta in rda function}
  \item{SCRDAmethod}{ one of ASCRDAk, ASCRDA1, SCRDA. Default is ASCRDAk}
  \item{\dots}{ optional arguments passed to rda}
}


\value{Vector of length 2, misclassification error rate and K, number of knn neighbours selected}

\references{
Guo Y, Hastie T, Tibshirani R (2007). 
Regularized linear discriminant analysis and its application in microarrays. 
Biostatistics 8, 86-100. 

Islam, M. S. (2008). Periodicity, Change Detection and Prediction in Microarrays. Ph.D. thesis, 
The University of Western Ontario, Canada.
}

\author{M. S. Islam and A. I. McLeod, email: aimcleod@uwo.ca}

\examples{
#Example 1. Tiny test example - 1.4 seconds
set.seed(132561)#for exact reproducibility
M<-synma(n=c(10,10), nt=c(10,10), B=c(5,5), m=c(5,5))
startTime <- proc.time()[3]
ans<-ascrda(M$X, M$y, M$Xt, M$yt)
endTime <- proc.time()[3]
TotalTime <- endTime-startTime
out <- c(ans[[1]], TotalTime)
names(out)<-c("MisclassificationRate", "TotalTime")
out
#> out
#MisclassificationRate             TotalTime
#                 0.45                  1.51

\dontrun{
#Example 2.  large example. About 78 sec
#training sample size = 200
#test sample size = 1000
#number of genes = 1000
set.seed(132561)#for exact reproducibility
BScenario <- list(c(50, 10),c(10, 50), c(10, 10),  c(50, 50))
mScenario <- list(c(20,100),c(100, 20),c(100, 100),c(20, 20))
M <- synma(n=c(100,100), nt=c(500,500), rho = c(0.9, 0.9),
  B=BScenario[[1]], m=mScenario[[1]], fE=0.05)
StartTime <- proc.time()[3]
ANS1 <- ascrda(M$X, M$y, M$Xt, M$yt)
EndTime1 <- proc.time()[3]
ANS2 <- ascrda(M$X, M$y, M$Xt, M$yt, SCRDAmethod = "ASCRDA1")
EndTime2 <- proc.time()[3]
ANS3 <- ascrda(M$X, M$y, M$Xt, M$yt, SCRDAmethod = "SCRDA")
EndTime3 <- proc.time()[3]
TotalTimes <-c(EndTime1-StartTime,EndTime2-EndTime1,EndTime3-EndTime2)
Errs <- c(ANS1[[1]],ANS2[[1]],ANS3[[1]])
out <- matrix(c(Errs, TotalTimes),ncol=2)
rownames(out) <- c("ASCRDAk", "ASCRDA1", "SCRDA")
colnames(out) <- c("ErrorRate", "Time")
out
#> out
#        ErrorRate  Time
#ASCRDAk     0.076 43.39
#ASCRDA1     0.185 22.02
#SCRDA       0.300 12.22
}

#Example 3. With 3 classes.  About 2.5 sec
#Spira Dataset
y <- as.numeric(Spira$Patients$STATUS)
X <- Spira$X
genes <- Spira$genes
X <- X[,genes]
set.seed(132561)#for exact reproducibility
M <-  getTrainTest(X, y, d=0.25)
StartTime <- proc.time()[3]
out <- ascrda(M$X, M$y, M$Xt, M$yt)
EndTime <- proc.time()[3]
TotalTime <-EndTime-StartTime
out
TotalTime
}

\keyword{ classif }

