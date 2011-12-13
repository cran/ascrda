\name{FitSvm}
\alias{FitSvm}

\title{Classification using Support Vector Machine }
\description{
Fits Support Vector Machine, computes mis-classification rates on training and test data.
}
\usage{
FitSvm(X,Y,Xt,Yt, tuningQ=FALSE, cross=10, ...)
}

\arguments{
  \item{X}{training, data matrix, input }
  \item{Y}{training output, factor}
  \item{Xt}{test, data matrix, input }
  \item{Yt}{test output, factor}
  \item{tuningQ}{use cross-validation to choose tuning parameter}
  \item{cross}{k for k-fold cv}
  \item{\ldots}{optional arguments passed to svm()}
  }

\value{
Error rates for training and test data
}


\references{ 
R package: e1071.
sfsmisc: Utilities from Seminar fuer Statistik ETH Zurich.
http://CRAN.R-project.org/package=e1071.
Evgenia Dimitriadou, Kurt Hornik, Friedrich Leisch, David Meyer and Andreas Weingessel. 2010.

Support Vector Machines in R.
Journal of Statistical Software, V.15, Number 9.
Karatzoglou, A., Meyer, D. and Hornik, K. 2006.

Hastie T, Tibshirani R, Frienman J: The Elements of Statistical Learning: Data Mining, Inference, and
Prediction. Springer-Verlag 2009.
}

\author{M. S. Islam and A. I. McLeod, email: aimcleod@uwo.ca}

\examples{
#Example 1.
#Simple fast example - about 0.25 sec.
#Small example. p=400 genes; training sample size = 100; test sample size = 200.
set.seed(132561)#for exact reproducibility
M<-synma(n=c(50,50), nt=c(100,100), B=c(20,4), m=c(20,100), fE=0.05)
FitSvm(M$X, M$y, M$Xt, M$yt)

#Examle 2.
#Small example. Compare tuning and no tuning. 
#p=400 genes; training sample size = 100; test sample size = 200.
\dontrun{
set.seed(132561)#for exact reproducibility
M<-synma(n=c(50,50), nt=c(100,100), B=c(20,4), m=c(20,100), fE=0.05)
StartTime <- proc.time()[3]
ANS1 <- FitSvm(M$X, M$y, M$Xt, M$yt)
EndTime1 <- proc.time()[3]
ANS2 <- FitSvm(M$X, M$y, M$Xt, M$yt, tuningQ=TRUE)
EndTime2 <- proc.time()[3]
TotalTime1 <-EndTime1-StartTime
TotalTime2 <-EndTime2-EndTime1
TotalTime <- c(TotalTime1, TotalTime2)
names(TotalTime) <- c("default","tuning")
Err <- c(ANS1$r[2],ANS2$r[2])
names(Err) <- c("default","tuning")
m<-matrix(c(Err, TotalTime), ncol=2)
dimnames(m) <- list(names(Err), c("MisclassificationRate", "TotalTime"))
m
#
#OUTPUT
##> m
##        MisclassificationRate TotalTime
##default                  0.26      0.25
##tuning                   0.39     35.48
#
#Note: tuning did not improve accuracy and timing increased dramatically!
}

#Example 3. Large Example
\dontrun{
#training sample size = 200
#test sample size = 1000
#number of genes = 1000
set.seed(132561)#for exact reproducibility
BScenario <- list(c(50, 10),c(10, 50), c(10, 10),  c(50, 50))
mScenario <- list(c(20,100),c(100, 20),c(100, 100),c(20, 20))
M <- synma(n=c(100,100), nt=c(500,500), rho = c(0.9, 0.9),
  B=BScenario[[1]], m=mScenario[[1]], fE=0.05)
StartTime <- proc.time()[3]
ANS <- FitSvm(M$X, M$y, M$Xt, M$yt)
EndTime <- proc.time()[3]
TotalTime <-EndTime-StartTime
ans <- c(ANS$r[2], TotalTime)
names(ans) <- c("MisclassificationRate", "TotalTime")
ans
#> ans
#MisclassificationRate             TotalTime
#                0.219                 1.610
}

#Example 4. With 3 classes
#
#Spira Dataset
#takes about 0.25 seconds.
#So full cross-validation with 1000 replications ~8.6 hours
#Here p=9968; n=75; 3 classes
\dontrun{
#Spira Dataset
y <- as.numeric(Spira$Patients$STATUS)
X <- Spira$X
genes <- Spira$genes
X <- X[,genes]
set.seed(132561)#for exact reproducibility
M <-  getTrainTest(X, y, d=0.25)
StartTime <- proc.time()[3]
ANS <- FitSvm(M$X, M$y, M$Xt, M$yt)
EndTime <- proc.time()[3]
TotalTime <-EndTime-StartTime
out <- c(ANS$r[2], TotalTime)
names(out) <- c("MisclassificationRate", "Time")
out
TotalTime
}
}

\keyword{ classif }
