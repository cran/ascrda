getTrainTest <-
function(X, y, d=0.25){
  uy <- sort(unique(y))
  ok <- all(uy==seq(1,length(uy)))
  if (!ok) stop("getTrainTest: error - classes must be consecutive integers: 1, 2,...")
  if (d<=0 || d>=1) stop("getTrainTest: error - d must be in (0,1)")
#
  ny <- length(y)
  tabY <- table(y)
  fn<- (1-d)*tabY # looking for the training set
  fn1 <- ifelse(fn<1, 1, fn)  # we want at least one observation from each class
  sampN<-round(fn1)
  IndY <- split(1:ny, y)
  newInd<-list()
  for ( i in 1:length(tabY)){
    if (length(IndY[[i]])>1)
       newInd[[i]]<-sample(IndY[[i]], sampN[i]) else
          newInd[[i]]<-IndY[[i]]
    }
  itrain <- unlist(newInd)
  itest <- setdiff(1:length(y), itrain)
  list(X=X[itrain,], y=y[itrain], Xt=X[itest,], yt=y[itest])
}
