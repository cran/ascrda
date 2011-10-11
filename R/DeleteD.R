DeleteD <-
function(y, d=0.25){
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
    unlist(newInd)
}

