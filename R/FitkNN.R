`FitkNN` <-
function(X, Y, Xt, Yt, k=1){
##requires library(class)
##Note: knn() requires factor for response variable
##      but we use classes 1, 2, ... for compatibilty with DeleteD
#
#input validation
ok1 <- is.matrix(X)&&is.matrix(Xt)
ok2 <- nrow(X)==length(Y)&&nrow(Xt)==length(Yt)
uy <- sort(unique(Y)) #sorted unique y's
uyt <- sort(unique(Yt)) #sorted unique y's
ok3 <- all(uy==seq(length(uy)))
ok4 <- all(uyt==seq(length(uyt)))
ok <- ok1&&ok2&&ok3&&ok4
if (!ok) stop("FitkNN - Error code:", c("ok1","ok2","ok3","ok4")[!c(ok1,ok2,ok3,ok4)])
#
y <- as.factor(Y)
yt <- as.factor(Yt)
#training
    yp<-knn(X, X, y, k=k)
    cmat<-table(y,yp,dnn=c("Observed","Predicted"))
    r<-1-sum(diag(cmat))/sum(cmat)
#test data
    ytp<-knn(X, Xt, y, k=k)
    cmatTest<-table(yt,ytp,dnn=c("Observed","Predicted"))
    rTest<-1-sum(diag(cmatTest))/sum(cmatTest)
    r<-c(r,rTest)
    names(r)<-c("Train","Test")
    list(r=r, CMats=list(cmat,cmatTest))
}

