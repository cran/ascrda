FitSvm <-
function(X,Y,Xt,Yt, tuningQ=FALSE, cross=10, ...){
y <- Y
yt <- Yt
#  svm, lssvm require factors. But for compatibility with DeleteD cross-validation
#   we use consecutive integers and then convert to factors
ok1  <-nrow(X)==length(y)&&nrow(Xt)==length(yt)&&ncol(Xt)==ncol(X)
uy   <- sort(unique(y)) #sorted unique y's
uyt  <- sort(unique(yt)) #sorted unique y's
ok3  <- all(uy==seq(length(uy)))
ok4  <- all(uyt==seq(length(uyt)))
ok   <- ok1&&ok3&&ok4
if (!ok) stop("FitSvm - Error code:", c("ok1","ok3","ok4")[!c(ok1,ok3,ok4)])
#
zY <- as.factor(Y)
zYt <- as.factor(Yt)
# training
   xy.df<-data.frame(X=X, Y=zY)
#
   if (tuningQ) {
   	tobj <- tune.svm(Y~., data=xy.df, gamma=10^(-6:-3), cost=10^(1:2))
   	bestG <- tobj$best.parameters[[1]]
   	bestC <- tobj$best.parameters[[2]]
   	model <- svm(Y~., data=xy.df, cost=bestC, gamma=bestG, cross=cross, ...)
   } else model <- svm(Y~., data=xy.df)
#training data, confusion matrix
   yp <- predict(model)
   cmat<-table(zY,yp,dnn=c("Observed","Predicted"))
   r<-1-sum(diag(cmat))/sum(cmat)
# test data
    xt.df<-data.frame(X=Xt)
    ytp<-predict(model, xt.df)
#confusion matrix
    cmatTest<-table(zYt,ytp,dnn=c("Observed","Predicted"))
    rTest<-1-sum(diag(cmatTest))/sum(cmatTest)
    r<-c(r,rTest)
    names(r)<-c("Train","Test")
#output
    out <- list(r=r,  model=model, cmatTrain=cmat, cmatTest=cmatTest)
    out
}
