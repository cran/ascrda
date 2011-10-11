FitSVM <-
function(X,y,Xt,yt){
 stopifnot(is.vector(y), is.vector(yt),
                nrow(X)==length(y),
                nrow(Xt)==length(yt),
                ncol(Xt)==ncol(X))
# training
   y<-as.factor(y)
   yt<-as.factor(yt)
    xy.df<-data.frame(X,Y=y)
    ans<-svm(Y~., data=xy.df)
    yp <- predict(ans)
    cmat<-table(y,yp,dnn=c("Observed","Predicted"))
    r<-1-sum(diag(cmat))/sum(cmat)
# test data
    xt.df<-data.frame(Xt)
    ytp<-predict(ans, xt.df)
    cmatTest<-table(yt,ytp,dnn=c("Observed","Predicted"))
    rTest<-1-sum(diag(cmatTest))/sum(cmatTest)
    r<-c(r,rTest)
    names(r)<-c("Train","Test")
    r
}

