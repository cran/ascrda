`FitDLDA` <-
function(X, Y, Xt, Yt, pool=TRUE){
#for consistency
y <- Y
yt <- Yt
#requires sfsmisc
#y and yt class labels: must be consecutive integers 1, 2,...
#Assume X is design matrix form
#input validation
ok1 <- is.matrix(X)&&is.matrix(Xt)&&is.vector(y)&&is.vector(yt)&&is.logical(pool)
ok2 <- nrow(X)==length(y)&&nrow(Xt)==length(yt)&&ncol(Xt)==ncol(X)
uy <- sort(unique(y)) #sorted unique y's
uyt <- sort(unique(yt)) #sorted unique y's
ok3 <- all(uy==seq(length(uy)))
ok4 <- all(uyt==seq(length(uyt)))
ok <- ok1&&ok2&&ok3&&ok4
if (!ok) stop("FitDLDA - Error code:", c("ok1","ok2","ok3","ok4")[!c(ok1,ok2,ok3,ok4)])
#
ytp <- diagDA(X, y, Xt, pool=pool)
cmatTest <- table(yt, ytp, dnn = c("Observed", "Predicted"))
Err <- 1 - sum(diag(cmatTest))/sum(cmatTest)
list(Err=Err, confusionMatrix=cmatTest)
}

