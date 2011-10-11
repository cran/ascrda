khat <-
function(X, Y, kmax=ceiling(length(Y)*0.5), plot=TRUE){
d<-numeric(kmax)
Q <- length(unique(Y))
if (Q<2)
    stop("error - fewer than 2 levels in Y")
y <- factor(Y)
for (k in 1:kmax){
    z <- nnc(X = X, Y = y, k=k)
    if (Q == 2)
        d[k]<- deviance(glm.fit(x = z, y = y, family = binomial(link = "logit")))
    else {
        sink("junk.txt")
        d[k] <- multinom(y~z)$deviance
        sink()
    }
}
unlink("junk.txt")        
khat <- which.min(d)
#plot and/or return
if (plot) {
    plot(d, xlab="k", ylab="deviance")
    points(khat, d[khat], col="red", pch=16, cex=1.4)
    title(sub=bquote(hat(k) == .(khat)))
    return(invisible(khat))
    }
else return(khat)
}
