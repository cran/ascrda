nncTest<-function (X, Y, XTest, k) { 
    n <- length(Y)
    nt<-nrow(XTest)
    if (nrow(X) != n) 
        stop("nrow(X) != length(Y)")
    classes <- unique(Y)
    Q <- length(classes)
    if (Q < 2) 
        stop("invalid, need Q >= 2, Q = number of classes")
    y <- numeric(n)
    if (Q == 2) {
        if (!all(y %in% c(-1, 1))) {
            ind <- Y == classes[1]
            y[ind] <- -1
            y[!ind] <- 1
        }
        ans<- knn(train = X, test = XTest, cl = as.factor(y), k = k, prob = TRUE, use.all = FALSE)
        pr<- attr(ans, "prob") #proportion of votes for winning class
        yfit<-as.numeric(as.character(ans)) #predicted classes as numeric -1 or +1
        NumberOfVotesForWinner <- pr*k
        NumberOfVotesForLoser <- k - NumberOfVotesForWinner
        ConsensusProportion <- (NumberOfVotesForWinner - NumberOfVotesForLoser)/k
        z <- yfit * ConsensusProportion
    }
    else {
        ind1 <- Y == classes[1]
        y[ind1] <- -1
        z <- matrix(numeric(nt * (Q - 1)), nrow = nt)
        for (j in 2:Q) {
            indk <- Y == classes[j]
            indOther <- !(ind1 | indk)
            y[indk] <- 1
            y[indOther] <- -1
            zA <- nncTest(X, y, XTest, k)
            y[indOther] <- 1
            zB <- nncTest(X, y, XTest, k)
            z[, j - 1] <- (zA + zB)/2
        }
    }
    z
}
