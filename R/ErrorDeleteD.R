ErrorDeleteD <-
function(Data, cl, FUN=FitRda,  d=0.33, ...){
#
# It calculates error rate using delete-d method. Number of variables selected beforehand.
#
# Data    : Combination of training and test data
# cl      : Combination of training and test class
# FUN     : Methods to be implemented, default is rda
# d       : Proportion to hold out for delete-d method
# nsim    : Number of times the CV should be performed
# lengthRes : length of output from FUN
#

    IndSample <- seq(cl)[-DeleteD(cl,d)]
    ErrorRow( IndSample, Data=Data, cl=cl,  FUN=FUN, ... )
}

