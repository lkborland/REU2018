library(ape)
library(data.table)

## Acceleration function

accFunction <- function(index, tagTrialList, accelerationByCell,
                        firstOrder, allCells){

    ## Extract out fish
    ac1 <- accelerationByCell[ TagCodeTrial == tagTrialList[ index ,
                                                            TagCodeTrial] &
                               Period2 == tagTrialList[ index,
                                                       Period2] , ][
        order(xCell, yCell), ]
    
    setkey(ac1, xCell, yCell)
    ## Create data.table with all cells 
    acUse <- ac1[allCells]
    ## set missing cells to be zero 
    acUse[ is.na(acc), acc := 0]
    
    ## Create Markov Matrix    
    inputVector = firstOrder[ TagCodeTrial == tagTrialList[ index ,
                                                           TagCodeTrial] &
                              Period2 == tagTrialList[ index, Period2] ]
    
    allCells[ , cellCoord := paste( xCell, yCell, sep = "-")]
    
    allTransitionsDT <-
        data.table(allCells[ , expand.grid(first = cellCoord,
                                           second = cellCoord)])
    allTransitionsDT[ , firstTransition := paste(first, second, sep = ":")]
    
    setkey(allTransitionsDT, "firstTransition")
    setkey(inputVector, "firstTransition")
    
    allTransitionsDT2 <-
        inputVector[ , .(firstTransition, prob)][allTransitionsDT]
    allTransitionsDT2[ is.na(prob), prob := 0]
    
    nRows <- sqrt(nrow(allTransitionsDT2))
    rowColNames <- as.character(allTransitionsDT2[ , unique(first)])
    
    firstOrderMatrix <- matrix(data = allTransitionsDT2[ , prob],
                               nrow = nRows,
                               byrow = TRUE,
                               dimnames = list(rowColNames,
                                               rowColNames))
    
    out <- as.data.frame(ape::Moran.I(acUse[ , acc], firstOrderMatrix))
    out$TagCodeTrial <- tagTrialList[ index , TagCodeTrial]
    out$Period <-  tagTrialList[ index, Period2]   
    return(out)
}

## First order observations
firstFunction <- function(index, tagTrialList,
                          trialCellAve,
                          firstOrder, allCells,
                          varUse = "dist"
                          ){

    ## Extract out fish
    ac1 <- trialCellAve[ TagCodeTrial == tagTrialList[ index ,
                                                            TagCodeTrial] &
                               Period2 == tagTrialList[ index,
                                                       Period2] , ][
        order(x.zone, y.zone), ]
    
    setkey(ac1, x.zone, y.zone)
    ## Create data.table with all cells 
    acUse <- ac1[allCells]
    ## set missing cells to be zero 
    acUse <-
        acUse[  is.na(eval(as.symbol(varUse))), (varUse) := 0]

    ## Create Markov Matrix    
    inputVector = firstOrder[ TagCodeTrial == tagTrialList[ index ,
                                                            TagCodeTrial] &
                               Period2 == tagTrialList[ index, Period2] ]
    
    allCells[ , cellCoord := paste( xCell, yCell, sep = "-")]
    
    allTransitionsDT <-
        data.table(allCells[ , expand.grid(first = cellCoord,
                                           second = cellCoord)])
    allTransitionsDT[ , firstTransition := paste(first, second, sep = ":")]
    
    setkey(allTransitionsDT, "firstTransition")
    setkey(inputVector, "firstTransition")
    
    allTransitionsDT2 <-
        inputVector[ , .(firstTransition, prob)][allTransitionsDT]
    allTransitionsDT2[ is.na(prob), prob := 0]
    
    nRows <- sqrt(nrow(allTransitionsDT2))
    rowColNames <- as.character(allTransitionsDT2[ , unique(first)])
    
    firstOrderMatrix <- matrix(data = allTransitionsDT2[ , prob],
                               nrow = nRows,
                               byrow = TRUE,
                               dimnames = list(rowColNames,
                                               rowColNames))
    
    out <- as.data.frame(ape::Moran.I(acUse[ , eval(as.symbol(varUse)) ], 
                                      firstOrderMatrix))
    out$TagCodeTrial <- tagTrialList[ index , TagCodeTrial]
    out$Period <-  tagTrialList[ index, Period2]   
    return(out)
}
