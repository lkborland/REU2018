## Define custom functions for Markov chains
library(data.table)


firstOrderMC <- function(tag    = '2059.535-1',
                         period = 'PreCO2',
                         species = "BHC",
                         n.iter = 100,
                         firstOrderDT = firstOrder){

    transition.vector <-
        firstOrderDT[ TagCodeTrial == tag & Period2 == period, ]
    fish.cell <- character(n.iter)

    fish.cell[1] <- transition.vector[ , sample(first.positions, 1)]
    for( index in 2:n.iter){
        fish.cell[index] <- transition.vector[ first.positions == fish.cell[index - 1],
                                              sample(second.positions, size = 1, prob =  prob)]
    }

    fish.cell.DT <- data.table::data.table(fish.cell)
    fish.summary <- fish.cell.DT[ , .(.N, prob = .N/nrow(fish.cell.DT)), by = .(fish.cell)]
    fish.summary[ , x := gsub("(\\d+)-(\\d+)", "\\1", fish.cell)]
    fish.summary[ , y := gsub("(\\d+)-(\\d+)", "\\2", fish.cell)]
    fish.summary[ , tag := tag]
    fish.summary[ , period := period]
    fish.summary[ , species := species]
    return(fish.summary)
}
