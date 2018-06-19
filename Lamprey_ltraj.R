#make variable where duplicate observations of fishId and seconds are removed
lampt1 <- Lamprey[!duplicated(Lamprey[c("fishID","seconds")]), ]

#change date and time to class POSIXct
lampdate <- as.POSIXct(lampt1$seconds,"1960-01-01", tz = "GMT")

#create object of class ltraj
lamprey1 <- as.ltraj(xy = lampt1[,c("x","y")], date = lampdate,
                     id = lampt1$fishID, typeII = TRUE,
                     infolocs = lampt1[,5:9])

# Creating net squared displacement histograms for Lampreys
lampreydf <- ld(lamprey1)
lamprey_ltraj <- dl(lampreydf)

lamp_split <- split(lampreydf,lampreydf$fishID)

plotlamp <- function(x) {
  hist(x$R2n, breaks=15, xlab="NSD", main=paste0("NSD ", unique(x$id), " n=", nrow(x[!is.na(x$R2n),])))
}

par(mfrow=c(2,3))
lapply(lamp_split, plotlamp)