#! /usr/bin/env Rscript
library(parallel)
library(vegan)
ot <- read.csv("soil_otu_reduce.csv",header = TRUE,row.names = 1)

func1 <- function(x) {
        require(vegan)
        ot <- read.csv("soil_otu_reduce.csv",header = TRUE,row.names = 1)
        bray.boot <- vegdist(t(ot)[sample(ncol(ot),replace = TRUE),])
        return(bray.boot)
}

bray.mean.boot <- matrix(0, nrow=ncol(ot), ncol=ncol(ot))
bray.sd.boot <- matrix(0, nrow=ncol(ot), ncol=ncol(ot))

x=1:1000
cl <- makeCluster(32)
results <- parLapply(cl,x,func1) 
stopCluster(cl) 

for(j in 1:1000){
    bray.mean.boot  <- bray.mean.boot + as.matrix(results[[j]])  
}
for(k in 1:1000){
    bray.sd.boot  <- (as.matrix(results[[k]])- bray.mean.boot/1000)^2
}

write.csv(round(bray.mean.boot/1000,4),"bray_mean_boot.csv")
write.csv(round(sqrt(bray.sd.boot/1000),4),"bray_sd_boot.csv")
