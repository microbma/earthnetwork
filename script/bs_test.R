#! /usr/bin/env Rscript
library(boot)
ot <- read.csv("soil_otu_reduce.csv",header = TRUE,row.names = 1)

mean.bs <- matrix(0, nrow=ncol(ot), ncol=ncol(ot))
sd.bs <- matrix(0, nrow=ncol(ot), ncol=ncol(ot))
cor.boot <- function(data, k) cor(data[k,],method = "spearman")[1,2]
for(i in 1:(ncol(ot)-1)){
        for (j in (i+1):ncol(ot)){
                cor.res <- boot(data=ot[,c(i,j)], 
                                statistic=cor.boot, 
                                R=999,parallel = "snow",
                                ncpus = 30)
                mean.bs[i,j] <- mean(cor.res$t)
                sd.bs[i,j] <- sd(cor.res$t)
        }
        print(paste("bs",i,sep = ""))
}

write.csv(mean.bs,"mean_bs.csv")
write.csv(sd.bs,"sd_bs.csv")

