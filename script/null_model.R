# packages
library(vegan)
otutable.name <- list.files("otutable/")
for(i in c(1,3,4,6:10,12,14)){
        otu.tb <- read.csv(paste("otutable/",otutable.name[i],sep=""),
                           header = TRUE, row.names = 1)
        null.tb <- permatfull(otu.tb, times = 1,fixedmar = "col")$perm
        write.csv(null.tb,paste("null_1/null_",otutable.name[i],sep=""),row.names = TRUE)
        print(i)
}


# Calculate Spearman correlation matrix
require(WGCNA)
null.name <- list.files("null_1/")
for(j in 1:10){
        null.table <- read.csv(paste("null_1/",null.name[j],sep=""),
                       header = TRUE, row.names = 1)
        null.cor <- corAndPvalue(null.table,method="spearman")
        null.mat <- matrix(0,
                           nrow=dim(null.table)[2],
                           ncol=dim(null.table)[2])
        null.p <- matrix(p.adjust(null.cor$p, method = "BY"),
                         nrow=dim(null.cor$p)[1])
        null.mat[null.cor$cor>.5 & null.cor$p<0.05] <- 1
        fdr <- (table(null.mat)[2]-dim(null.cor$p)[1])/(dim(null.cor$p)[1]^2)
        print(fdr)
}
        
# Calculate Bray-Curtis matrix



# Generate 
