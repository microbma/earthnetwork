spearman_p_matrix <- function(otutable="otutable/animal_corpus.csv"){
        
        ot <- read.csv(otutable,row.names = 1,header = TRUE)
        require(parallel)
        func3 <- function(x) {
                ot <- read.csv(otutable,row.names = 1,header = TRUE)
                corr.boot <- cor(ot[sample(nrow(ot),replace = TRUE),],
                                 method = "spearman")
                return(corr.boot)
        }
        

        
        corr.mean.boot <- matrix(0, nrow=ncol(ot), ncol=ncol(ot))
        corr.sd.boot <- matrix(0, nrow=ncol(ot), ncol=ncol(ot))
        
        x=1:960
        cl <- makeCluster(80)
        results3 <- parLapply(cl,x,func3) 
        stopCluster(cl) 
        
        for(j in 1:960){
                corr.mean.boot  <- corr.mean.boot + as.matrix(results3[[j]])  
        }
        for(k in 1:960){
                corr.sd.boot  <- corr.sd.boot+(as.matrix(results3[[k]])- corr.mean.boot/960)^2
        }
        print("bootstrap")
        corr.mean.boot <- round(corr.mean.boot/960,4)
        corr.sd.boot <- round(sqrt(corr.sd.boot/960),4)
        return(list(mean=corr.mean.boot,
                    sd=corr.sd.boot))
}

otutable.name <- list.files("null/")

for(i in c(1:14)){
        bray.p <- spearman_p_matrix(paste("null/", otutable.name[i],sep = ""))
        write.csv(bray.p$mean,paste("null_corrP/boot_mean_",otutable.name[i],sep=""))
        write.csv(bray.p$sd,paste("null_corrP/boot_sd_",otutable.name[i],sep=""))
        print(otutable.name[i])
}






