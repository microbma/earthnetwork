spearman_p_matrix <- function(otutable="otutable/rhizosphere.csv"){
        
        ot <- read.csv(otutable,row.names = 1,header = TRUE)
        
        require(parallel)
        func3 <- function(x){
                require(vegan)
                ot <- read.csv(otutable, row.names = 1, header = TRUE)
                bray.boot <- vegdist(t(ot[sample(nrow(ot),replace = TRUE),]),method = "bray")
                return(1-as.matrix(bray.boot))
        }
        
        bray.mean.boot <- matrix(0, nrow=ncol(ot), ncol=ncol(ot))
        bray.sd.boot <- matrix(0, nrow=ncol(ot), ncol=ncol(ot))
        
        x=1:960
        cl <- makeCluster(80)
        results3 <- parLapply(cl,x,func3) 
        stopCluster(cl) 
        
        for(j in 1:960){
                bray.mean.boot  <- bray.mean.boot + as.matrix(results3[[j]])  
        }
        for(k in 1:960){
                bray.sd.boot  <- bray.sd.boot+(results3[[k]]- bray.mean.boot/960)^2
        }
        print("bootstrap")
        bray.mean.boot <- round(bray.mean.boot/960,4)
        bray.sd.boot <- round(sqrt(bray.sd.boot/960),4)
        return(list(mean = bray.mean.boot,
                    sd = bray.sd.boot))
}

otutable.name <- list.files("null/")

for(i in c(1:14)){
        bray.p <- spearman_p_matrix(otutable=paste("null/", otutable.name[i],sep = ""))
        write.csv(bray.p$mean, paste("null_brayP/boot_mean_",otutable.name[i],sep=""))
        write.csv(bray.p$sd, paste("null_brayP/boot_sd_",otutable.name[i],sep=""))
        print(otutable.name[i])
}










