bray_p_matrix <- function(otutable="otutable/plant_corpus.csv"){
        require(parallel)
        ot <- read.csv(otutable,row.names = 1,header = TRUE)
        ## For Bray-curtis
        func2 <- function(x) {
                require(vegan)
                perm <- function(x){return(sample(x))}
                ot <- read.csv(otutable,row.names = 1,header = TRUE)
                ot.perm <- apply(ot, MARGIN = 2, FUN = perm)
                bray.perm <- vegdist(t(ot.perm))
                return(bray.perm)
        }
        
        x=1:1000
        cl <- makeCluster(80)
        results2 <- parLapply(cl,x,func2) 
        stopCluster(cl) 
        
        bray.mean.perm <- matrix(0, nrow=ncol(ot), ncol=ncol(ot))
        
        for(j in 1:1000){
                bray.mean.perm  <- bray.mean.perm + as.matrix(results2[[j]])  
        }
        
        bray.mean.perm <- round(bray.mean.perm/1000,4)
        
        # Calculate Mean and SD with bootstrapping
        
        ## Bray-curtis 
        
        func4 <- function(x) {
                require(vegan)
                ot <- read.csv(otutable,row.names = 1,header = TRUE)
                bray.boot <- vegdist(t(ot[sample(nrow(ot),replace = TRUE),]))
                return(bray.boot)
        }
        
        bray.mean.boot <- matrix(0, nrow=ncol(ot), ncol=ncol(ot))
        bray.sd.boot <- matrix(0, nrow=ncol(ot), ncol=ncol(ot))
        
        x=1:1000
        cl <- makeCluster(80)
        results4 <- parLapply(cl,x,func4) 
        stopCluster(cl) 
        
        for(j in 1:1000){
                bray.mean.boot  <- bray.mean.boot + as.matrix(results4[[j]])  
        }
        for(k in 1:1000){
                bray.sd.boot  <- (as.matrix(results4[[k]])- bray.mean.boot/1000)^2
        }
        
        bray.mean.boot <- round(bray.mean.boot/1000,4)
        bray.sd.boot <- round(sqrt(bray.sd.boot/1000),4)
        
        # Calculate P matrix
        
        cor.p <- matrix(0, nrow=ncol(ot), ncol=ncol(ot))
        for(i in 1:ncol(bray.mean.perm)) {
                for (j in 1:ncol(bray.mean.perm)){        
                        p <- pnorm(q = abs(bray.mean.perm[i,j]),
                                   mean = abs(bray.mean.boot[i,j]),
                                   sd = bray.sd.boot[i,j],
                                   lower.tail = TRUE)
                        cor.p[i,j] <- p
                }
        }
        return(cor.p)
}

otutable.name <- list.files("otutable/")

for(i in c(1:14)){
        bray.p <- bray_p_matrix(otutable = paste("otutable/", otutable.name[i],sep = ""))
        write.csv(bray.p,paste("brayP/",otutable.name[i],sep=""))
        print(otutable.name[i])
}

