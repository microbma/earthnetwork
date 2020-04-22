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
        
        x=1:100
        cl <- makeCluster(4)
        results2 <- parLapply(cl,x,func2) 
        stopCluster(cl) 
        
        bray.mean.perm <- matrix(0, nrow=ncol(ot), ncol=ncol(ot))
        
        for(j in 1:100){
                bray.mean.perm  <- bray.mean.perm + as.matrix(results2[[j]])  
        }
        
        bray.mean.perm <- round(bray.mean.perm/100,4)
        return(1-bray.mean.perm)
}

otutable.name <- list.files("null/")

for(i in c(1:14)){
        bray.p <- bray_p_matrix(otutable = paste("null/", otutable.name[i],sep = ""))
        write.csv(bray.p,paste("null_brayP/perm_",otutable.name[i],sep=""))
        print(otutable.name[i])
}
