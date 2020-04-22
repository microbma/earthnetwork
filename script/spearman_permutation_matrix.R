spearman_p_matrix <- function(otutable="otutable/rhizosphere.csv"){
        
        ot <- read.csv(otutable,row.names = 1,header = TRUE)
        func1 <- function(x) {
                require(vegan)
                perm <- function(x){return(sample(x))}
                ot <- read.csv(otutable,row.names = 1,header = TRUE)
                ot.perm <- apply(ot, MARGIN = 2, FUN = perm)
                cor.perm <- cor(decostand(ot.perm, method = "normalize",MARGIN = 2),method = "spearman")
                return(cor.perm)
        }
        
        require(parallel)
        
        x=1:20
        cl <- makeCluster(4)
        results1 <- parLapply(cl,x,func1) 
        stopCluster(cl) 
        
        corr.mean.perm <- matrix(0, nrow=ncol(ot), ncol=ncol(ot))
        
        for(j in 1:1000){
                corr.mean.perm  <- corr.mean.perm + as.matrix(results1[[j]])  
        }
        
        corr.mean.perm <- round(corr.mean.perm/1000,4)
        print("permutation")
        return(corr.mean.perm)
}

otutable.name <- list.files("otutable/")

for(i in c(4:14)){
        bray.p <- spearman_p_matrix(paste("otutable/", otutable.name[i],sep = ""))
        write.csv(bray.p,paste("corrP/perm_",otutable.name[i],sep=""))
        print(otutable.name[i])
}









