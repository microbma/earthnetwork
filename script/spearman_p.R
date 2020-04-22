otutable.name <- list.files("otutable/")

for(k in 1:14){
        
        bray.mean.perm <- read.csv(paste("corrP/perm_",otutable.name[k],sep=""),row.names = 1)
        bray.mean.boot <- read.csv(paste("corrP/boot_mean_",otutable.name[k],sep=""),row.names = 1)
        bray.sd.boot <- read.csv(paste("corrP/boot_sd_",otutable.name[k],sep=""),row.names = 1)
        if (bray.mean.perm[1,1]==0){bray.mean.perm <- 1- bray.mean.perm}
        if (bray.mean.boot[1,1]==0){bray.mean.boot <- 1- bray.mean.boot}
        
        cor.p <- matrix(0, nrow=ncol(bray.mean.perm), ncol=ncol(bray.mean.perm))
        for(i in 1:ncol(bray.mean.perm)) {
                for (j in 1:ncol(bray.mean.perm)){        
                        p <- ks.test(x = abs(bray.mean.perm[i,j]),
                                     "pnorm",
                                     mean = abs(bray.mean.boot[i,j]),
                                     sd = bray.sd.boot[i,j],
                                     alternative = "greater")
                        cor.p[i,j] <- p$p.value
                }
        }
        write.csv(cor.p,paste("corrP/p_",otutable.name[k],sep=""))
        print(otutable.name[k])
}


