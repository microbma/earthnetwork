# Calculate P matrix
otutable.name <- list.files("otutable/")

for(k in 1:14){
        
        bray.mean.perm <- read.csv(paste("brayP/perm_",otutable.name[k],sep=""),row.names = 1)
        bray.mean.boot <- read.csv(paste("brayP/boot_mean_",otutable.name[k],sep=""),row.names = 1)
        bray.sd.boot <- read.csv(paste("brayP/boot_sd_",otutable.name[k],sep=""),row.names = 1)
        
        cor.p <- matrix(0, nrow=ncol(bray.mean.perm), ncol=ncol(bray.mean.perm))
        for(i in 1:ncol(bray.mean.perm)) {
                for (j in 1:ncol(bray.mean.perm)){        
                        p <- ks.test(x = abs(1-bray.mean.perm[i,j]),
                                     "pnorm",
                                     mean = abs(bray.mean.boot[i,j]),
                                     sd = bray.sd.boot[i,j],
                                     alternative="two.side")
                        cor.p[i,j] <- p$p.value
                }
        }
        write.csv(cor.p,paste("brayP/p_",otutable.name[k],sep=""))
        print(otutable.name[k])
}

