#! /usr/bin/env Rscript

# Calculate P value
mean0 <- read.csv("brayP/perm_rhizosphere.csv",row.names = 1)
mean <- read.csv("brayP/boot_mean_rhizosphere.csv",row.names = 1)
sd <- read.csv("brayP/boot_sd_rhizosphere.csv",row.names = 1)

cor.p <- matrix(0, nrow=ncol(mean), ncol=ncol(mean))
for(i in 1:ncol(mean)) {
    for (j in 1:ncol(mean)){        
        p <- pnorm(q = 1-abs(mean0[i,j]),
                   mean = abs(mean[i,j]),
                   sd = sd[i,j], lower.tail = TRUE)
        cor.p[i,j] <- p
    }
        print(i)
}

p.fdr <- matrix(p.adjust(cor.p),nrow = ncol(mean0))

score <- mean
score[p.fdr>0.05] <- 0
diag(score) <- 0

g1 <- graph_from_adjacency_matrix(as.matrix(score),mode = "undirected",weighted = TRUE,diag = F)

write.graph(g1,file = "g1.graphml",format = "graphml")
