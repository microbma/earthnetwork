# overrepression of vertex envronmental types in modules

require(ggplot2)
require(igraph)
require(corrplot)

# phyper(x, m, n, k)

mod.number <- table(modular.df1$modularity_class)
mod.number1 <- mod.number[pmatch(mod.name, names(mod.number))]

overrep.mod <- function(mod = mod.inhab,i = 2, j = 1){
        x =  mod[j,i] #the number of vertex from env_i in module_j
        m =  sum(mod[,i]) #vertex from env_i
        n =  2928 - m #total vertex - m
        k =  mod.number1[j] #vertex number in module_j
        p = phyper(x, m, n, k)
        return(p)
        }

# cycle to get the p values of overrepression in modules

p.over.mod <- NULL

for(i in 1:14){
        for (j in 1:8){
                p = 1 - overrep.mod(mod = mod.inhab, i, j)
                p.over.mod <- c(p.over.mod, p)
        }
}

p.over.mod.adj <- matrix(p.adjust(p.over.mod, method = "fdr"),nrow = 8,ncol = 14)

row.names(p.over.mod.adj) <- paste("M",c(2,3,1,3,5,6,7,8),sep="")
colnames(p.over.mod.adj) <- colnames(mod)

# make a heatmap
require(corrplot)

corrplot(1-p.over.mod.adj,
         p.mat = p.over.mod.adj,
         method = "square",
         sig.level = .05,
         insig = "blank")

over.mod.ht <- 1-p.over.mod.adj
over.mod.ht[p.over.mod.adj>0.05] <- 0

pdf("overexp_mod.pdf",width = 5,height = 5)
pheatmap(over.mod.ht[order(c(2,4,1,3,5,6,7,8)),order(colnames(over.mod.ht))],
         color = c("black","#006658"),
         cluster_rows = FALSE,
         cluster_cols = FALSE,
         labels_row = paste("M",1:8,sep=""),
         legend = FALSE,
         cellwidth = 20,
         cellheight = 20)
dev.off()



