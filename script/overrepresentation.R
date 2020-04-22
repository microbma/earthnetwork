# Over-representation
require(ggplot2)
require(igraph)
require(corrplot)

g1 <- g.link

overrep <- function(p1 = "Firmicutes", 
                    p2 = "Firmicutes",
                    g1 = g.link){
        vertex.phylum <- V(g1)$class
        V(g1)$name <- vertex.phylum
        
        if (p1 != p2){
                v1 = table(vertex.phylum)[p1]
                v2 = table(vertex.phylum)[p2]
                m <- v1*v2*(vcount(g1)-1)/vcount(g1)
                n <- vcount(g1)*(vcount(g1)-1) - m
                k <- ecount(g1)
                x <- length(E(g1)[grep(p1, vertex.phylum) 
                                  %--% grep(p2, vertex.phylum)])
                p = phyper(x, m, n, k)}
        else {
                v1 = table(vertex.phylum)[p1]
                m <- v1*c(v1-1)*(vcount(g1)-1)/vcount(g1)
                n <- vcount(g1)*(vcount(g1)-1) - m
                k <- ecount(g1)
                x <- length(E(g1)[grep(p1, vertex.phylum) 
                                  %--% grep(p2, vertex.phylum)])
                p = phyper(x, m, n, k)
        }
}
vertex.phylum <- V(g.link)$class
overrep <- function(p1 = "Firmicutes", 
                    p2 = "Firmicutes",
                    g1 = g.link){
        vertex.phylum <- V(g1)$class
        V(g1)$name <- vertex.phylum
        
        if (p1 != p2){
                v1 = table(vertex.phylum)[p1]
                v2 = table(vertex.phylum)[p2]
                prob <- v1*v2/(vcount(g1)*(vcount(g1)-1))
                size <- ecount(g1)
                q <- length(E(g1)[grep(p1, vertex.phylum) 
                                  %--% grep(p2, vertex.phylum)])
                p = pbinom(q,size,prob,lower.tail = FALSE)}
        else {
                v1 = table(vertex.phylum)[p1]
                prob <- v1*(v1-1)/(vcount(g1)*(vcount(g1)-1))
                size <- ecount(g1)
                q <- length(E(g1)[grep(p1, vertex.phylum) 
                                  %--% grep(p2, vertex.phylum)])
                p = pbinom(q,size,prob,lower.tail = FALSE)}
        }



p.over <- matrix(nrow = 85,ncol = 85)
for(i in 1:85){
        for(j in 1:85){
                vertex.phylum <- V(g.link)$class
                po <- overrep(p1 = names(sort(table(vertex.phylum),decreasing = TRUE))[i], 
                              p2 = names(sort(table(vertex.phylum),decreasing = TRUE))[j])
                p.over[i,j] <- po
                }
        print(i)
        }

p.over.adj <- matrix(p.adjust(p.over,method = "fdr"),nrow = 85)
require(corrplot)

row.names(p.over.adj) <- names(sort(table(vertex.phylum),decreasing = TRUE))
colnames(p.over.adj) <- names(sort(table(vertex.phylum),decreasing = TRUE))
tiff(filename = "manuscript/overrep_dom.tiff",
     type = "cairo",
     compression = "lzw",
     width = 2000,
     height = 4000,
     res = 500)
as.character(dom.class) -> dom.class
corrplot(corr = 1-p.over.adj[dom.class,dom.class],
         #order = "hclust",
         #hclust.method = "single",
         col = getPalette(1),
         tl.col = "black",
         tl.cex=.6,
         #tl.pos = "ld",
         cl.pos = "n",
         method = "square",
         insig = "blank",
         sig.level = 0.05, 
         p.mat = p.over.adj[dom.class,dom.class],
         type = "lower"
         )
dev.off()

# Provence 

dom.class <- df1$region
df.tax.link <- function(g = g.p1){
        el <- get.edgelist(g)
        a <- V(g.link)$class
        names(a) <- V(g.link)$name
        a.rank <- a[el[,1]]
        b.rank <- a[el[,2]]
        ab.edge <- paste(a.rank,b.rank,sep=" ")
        ab.edge <- ab.edge[a.rank %in% dom.class & b.rank %in% dom.class]
        tax.df <- NULL
        for(i in 1:15){
            edge.spe <- ab.edge[grep(dom.class[i],ab.edge)]
            tax.link <- 1:15
            for(j in c(1:15)[-i]){
                    tax.link[j] <-  length(grep(dom.class[j], edge.spe))
            }
            tax.link[i] <- length(grep(paste(dom.class[i],dom.class[i]),ab.edge))
            tax.df <- rbind(tax.df,tax.link)
        }
       return(tax.df)
} 

df.g.entire <- df.tax.link(g.link)
file.list <- list.files("graph_1/")

tiff(filename = "manuscript/overrep_subnet.tiff",
     type = "cairo",
     compression = "lzw",
     width = 4000,
     height = 2000,
     res = 500)
par(mfrow = c(2,7))
for(i in 1:14){
        g.p1 <- read.graph(paste("graph/",file.list[i],sep=""),"graphml")
        df.p1 <- df.tax.link(g.p1)
        overrep.p.mat <- matrix(nrow = 15,ncol = 15)
        for(j in 1:15){
                for(k in 1:15){
                    p <- pbinom(q = df.p1[j,k], 
                                size = df.g.entire[i,j], 
                                prob = ecount(g.p1)/ecount(g.link),
                                lower.tail = FALSE)
                    overrep.p.mat[j,k] <- p
                }}
        write.csv(x = overrep.p.mat,
                  file = paste("overrepresentation/",
                               file.list[i],
                               sep=""))
        overrep.p.adj <- matrix(p.adjust(overrep.p.mat,method = "fdr"),nrow = 15)
        corrplot(corr = 1-overrep.p.adj,
                 col = getPalette(1),
                 tl.pos = "n",
                 cl.pos = "n",
                 #method = "square",
                 insig = "blank",
                 sig.level = 0.05,
                 type = "lower",
                 p.mat = overrep.p.adj
                )
        
        print(i)
        }

dev.off()

