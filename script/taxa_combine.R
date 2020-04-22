# Find genus association combination

V(g.link)$species <- otu.df$taxa.df[V(g.link)$name,7]
el.df <- get.edgelist(g.link)

combine.1 <- paste(otu.df$taxa.df[el.df[,1],7],otu.df$taxa.df[el.df[,2],7],sep = "-")
combine.1c <- combine.1[otu.df$taxa.df[el.df[,1],7] != otu.df$taxa.df[el.df[,2],7]]

combine <- sort(table(c(combine.1c)),decreasing = TRUE)
combine.c <- combine[-grep("Unknown",names(combine))]

# Overrepresentations
overrep.sp <- function(p1 = "Gemmatimonas_aurantiaca", 
                       p2 = "Elusimicrobium_minutum",
                       g1 = g.link){
        vertex.sp <- V(g1)$species
        V(g1)$name <- vertex.sp
        v1 = table(vertex.sp)[p1]
        v2 = table(vertex.sp)[p2]
        prob <- v1*v2/(vcount(g1)*(vcount(g1)-1))
        size <- ecount(g1)
        q <- length(E(g1)[grep(p1, vertex.sp) 
                          %--% grep(p2, vertex.sp)])
        p <- pbinom(q,size,prob,lower.tail = FALSE)
        return(p)
        }


p.sp.combine <- NULL

p1.sp <- unlist(base::strsplit(names(combine.c[1:830]),split = "-"))[seq(1,1660,2)]
p2.sp <- unlist(base::strsplit(names(combine.c[1:830]),split = "-"))[seq(2,1660,2)]

for(i in 1:830){
        p.sp <- overrep.sp(p1 = p1.sp[i],
                           p2 = p2.sp[i])
        p.sp.combine <- c(p.sp.combine,p.sp)
        print(i)
}
p.sp.adj <- p.adjust(p.sp.combine,method = "fdr")

sp.combine <- names(combine.c[1:830])
sp.combine.c <- sp.combine[p.sp.adj<0.05]

sp.df <- data.frame(species1 = unlist(base::strsplit(sp.combine.c,split = "-"))[seq(1,1624,2)],
                    species2 = unlist(base::strsplit(sp.combine.c,split = "-"))[seq(2,1624,2)],
                    time = combine.c[1:830][p.sp.adj<0.05],
                    over.expression = p.sp.adj[p.sp.adj<0.05])

write.csv(sp.df, "sp_combine.csv")
