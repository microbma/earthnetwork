# BA model merge delete node with smallest degree
library(igraph)

g0 = conet
for(i in 1:2927){
        deg <- degree(g0)
        pos <- sample(c(1:(2929-i))[deg == min(deg)],size = 1)
        g0.ba <- delete_vertices(graph = g0,v = pos) 
        write.graph(g0.ba,file = paste("graph_growth_ba/g", i,".graphml", sep = ""),format = "graphml")
        g0 <- g0.ba
        print(i)
}

