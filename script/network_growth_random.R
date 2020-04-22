# track the growth of tree
require(igraph)
require(ape)
library(ggtree)

# Read tree file
tree <- read.tree("tree/mc/otu.tree")

# molecular dating by penalised likelihood
chr <- chronos(tree)

# Give the edge length in chr to tree
tree1 <- tree
tree1$edge.length <- chr$edge.length

# get the value of nodes
gdata <- ggtree(tree1)
node.df <-  gdata$data[-c(1:2929),]
node.df.rank <- node.df[order(node.df$x,decreasing = TRUE),]
tip.df <- gdata$data[1:2928,]

### define a function for the m value (inital degree) and t values for each nodes ###

# define a function for merge the nodes and edges

MergeNode <- function(g = g0, p1 = 1119, p2 = 1672){
        a <- E(g)[from(p1)]
        b <- E(g)[from(p2)]
        g1 <- g-a[!a %in% b]
        c <- E(g1)[from(p2)]
        g2 <- g1-c[!c %in% a]
        return(g2)
}

# network growth
m0 <- NULL # initial links
g0 <- conet # assign network to g0

for(i in 2925:1){  # a cycle for all nodes in the tree
        ran.order <- sample(x = 1:i,size = 2,replace = FALSE)
        p1 = ran.order[1]
        p2 = ran.order[2]
        # contract the network
        mapping <- 1:(length(V(g0)$name)) 
        
        # define a numeric vector for mapping 
        # contract the network
        if(p1 < p2){
                mapping[p2] <- p1
                mapping1 <- mapping
                mapping1[mapping > p2] <- mapping[mapping > p2] - 1
                g1 <- contract(MergeNode(g0,p1=p1,p2=p2),mapping1,vertex.attr.comb=toString)
                g1 <- simplify(g1)
                V(g1)$name[p1] <- paste("node",i,sep="")
                m <- neighbors(g0,V(g0)$name[p1]) %in% neighbors(g0,V(g0)$name[p2])
                m0 <- c(m0,length(m[m==TRUE]))
                names(m0) <- paste("node",1:c(2926-i),sep="")
        }else{
                mapping[p1] <- p2
                mapping1 <- mapping
                mapping1[mapping > p1] <- mapping[mapping > p1] - 1
                g1 <- contract(MergeNode(g0,p1=p1,p2=p2),mapping1,vertex.attr.comb=toString)
                g1 <- simplify(g1)
                V(g1)$name[p2] <- paste("node",i,sep="")
                m <- neighbors(g0,V(g0)$name[p1]) %in% neighbors(g0,V(g0)$name[p2])
                m0 <- c(m0,length(m[m==TRUE]))
                names(m0) <- paste("node",1:c(2926-i),sep="")
        }
              
        g0 <- g1
        write_graph(g0,paste("graph_growth_random/g",i,sep=""),"graphml")
        print(i)
}
