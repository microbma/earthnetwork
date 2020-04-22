# toy example
require(ape)
require(igraph)
require(ggtree)

## tree
tree <- read.tree("toy_tree")
gtre <- ggtree(tree,right = TRUE,size=1,color=grey(.3))
arrow.df <- data.frame(x=4:0,
                       yend=rep(.5,5),
                       y=rep(5.5,5),
                       xend=4:0)

gtre+geom_nodepoint(size=3)+
        geom_tiplab(hjust=2)+
        geom_tippoint()+
        geom_text2(aes(label=label,subset=!isTip),hjust=1.2)+
        geom_text(aes(x=x,y=y,label=label),fontface="italic",size=6,
                  color= qualitative_hcl(5,palette = "Dark2",alpha = 1),
                  data=data.frame(x=4:0,y=5.5,label=paste("t",4:0,sep="")))+
        geom_segment(data=arrow.df,aes(x=x,y=y,xend=xend,yend=yend),
                     inherit.aes = FALSE,arrow = arrow(length = unit(0.1, "inches"),type = "closed"),
                   linetype=1, size=1,
                   color = qualitative_hcl(5,palette = "Dark2",alpha = .3))+
        scale_x_reverse()+
        ylim(c(0,6))

ggsave("toy_tree.pdf",height = 2,width=6)

## network

### T1
pdf("toy_net.pdf",width = 12,height = 3)
par(mar=rep(2,4),mfrow=c(1,4))
g1 <- matrix(c("A","B","A","C","A","D","A","E","B","C","B","E","C","D","C","E","D","E"),
             ncol = 2,byrow = TRUE)
g.toy1 <- graph_from_edgelist(g1,directed = FALSE)
g.ly <- layout.circle(g.toy1)
plot(g.toy1,layout=g.ly,
     vertex.label.cex=2,
     vertex.size=50,
     vertex.label.dist=1,
     xlim=c(-1,1),ylim=c(-1,1),
     vertex.frame.color=NA,
     vertex.color= c(qualitative_hcl(3)[c(1,2)],rep("grey",3)),
     edge.lty = c(1,1,2,rep(1,6)),
     edge.width=c(1,5,5,5,5,5,1,1,1),
     edge.color=c("grey",qualitative_hcl(3)[c(1,1,1,2,2)],rep("grey",3)))

### T2
g2 <- matrix(c("Ancestor_1","C","Ancestor_1","E","C","D","C","E","D","E"),
             ncol = 2,byrow = TRUE)
g.toy2 <- graph_from_edgelist(g2,directed = FALSE)
plot(g.toy2,layout=g.ly[c(2,3,5,4),],
     xlim=c(-1,1),ylim=c(-1,1),
     vertex.label.cex=2,
     vertex.size=50,
     vertex.color= c(qualitative_hcl(3)[c(1,2)],rep("grey",3)),
     vertex.frame.color=NA,
     edge.lty = c(1,1,2,1,1),
     edge.width=c(1,5,5,5,1),
     edge.color=c("grey",qualitative_hcl(3)[c(1,2,2)],rep("grey",2)))

### T3
g3 <- matrix(c("Ancestor_2","E","D","E"),
             ncol = 2,byrow = TRUE)
g.toy3 <- graph_from_edgelist(g3,directed = FALSE)
plot(g.toy3,
     layout=g.ly[c(3,5,4),],
     vertex.frame.color=NA,
     vertex.label.cex=2,
     vertex.size=50,
     vertex.color= c(qualitative_hcl(3)[c(1,2)],rep("grey",3)),
     edge.lty = c(1,2),
     edge.width=c(1,5),
     edge.color=c("grey",qualitative_hcl(3)[c(1)]))

### T4
g4 <- matrix(c(0,0,0,0),
             ncol = 2,byrow = TRUE)
g.toy4 <- graph_from_adjacency_matrix(g4)
plot(g.toy4,layout = g.ly[c(4,5),],
     vertex.frame.color=NA,
     vertex.label=c("Ancestor_3","E"),
     vertex.color= c(qualitative_hcl(3)[c(1,2)],rep("grey",3)),
     edge.lty = c(1,2),
     vertex.label.cex=2,
     vertex.size=50,
     edge.width=c(1,3),
     edge.color=c("grey",qualitative_hcl(3)[c(1)]))
dev.off()
