#connection between modulars

## modular #1

cluster.group <- modular.df1$modularity_class
names(cluster.group) <- modular.df1$v_name

cluster.connect.matrix <- NULL
for(i in 1:8){
        cluster.connect <- rep(0,8)
        names(cluster.connect) <- mod.name[1:8]
        cluster.node <- modular.df1$v_name[modular.df1$modularity_class == mod.name[i]]
        edgelist.glink <- as_edgelist(g.link)
        edgelist.cluster <- edgelist.glink[edgelist.glink[,1] %in% cluster.node|edgelist.glink[,2] %in% cluster.node,]
        
        cluster.cluster.group <- cluster.group[c(edgelist.cluster)]
        cluster.table <- table(cluster.cluster.group)
        cluster.table1 <- cluster.table[names(cluster.table) %in% mod.name[1:8]]
        cluster.connect[names(cluster.table1)] <- cluster.table1
        cluster.connect.matrix <- rbind(cluster.connect.matrix,cluster.connect)
}

cluster.sequence <- c(2, 4, 1, 3, 5, 6, 7, 8)
#cluster.connect.matrix <- as.matrix(cluster.connect.matrix[cluster.sequence,cluster.sequence])

row.names(cluster.connect.matrix) <- paste("M",cluster.sequence,sep=" ")
colnames(cluster.connect.matrix) <- paste("M",cluster.sequence,sep=" ")

cluster.connect.matrix <- as.matrix(cluster.connect.matrix[order(cluster.sequence),order(cluster.sequence)])

#cluster.connect.matrix[lower.tri(cluster.connect.matrix)] <- 0
#diag(cluster.connect.matrix) <- 0
ggplot(melt(cluster.connect.matrix),aes(x=Var1,y=Var2,size=log(value,10)))+
        geom_point(color="#339900")+
        scale_size_continuous(range = c(0,10),
                              name = "Edge\nnumber",
                              breaks = c(1,2,3,4),
                              labels = expression(10^1,10^2,10^3,10^4))+
        theme_bw()+xlab("Network module")+ylab("Network module")+
        #scale_color_continuous(low = "black",high = "green")+
        theme(panel.background = element_rect(fill="black"),
              panel.grid = element_blank(),
              #axis.text.x = element_text(angle = -45, hjust = 0),
              aspect.ratio = 1)

ggsave("cluster_links.tiff",height = 4,width = 4,dpi = 600)   
