topo_sub <- NULL

for (i in 1:14){
        a <- read.graph(paste("graph_1/",graph.ls[i],sep=""),format = "graphml")
        a <- subgraph(a, v = igraph::degree(a)>0)
        topo <- c(vnum = length(V(a)$name),
                  vedge = length(E(a)),
                  d = diameter(a),
                  c = transitivity(a),
                  s = mean_distance(a),
                  b = mean(betweenness(a)),
                  m = modularity(a,membership(cluster_fast_greedy(a))),
                  cnum <- length(unique(membership(cluster_fast_greedy(a))))
        )
        topo_sub <- rbind(topo_sub,topo)
}

row.names(topo_sub) <- env.name1

write.csv(topo_sub,"topology_subnetworks.csv")

vertic.edge.df <- data.frame(vertex = topo_sub[,1],
                             edge = topo_sub[,2])

row.names(vertic.edge.df) <- env.name1

ggplot(vertic.edge.df, aes(x=vertex,y=edge,size=edge/vertex))+
        geom_point(aes(color=factor(1:14)))+
        theme_bw()+labs(x="Vertex number",y="Edge number")+
        geom_smooth(method = "glm",se=FALSE,color=grey(.6,.8),size=.5,linetype=2)+
        scale_color_d3(palette = "category20",alpha = .7)+
        theme(panel.grid = element_blank(),
              legend.position = c(.8,.2))+
        geom_text_repel(aes(label=env.name1),size=3)+
        scale_x_log10(breaks = c(10,100,1000,10000,10000),
                      labels = expression(10^1,10^2,10^3,10^4,10^5),
                      limits= c(10,10000)
                      )+
        scale_y_log10(breaks = c(10,100,1000,10000,100000),
                      limits= c(100,100000),
                      labels = expression(10^1,10^2,10^3,10^4,10^5)
                      )+
        scale_size_continuous(name="Mean\ndegree",
                              breaks = c(1,10,20,30),
                              limits=c(5,50))+
        guides(color=FALSE)

ggsave("manuscript/ver_edge.pdf", height = 4, width = 4)

clust.mod.df <- data.frame(vertex = topo_sub[,4],
                             edge = topo_sub[,7])

ggplot(clust.mod.df, aes(x=vertex,y=edge))+
        geom_point(aes(color=factor(1:14)))+
        theme_bw()+labs(x="Clusterring coefficient",y="Modularity")+
        geom_smooth(method = "glm",se=FALSE,color="#983421")+
        scale_color_d3(palette = "category20",alpha = .7)+
        theme(panel.grid = element_blank(),
              legend.position = "none")+
        geom_text_repel(aes(label=inhabit.name),size=3)
ggsave("clust_mod.pdf", height = 4, width = 4)

                
write.csv(topo_sub,file = "topo_sub.csv")

topo.name <- c("Vertex number",	"Edge number",	"Diameter",
               "Clustering coefficient","Average separation",
               "Average betweenness","Modularity","Number of modules")

colnames(topo_sub) <- topo.name

ggplot(melt(topo_sub),aes(x=Var1,y=value,fill=Var1))+
        geom_bar(stat = "identity")+
        facet_wrap(~Var2,scales = "free_y",nrow = 4)+
        theme_bw()+labs(x="",y="Topological property")+
        scale_fill_d3("category20")+
        theme(legend.position = "none",
              panel.grid = element_blank(),
              axis.text.x = element_text(angle = -80,
                                         hjust = 0))

ggsave("tope_sub.pdf",height = 6,width = 4)

t.test(topo_sub[1:8,2],
       topo_sub[9:14,2])

