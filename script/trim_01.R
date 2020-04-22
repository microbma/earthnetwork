# Keystone nodes

hub.finder <- function(g.cor){
        d <- degree(g.cor)
        bc <- betweenness(g.cor,normalized = TRUE)
        return(V(g.cor)$name[order(d,decreasing = TRUE)[1:10]])
}

hub.list <- list(hub.finder(read_graph(paste("trim_otu/trim_graph/",
                                             gsub("csv","graphml",otutb.name[2]),sep = ""),
                                       format = "graphml")),
                 hub.finder(read_graph(paste("trim_otu/trim_graph/",
                                             gsub("csv","graphml",otutb.name[3]),sep = ""),
                                       format = "graphml")),
                 hub.finder(read_graph(paste("trim_otu/trim_graph/",
                                             gsub("csv","graphml",otutb.name[4]),sep = ""),
                                       format = "graphml")),
                 hub.finder(read_graph(paste("trim_otu/trim_graph/",
                                             gsub("csv","graphml",otutb.name[5]),sep = ""),
                                       format = "graphml")),
                 hub.finder(read_graph(paste("trim_otu/trim_graph/",
                                             gsub("csv","graphml",otutb.name[7]),sep = ""),
                                       format = "graphml")),
                 hub.finder(read_graph(paste("trim_otu/trim_graph/",
                                             gsub("csv","graphml",otutb.name[8]),sep = ""),
                                       format = "graphml")),
                 hub.finder(read_graph(paste("trim_otu/trim_graph/",
                                             gsub("csv","graphml",otutb.name[9]),sep = ""),
                                       format = "graphml")),
                 hub.finder(read_graph(paste("trim_otu/trim_graph/",
                                             gsub("csv","graphml",otutb.name[10]),sep = ""),
                                       format = "graphml")),
                 hub.finder(read_graph(paste("trim_otu/trim_graph/",
                                             gsub("csv","graphml",otutb.name[11]),sep = ""),
                                       format = "graphml")),
                 hub.finder(read_graph(paste("trim_otu/trim_graph/",
                                             gsub("csv","graphml",otutb.name[12]),sep = ""),
                                       format = "graphml")),
                 hub.finder(read_graph(paste("trim_otu/trim_graph/",
                                             gsub("csv","graphml",otutb.name[13]),sep = ""),
                                       format = "graphml")),
                 hub.finder(read_graph(paste("trim_otu/trim_graph/",
                                             gsub("csv","graphml",otutb.name[14]),sep = ""),
                                       format = "graphml")))

hub.dis.df <- matrix(0,nrow = 60,ncol=12)
for(i in 1:12){
        hub.dis.df[unique(unlist(hub.list)) %in% hub.list[[i]],i] <- 1
}

trim.otu.name <- as.numeric(gsub("otu","",vertex.unique))
trim.otu.net <- DNAStringSet(otu.id[trim.otu.name])
names(trim.otu.net) <- vertex.unique
writeXStringSet(trim.otu.net,file = "trim_otu_net.fasta")

trim.taxonomy <- read.table("trim_otu.sintax",sep="\t")

source("sintax.R")
otu.df <- taxa.df(trim.taxonomy)

row.names(hub.dis.df) <- gsub("_genera_incertae_sedis","",otu.df$taxa.df[unique(unlist(hub.list)),6])
colnames(hub.dis.df) <- env.name1[c(2:5,7:14)]

paste(paste("italic(",row.names(hub.dis.df),")",sep=""),collapse = ",")

hub.dis1 <- hub.dis.df[rowSums(hub.dis.df)>1,
           order(colSums(hub.dis.df),decreasing = TRUE)]

pdf("hub_heatmap1.pdf",width = 4,height = 9)
pheatmap(hub.dis.df, 
         cellwidth = 4,
         cellheight = 4,
         legend = FALSE,         
         color = c("black","darkgreen"),
         treeheight_row = 20,
         treeheight_col = 20,
         cutree_rows = 3, 
         cutree_cols = 2,
         cluster_rows = TRUE, 
         cluster_cols = TRUE, 
         fontsize = 5)
dev.off()

d <- sort(table(otu.df$taxa.df[unique(unlist(hub.list)),3]),decreasing = TRUE)
d1 <- c(d[1:11], Others = sum(d[-c(1:11)]))
set.seed(1234)
pdf("hub_class.pdf",height = 5,width = 6)
pie(d1,border = .2,
    col = qualitative_hcl(n=71,
                          palette = "Dark3",
                          alpha = .9)[sample(71)][c(factor(names(d1),
                                                           levels = class.level))],
    radius = .6,init.angle = -45)
draw.circle(0,0,radius = .3,col = "white",border = F)
dev.off()

hub.time <- data.frame(times = names(table(rowSums(hub.dis.df))),
                       abundance = table(rowSums(hub.dis.df))
                       )
set.seed(1234)
ggplot(hub.time,aes(x = times, y = log(abundance.Freq)+1, fill = times))+
    geom_bar(stat = "identity")+
    geom_text(aes(label = paste(hub.time$abundance.Freq,"\n(",hub.time$times," networks)",sep=""),
                  y=log(abundance.Freq)+1.2,
                  angle = 360-c(360/14*seq(1,14,2))),
              fontface = 2)+
    geom_text(x=0,y=-2,label="Prevelance\nof hubs",
              size=4,fontface=2,lineheight=.8,)+
    scale_fill_manual(values = morandi.color[sample(1:36,12)])+
    ylim(c(-2,5))+
    #geom_hline(yintercept = c(.1),linetype=2,color="red")+
    coord_polar()+guides(fill=FALSE)+
    theme(axis.line = element_blank(),
          text = element_text(lineheight = .2),
          axis.text = element_blank(),
          axis.ticks = element_blank(),
          axis.title = element_blank(),
          panel.background = element_blank())    

