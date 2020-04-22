# Loading package
require(igraph)
getPalette = colorRampPalette(brewer.pal(8,"Dark2"))

fisher.comb <- function (pvalues){
        df=length(pvalues)
        ch2=(-2*sum(log(pvalues)))
        pchisq(ch2, df=df, lower.tail=FALSE)
} 

Net_Infer <- function(p.mat.corr = p.mat.corr, 
                      p.mat.bray = p.mat.bray,
                      cor.mat = cor.mat, 
                      bray.mat = bray.mat,
                      otu.tb = otu.tb,
                      graph.name = graph.name,
                      cor.threshold = cor.threshold,
                      bray.threshold = bray.threshold
                      ){
        # p merge
        q = -2*log(p.mat.corr)-2*log(p.mat.bray)
        p.mat <- matrix(nrow = dim(p.mat.corr)[1],ncol = dim(p.mat.corr)[1])
        for(j in 1:dim(p.mat.corr)[1]){
                for(k in 1:dim(p.mat.corr)[1]){
                p.mat[j,k] <- pchisq(q[j,k],4,lower.tail = FALSE)
        }}
        # p adjust
        p.adj <- matrix(p.adjust(p.mat),nrow = dim(p.mat))
        diag(p.adj) <- 1
        # correlation matrix
        adj.mat <- matrix(nrow = dim(p.mat)[1],ncol = dim(p.mat)[1])
        adj.mat[] <- 1
        adj.mat[p.mat>0.05] <- 0
        adj.mat[abs(cor.mat) < cor.threshold] <- 0
        adj.mat[abs(bray.mat)< bray.threshold] <- 0
        g.adj <- graph_from_adjacency_matrix(adj.mat,diag = FALSE,mode = "undirected")
        # annotation the nodes information with OTU number
        V(g.adj)$name <- colnames(otu.tb)
        E(g.adj)$direction <- cor(otu.tb, method = "spearman")[lower.tri(adj.mat)][adj.mat[lower.tri(adj.mat)]==1]
        write.graph(g.adj,
                    file = paste("graph_cor/", 
                                  graph.name,
                                 sep = ""),
                    format = "graphml")
}

plot(g.adj)

# Reading a P matrix
p.corr.name <- list.files("corrP/")[1:14+28]
p.bray.name <- list.files("brayP/")[1:14+28]
cor.name <- list.files("corrP/")[1:14]
bray.name <- list.files("brayP/")[1:14]
otutable.name <- list.files("otutable/")

# RMT
require(RMThreshold)
for(i in 1:14){
cor.mat <- read.csv(paste("corrP/",cor.name[i],sep=""),
                    header = TRUE, row.names = 1)
rm.get.threshold(as.matrix(cor.mat),interactive = F,interval = c(.2,.8))
}

cor.threshold.list <- c(0.4,.53,0.47,.44,.51,.56,.58,.44,.52,.56,.53,.54,.52,.51)

for(i in 1:14){
        i=i+1
        bray.mat <- read.csv(paste("brayP/",cor.name[i],sep=""),
                            header = TRUE, row.names = 1)
        rm.get.threshold(as.matrix(bray.mat),interactive = F)
}

bray.threshold.list <- c(.27,.45,.44,.36,.38,.52,.34,.38,.47,.5,.36,.32,.35,.39) - 0.1

# Network
for(i in 1:14) {
        p.mat.corr <- read.csv(paste("corrP/",p.corr.name[i],sep=""),
                               header = TRUE,row.names = 1)
        p.mat.bray <- read.csv(paste("brayP/",p.bray.name[i],sep=""),
                               header = TRUE,row.names = 1)
        cor.mat <- read.csv(paste("corrP/",cor.name[i],sep=""),
                            header = TRUE, row.names = 1)
        bray.mat <- read.csv(paste("brayP/",bray.name[i],sep=""),
                             header = TRUE, row.names = 1)
        otu.tb <- read.csv(paste("otutable/",otutable.name[i],sep=""),
                           header = TRUE, row.names = 1)
        graph.name <- gsub("csv","graphml",otutable.name[i])
        bray.threshold <- bray.threshold.list[i]
        cor.threshold <- cor.threshold.list[i]
        
        Net_Infer(p.mat.corr = p.mat.corr, 
                  p.mat.bray = p.mat.bray,
                  cor.mat = cor.mat, 
                  bray.mat = bray.mat,
                  otu.tb = otu.tb,
                  cor.threshold = cor.threshold,
                  bray.threshold = bray.threshold,
                  graph.name = graph.name
                  )
        print(i)
}

# Check the direction
par(mfcol=c(4,4))
for(i in 1:14){
        g.cor <- read_graph(paste("graph_cor/",graph.ls[i],sep = ""),format = "graphml")
        E(g.cor)$direction->a
        hist(a)
}

# Entire network
graph.ls <- list.files("graph_1/")
g.entire <- read.graph(paste("graph_1/",graph.ls,sep=""),format = "graphml")
for (i in 2:14){
        g.entire <- g.entire + read.graph(paste("graph_1/",graph.ls[i],sep=""),format = "graphml")
}

for (i in 2:14){
        a <-read.graph(paste("graph_1/",graph.ls[i],sep=""),format = "graphml")
        b <- E(a)$direction
        print(b[b<0])
}

V(g.entire)$group.number <- vertex.presence.sum
E(g.entire)$group.number <- edge.presence.sum

vertex.presence.pos <- c(1:4591)[vertex.presence[,1]] 
vertex.sample <- rep(0,4591)
for(i in vertex.presence.pos){
        vertex.sample[i] <-  c(1:14)[vertex.presence[i,]]
}

V(g.entire)$s1 <- rep(0,4591)
V(g.entire)$s1[vertex.presence[,1]]<-1
V(g.entire)$s2 <- rep(0,4591)
V(g.entire)$s2[vertex.presence[,2]]<-1
V(g.entire)$s3 <- rep(0,4591)
V(g.entire)$s3[vertex.presence[,3]]<-1
V(g.entire)$s4 <- rep(0,4591)
V(g.entire)$s4[vertex.presence[,4]]<-1
V(g.entire)$s5 <- rep(0,4591)
V(g.entire)$s5[vertex.presence[,5]]<-1
V(g.entire)$s6 <- rep(0,4591)
V(g.entire)$s6[vertex.presence[,6]]<-1
V(g.entire)$s7 <- rep(0,4591)
V(g.entire)$s7[vertex.presence[,7]]<-1
V(g.entire)$s8 <- rep(0,4591)
V(g.entire)$s8[vertex.presence[,8]]<-1
V(g.entire)$s9 <- rep(0,4591)
V(g.entire)$s9[vertex.presence[,9]]<-1
V(g.entire)$s10 <- rep(0,4591)
V(g.entire)$s10[vertex.presence[,10]]<-1
V(g.entire)$s11 <- rep(0,4591)
V(g.entire)$s11[vertex.presence[,11]]<-1
V(g.entire)$s12 <- rep(0,4591)
V(g.entire)$s12[vertex.presence[,12]]<-1
V(g.entire)$s13 <- rep(0,4591)
V(g.entire)$s13[vertex.presence[,13]]<-1
V(g.entire)$s14 <- rep(0,4591)
V(g.entire)$s14[vertex.presence[,14]]<-1

g.link <- subgraph(g.entire, degree(g.entire)>0)
V(g.link)$phylum <- otu.df$taxa.df[,2]
 
write.graph(g.entire,file = "entire.graphml",format = "graphml")
write.graph(g.link,file = "glink.graphml",format = "graphml")

otu.name <- as.numeric(gsub("otu","",V(g.link)$name))
otu.net <- DNAStringSet(otu.id[otu.name])
names(otu.net) <- V(g.link)$name

writeXStringSet(otu.net,file = "otu_net.fasta")

taxonomy <- read.table("otu.sintax",sep="\t")

source("../Wangmengchen/rice_endo/sintax.R")
otu.df <- taxa.df(taxonomy)

edge.df <- data.frame(edge = degree(g.link),
                      phylum = V(g.link)$phylum,
                      class = V(g.link)$class)

require(tidyverse)
edge.mean <- edge.df %>%
        group_by(class) %>%
        summarise(md = mean(edge)) %>%
        arrange(desc(md))

edge.df$class <- factor(edge.df$class, 
                         levels = as.character(edge.mean$class))

ggplot(edge.df, aes(x = class, y = edge, fill = class))+
        geom_boxplot()+
        #geom_jitter()+
        scale_y_log10()+
        xlab("Class")+
        ylab("Edge number")+
        theme_bw()+guides(fill=FALSE,color=FALSE)+
        scale_fill_discrete_divergingx()+
        theme(axis.text.x = element_text(angle=-90,hjust = 0))
ggsave("edge_class.tiff",width = 6,height = 4)


# Pie

modular.df <- read.csv("network/entire_graph_store.csv", header = TRUE)
mod.name <- names(sort(table(modular.df$modularity_class),decreasing = TRUE)[1:9])
class.level <- factor(unique(otu.df$taxa.df[,3]))
color <- getPalette(85)
c(class.level[class.level %in% names(mod.tax.per1)[1:15]])

par(mfcol=c(3,3),mar=rep(0,4))
for(i in 1:9){
        tiff(paste("network/pie_mod_",i,".tiff",sep=""),
             width = 500,
             height = 500,
             compression = "lzw",
             type = "cairo",
             bg = "black")
        mod.tax <- otu.df$taxa.df[,3][modular.df$modularity_class == mod.name[i]]
        mod.tax.per <- sort(table(mod.tax), decreasing = TRUE)
        mod.tax.per1 <- c(mod.tax.per[1:15], 
                          Others = sum(mod.tax.per[-c(1:15)]))
        col.level <- class.level[pmatch(names(mod.tax.per1)[1:15], class.level)]
        print(col.level)
        pie(mod.tax.per1,
            labels = NA,
            radius = .5,
            border = FALSE, 
            col =  c(color[c(col.level)],"grey"),
            init.angle = -90,
            cex=1)
        pie.labels(angles = c(-1,-.3,.5,.8,1,seq(2,3.7,.18),4.5), 
                   minangle = .1,
                   col= c("white"),
                   labels = names( mod.tax.per1),
                   radius = .6)
        draw.circle(0,0,radius = .2, border = "black", col = "black")
        dev.off()
}

# modular V.S. inhabit
modular.df1 <- modular.df[,order(colnames(modular.df))]
mod.inhab <- NULL
for(i in 1:8){
modular.df2 <- colSums(modular.df1[modular.df1$modularity_class == as.numeric(mod.name[i]),10:23])
mod.inhab <- rbind(mod.inhab,modular.df2)
}

row.names(mod.inhab) <- paste("Cluster",c(2, 4, 1, 3, 5, 6, 7, 8),sep="")
colnames(mod.inhab) <- env.name1[c(1,10:14,2:9)]

mod.inhab.df <- melt(mod.inhab)

mod.inhab.df$Var2 <- factor(mod.inhab.df$Var2,levels = env.name1)

mod.inhab.df$Var1 <- factor(mod.inhab.df$Var1, 
                            labels = paste("M",1:8,sep=""),
                            levels = paste("Cluster",1:8,sep=""))

ggplot(mod.inhab.df, aes(x=Var2, y=Var1, size = value, color=value))+
        geom_point(color="#339900")+
        scale_size_continuous(range = c(0,10),name = "Vertex\nnumber",
                              breaks = c(10,50,100,200),
                              labels = expression(10,50,100,200))+
        theme_bw()+xlab("Environmental type")+ylab("Network module")+
       # scale_color_continuous(low = "#339900",high = "#339900")+
        theme(panel.background = element_rect(fill="black"),
              panel.grid = element_blank(),
              axis.text.x = element_text(angle = -45, hjust = 0),
              aspect.ratio = 2/3)

ggsave("cluster_source.pdf",
       height = 5,
       width = 5)   

ggsave("cluster_source.tiff",
       height = 3.5,
       width = 6,
       dpi = 600)   


# cluster the inhibitation

pdf("manuscript/inh_mod_hc.pdf",width = 4,height = 6)
plot(agnes(daisy(t(mod.inhab))),which.plots=2,fromLeft=FALSE,
       main = "Dendrogram of environmental types",
       xlab="",ylab="",
       mar=c(0,0,0,0))
dev.off()

# Network

cor.inh <- vegdist(t(mod.inhab),method = "jaccard")
cor.inh.d <- Network_Enhancement(as.matrix(cor.inh))

row.names(cor.inh.d) <- colnames(mod.inhab)
colnames(cor.inh.d) <- colnames(mod.inhab)

cor.inh.d1 <- as.matrix(cor.inh.d)
cor.inh.d1[cor.inh.d<0.5] <- 0

g.inh.mod <- graph_from_adjacency_matrix(cor.inh.d1, 
                                         diag = F,
                                         mode = "undirected", 
                                         weighted = TRUE)

plot(g.inh.mod)

write.graph(g.inh.mod,"inh_mod.graphml",format = "graphml")

