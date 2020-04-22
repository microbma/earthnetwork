# OS pattern
edgelist <- get.edgelist(g.link)
edge <- E(g.link)

# make dataset
## find node from all otu_table
otutable.name <- list.files("otutable/")
otu.tb <- list(read.csv(paste("otutable/",otutable.name[1],sep=""),
                   header = TRUE, row.names = 1),
               read.csv(paste("otutable/",otutable.name[2],sep=""),
                        header = TRUE, row.names = 1),
               read.csv(paste("otutable/",otutable.name[3],sep=""),
                        header = TRUE, row.names = 1),
               read.csv(paste("otutable/",otutable.name[4],sep=""),
                        header = TRUE, row.names = 1),
               read.csv(paste("otutable/",otutable.name[5],sep=""),
                        header = TRUE, row.names = 1),
               read.csv(paste("otutable/",otutable.name[6],sep=""),
                        header = TRUE, row.names = 1),
               read.csv(paste("otutable/",otutable.name[7],sep=""),
                        header = TRUE, row.names = 1),
               read.csv(paste("otutable/",otutable.name[8],sep=""),
                        header = TRUE, row.names = 1),
               read.csv(paste("otutable/",otutable.name[9],sep=""),
                        header = TRUE, row.names = 1),
               read.csv(paste("otutable/",otutable.name[10],sep=""),
                        header = TRUE, row.names = 1),
               read.csv(paste("otutable/",otutable.name[11],sep=""),
                        header = TRUE, row.names = 1),
               read.csv(paste("otutable/",otutable.name[12],sep=""),
                        header = TRUE, row.names = 1),
               read.csv(paste("otutable/",otutable.name[13],sep=""),
                        header = TRUE, row.names = 1),
               read.csv(paste("otutable/",otutable.name[14],sep=""),
                        header = TRUE, row.names = 1)
               )

vertex.presence <- NULL

for(k in 1:2928){
        v.name <- V(g.link)$name[k]
        presence <- NULL
        for(i in 1:14){
                presence <- c(presence, v.name %in% colnames(otu.tb[[i]]))
        }
        vertex.presence <- rbind(vertex.presence,presence)
        print(k)
        }
        
row.names(vertex.presence) <- V(g.link)$name
vertex.mat <- matrix(nrow = 2928, ncol = 14)
vertex.mat[vertex.presence] <- 1
vertex.mat[!vertex.presence] <- 0
vertex.presence.sum <- rowSums(vertex.mat)

vertex.mat.1 <- vertex.mat[rowSums(vertex.mat) == 1,]
colSums(vertex.mat.1)

vertex.mat.2 <- vertex.mat[rowSums(vertex.mat) == 2,]
pheatmap(t(vertex.mat.2))

vertex.mat.8 <- vertex.mat[rowSums(vertex.mat) == 4,]
#pheatmap(t(vertex.mat.8))
g <- graph_from_incidence_matrix(vertex.mat.8)
plot(g,vertex.size=1,layout=layout_nicely,vertex.label = NA)

# Presence of edges
edgelist.sep <- list(get.edgelist(read.graph(paste("graph_1/",graph.ls[1],sep=""),format = "graphml")),
                     get.edgelist(read.graph(paste("graph_1/",graph.ls[2],sep=""),format = "graphml")),
                     get.edgelist(read.graph(paste("graph_1/",graph.ls[3],sep=""),format = "graphml")),
                     get.edgelist(read.graph(paste("graph_1/",graph.ls[4],sep=""),format = "graphml")),
                     get.edgelist(read.graph(paste("graph_1/",graph.ls[5],sep=""),format = "graphml")),
                     get.edgelist(read.graph(paste("graph_1/",graph.ls[6],sep=""),format = "graphml")),
                     get.edgelist(read.graph(paste("graph_1/",graph.ls[7],sep=""),format = "graphml")),
                     get.edgelist(read.graph(paste("graph_1/",graph.ls[8],sep=""),format = "graphml")),
                     get.edgelist(read.graph(paste("graph_1/",graph.ls[9],sep=""),format = "graphml")),
                     get.edgelist(read.graph(paste("graph_1/",graph.ls[10],sep=""),format = "graphml")),
                     get.edgelist(read.graph(paste("graph_1/",graph.ls[11],sep=""),format = "graphml")),
                     get.edgelist(read.graph(paste("graph_1/",graph.ls[12],sep=""),format = "graphml")),
                     get.edgelist(read.graph(paste("graph_1/",graph.ls[13],sep=""),format = "graphml")),
                     get.edgelist(read.graph(paste("graph_1/",graph.ls[14],sep=""),format = "graphml"))
                     )

edgelist.sep.p <- list(paste(edgelist.sep[[1]][,1],edgelist.sep[[1]][,2],sep="-"),
                       paste(edgelist.sep[[2]][,1],edgelist.sep[[2]][,2],sep="-"),
                       paste(edgelist.sep[[3]][,1],edgelist.sep[[3]][,2],sep="-"),
                       paste(edgelist.sep[[4]][,1],edgelist.sep[[4]][,2],sep="-"),
                       paste(edgelist.sep[[5]][,1],edgelist.sep[[5]][,2],sep="-"),
                       paste(edgelist.sep[[6]][,1],edgelist.sep[[6]][,2],sep="-"),
                       paste(edgelist.sep[[7]][,1],edgelist.sep[[7]][,2],sep="-"),
                       paste(edgelist.sep[[8]][,1],edgelist.sep[[8]][,2],sep="-"),
                       paste(edgelist.sep[[9]][,1],edgelist.sep[[9]][,2],sep="-"),
                       paste(edgelist.sep[[10]][,1],edgelist.sep[[10]][,2],sep="-"),
                       paste(edgelist.sep[[11]][,1],edgelist.sep[[11]][,2],sep="-"),
                       paste(edgelist.sep[[12]][,1],edgelist.sep[[12]][,2],sep="-"),
                       paste(edgelist.sep[[13]][,1],edgelist.sep[[13]][,2],sep="-"),
                       paste(edgelist.sep[[14]][,1],edgelist.sep[[14]][,2],sep="-"))
                       

edge.presence1 <- NULL
for(i in 1:54299){
        e.name <- paste(edgelist[i,1],edgelist[i,2], sep = "-")
        presence <- NULL
        for(j in 1:14){
                presence <- c(presence, e.name %in% edgelist.sep.p[[j]])
        }
        edge.presence1 <- rbind(edge.presence1,presence)
        print(i)
}

row.names(edge.presence1) <- 1:54299
edge.mat1 <- matrix(nrow = 54299,ncol = 14)
edge.mat1[edge.presence1] <- 1
edge.mat1[!edge.presence1] <- 0
edge.presence.sum1 <- rowSums(edge.mat1)

edge.presence2 <- NULL
for(i in 1:54299){
        e.name <- paste(edgelist[i,2],edgelist[i,1], sep = "-")
        presence <- NULL
        for(j in 1:14){
                presence <- c(presence, e.name %in% edgelist.sep.p[[j]])
        }
        edge.presence2 <- rbind(edge.presence2,presence)
        print(i)
}

row.names(edge.presence2) <- 1:54299
edge.mat2 <- matrix(nrow = 54299,ncol = 14)
edge.mat2[edge.presence2] <- 1
edge.mat2[!edge.presence2] <- 0
edge.presence.sum2 <- rowSums(edge.mat2)

edge.presence.sum <- edge.presence.sum1 + edge.presence.sum2
names(edge.presence.sum) <- 1:54299



ggplot(data.frame(x=1:14,y = table(vertex.presence.sum)),
        aes(x,y.Freq,fill=factor(x)))+
        geom_bar(stat = "identity")+
        geom_text(aes(y=y.Freq+50,label=y.Freq),
                  color = getPalette(14))+
        scale_fill_manual(values = alpha(getPalette(14),.8))+
        theme_bw()+
        guides(fill=FALSE)+
        theme(panel.grid = element_blank(),
              aspect.ratio = 5/14)+
        scale_x_continuous(breaks = 1:14)+
        labs(x="Number of province",
             y="Frequncy")
ggsave("vertex_occurrence.pdf",height = 3,width = 6)

ggplot(data.frame(x=0:12,y=table(edge.presence.sum)),
        aes(x,y.Freq,fill=factor(x)))+
        geom_bar(stat = "identity")+
        geom_text(aes(y=y.Freq+800,label=y.Freq),
                  color = getPalette(13))+
        scale_fill_manual(values = alpha(getPalette(13),.8))+
        theme_bw()+
        guides(fill=FALSE)+
        theme(panel.grid = element_blank(),
              aspect.ratio = 5/11)+
        scale_x_continuous(breaks = 0:12)+
        labs(x="Number of province",
             y="Frequncy")

ggsave("edge_occurrence.pdf",height = 3,width = 6)


# omission score
## edges presented in more than one graph

edge.presence <- edge.presence1 + edge.presence2

for(n in 2:12){
multi.edge <- c(1:54299)[edge.presence.sum == n]

p <- NULL
os <- NULL
for(j in 1:length(multi.edge)){
        graph.p2 <- c(1:14)[edge.presence[multi.edge[j],]==1]
        ot.os <- NULL
        ot <- NULL
        for(k in 1:n){
        ot.os <- c(ot.os,list(otu.tb[[graph.p2[k]]][,c(edgelist[multi.edge[j],])]))
        ot <- rbind(ot,ot.os[[k]])
        }
        
        cor.org <- cor(ot[,1],ot[,2],method = "spearman")
        
        cor.os <- NULL
        for(k in 1:n){
                os1 <- lapply(ot.os, FUN = function(x){x[,1]})
                os2 <- lapply(ot.os, FUN = function(x){x[,2]})
                cor.os <- c(cor.os, 
                            cor(unlist(os1), 
                                unlist(os2),
                                method = "spearman"))
        }
        
        p.os <- NULL
        for(k in 1:n){
                p.os1 <- 0
                for(i in 1:500){
                        cor.random <- cor(ot[-sample(1:dim(ot)[1],dim(ot.os[[k]])[1]),])[1,2]/cor.os[k]
                        if (cor.random < 1){ p.os1 <- p.os1+1 }
                }
                p.os <- c(p.os, p.os1)
        }
        
        p <- rbind(p,p.os/500)
        os <- rbind(os, c(cor.os/cor.org))
}

write.csv(p,file = paste("os_test/p_",n,sep=""),row.names = TRUE,col.names = TRUE)
write.csv(os,file = paste("os_test/os_",n,sep=""),row.names = TRUE,col.names = TRUE)
print(n)
}

local <- NULL
local1 <- NULL
for(i in c(2:6,8:10)){
        p.file <- read.csv(file = paste("os_test/p_",i,sep=""))
        p.file.d <- p.file[,-1]
        p.file.d[p.file[,-1] > 0.05] <- 1
        p.file.d[p.file[,-1] <= 0.05] <- 0
        num.local <- p.file[rowSums(p.file.d)!=i,-1]
        edge.mat.local <- edge.mat[rowSums(edge.mat) == i,][rowSums(p.file.d)!=i,]
        edge.mat.local.t <- t(edge.mat.local)
        edge.mat.local.t[edge.mat.local.t==1] <- unlist(c(t(num.local)))
        edge.mat.local.t2 <- t(edge.mat.local.t)
        edge.mat.local.t2[edge.mat.local.t2 > 0.05] <- 0
        edge.mat.local.t2[edge.mat.local.t2 != 0] <- 1
        local <- rbind(local,colSums(edge.mat.local.t2))
        local1 <- c(local1, dim(edge.mat.local.t2)[1])
}

edge.mat <- edge.mat1 + edge.mat2

spe.per <- NULL
for(i in 1:14){
        dis <- table(edge.presence.sum[edge.mat[,i]==1])
        spe <- c(dis[1],sum(dis[-1]))
        spe.per <- rbind(spe.per, spe)
}

node.per <- NULL
for(i in 1:14){
        dis <- table(vertex.presence.sum[vertex.mat[,i]==1])
        spe <- c(dis[1],sum(dis[-1]))
        node.per <- rbind(node.per, spe)
}

row.names(node.per) <- env.name1
row.names(spe.per) <- env.name1
colnames(node.per) <- c("s","u")
colnames(spe.per) <- c("s","u")

node.per1 <- node.per[order(rowSums(spe.per),decreasing = TRUE),]
spe.per1 <- spe.per[order(rowSums(spe.per),decreasing = TRUE),]

require(reshape2)

spe.df <- data.frame(melt(rbind(spe.per1,node.per1)),
                     type=rep(c("Edge","Vertex"),each=14))

spe.df$type <- factor(spe.df$type,levels = c("Vertex","Edge"))

percent <- round(rbind(spe.per1,node.per1)[,1]/rowSums(rbind(spe.per1,node.per1)),2)*100

spe.100 <- data.frame(percent = paste(percent,"%",sep=""),
                      x = rownames(rbind(spe.per1,node.per1)),
                      type=rep(c("Edge","Vertex"),each=14),
                      y=rowSums(rbind(spe.per1+2000,node.per1+100)))

ggplot(spe.df,aes(x=Var1,
                         y=value,
                         fill=Var2))+
               geom_bar(stat = "identity",
                        width = .6,
                        position = "stack")+
        facet_grid(type~., scales = "free_y")+
        theme_bw()+
        scale_fill_discrete_qualitative(labels=c("Local","Global"),
                                        name="","Set 2",alpha=1)+
        xlab("Environmental type")+ylab("Count number")+
        theme(aspect.ratio = 1/3.5,plot.margin = unit(c(4,4,1,4),units = "mm"),
              legend.position = c(.8,.85),
              legend.key.size = unit(3,"mm") ,
              panel.grid = element_blank(),
              axis.text.x = element_text(angle=-65,
                                         hjust = 0,
                                         vjust = 0))+
        geom_text(aes(label = spe.100$percent,x=x,y=y),size=3,
                  data = spe.100,inherit.aes = FALSE)

ggsave("manuscript/edge_type.eps",height = 6,width = 5)

env.name1 <- c("Animal corpus", "Animal distal gut", "Animal proximal gut",
               "Animal secretion","Animal surface", "Plant corpus" ,
               "Plant surface","Rhizosphere","Sediment (nonsaline)",
               "Sediment (saline)","Soil","Surface (nonsaline)" ,
               "Water (nonsaline)","Water (saline)")

node.pie <- sort(node.per[,1],decreasing = TRUE)
node.pie1 <- c(node.pie[1:8],Others=sum(node.pie[9:14]))

spe.pie <- sort(spe.per[,1],decreasing = TRUE)
spe.pie1 <- c(spe.pie[1:8],Others=sum(spe.pie[9:14]))

require(plotrix)
set.seed(1234)
postscript("manuscript/pie.eps",height = 6,width = 8)
par(mfcol=c(1,2),mar=rep(0,3,0,3))
set.seed(1234)
pie(node.pie1,init.angle = 90,radius = .5,
    col = qualitative_hcl(14,"Dark3")[sample(1:14,14)], 
    border = NA,
    cex=.7) 
draw.circle(x=0,y=0,
            radius = .3,
            border = NA,
            col = "white")

set.seed(1234)
pie(spe.pie1,init.angle = 45,radius = .5,
    col = qualitative_hcl(14,"Dark3")[sample(1:14,14)], 
    border = NA,
    cex=.7) 
draw.circle(x=0,y=0,
            radius = .3,
            border = NA,
            col = "white")
dev.off()


## RiverPlot

library(networkD3)
library(riverplot)

nodes <- data.frame(ID = c("n101",paste("n2",1:14,sep="0"),paste("n3",1:28,sep="0")),
                    x = c(1,rep(2,14),rep(3,28)),
                    col = c("grey", qualitative_hcl(14,"Dark3"),rep(c("black","grey"),14)),
                    labels=rep("",43),
                    stringsAsFactors= FALSE
)

row.names(spe.per) <- env.name
spe.per <- spe.per[order(rowSums(spe.per),decreasing = TRUE),]

node2 <- rowSums(spe.per)
colnames(spe.per) <- c("u","g")
names(node2) <- paste("n2",1:14,sep="0")

edge <- list(n101 = as.list(node2),
             n201  = list(n301=spe.per[1,1],n302=spe.per[1,2]),
             n202  = list(n303=spe.per[2,1],n304=spe.per[2,2]),
             n203  = list(n305=spe.per[3,1],n306=spe.per[3,2]),
             n204  = list(n307=spe.per[4,1],n308=spe.per[4,2]),
             n205  = list(n309=spe.per[5,1],n3010=spe.per[5,2]),
             n206  = list(n3011=spe.per[6,1],n3012=spe.per[6,2]),
             n207  = list(n3013=spe.per[7,1],n3014=spe.per[7,2]),
             n208  = list(n3015=spe.per[8,1],n3016=spe.per[8,2]),
             n209  = list(n3017=spe.per[9,1],n3018=spe.per[9,2]),
             n2010 = list(n3019=spe.per[10,1],n3020=spe.per[10,2]),
             n2011 = list(n3021=spe.per[11,1],n3022=spe.per[11,2]),
             n2012 = list(n3023=spe.per[12,1],n3024=spe.per[12,2]),
             n2013 = list(n3025=spe.per[13,1],n3026=spe.per[13,2]),
             n2014 = list(n3027=spe.per[14,1],n3028=spe.per[14,2])
)

river <- makeRiver(nodes,edge)
riverplot(river, node_margin = .03,nsteps = 500)

# Sample number and EVS number 

sam.evs <- mapply(otu.tb, FUN = function(x) {dim(x)},USE.NAMES = TRUE, SIMPLIFY = TRUE)
sam.evs.df <- data.frame(sample = sam.evs[1,],
                         evs = sam.evs[2,],
                         env = env.name1
                         )

row.names(sam.evs.df) <- env.name1

sam.evs.df1 <- cbind(sam.evs.df[as.character(spe.df[1:14,1]),],
                     edge = spe.df[1:14,3],
                     vertex = spe.df[15:28,3],
                     all.edge = spe.df[1:14,3]+spe.df[1:14+28,3],
                     all.vertex = spe.df[15:28,3]+spe.df[15:28+28,3])

ggplot(sam.evs.df1, aes(x = sample, y = evs, color=env))+
    geom_point()+
    theme_bw()+
    labs(x="Sample numbers",y="EVS number")+
    geom_smooth(method = "glm",se=TRUE,color=grey(.1,.8),size=.5,linetype=2)+
    scale_color_d3(palette = "category20",alpha = 1)+
    theme(panel.grid = element_blank(),
          legend.position = c(.2,.8))+
    geom_text_repel(aes(label=env),size=3)+
    guides(color=FALSE)+
    scale_x_log10(breaks = c(200,500,1000,2000,5000))+scale_y_log10()

ggsave("manuscript/sample_evs.pdf",height = 4,width = 4)

ggplot(sam.evs.df1, aes(x = sample, y = all.vertex, color=env,size=vertex))+
    geom_point()+
    theme_bw()+
    labs(x="Sample numbers",y="Subnetwork vertex number")+
    geom_smooth(method = "glm",se=TRUE, color=grey(.1,.8),size=.5,linetype=2)+
    scale_color_d3(palette = "category20",alpha = 1)+
    theme(panel.grid = element_blank(),
          legend.position = c(.3,.8),legend.key.height = unit(x = 1,units = "mm"))+
    geom_text_repel(aes(label=env),size=3)+
    scale_size_continuous(breaks = c(10,50,100,200),name="Local vertex number")+
    guides(color=FALSE)+
    scale_x_log10(breaks = c(200,500,1000,2000,5000))+
    scale_y_log10(breaks = c(100,500,1000,2000))

ggsave("manuscript/sample_vertex.pdf",height = 4,width = 4)

ggplot(sam.evs.df1, aes(x = sample, y = all.edge, color=env,size=edge))+
    geom_point()+
    theme_bw()+
    labs(x="Sample numbers",y="Subnetwork edge number")+
    geom_smooth(method = "glm",se=TRUE, color=grey(.1,.8),size=.5,linetype=2)+
    scale_color_d3(palette = "category20",alpha = 1)+
    theme(panel.grid = element_blank(),
          legend.position = c(.3,.8),legend.key.height = unit(x = 1,units = "mm"))+
    geom_text_repel(aes(label=env),size=3)+
    scale_size_continuous(breaks = c(500,1000,5000),name="Local edge number")+
    guides(color=FALSE)+
    scale_x_log10(breaks = c(200,500,1000,2000,5000))+scale_y_log10()

ggsave("manuscript/sample_edge.pdf",height = 4, width = 4)


write.csv(sam.evs.df1[,c(3,1,2,4,5)],file = "manuscript/table 3.csv",row.names = FALSE,col.names = TRUE)

summary(lm(edge~sample,data = sam.evs.df1))

