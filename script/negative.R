# Negative edge

edgelist.trim.neg.df <- matrix(0, nrow = 14201, ncol = 12)
for(i in 1:12){
        g.cor <- read_graph(paste("trim_otu/trim_graph/",
                                  gsub("csv","graphml",otutb.name[c(2:5,7:14)][i]),
                                  sep = ""),format = "graphml")
        pos <- 1:14201
        pos1 <- pos[c(edgelist.unique %in% edgelist.trim.p[[i]])][E(g.cor)$dir <0]
        edgelist.trim.neg.df[pos1, i] <- 1
}

neg.num <- NULL
for(i in 1:12){
        edgelist.class[edgelist.trim.neg.df[,i] == 1] -> a
        neg.num <- c(neg.num, length(a))
        print(length(a))
}

sort(round(neg.num/topo_trim[,2],3))

neg.df1 <- data.frame(neg = round(neg.num/topo_trim[,2],3),
                      neg1 = neg.num,
                      env = names(topo_trim[,2]))

neg.df1$env <- factor(neg.df1$env,levels = neg.df1$env[order(neg.df1$neg)])

set.seed(1234)
ggplot(neg.df1[order(neg.df1$neg),], aes(x = env, y = neg, fill=env))+
        geom_bar(stat = "identity")+
        geom_text(aes(y=c(.2,rep(.3,5),neg[7:12]+.1),
                      label = paste(neg*100,"%\n(",neg1,")",sep="")),
                  lineheight=.8,
                  angle = c(-30*1:3+15,
                            -30*4:6-165,
                            -30*7:9-165,
                            -30*10:12+15),
                  fontface = 2)+
        geom_text(aes(label = gsub(" ","\n",env),
                      y=0),
                  hjust=0,lineheight=.8,
                  size=3,
                  angle = c(-30*1:12+105),
                  fontface = 1)+
        geom_text(x=0,y=-.3,label="Negative\nedges",
                  size=4,fontface=2,lineheight=.8,)+
        scale_fill_manual(values = morandi.color[sample(1:36,12)])+
        ylim(c(-.3,.6))+
        #geom_hline(yintercept = c(.1),linetype=2,color="red")+
        coord_polar()+guides(fill=FALSE)+
        theme(axis.line = element_blank(),
              axis.text = element_blank(),
              axis.ticks = element_blank(),
              axis.title = element_blank(),
              panel.background = element_blank())

ggsave("neg_rose.pdf",height = 8,width = 6)

# negative edge list

neg.list <- NULL
for(i in 1:12){
        a <- sort(table(edgelist.genus[edgelist.trim.neg.df[,i] == 1]),decreasing = TRUE)[1:10]
        b <- data.frame(edge = a,env = env.name1[-c(1,6)][i])
        neg.list <- rbind(neg.list, b)
}

neg.class <- NULL
for(i in 1:12){
        a <- edgelist.class[edgelist.trim.neg.df[,i] == 1]
        b <- sort(table(unlist(base::strsplit(a,"-"))),decreasing = TRUE)
        neg.class <- c(neg.class, names(b[1:9]))}

neg.class.level <- c(unique(neg.class),"Others")


ly.mat <- matrix(0,nrow = 6,ncol = 4)
env.name2 <- factor(env.name1[-c(1,6)],levels = neg.df1$env)

#ly.mat[1,] <- c(11,12,1,2)
#ly.mat[2:4,4] <- 3:5
#ly.mat[4,1:3] <- 8:6
#ly.mat[2:3,1] <- 10:9

ly.mat[2:5,1] <- 11:8
ly.mat[2:5,4] <- 2:5
ly.mat[1,2:3] <- c(12,1)
ly.mat[6,2:3] <- c(7,6)

pdf("neg_class.pdf",width = 16,height = 8)
layout(ly.mat)
par(mar=c(0,4,0,4))
set.seed(1234)
for(i in 1:12){
        a <- edgelist.class[edgelist.trim.neg.df[,i] == 1]
        b <- sort(table(unlist(base::strsplit(a,"-"))),decreasing = TRUE)
        d <- c(b[1:10], Others = sum(b[-(1:10)]))
        set.seed(1234)
        pie(d,border = FALSE, cex=1.2,
            radius = .9, init.angle = -45,
            col = morandi.color[sample(x = 11:35,25)][c(factor(names(d),
                                                               levels = neg.class.level))])
        draw.circle(0,0,radius = .5,col = "white",border = F)
        set.seed(1234)
        text(0,0,label=gsub(" ","\n",neg.df1$env[i]),
             font = 2,cex=1.2,
             col = "black")#morandi.color[sample(1:36,12)][i])
}
dev.off()

## Negative V.S. dissimilarity
### PD values
require(PhyloMeasures)
tree.trim <- read.tree("tree/trim_tree.tre")
pd.trim <- NULL
for(i in c(2:5,7:14)){
        ot <- read.csv(paste("trim_otu/",otutb.name[i],sep=""),row.names = 1)
        tr <- keep.tip(tree.trim,row.names(ot))
        pd.v <- pd.query(tree = tr,t(ot))
        pd.trim <- c(pd.trim, max(pd.v))
}

neg.per.trim <- NULL
for(i in c(2:5,7:14)){
        g.cor <- read_graph(paste("trim_otu/trim_graph/",
                                  gsub("csv","graphml",otutb.name[i]),
                                  sep = ""),format = "graphml")
        ot <- read.csv(paste("trim_otu/",otutb.name[i],sep=""),row.names = 1)
        neg.per <- NULL
        for (j in 1:360){
                gsub.cor <- subgraph(graph = g.cor, v = ot[,j] != 0)
                neg.pro <- table(E(gsub.cor)$dir < 0)
                neg.per <- c(neg.per,neg.pro[2]/sum(neg.pro))
        }
        neg.per.trim <- c(neg.per.trim, neg.per)
        print(i)
}

neg.pd <- data.frame(neg = neg.num/topo_trim[,2], # neg.per.trim, 
                     pd = pd.trim,
                     env = rep(row.names(topo_trim),each = 1))

ggplot(neg.pd, aes(y = neg*100,x = pd, color = env))+
        geom_point(size = 3,alpha=.6)+ #color = rainbow(1,alpha = .5),size=3)+
        #geom_smooth(method = "lm",color="black",linetype = 2)+
        geom_text_repel(aes(label = env))+
        ylab("Proportion of negative edge (%)")+
        xlab("Unrooted phylogenetic diversity")+
        guides(color=FALSE)+
        theme_bw()+ #ylim(c(0,2e+05))+#xlim(57, 62)+
        theme(panel.grid = element_blank(),aspect.ratio = 1)

ggsave("neg_pd.pdf",width = 4,height = 4)
