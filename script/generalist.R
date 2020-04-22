# Specific edges

edgelist.unique <- unique(unlist(edgelist.trim.p)) # find unique edges
edgelist.trim.df <- matrix(0, nrow = 14201, ncol = 12) # make a blank matrix

# fill the matrix showing edge presentation
for(i in 1:12){
        edgelist.trim.df[edgelist.unique %in% edgelist.trim.p[[i]],i] <- 1
}

# two ESV names
left <- unlist(base::strsplit(edgelist.unique,split = "-"))[seq(1,14201*2,2)]
right <- unlist(base::strsplit(edgelist.unique,split = "-"))[seq(2,14201*2,2)]

# ESVs presenting in more than 1 environments
vertex.unique.unspec <- vertex.unique[rowSums(vertex.trim.df) > 1]

# the presenting times of the EVSs combination in a edge 
edge.combine.present <- NULL
for (i in 1:14201) {
        ver.comb <- unlist(strsplit(edgelist.unique[i],"-"))
        ver.comb.df <- colSums(vertex.trim.df[vertex.unique %in% ver.comb,])
        edge.combine.present <- c(edge.combine.present,
                                  length(ver.comb.df[ver.comb.df == 2]))
}

# edges with class names
edgelist.class <- gsub("_genera_incertae_sedis","",
                       paste(otu.df$taxa.df[left[],3],
                             otu.df$taxa.df[right[],3],
                             sep = "-"))

edgelist.genus <- gsub("_genera_incertae_sedis","",
                       paste(otu.df$taxa.df[left[],6],
                             otu.df$taxa.df[right[],6],
                             sep = "-"))

# the number of specific edge in each environments
## edge presented once but the combination of ESVs presented more than once
specific.edge.number <- NULL 
for(i in 1:12){
        edgelist.class[rowSums(edgelist.trim.df) == 1 & edgelist.trim.df[,i] == 1 & edge.combine.present > 1] -> a
        b <- sort(table(a),decreasing = TRUE)
        specific.edge.number <- c(specific.edge.number, length(a))
        print(length(a))
        #print(b[1:5])
}

sp.df <- data.frame(sp = sort(round(specific.edge.number/topo_trim[,2],3)),
                    env = env.name1[-c(1,6)][order(round(specific.edge.number/topo_trim[,2],3))])

sp.df$env <- factor(sp.df$env,sp.df$env)
sp.df2 <- sp.df %>%
        arrange(env)

set.seed(1234)
ggplot(sp.df2, aes(x = env, y = sp,fill=env))+
        geom_bar(stat = "identity")+
        geom_text(aes(y=sp+.02,label = paste(sp*100,"%",sep="")),
                  angle = c(-30*1:3+15,
                            -30*4:6-165,
                            -30*7:9-165,
                            -30*10:12+15),
                  fontface = 2)+
        geom_text(aes(label = gsub(" ","\n",env),y=.05),hjust=.5,
                  size=3,lineheight=.8,
                  angle = c(-30*1:6+105,-30*7:12-75),
                  fontface = 2)+
        geom_text(x=0,y=-.1,label="Specialist edges\n(general vertex pair)",
                  size=4,fontface=2,lineheight=.8,)+
        scale_fill_manual(values = qualitative_hcl(12)[sample(12)])+
        ylim(c(-.1,.32))+
        geom_hline(yintercept = .2,linetype=2,color="black",size=1)+
        coord_polar()+guides(fill=FALSE)+
        theme(axis.line = element_blank(),
              axis.text = element_blank(),
              axis.ticks = element_blank(),
              axis.title = element_blank(),
              panel.background = element_blank())
ggsave("sp_rose.pdf",width = 8,height = 8)

# specialist edge special vertex

specific.edge1.number <- NULL 
for(i in 1:12){
        a <- edgelist.class[rowSums(edgelist.trim.df) == 1 & edgelist.trim.df[,i] == 1 & edge.combine.present == 1]
        b <- sort(table(a),decreasing = TRUE)
        specific.edge1.number <- c(specific.edge1.number, length(a))
        print(length(a))
}

sp.df1 <- data.frame(sp = sort(round(specific.edge1.number/topo_trim[,2],3)),
                     env = env.name1[-c(1,6)][order(round(specific.edge1.number/topo_trim[,2],3))])

sp.df3$env <- factor(sp.df1$env,sp.df1$env)
sp.df3 <- sp.df1 %>%
        arrange(env)

set.seed(1234)
ggplot(sp.df3, aes(x = env, y = sp,fill=env))+
        geom_bar(stat = "identity")+
        geom_text(aes(y=sp+.02,label = paste(sp*100,"%",sep="")),
                  angle = c(-30*1:3+15,
                            -30*4:6-165,
                            -30*7:9-165,
                            -30*10:12+15),
                  fontface = 2)+
        geom_text(aes(label = gsub(" ","\n",env),y=.075),hjust=.5,
                  size=3,lineheight=.8,
                  angle = c(-30*1:6+105,-30*7:12-75),
                  fontface = 2)+
        geom_text(x=0,y=-.2,label="Specialist edges\n(special vertex pair)",
                  size=4,fontface=2,lineheight=.8,)+
        scale_fill_manual(values = qualitative_hcl(12)[sample(12)])+
        ylim(c(-.2,.6))+
        geom_hline(yintercept = .2,linetype=2,color="black",size=1)+
        coord_polar()+guides(fill=FALSE)+
        theme(axis.line = element_blank(),
              axis.text = element_blank(),
              axis.ticks = element_blank(),
              axis.title = element_blank(),
              panel.background = element_blank())
ggsave("sp_rose1.pdf",width = 8,height = 8)

# edge table
for(i in 1:12){
        edgelist.class[rowSums(edgelist.trim.df) == 1 & edgelist.trim.df[,i] == 1 & edge.combine.present > 1] -> a
        print(sort(table(a),decreasing = TRUE)[1:10])
}

# class table
class.df <- NULL
a <- edgelist.class        
c <- base::strsplit(x = as.character(a), split = "-")
d <- sort(table(unlist(c)),decreasing = TRUE)
#d1 <- c(d[1:9], Others = sum(d[-c(1:9)]))
class.df <- c(class.df, names(d))

class.level <- c(unique(class.df),"Others")

pdf(file = "sp_edge.pdf",height = 8,width = 6)
par(mfcol=c(6,2),mar=rep(0,4))
for(i in 1:12){
        a <- edgelist.class[rowSums(edgelist.trim.df) == 1 & edgelist.trim.df[,i] == 1 & edge.combine.present > 1]
        c <- base::strsplit(x = as.character(a), split = "-")
        d <- sort(table(unlist(c)),decreasing = TRUE)
        d1 <- c(d[1:9], Others = sum(d[-c(1:9)]))
        set.seed(1234)
        pie(d1,border = .2,
            col = qualitative_hcl(n=71,
                                  palette = "Dark3",
                                  alpha = .9)[sample(71)][c(factor(names(d1),
                                                                   levels = class.level))],
            radius = .6,
            init.angle = -45)
        draw.circle(0,0,radius = .3,col = "white",border = F)
        text(0,0,label=gsub(" ","\n",env.name1[-c(1,6)])[i],font = 2)
}
dev.off()

pdf(file = "sp_edge1.pdf",height = 12,width = 8)
par(mfcol=c(6,2),mar=rep(0,4))
for(i in 1:12){
        edgelist.class[rowSums(edgelist.trim.df) == 1 & edgelist.trim.df[,i] == 1 & edge.combine.present == 1] -> a
        c <- base::strsplit(x = as.character(a), split = "-")
        d <- sort(table(unlist(c)),decreasing = TRUE)
        d1 <- c(d[1:9], Others = sum(d[-c(1:9)]))
        set.seed(1234)
        pie(d1,border = .2,
            col = qualitative_hcl(n=71,
                                  palette = "Dark3",
                                  alpha = .9)[sample(71)][c(factor(names(d1),
                                                                   levels = class.level))],
            radius = .6,init.angle = -45)
        draw.circle(0,0,radius = .3,col = "white",border = F)
        text(0,0,label=gsub(" ","\n",env.name1[-c(1,6)][i]),font=2)
}
dev.off()

sp.edge.list <- NULL
for(i in 1:12){
        edgelist.genus[rowSums(edgelist.trim.df) == 1 & edgelist.trim.df[,i] == 1 & edge.combine.present > 1] -> a
        ge.list <- data.frame(edge = sort(table(a),decreasing = TRUE), 
                              environment=env.name1[-c(1,6)][i])
        sp.edge.list <- rbind(sp.edge.list, ge.list[ge.list[,2]>0,])
}

a <- sp.edge.list[grep("Microgenomates",sp.edge.list$edge.a),]
b <- sp.edge.list[grep("Armatimonadetes",sp.edge.list$edge.a),]

sort(tapply(sp.edge.list$edge.Freq, sp.edge.list$edge.a, sum),decreasing = TRUE)[1:5]

sp1.edge.list <- NULL
for(i in 1:12){
        edgelist.genus[rowSums(edgelist.trim.df) == 1 & edgelist.trim.df[,i] == 1 & edge.combine.present == 1] -> a
        ge.list <- data.frame(edge = sort(table(a),decreasing = TRUE), 
                              environment=env.name1[-c(1,6)][i])
        sp1.edge.list <- rbind(sp1.edge.list, ge.list[ge.list[,2]>1,])
}


# Generalist
generalist.edge.number <- NULL
for(i in 1:12){
        edgelist.class[rowSums(edgelist.trim.df) > 1 & edgelist.trim.df[,i] == 1] -> a
        generalist.edge.number <- c(generalist.edge.number, length(a))
        print(length(a))
}

generalist.edge.list <- NULL
for(i in 1:12){
        edgelist.genus[rowSums(edgelist.trim.df) > 1 & edgelist.trim.df[,i] == 1] -> a
        ge.list <- data.frame(edge = sort(table(a),decreasing = TRUE), 
                              environment=env.name1[-c(1,6)][i])
        generalist.edge.list <- rbind(generalist.edge.list, ge.list[ge.list[,2]>1,])
}


ge.df <- data.frame(ge = sort(round(generalist.edge.number/topo_trim[,2],3)),
                    env = env.name1[-c(1,6)][order(round(generalist.edge.number/topo_trim[,2],3))])

ge.df$env <- factor(ge.df$env, levels = env.name1[-c(1,6)])

ge.df1 <- ge.df %>%
        arrange(env)

set.seed(1234)
ggplot(ge.df1, aes(x = env, y = ge, fill=env))+
        geom_bar(stat = "identity")+
        geom_text(aes(y=ge+.02,label = paste(ge*100,"%",sep="")),
                  angle = c(-30*1:3+15,
                            -30*4:6-165,
                            -30*7:9-165,
                            -30*10:12+15),
                  fontface = 2)+
        geom_text(aes(label = gsub(" ","\n",env),y=.15),hjust=.5,
                  size=3,lineheight=.8,
                  angle = c(-30*1:6+105,-30*7:12-75),
                  fontface = 2)+
        geom_text(x=0,y=-.2,label="Generalist edge",size=4,fontface=2)+
        scale_fill_manual(values = qualitative_hcl(12)[sample(12)])+
        ylim(c(-.2,.7))+
        geom_hline(yintercept = .4,linetype=2,color="black",size=1)+
        coord_polar()+guides(fill=FALSE)+
        theme(axis.line = element_blank(),
              axis.text = element_blank(),
              axis.ticks = element_blank(),
              axis.title = element_blank(),
              panel.background = element_blank())
ggsave("ge_rose.pdf",width = 8,height = 8)

all.edge.df <- data.frame(edge = c(ge.df1$ge,sp.df3$sp,sp.df2$sp),
                          type = rep(c("g","sp1","sp2"),each=12),
                          env = c(as.character(ge.df1$env),
                                  as.character(sp.df3$env),
                                  as.character(sp.df2$env)))

all.edge.df$env <- factor(all.edge.df$env, 
                          levels = ge.df$env)

ggplot(all.edge.df, aes(x = env, y = edge*100, fill=type))+
        geom_bar(stat = "identity")+
        geom_hline(yintercept = 50,
                   linetype=2,
                   color=grey(level = .1,
                              alpha = .5))+
        scale_fill_manual(values = qualitative_hcl(3,alpha = .9),
                          name="Edge type",label=c("Generalist edge",
                                                   "Specialist edges (specialist vertex pair)",
                                                   "Specialist edges (generalist vertex pair)"))+
        theme_bw()+xlab("")+ylab("%")+
        theme(axis.text.x = element_text(angle = -45,hjust = 0),
              panel.grid = element_blank())
ggsave("edge_percent.pdf",width = 8,height = 3)

class.df1 <- NULL
for(i in 1:12){
        edgelist.class[rowSums(edgelist.trim.df) > 1 & edgelist.trim.df[,i] == 1] -> a
        c <- base::strsplit(x = as.character(a), split = "-")
        d <- sort(table(unlist(c)),decreasing = TRUE)
        d1 <- c(d[1:9], Others = sum(d[-c(1:9)]))
        class.df1 <- c(class.df1, names(d1))
}

pdf(file = "ge_edge.pdf",height = 12,width = 8)
par(mfcol=c(6,2),mar=rep(0,4))
for(i in 1:12){
        edgelist.class[rowSums(edgelist.trim.df) > 1 & edgelist.trim.df[,i] == 1] -> a
        c <- base::strsplit(x = as.character(a), split = "-")
        d <- sort(table(unlist(c)),decreasing = TRUE)
        d1 <- c(d[1:19], Others = sum(d[-c(1:19)]))
        set.seed(1234)
        pie(d1,border = .2,
            col = qualitative_hcl(n=71,
                                  palette = "Dark3",
                                  alpha = .9)[sample(71)][c(factor(names(d1),
                                                                   levels = class.level))],
            radius = .6,init.angle = -45)
        draw.circle(0,0,radius = .3,col = "white",border = F)
        text(0,0,label=gsub(" ","\n",env.name1[-c(1,6)])[i],font = 2)
}
dev.off()

par(mfcol=c(3,1),mar=c(0,0,0,0))
a <- edgelist.genus[rowSums(edgelist.trim.df) > 1]
c <- base::strsplit(x = as.character(a), split = "-")
d <- sort(table(unlist(c)),decreasing = TRUE)
d1 <- c(d[1:50], Others = sum(d[-c(1:50)]))
set.seed(1234)
pie(d1,border = .2,cex=.6,
    col = qualitative_hcl(31),
    radius = .6,init.angle = -45)
draw.circle(0,0,radius = .3,col = "white",border = F)
text(x=0,y=0,label="Generalist\nedge",font=2)

a <- edgelist.genus[rowSums(edgelist.trim.df) == 1 & edge.combine.present > 1]
c <- base::strsplit(x = as.character(a), split = "-")
d <- sort(table(unlist(c)),decreasing = TRUE)
d2 <- c(d[1:50], Others = sum(d[-c(1:50)]))
set.seed(1234)
pie(d2,border = .2,
    col = qualitative_hcl(31),
    radius = .6,init.angle = -45)
draw.circle(0,0,radius = .3,col = "white",border = F)
text(x=0,y=0,label="Specialist edge\ngroup 1",font=2)

a <- edgelist.genus[rowSums(edgelist.trim.df) == 1 & edge.combine.present == 1]
c <- base::strsplit(x = as.character(a), split = "-")
d <- sort(table(unlist(c)),decreasing = TRUE)
d3 <- c(d[1:50], Others = sum(d[-c(1:50)]))
set.seed(1234)
pie(d3,border = FALSE,
    col = qualitative_hcl(51)[sample(51)],
    radius = .6,init.angle = -45)
draw.circle(0,0,radius = .3,col = "white",border = F)
text(x=0,y=0,label="Specialist edge\ngroup 2",font=2)


edge.genus <- data.frame(taxa = c(names(d1[1:50]),
                                  names(d2[1:50]),
                                  names(d3[1:50])),
                         percent = c(d1[1:50]/sum(d1),
                                     d2[1:50]/sum(d2),
                                     d3[1:50]/sum(d3)),
                         type = rep(c("g","sp1","sp2"),each=50))
edge.mean <- tapply(edge.genus$percent,INDEX = edge.genus$taxa,FUN = mean)
edge.genus$type <- factor(edge.genus$type, levels = c("sp2","sp1","g"))



edge.genus$taxa <- gsub("1","Gp1",edge.genus$taxa,fixed = TRUE)
edge.genus$taxa <- gsub("_sensu_stricto","",edge.genus$taxa,fixed = TRUE)
edge.genus$taxa <- gsub("_Incertae_Sedis","",edge.genus$taxa,fixed = TRUE)

edge.genus$taxa <- factor(edge.genus$taxa, 
                          levels = names(sort(edge.mean,decreasing = TRUE)))

ggplot(edge.genus, aes(y = type, size = percent, x= taxa, color=taxa))+
        geom_point(shape=16)+
        scale_colour_manual(values = qualitative_hcl(78,alpha = .5)[sample(78)])+
        scale_size_continuous(breaks = c(0.01,.02,.03),
                              labels = c("1%","2%","3%"),
                              name = "Relative abundance")+
        scale_y_discrete(labels = c("Generalist edge",
                                    "Specialist edge (specialist vertex pair)",
                                    "Specialist edge (generalist vertex pair)")[3:1])+
        guides(color=FALSE)+
        theme_bw()+xlab("")+ylab("")+
        theme(axis.text.x = element_text(angle = -80,hjust = 0,size=6),
              legend.position = "top")

ggsave("edge_class_genus.pdf",width = 8,height = 2.5)    

# Compare edge types
edge.compare.df <- data.frame(g = ge.df1$ge*100,
                              sp1 = sp.df3$sp*100,
                              sp2 = sp.df2$sp*100,
                              env = ge.df1$env)



edgetype.g1 <- ggplot(edge.compare.df, aes(x = g, y = sp1, color = env))+
        theme_bw()+xlab("Generalist edge (%)")+
        ylab("Specialist edge linking\ngeneralist vertex pair (%)")+
        geom_point(size = 3,alpha=.6)+
        #geom_smooth(method = "lm",aes(x = g, y = sp1),inherit.aes = FALSE,color="black")+
        geom_text_repel(aes(label = env),color="black",size=2.5)+
        scale_color_discrete_qualitative()+
        ggtitle("A")+
        guides(color= FALSE)+
        theme(aspect.ratio = 1,
              panel.grid  = element_blank())


edgetype.g2 <- ggplot(edge.compare.df, aes(x = g, y = sp2,color = env))+
        theme_bw()+xlab("Generalist edge (%)")+
        ylab("Specialist edge linking\ngeneralist vertex pair (%)")+
        geom_point(size = 3,alpha=.6)+
        geom_smooth(method = "lm")+
        geom_text_repel(aes(label = env),color="black",size=2.5)+
        scale_color_discrete_qualitative()+
        ggtitle("B")+
        guides(color= FALSE)+
        theme(aspect.ratio = 1,
              panel.grid  = element_blank())

edgetype.g3 <- ggplot(edge.compare.df, aes(x = sp1, y = sp2,color = env))+
        geom_point(size = 3,alpha=.6)+
        geom_smooth(method = "lm")+
        geom_text_repel(aes(label = env),color="black",size=2.5)+
        scale_color_discrete_qualitative()+
        ggtitle("C")+
        guides(color= FALSE)+
        theme_bw()+
        theme(aspect.ratio = 1,
              panel.grid  = element_blank())+
        xlab("Specialist edge linking\nspecialist vertex pair (%)")+
        ylab("Specialist edge linking\nspecialist vertex pair (%)")

pdf("edgetype.pdf",width = 4,height = 8)
grid.arrange(edgetype.g1,
             edgetype.g2,
             edgetype.g3)
dev.off()

edge.compare.df2 <- data.frame(all.edge.df,edge = topo_sub[-c(1,6),2])

edge.number.g1 <- ggplot(edge.compare.df2[edge.compare.df2$type=="g",], aes(x = edge.1, y = edge*100,color = env))+
        geom_point(size = 3,alpha=.6)+
        geom_smooth(method = "lm")+
        geom_text_repel(aes(label = env),color="black",size=2.5)+
        theme_bw()+xlab("Edge number")+
        ylab("Generalist edge (%)")+
        scale_color_discrete_qualitative()+
        ggtitle("A")+
        guides(color= FALSE)+
        theme(aspect.ratio = 1,
              panel.grid  = element_blank())

edge.number.g2 <- ggplot(edge.compare.df2[edge.compare.df2$type=="sp1",], aes(x = edge.1, y = edge*100,color = env))+
        geom_point(size = 3,alpha=.6)+
        geom_smooth(method = "lm")+
        geom_text_repel(aes(label = env),color="black",size=2.5)+
        theme_bw()+xlab("Edge number")+
        ylab("Specialist edge linking\ngeneralist vertex pair (%)")+
        scale_color_discrete_qualitative()+
        ggtitle("B")+
        guides(color= FALSE)+
        theme(aspect.ratio = 1,
              panel.grid  = element_blank())

edge.number.g3 <- ggplot(edge.compare.df2[edge.compare.df2$type=="sp2",], aes(x = edge.1, y = edge*100,color = env))+
        geom_point(size = 3,alpha=.6)+
        geom_smooth(method = "lm")+
        geom_text_repel(aes(label = env),color="black",size=2.5)+
        theme_bw()+xlab("Edge number")+
        ylab("Specialist edge linking\nspecialist vertex pair (%)")+
        scale_color_discrete_qualitative()+
        ggtitle("C")+
        guides(color= FALSE)+
        theme(aspect.ratio = 1,
              panel.grid  = element_blank())

pdf("edgenumber.pdf",width = 4,height = 8)
grid.arrange(edge.number.g1,
             edge.number.g2,
             edge.number.g3)
dev.off()


ggsave("edgenumebr_type.pdf",width = 4,height = 8)
