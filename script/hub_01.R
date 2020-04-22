# Find the hubs in the network
g.hub <- conet
cen.bet.meta <- centr_clo(g.hub,mode = "total")[[1]]
names(cen.bet.meta) <- V(g.hub)$name
#cen.bet.meta.norm <- (cen.bet.meta)/max(cen.bet.meta)

# Local network

gh.file <- list.files("graph_1/")

cen.bet.df <- data.frame(local=NULL,
                         metacommunity=NULL,
                         env=NULL,
                         name = NULL)
for(i in 1:14){
        local.gh <- read_graph(file = paste("graph_1/", 
                                            gh.file[i],sep=""),
                               format = "graphml")
        cen.bet.local <- centr_clo(local.gh,
                                    normalized = FALSE)[[1]]
        names(cen.bet.local) <- V(local.gh)$name
        df.tem <- data.frame(local = cen.bet.local, 
                             metacommunity = cen.bet.meta[V(local.gh)$name],
                             env = gsub(".graphml","",gh.file[i]),
                             name = V(local.gh)$name)
        cen.bet.df <- rbind(cen.bet.df, df.tem)
}

hub <- rep("Peripheral",dim(cen.bet.df)[1])
hub[cen.bet.df$local>.4&cen.bet.df$metacommunity>.4] <- "Global hub"
hub[cen.bet.df$local<.4&cen.bet.df$metacommunity>.4] <- "Metacommunity hub"
hub[cen.bet.df$local>.4&cen.bet.df$metacommunity<.4] <- "Local hub"

cen.bet.df <- cbind(cen.bet.df[,1:4],hub)
cen.bet.df.r <- cen.bet.df[cen.bet.df$name %in% V(g.hub)$name,]
#cen.bet.df.r <- cen.bet.df.r[cen.bet.df.r$hub != "Global hub",]

ggplot(cen.bet.df.r,aes(x=local,
                      y=metacommunity,
                      color=hub))+
        geom_point(size=2)+
        facet_wrap(~env,nrow = 7,scales = "free")+
        theme_bw()+xlab("Local betweenness")+
        ylab("Metacommunity betweenness")+
        #geom_vline(xintercept = .4,linetype=2,color = pal_npg()(4)[4])+
        #geom_hline(yintercept = .4,linetype=2,color = pal_npg()(4)[4])+
        scale_color_manual(values = c(pal_tron()(5)[c(1,4,3)],grey(.5,alpha = .8)),name="Hub type")

ggsave("manuscript/hub_local.tiff",height = 10,width = 7,dpi = 500)


# statistic test
saline <- grep("(saline)", cen.bet.df.r$env)
nonsaline <- grep("(nonsaline)", cen.bet.df.r$env)

t.test(cen.bet.df.r$local[saline],cen.bet.df.r$local[nonsaline])

animal <- grep("gut", cen.bet.df.r$env)
plant <- grep("plant", cen.bet.df.r$env)

wilcox.test(cen.bet.df.r$metacommunity[animal],cen.bet.df.r$metacommunity[plant])


hub.df <- cen.bet.df.r[cen.bet.df.r$hub!="Peripheral",]
hub.tax <- otu.df$taxa.df[as.character(hub.df$name),]
hub.table <- data.frame(Env.type = hub.df$env,
                        Hub.type = hub.df$hub,
                        Local.betweenness = hub.df$local,
                        Meta.betweenness = hub.df$metacommunity,
                        OTU = hub.df$name,
                        hub.tax
                        ) 
hub.table <- hub.table[hub.table$OTU %in% V(g.hub)$name, ]

write.csv(x = hub.table,file = "manuscript/hub_table.csv")

g.hub <- hub.table %>%
             filter(Hub.type == "Global hub")

# Taxa profile

meta.hub <- hub.table[hub.table$Hub.type != "Local hub",]

meta.hub.tb <- sort(table(meta.hub$Genus))
meta.hub.tb <- meta.hub.tb[meta.hub.tb!=0]
pie(meta.hub.tb, cex=.5, 
    col = hcl_palettes(n = 10))

# Network Graph
row.names(layout.df) <- layout.df$name
ly.new <- as.matrix(layout.df[V(g.hub)$name,1:2]) # the layout for current graph

V(g.hub)$hub <- "Peripheral"
hub1.node.df <- cen.bet.df[grep(pattern = "Metacommunity hub", hub),]
V(g.hub)$hub[V(g.hub)$name %in% as.character(hub1.node.df$name)] <- "Metacommunity hub"
hub2.node.df <- cen.bet.df[grep(pattern = "Local hub", hub),]
V(g.hub)$hub[V(g.hub)$name %in% as.character(hub2.node.df$name)] <- "Local hub"
hub3.node.df <- cen.bet.df[grep(pattern = "Global hub",hub),]
V(g.hub)$hub[V(g.hub)$name %in% as.character(hub3.node.df$name)] <- "Global hub"
shape <- rep("circle",2928)

shape.meta <- shape
shape.meta[V(g.hub)$hub == "Metacommunity hub"] <- "sphere"

shape.local <- shape
shape.local[V(g.hub)$hub == "Local hub"] <- "sphere"

shape.global <- shape
shape.global[V(g.hub)$hub == "Global hub"] <- "sphere"

shape.all <- shape
shape.all[V(g.hub)$hub != "Peripheral"] <- "sphere"

E(g.hub)$color <- grey(.8,alpha = .01)
g1.edgelist <- get.edgelist(g.hub)
hub.meta <- as.character(unique(cen.bet.df.r[cen.bet.df.r$hub == "Metacommunity hub","name"]))
E(g.hub)$color[g1.edgelist[,1] %in%  hub.meta|g1.edgelist[,2] %in%  hub.meta] <- pal_tron(alpha=.4)(5)[c(4)]
hub.local <- as.character(unique(cen.bet.df.r[cen.bet.df.r$hub == "Local hub","name"]))
E(g.hub)$color[g1.edgelist[,1] %in%  hub.local|g1.edgelist[,2] %in%  hub.local] <- pal_tron(alpha=.4)(5)[c(3)]
hub.global <- as.character(unique(cen.bet.df.r[cen.bet.df.r$hub == "Global hub","name"]))
E(g.hub)$color[g1.edgelist[,1] %in%  hub.global|g1.edgelist[,2] %in%  hub.global] <- pal_tron(alpha=.4)(5)[c(1)]

hub.meta <- hub.meta[!hub.meta %in% hub.global]
E.color.meta <- rep(grey(.8,alpha = .01),54299)
E.color.meta[g1.edgelist[,1] %in%  hub.meta|g1.edgelist[,2] %in%  hub.meta] <- pal_tron(alpha=.4)(5)[c(4)]
E.color.global <- rep(grey(.8,alpha = .01),54299)
E.color.global[g1.edgelist[,1] %in%  hub.global|g1.edgelist[,2] %in%  hub.global] <- pal_tron(alpha=.4)(5)[c(1)]
E.color.local <- rep(grey(.8,alpha = .01),54299)
E.color.local[g1.edgelist[,1] %in%  hub.local|g1.edgelist[,2] %in%  hub.local] <- pal_tron(alpha=.4)(5)[c(3)]


tiff("manuscript/network_hub_meta.tiff",
     height = 2500,
     width = 2500,
     res = 400,
     compression = "lzw",
     type = "cairo")
par(bg="black",mar=rep(0,4))
plot(g.hub,
     layout=ly.new,
     vertex.shape = shape.meta,
     vertex.label = NA,
     vertex.size=log(degree(g1)),
     vertex.color=c(grey(.5,alpha = .2),pal_tron(alpha=.8)(5)[c(4)])[factor(shape.meta)],
     vertex.frame.color=NA,
     edge.color= E.color.meta,   
     edge.width=.3,
     edge.curved = TRUE)
dev.off()

plot(g.hub,
     layout=ly.new,
     vertex.shape = shape.global,
     vertex.label = NA,
     vertex.size=log(degree(g1)),
     vertex.color=c(grey(.5,alpha = .2),pal_tron(alpha=.8)(5)[c(1)])[factor(shape.global)],
     vertex.frame.color=NA,
     edge.color= E.color.global,   
     edge.width=.3,
     edge.curved = TRUE)

plot(g.hub,
     layout=ly.new,
     vertex.shape = shape.local,
     vertex.label = NA,
     vertex.size=log(degree(g1)),
     vertex.color=c(grey(.5,alpha = .2),pal_tron(alpha=.8)(5)[c(3)])[factor(shape.local)],
     vertex.frame.color=NA,
     edge.color= E.color.local,   
     edge.width=.3,
     edge.curved = TRUE)

## Seperate environments
env.name <- gsub(".graphml","",gh.file)
local <- data.frame(V(g.hub)$s1,
                    V(g.hub)$s2,
                    V(g.hub)$s3,
                    V(g.hub)$s4,
                    V(g.hub)$s5,
                    V(g.hub)$s6,
                    V(g.hub)$s7,
                    V(g.hub)$s8,
                    V(g.hub)$s9,
                    V(g.hub)$s10,
                    V(g.hub)$s11,
                    V(g.hub)$s12,
                    V(g.hub)$s13,
                    V(g.hub)$s14)
for(i in 1:14){
        local.hub <- unique(hub.df[hub.df$hub == "Local hub" & hub.df$env == env.name[i], 4])
        meta.hub <- unique(hub.df[hub.df$hub == "Metacommunity hub" & hub.df$env == env.name[i], 4])
        global.hub <- unique(hub.df[hub.df$hub == "Global hub" & hub.df$env == env.name[i], 4])
        
        v.shape <- local[,i]
        names(v.shape) <- V(g.hub)$name
        v.shape[V(g.hub)$name[] %in% local.hub] <- 2
        v.shape[V(g.hub)$name[] %in% meta.hub] <- 3
        v.shape[V(g.hub)$name[] %in% global.hub] <- 4
        
        E(g.hub)$color <- grey(.8,alpha = .01)
        E(g.hub)$color[g1.edgelist[,1] %in%  local.hub|g1.edgelist[,2] %in%  local.hub] <- pal_tron(alpha=.6)(5)[c(3)]
        E(g.hub)$color[g1.edgelist[,1] %in%  meta.hub|g1.edgelist[,2] %in%  meta.hub] <- pal_tron(alpha=.6)(5)[c(4)]
        E(g.hub)$color[g1.edgelist[,1] %in%  global.hub|g1.edgelist[,2] %in%  global.hub] <- pal_tron(alpha=.6)(5)[c(1)]
        
        tiff(paste("manuscript/local_hub/",env.name[i],".tiff",sep=""),
             height = 2500,
             width = 2500,
             res = 400,
             compression = "lzw",
             type = "cairo")
        par(bg="black",mar=rep(0,4))
        plot(g.hub,
             layout=ly.new,
             vertex.shape = c("circle","circle","sphere","sphere","sphere")[v.shape+1],
             vertex.label = NA,
             vertex.size = log(degree(g1))/2,
             vertex.color =c(grey(.7,alpha = .3),pal_tron(alpha=.3)(5)[c(2,3,4,1)])[v.shape+1],
             vertex.frame.color=NA,
             #edge.color = NA,   
             edge.width =.5,
             edge.curved = TRUE)
        dev.off()
        print(i)
}

# compare betweenness centrility among different environmental types

## local networks
require(tidyverse)
require(ggplot2)
require(ggbeeswarm)
require(colorspace)

mean.ord <- cen.bet.df %>%
    group_by(env) %>%
    filter(local != 0) %>%
    summarise(mean = mean(log10(local)))
mean.ord1 <- order(mean.ord$mean,decreasing = TRUE)

mean.ord2 <- cen.bet.df %>%
    group_by(env) %>%
    filter(metacommunity != 0) %>%
    summarise(mean = quantile(log10(metacommunity))[3])
mean.ord3 <- order(mean.ord2$mean,decreasing = TRUE)

local.cen.aov <- aov(local~env,data = cen.bet.df)
local.cen.aov.tuk <- TukeyHSD(local.cen.aov,ordered = TRUE,conf.level = 0.001)

cen.bet.df$env1 <- factor(cen.bet.df$env, 
                          levels = levels(cen.bet.df$env)[mean.ord1[c(c(1,4,6,7,8,9,10,12,2,3,5,11,13,14))]])
                         

ggplot(cen.bet.df, aes(x = env1, y = local))+
    geom_quasirandom(size=.6,width = .3,aes(color=env1))+
    geom_violin(width = 1,draw_quantiles = c(.5),size=.3,fill=NA)+
    scale_y_log10(breaks = c(0.000001,0.0001,0.01,1),
                  labels=expression(10^-6,10^-4,10^-2,10^0))+
    scale_color_discrete_qualitative(alpha = .6)+
    theme_bw()+
    ylab("Local betweenness centrility")+
    xlab("Environmental type")+
    theme(panel.grid = element_blank(),
          axis.text.x = element_text(angle=-70,hjust = 0)) +
    scale_x_discrete(labels =  env.name1[mean.ord1[c(1,4,6,7,8,9,10,12,2,3,5,11,13,14)]]) +
    guides(color=FALSE)+geom_vline(xintercept = 8.5)
ggsave(filename = "local_bet_cen.pdf",height = 4,width = 6)

## Meta-network

cen.bet.df$env2 <- factor(cen.bet.df$env, 
                          levels = levels(cen.bet.df$env)[mean.ord3[c(1,4,6,8,9,10,11,14,2:3,5,7,12,13)]])

meta.cen.aov <- aov(metacommunity~env,data = cen.bet.df)
meta.cen.aov.tuk <- TukeyHSD(meta.cen.aov)

ggplot(cen.bet.df, aes(x = env2, y = metacommunity))+
    geom_quasirandom(size=.6,width = .3,aes(color=env2))+
    geom_violin(width = 1,draw_quantiles = c(.5),size=.3,fill=NA)+
    scale_y_log10(breaks = c(0.000001,0.0001,0.01,1),labels=expression(10^-6,10^-4,10^-2,10^0))+
    scale_color_discrete_qualitative(alpha = .6)+
    theme_bw()+
    ylab("metacommunity betweenness centrility")+
    xlab("Environmental type")+
    theme(panel.grid = element_blank(),
          plot.margin = unit(x=c(5,3,1,3),units = "mm"),
          axis.text.x = element_text(angle=-70,hjust = 0)) +
    scale_x_discrete(labels =  env.name1[mean.ord3[c(1,4,6,8,9,10,11,14,2,3,5,7,12,13)]]) +
    guides(color=FALSE) + geom_vline(xintercept = 8.5, linetype=2)
ggsave(filename = "meta_bet_cen.pdf",height = 4,width = 6)
    

