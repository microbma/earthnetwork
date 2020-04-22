# Change of lamda

lamda_r <- NULL
for(i in c(2:2925)){
        g <- read_graph(file = paste("graph_growth_random/",gfile[i],sep=""),format = "graphml")
        g <- subgraph(g,degree(g)>0)
        deg.initial <- degree(g)
        deg.initial <- deg.initial[deg.initial!=0]
        deg.dis <- sort(table(deg.initial+10))
        deg.dis.df <- data.frame(degree = log(as.numeric(names(deg.dis)),10),
                                 pk = log((deg.dis)/sum(deg.dis),10))
        if(dim(deg.dis.df)[1]>5){
                ab.model <- lm(pk.Freq ~ degree, 
                               data = deg.dis.df,
                               na.action = na.omit)
                lamda_r <- c(lamda_r, -summary(ab.model)$coefficients[2,1])
        }else{
                lamda_r <- c(lamda_r, NA)}
        print(i)
}

lamda.df <- data.frame(x = c(1:2924,1:2924), 
                       lamda = c(lamda[2925:2],lamda_r),
                       sam = rep(c("evo","random"),each=2924))

ggplot(lamda.df, aes(x= x,y=lamda,color=sam))+
        geom_point(size = .5)+
        scale_color_discrete_diverging(alpha = .5,
                                       name="",
                                       labels=c("Evolution","Random"),)+
        theme_bw()+xlab("t") + ylab(expression(lambda))+
        theme(axis.title = element_text(face="italic",family = "Times"),
              aspect.ratio = 1/3,
              panel.grid = element_blank())

ggsave("manuscript/lambda_subset.tiff",height = 2,width = 6,dpi = 500)

# Average separeation

separation_r <- NULL
for(i in 2:2925){
        g <- read_graph(file = paste("graph_growth_random/",gfile[i],sep=""),format = "graphml")
        g <- subgraph(g,degree(g)>0)
        apl <- average.path.length(g)
        separation_r <- c(separation_r, apl)
        print(i)
}

sep.df <- data.frame(x = c(1:2900,1:2900), 
                     separation = c(separation[2900:1],separation_r[1:2900]),
                     sam = rep(c("evo","random"),each=2900))
                     

ggplot(sep.df, aes(x= x,y=separation,color=sam))+
        scale_color_discrete_diverging(alpha = .5,
                                       name="",
                                       labels=c("Evolution","Random"))+
        #geom_rect(xmin=1400,xmax=1800,ymin=.8,ymax=1.6,fill = "grey80")+
        geom_point(size =.5)+
        theme_bw()+xlab("t")+ylab("d")+
        xlim(c(120,3000))+#ylim(c(3,4.5))+
        theme(axis.title = element_text(face="italic",family = "Times"),
              aspect.ratio = 1/3,
              panel.grid = element_blank())

ggsave("manuscript/separate.tiff",height = 2,width = 6,dpi = 500)

# clustering coefficient decays
cluster_r <- NULL
for(i in 2:2928){
        g <- read_graph(file = paste("graph_growth_random/",gfile[i],sep=""),format = "graphml")
        subset <- subgraph(g,degree(g)>0)
        apl <- transitivity(subset)
        cluster_r <- c(cluster_r, apl)
        print(i)
}

cluster.df <- data.frame(cluster = c(cluster[2870:1],cluster_r[400:2924]),
                         x = c(1:2870,400:2924),
                         sam = c(rep("evo",2870),rep("random",2525)))

ggplot(cluster.df, aes(x= x,y=cluster,color=sam))+
        scale_color_discrete_diverging(alpha = .5,
                                       name="",
                                       labels=c("Evolution","Random"),)+
        geom_point(size =.5)+
        theme_bw()+xlab("t")+ylab("C")+
        xlim(c(100,3000))+ #ylim(c(.35,.5))+
        theme(axis.title = element_text(face="italic",family = "Times"),
              aspect.ratio = 1/3,
              panel.grid = element_blank())

ggsave("manuscript/cluster.tiff",height = 2,width = 6,dpi = 500)

# cluster fraction

fraction.r <- NULL
for(i in 2:2925){
        g <- read_graph(file = paste("graph_growth/",gfile[i],sep=""),format = "graphml")
        subset <- subgraph(g,degree(g)>0)
        mod <- cluster_walktrap(subset)
        apl <- max(table(mod$membership))/length(V(subset))
        fraction.r <- c(fraction.r, apl)
        print(i)
}

frac.df <- data.frame(x = 1:2925, 
                      kd = fraction[2925:1])

ggplot(frac.df, aes(x= x, y=kd, color = level1))+
        geom_point(color = pal_aaas(alpha = .3)(5)[1])+
        theme_bw()+xlab("t")+ylab("r")+
        xlim(c(10,3000))+ #ylim(c(0,.4))+
        theme(axis.title = element_text(face="italic"),
              aspect.ratio = 1/3,
              panel.grid = element_blank())

ggsave("manuscript/fraction.tiff",height = 2,width = 6,dpi = 500)


