# Package
require(ggplot2)
require(ggsci)
require(reshape2)

# Fitness tree growth
m1 <- m0
names(m1) <- node.df.rank$node
tip.initial.degree <- m1[as.character(tip.df$parent)]
ki <- degree(conet)[tip.df$label] 

t = node.df.rank$x
names(t) <- node.df.rank$node
t0 <- order(t[as.character(tip.df$parent)])

a <- log(c(ki+1)/c(tip.initial.degree+1))
b <- log(2929/(2929-t0))
fitness <- a/b

fit.dis <- sort(table(round(fitness[1:2928],1)))
fit.dis.df <- data.frame(bin = sort(as.numeric(names(fit.dis)),
                                    decreasing = F),
                         cum = fit.dis[order(as.numeric(names(fit.dis)),
                                             decreasing = F)])

ggplot(fit.dis.df[-c(1:2),], aes(x=bin, y = cum.Freq/2923))+
        geom_point(size=1,color = pal_aaas(alpha = .5)(6)[1])+
        theme_bw()+
        ylab(expression(italic(rho(beta))))+xlim(c(.2,30))+
        xlab(expression(beta))+
        scale_x_log10(breaks = c(.1,.5,1,5,10,50),limits = c(.1,60))+
        geom_segment(aes(x=30,y=0.001,xend=1,yend=.024),
                     size=.8,linetype = 1,
                     color = pal_aaas(alpha = .4)(6)[2])+
        geom_segment(aes(x=1,y=.024,xend=.2,yend=.1),
                     size=.8,linetype = 3,
                     color = pal_aaas(alpha = .4)(6)[2])+
        geom_vline(xintercept = 1,
                   linetype=2,
                   color=pal_npg()(6)[4])+
        theme(aspect.ratio = .66,
              panel.grid = element_blank(),
              axis.title = element_text(face = "italic"))+
        scale_y_log10(limits=c(.0005,.1),
                      breaks=c(.001,.01,.1),
                      labels=expression(10^-3,10^-2,10^-1))

ggsave("manuscript/fit_dis.tiff",height = 3,width = 4,dpi = 500)

fit.df <- data.frame(fitness=fitness[order(t0)],
                     t = 1:29)

fit.df <- data.frame(x=1-t,
                     fitness = fitness[1:2925],
                     sep = c(rep(.1,325),rep(.6,2600)))

ggplot(fit.df[,], aes(x=x,y=fitness))+
        geom_point(size=2,color=alpha(alpha = .3,pal_aaas()(6)[4]))+
        theme_bw()+
        scale_y_log10(breaks=c(.1,1,10,100,1000),limits=c(.1,2000),
                      labels=expression(10^-1,10^0,10^1,10^2,10^3))+
        xlab(expression(italic(t[age])))+
        ylab(expression(beta))+
        geom_smooth(method = "lm",color="black",linetype=1)+
        scale_alpha_continuous(range = c(.05,.3))+
        scale_x_log10(breaks=c(.001,.01,.1,1),
                      labels=expression(10^-3,10^-2,10^1,10^0))+
        theme(aspect.ratio = 1/5,
              panel.grid = element_blank(),
              axis.title = element_text(face="italic"))
        
ggsave("manuscript/fit_time.tiff",height = 3,width = 9,dpi = 500)

a=log(fit.df$x[])
b=log(fit.df$fitness[])

a <- lm(fitness~x,log(fit.df[1:2500,]),na.action = na.exclude)

fit.k.df <- data.frame(d = degree(conet),
                       fit = fitness[V(conet)$name],
                       t=1:2928)

ggplot(fit.k.df, aes(x=d,y=fit))+
        geom_point(size=.6,shape=16)+
        stat_summary_bin(size = .2, 
                         geom = "pointrange",
                         color = pal_aaas(alpha = .6)(6)[5])+
        theme_bw()+
        #geom_vline(xintercept = 10,linetype=2)+
        xlab(expression(italic(k[t])))+
        ylab(expression(italic(beta)))+
        geom_smooth(method = "glm",se=TRUE,color="red",linetype=2)+
        scale_x_log10(breaks=c(.1,1,10,100),limits=c(1,400),
                      labels=expression(10^-1,10^0,10^1,10^2))+
        theme(strip.text = element_text(face = "italic"),
              aspect.ratio = 2/3,
              panel.grid = element_blank())+
        scale_y_log10(breaks=c(.1,1,10,100,1000),limits=c(.1,2000),
                      labels=expression(10^-1,10^0,10^1,10^2,10^3))

ggsave("manuscript/fit_kt.tiff",height = 3,width = 4,dpi = 500)

ancestor.df <- data.frame(d0 = tip.initial.degree,
                          fitness)
ggplot(ancestor.df,aes(x=factor(d0), y =fitness,color=factor(d0)))+
        geom_boxplot(outlier.color = "white")+
        geom_quasirandom(size=.2,alpha=.5)+
        scale_y_log10(breaks = c(.1,1,10,100,1000),
                      labels=expression(10^-1,10^0,10^1,10^2,10^3),
                      limits=c(0.1,1000))+
        xlab("Ancestor degree (k')")+
        ylab(expression(beta))+
        guides(color=FALSE)+
        theme_bw()+
        theme(panel.grid = element_blank())

ggsave("manuscript/Figures/ancestor.pdf",height = 3,width = 4)

summary(aov(log(fitness+1,10)~factor(d0),data = ancestor.df[,]))

# Graph files

gfile <- list.files("graph_growth/")
gfile <- gfile[order(as.numeric(gsub("g","",gfile)))]

# Change of lamda

lamda_r <- NULL
for(i in c(2:2925)){
        g <- read_graph(file = paste("graph_growth_random//",gfile[i],sep=""),format = "graphml")
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

lamda.df <- data.frame(x = 2925:1, 
                       lamda = lamda,
                       lamda.r = lamda.r)

ggplot(lamda.df, aes(x= x,y=lamda))+
        geom_point(color = pal_aaas(alpha = .3)(5)[4])+
        theme_bw()+xlab("t") + ylab(expression(lambda))+
        theme(axis.title = element_text(face="italic",family = "Times"),
              aspect.ratio = 1/3,
              panel.grid = element_blank())

ggsave("manuscript/lambda_subset.tiff",height = 2,width = 6,dpi = 500)

# Average separeation

separation <- NULL
for(i in 1:2925){
        g <- read_graph(file = paste("graph_growth/",gfile[i],sep=""),format = "graphml")
        g <- subgraph(g,degree(g)>0)
        apl <- average.path.length(g)
        separation <- c(separation, apl)
        print(i)
}

sep.df <- data.frame(x = 1:2900, 
                     separation = separation[2900:1])

ggplot(sep.df, aes(x= x,y=separation))+
        #geom_rect(xmin=1400,xmax=1800,ymin=.8,ymax=1.6,fill = "grey80")+
        geom_point(color = pal_aaas(alpha = .3)(5)[3])+
        theme_bw()+xlab("t")+ylab("d")+
        xlim(c(120,3000))+#ylim(c(3,4.5))+
        theme(axis.title = element_text(face="italic",family = "Times"),
              aspect.ratio = 1/3,
              panel.grid = element_blank())
ggsave("manuscript/separate.tiff",height = 2,width = 6,dpi = 500)

# clustering coefficient decays
cluster <- NULL
for(i in 1:2928){
        g <- read_graph(file = paste("graph_growth/",gfile[i],sep=""),format = "graphml")
        subset <- subgraph(g,degree(g)>0)
        apl <- transitivity(subset)
        cluster <- c(cluster, apl)
        print(i)
}

cluster.df <- data.frame(x = 1:2925, 
                         cluster = cluster[2925:1])

ggplot(cluster.df[150:2925,], aes(x= x,y=cluster))+
        geom_point(color = pal_aaas(alpha = .3)(5)[2])+
        theme_bw()+xlab("t")+ylab("C")+
        xlim(c(100,3000))+ #ylim(c(.35,.5))+
        theme(axis.title = element_text(face="italic",family = "Times"),
              aspect.ratio = 1/3,
              panel.grid = element_blank())

ggsave("manuscript/cluster.tiff",height = 2,width = 6,dpi = 500)

# cluster fraction

fraction <- NULL
for(i in 1:2925){
        g <- read_graph(file = paste("graph_growth/",gfile[i],sep=""),format = "graphml")
        subset <- subgraph(g,degree(g)>0)
        mod <- cluster_walktrap(subset)
        apl <- max(table(mod$membership))/length(V(subset))
        fraction <- c(fraction, apl)
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


# average degree

degree.mean <- NULL
for(i in 1:2928){
        subnet <- subgraph(conet,names(node.seq)[(i):2928])
        apl <- length(E(subnet))/length(V(subnet))
        degree.mean <- c(degree.mean, apl)
        print(i)
}

km.df <- data.frame(x = 1:2928, 
                    kd = degree.mean[2928:1])

ggplot(km.df, aes(x= x, y=kd))+
        geom_point(color = pal_aaas(alpha = .3)(5)[5])+
        theme_bw()+xlab("t")+ylab("<k>")+
        xlim(c(10,3000))+ ylim(c(0,20))+
        theme(axis.title = element_text(face="italic"),
              aspect.ratio = 1/3,
              panel.grid = element_blank())

ggsave("kmean.tiff",height = 2,width = 6,dpi = 500)

# turnover

fitness_turnover <- data.frame(c1=fitness-.1/(1-.1),
                               c2=fitness-.2/(1-.2),
                               c3=fitness-.3/(1-.3),
                               c4=fitness-.4/(1-.4),
                               c5=fitness-.5/(1-.5),
                               c6=fitness-.6/(1-.6),
                               c7=fitness-.7/(1-.7),
                               c8=fitness-.8/(1-.8),
                               c9=fitness-.9/(1-.9),
                               c0=fitness)

turnover.df <- data.frame(fitness,melt(fitness_turnover),
                          c.value = rep(c(seq(.1,.9,.1),0),each=length(fitness)))

g.turnover <- ggplot(turnover.df,
       aes(x=fitness,y=value,color=c.value))+
        geom_abline(slope = 1,color="grey",linetype=2)+
        geom_point(size=1,shape=1,alpha=.2)+
        scale_x_log10(breaks = c(.0001,.01,1,100,10000),
                      limits=c(0.01,10000),
                      labels = expression(10^-4,10^-2,10^0,10^2,10^4))+
        scale_y_log10(breaks = c(.0001,.01,1,100,10000),
                      limits=c(0.01,10000),
                      labels = expression(10^-4,10^-2,10^0,10^2,10^4))+
        scale_color_continuous_sequential(breaks = seq(0,.9,.1),name = "c")+
        theme_bw()+
        theme(aspect.ratio = 1,
              panel.grid = element_blank(),
              legend.title.align = .5,
              legend.position = c(.8,.3),
              legend.key.height = unit(5,"mm")
              )+
        xlab(expression(italic(beta)))+
        ylab(expression(italic(beta[turnover])))

ggsave("manuscript/turnover.tiff",g.turnover,height = 3,width = 5,dpi = 500)

# decay

v=-1
ti=500
d.v1.t2 <- fitness + log(((ti)^(-v)))/log(2928/(2928-ti))
ti=1000
d.v1.t3 <- fitness + log(((ti)^(-v)))/log(2928/(2928-ti))
ti=1500
d.v1.t4 <- fitness + log(((ti)^(-v)))/log(2928/(2928-ti))
ti=2000
d.v1.t5 <- fitness + log(((ti)^(-v)))/log(2928/(2928-ti))
ti=2500
d.v1.t6 <- fitness + log(((ti)^(-v)))/log(2928/(2928-ti))

d.v0.df <- data.frame(d.v1.t2,
                      d.v1.t3,
                      d.v1.t4,
                      d.v1.t5,
                      d.v1.t6)

v=0
ti=500
d.v1.t2 <- fitness + log(((ti)^(-v)))/log(2928/(2928-ti))
ti=1000
d.v1.t3 <- fitness + log(((ti)^(-v)))/log(2928/(2928-ti))
ti=1500
d.v1.t4 <- fitness + log(((ti)^(-v)))/log(2928/(2928-ti))
ti=2000
d.v1.t5 <- fitness + log(((ti)^(-v)))/log(2928/(2928-ti))
ti=2500
d.v1.t6 <- fitness + log(((ti)^(-v)))/log(2928/(2928-ti))

d.v2.df <- data.frame(d.v1.t2,
                      d.v1.t3,
                      d.v1.t4,
                      d.v1.t5,
                      d.v1.t6)

v=1
ti=500
d.v1.t2 <- fitness + log(((ti)^(-v)))/log(2928/(2928-ti))
ti=1000
d.v1.t3 <- fitness + log(((ti)^(-v)))/log(2928/(2928-ti))
ti=1500
d.v1.t4 <- fitness + log(((ti)^(-v)))/log(2928/(2928-ti))
ti=2000
d.v1.t5 <- fitness + log(((ti)^(-v)))/log(2928/(2928-ti))
ti=2500
d.v1.t6 <- fitness + log(((ti)^(-v)))/log(2928/(2928-ti))

d.v3.df <- data.frame(d.v1.t2,
                      d.v1.t3,
                      d.v1.t4,
                      d.v1.t5,
                      d.v1.t6)

v=10
ti=500
d.v1.t2 <- fitness + log(((ti)^(-v)))/log(2928/(2928-ti))
ti=1000
d.v1.t3 <- fitness + log(((ti)^(-v)))/log(2928/(2928-ti))
ti=1500
d.v1.t4 <- fitness + log(((ti)^(-v)))/log(2928/(2928-ti))
ti=2000
d.v1.t5 <- fitness + log(((ti)^(-v)))/log(2928/(2928-ti))
ti=2500
d.v1.t6 <- fitness + log(((ti)^(-v)))/log(2928/(2928-ti))

d.v4.df <- data.frame(d.v1.t2,
                      d.v1.t3,
                      d.v1.t4,
                      d.v1.t5,
                      d.v1.t6)

decay.df <- data.frame(fitness,melt(data.frame(v0=melt(d.v0.df),
                                               v1=melt(d.v1.df),
                                               v2=melt(d.v2.df),
                                               v3=melt(d.v3.df),
                                               v4=melt(d.v4.df))))

decay.df$v <- factor(decay.df$variable,
                             labels = c("v = -10",
                                        "v = 0",
                                        "v = 1",
                                        "v = 10"))

decay.df$t <- factor(rep(rep(c(2500,2000,1500,1000,500),each=2928),4),
                     levels = c(500,1000,1500,2000,2500))

g.decay <- ggplot(decay.df,aes(x=fitness,y=value,color=t))+
        geom_point(size=1,shape=16)+
        scale_color_aaas(alpha = .9,name=expression(t[i]))+
        theme_bw()+
        theme(aspect.ratio = 1,panel.grid = element_blank(),
              legend.title.align = .5)+
        facet_grid(~v)+
        xlab(expression(beta))+
        ylab(expression(beta[decay]))+
        scale_x_log10(breaks = c(1,10,100,1000,10000),
                      labels = expression(10^0,10^1,10^2,10^3,10^4))+
        scale_y_log10(breaks = c(.0001,.01,1,100,10000),
                      labels = expression(10^-4,10^-2,10^0,10^2,10^4))
        
ggsave("manuscript/decay.tiff",plot = g.decay,height = 3,width = 12,dpi = 500)


v=c(-1,seq(0,1,.1))
ti=2927
d.df <- fitness
for(i in 1:13){
        a1 <- log(c(ki*(ti^-v[i])+1)/c(tip.initial.degree+1))
        b <- log(2929/(2929-t0))
        d.vseq <- a1/b
        d.df <- data.frame(d.df,d.vseq)
}

decay.df1 <- data.frame(f = fitness,
                        d = melt(d.df[,c(2,4,6,8,10,12)])$value,
                        a = factor(rep(paste("v=",c(-1,seq(.2,1,.2)),sep=""),each=2928)))

g.decay <- ggplot(decay.df1,aes(x=f,y=d,color=a))+
        geom_point()+
        scale_color_aaas(alpha = .1)+
        xlab(expression(beta))+
        ylab(expression(beta[decay]))+
        scale_x_log10(limits = c(0.0001,1000),
                      breaks=c(0.001,0.1,10,1000),
                      labels=expression(10^-3,10^-1,10^1,10^3))+
        scale_y_log10(limits = c(0.0001,1000),
                      breaks=c(0.001,0.1,10,1000),
                      labels=expression(10^-3,10^-1,10^1,10^3))+
        theme_bw()+
        facet_wrap(~a,ncol = 6)+
        guides(color=FALSE)+
        theme(aspect.ratio = 1.25,
              strip.text = element_text(face="italic"))+
        geom_abline(slope=1,color="grey",linetype=2)

ggsave("manuscript/decay.pdf",plot = g.decay,height = 3,width = 10)


