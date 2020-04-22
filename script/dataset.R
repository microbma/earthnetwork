otu.abundance <- rep(0,2928)
names(otu.abundance) <- V(g.link)$name

sample1 <- V(g.link)$name[V(g.link)$s2==1]
table1 <- colSums(read.csv(paste("otutable/",
                                 inhabit.name[2],
                                 ".csv",sep = ""),
                           header = TRUE,row.names = 1)[,sample1])
otu.abundance[sample1] <- table1

sample1 <- V(g.link)$name[V(g.link)$s3==1]
table1 <- colSums(read.csv(paste("otutable/",
                                 inhabit.name[3],
                                 ".csv",sep = ""),
                           header = TRUE,row.names = 1)[,sample1])
otu.abundance[sample1] <- table1

sample1 <- V(g.link)$name[V(g.link)$s4==1]
table1 <- colSums(read.csv(paste("otutable/",
                                 inhabit.name[4],
                                 ".csv",sep = ""),
                           header = TRUE,row.names = 1)[,sample1])
otu.abundance[sample1] <- table1

sample1 <- V(g.link)$name[V(g.link)$s5==1]
table1 <- colSums(read.csv(paste("otutable/",
                                 inhabit.name[5],
                                 ".csv",sep = ""),
                           header = TRUE,row.names = 1)[,sample1])
otu.abundance[sample1] <- table1

sample1 <- V(g.link)$name[V(g.link)$s6==1]
table1 <- colSums(read.csv(paste("otutable/",
                                 inhabit.name[6],
                                 ".csv",sep = ""),
                           header = TRUE,row.names = 1)[,sample1])
otu.abundance[sample1] <- table1

sample1 <- V(g.link)$name[V(g.link)$s7==1]
table1 <- colSums(read.csv(paste("otutable/",
                                 inhabit.name[7],
                                 ".csv",sep = ""),
                           header = TRUE,row.names = 1)[,sample1])
otu.abundance[sample1] <- table1

sample1 <- V(g.link)$name[V(g.link)$s8==1]
table1 <- colSums(read.csv(paste("otutable/",
                                 inhabit.name[8],
                                 ".csv",sep = ""),
                           header = TRUE,row.names = 1)[,sample1])
otu.abundance[sample1] <- table1

sample1 <- V(g.link)$name[V(g.link)$s9==1]
table1 <- colSums(read.csv(paste("otutable/",
                                 inhabit.name[9],
                                 ".csv",sep = ""),
                           header = TRUE,row.names = 1)[,sample1])
otu.abundance[sample1] <- table1

sample1 <- V(g.link)$name[V(g.link)$s10==1]
table1 <- colSums(read.csv(paste("otutable/",
                                 inhabit.name[10],
                                 ".csv",sep = ""),
                           header = TRUE,row.names = 1)[,sample1])
otu.abundance[sample1] <- table1

sample1 <- V(g.link)$name[V(g.link)$s11==1]
table1 <- colSums(read.csv(paste("otutable/",
                                 inhabit.name[11],
                                 ".csv",sep = ""),
                           header = TRUE,row.names = 1)[,sample1])
otu.abundance[sample1] <- table1

sample1 <- V(g.link)$name[V(g.link)$s12==1]
table1 <- colSums(read.csv(paste("otutable/",
                                 inhabit.name[12],
                                 ".csv",sep = ""),
                           header = TRUE,row.names = 1)[,sample1])
otu.abundance[sample1] <- table1

sample1 <- V(g.link)$name[V(g.link)$s13==1]
table1 <- colSums(read.csv(paste("otutable/",
                                 inhabit.name[13],
                                 ".csv",sep = ""),
                           header = TRUE,row.names = 1)[,sample1])
otu.abundance[sample1] <- table1

sample1 <- V(g.link)$name[V(g.link)$s14==1]
table1 <- colSums(read.csv(paste("otutable/",
                                 inhabit.name[14],
                                 ".csv",sep = ""),
                           header = TRUE,row.names = 1)[,sample1])
otu.abundance[sample1] <- table1

abu.deg.df <- data.frame(abundance = otu.abundance,
                         degree = igraph::degree(g.link))
ggplot(abu.deg.df, aes(y = abundance, x = degree))+
        geom_point(color= alpha("#a45656",.5))+
        scale_y_log10(breaks = c(100,1000,10000,100000,1000000),
                      labels=expression(10^2,10^3,10^4,10^5,10^6))+
        scale_x_log10(breaks = c(1,10,100),
                      labels=expression(10^0,10^1,10^2))+
        theme_bw()+ylab("Reads")+xlab("k")+geom_smooth(method = "lm")

ggsave("Abu_deg.tiff",
       height = 4,
       width = 4,
       dpi = 500)

cor.test(abu.deg.df$abundance,abu.deg.df$degree,method = "spearman")


ggplot(melt(table(abu.deg.df$degree)),aes(x=Var1,y=value))+
        geom_point(color= alpha("#aa5656",.6))+
        scale_x_log10(breaks = c(1,10,100),
                      labels=expression(10^0,10^1,10^2))+
        scale_y_log10(breaks = c(1,10,100,1000),
                      labels=expression(10^0,10^1,10^2,10^3))+
        theme_bw()+ylab("Count number")+xlab("k")+
        geom_smooth()

ggsave("Deg_dis.tiff",
       height = 4,
       width = 4,
       dpi = 500)

a <- melt(table(abu.deg.df$degree))
summary(lm(data = a, Var1~value))


mod.size <- sort(table(modular.df$modularity_class),decreasing = TRUE)

mod.size1 <- c(mod.size[1:8], others = sum(mod.size[-c(1:8)])) 

names(mod.size1) <- paste(c(paste("Module",c(2, 4, 1, 3, 5, 6, 7, 8)),
                            "Other modules")," (",
                          round(mod.size1/sum(mod.size1)*100,0),"%)",sep="")

mod.size2 <- mod.size1[order(c(2, 4, 1, 3, 5, 6, 7, 8,9))]

tiff(filename = "module_profile.tiff",
        width = 3000,
        height = 3000,
        res = 600,
        compression = "lzw",
        type = "cairo",font=4)

pdf(file = "module_profile.pdf")
par(mar=rep(0,4))
pie(mod.size2,
    #labels = c(paste("Module",1:8),"Other modules"),
    radius = .5,
    border = FALSE, 
    col =  pal_npg(alpha = .8)(9),
    init.angle = 90,
    cex=1.2,font=4)
draw.circle(0,0,radius = .3, 
            border = "white", 
            col = "white")
dev.off()

