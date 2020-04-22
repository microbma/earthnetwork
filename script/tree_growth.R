# track the growth of tree
require(igraph)
require(ape)
library(ggtree)

# Read tree file
tree <- read.tree("tree/mc/otu.tree")

# molecular dating by penalised likelihood
chr <- chronos(tree)

# Give the edge length in chr to tree
tree1 <- tree
tree1$edge.length <- chr$edge.length

# get the value of nodes
gdata <- ggtree(tree1,size=.1,color=grey(.1)) # layout = "circular"
node.df <-  gdata$data[-c(1:2929),]
node.df.rank <- node.df[order(node.df$x,decreasing = TRUE),]
tip.df <- gdata$data[1:2928,]

node.df.rank.label <- node.df.rank$x[2925-c(3,10,50,100,500,1000,2500)]

tip.df1 <- tip.df
tip.df1$fit <- fitness

gdata+scale_x_continuous(breaks = c(node.df.rank.label[-c(1:2)]),
                         name ="Simulating evoluiton process",
                         labels = c(50,100,500,1000,2500))+
        theme_tree2()+
        #geom_point2(size=.3,color=pal_aaas(alpha = .6)(8)[8])+
        geom_vline(xintercept = c(node.df.rank.label[-c(1:2)]),
                   color=pal_aaas(alpha = .8)(3)[3],
                   linetype=2)+
        theme(axis.text.x = element_text(angle=90),
              axis.line.x = element_line(color = "black"))

ggsave("tree.pdf",width = 5,height = 2)

gdata1 <- ggtree(tree1,size=.1,color=grey(.1),layout = "fan")

gdata1 + theme_tree2()+
         geom_point(data = tip.df1,
                    size=2,alpha=.1,
                    aes(x=x,y=y,
                        color=log(fit)),
                    inherit.aes = FALSE)+
                scale_color_continuous_qualitative()+
        theme(axis.line = element_blank(),axis.text = element_blank())

fit.taxa.df <- data.frame(taxa = V(g.link)$phylum,
                          fit = fitness[V(g.link)$name])

fit.taxa.rank <- fit.taxa.df %>%
        group_by(taxa) %>%
        summarise(ave = quantile(fit)[3],num = length(fit)) %>%
        filter(num >20) %>%
        arrange(ave)

fit.taxa.df$taxa <- factor(fit.taxa.df$taxa,levels = as.character(fit.taxa.rank$taxa[-c(6)]),exclude = TRUE)

fit.taxa.df <- fit.taxa.df[!is.na(fit.taxa.df$taxa),]

ggplot(fit.taxa.df, aes(x=taxa,y=fit,color=taxa))+
        geom_quasirandom(alpha=.3)+
        geom_violin(fill=NA,draw_quantiles = .5,size=.5)+
        scale_y_log10(limits = c(.1,8000),
                      breaks = c(.1,1,10,100,1000,10000),
                      labels= expression(10^-1,10^0,10^1,10^2,10^3,10^4))+
        scale_color_discrete_qualitative()+
        guides(color=FALSE)+
        theme_bw()+xlab("Phylum")+ylab(expression(beta))+
        theme(axis.text.x = element_text(angle = -75,hjust = 0),panel.grid = element_blank())
        
TukeyHSD(aov(fit~taxa,fit.taxa.df))

ggsave("taxa_fit.pdf",height = 4,width = 4)

### define a function for the m value (inital degree) and t values for each nodes ###

# define a function for merge the nodes and edges

MergeNode <- function(g = g0, p1 = 1119, p2 = 1672){
        a <- E(g)[from(p1)]
        b <- E(g)[from(p2)]
        g1 <- g-a[!a %in% b]
        c <- E(g1)[from(p2)]
        g2 <- g1-c[!c %in% a]
        return(g2)
}

# network growth
m0 <- NULL # initial links
g0 <- conet # assign network to g0

for(i in 1:dim(node.df.rank)[1]){  # a cycle for all nodes in the tree
        child.pos <- grep(as.character(node.df.rank$node[i]),
                          as.character(gdata$data$parent)) # the childrens of the nodes
        
        # set a if...else section to find out the corresponding vertex of both children nodes in network
        if (child.pos[1] <= 2928){
                p1 <- pmatch(tip.df$label[child.pos[1]],
                           unlist(V(g0)$name))
        }else{
                p1.name <- names(m0)[pmatch(gdata$data$node[child.pos[1]],
                                            node.df.rank$node)]
                p1 <- pmatch(p1.name,unlist(V(g0)$name))[1]
        } # P1
        
        if (child.pos[2] <= 2928){
                p2 <- pmatch(tip.df$label[child.pos[2]],
                             unlist(V(g0)$name))[1]
        }else{
                p2.name <- names(m0)[pmatch(gdata$data$node[child.pos[2]],
                                       node.df.rank$node)]
                p2 <- pmatch(p2.name, unlist(V(g0)$name))[1]
        } # P2
        
        mapping <- 1:(length(V(g0)$name)) # define a numeric vector for mapping 
        # contract the network
        if(p1 < p2){
                mapping[p2] <- p1
                mapping1 <- mapping
                mapping1[mapping > p2] <- mapping[mapping > p2] - 1
                g1 <- contract(MergeNode(g0,p1=p1,p2=p2),mapping1,vertex.attr.comb=toString)
                g1 <- simplify(g1)
                V(g1)$name[p1] <- paste("node",i,sep="")
                m <- neighbors(g0,V(g0)$name[p1]) %in% neighbors(g0,V(g0)$name[p2])
                m0 <- c(m0,length(m[m==TRUE]))
                names(m0) <- paste("node",1:i,sep="")
        }else{
                mapping[p1] <- p2
                mapping1 <- mapping
                mapping1[mapping > p1] <- mapping[mapping > p1] - 1
                g1 <- contract(MergeNode(g0,p1=p1,p2=p2),mapping1,vertex.attr.comb=toString)
                g1 <- simplify(g1)
                V(g1)$name[p2] <- paste("node",i,sep="")
                m <- neighbors(g0,V(g0)$name[p1]) %in% neighbors(g0,V(g0)$name[p2])
                m0 <- c(m0,length(m[m==TRUE]))
                names(m0) <- paste("node",1:i,sep="")
        }
        
        g0 <- g1
        write_graph(g0,paste("graph_growth/g",i,sep=""),"graphml")
        print(i)
        }

# Assign the initial degree to each vertex (tips in tree)

m1 <- m0
names(m1) <- node.df.rank$node
tip.initial.degree <- m1[as.character(tip.df$parent)]
ki <- degree(conet)[tip.df$label] 

plot(g,vertex.size=1,vertex.label=NA)


# Dynamic graph

g.layout <- read.graph("graph_growth/g_link.net",format = "pajek")

layout.df <- data.frame(x=V(g.layout)$x,
                        y=V(g.layout)$y,
                        name = V(conet)$name)

png(file="animation/output_%03d.png", width=1000,height=1000)

#Time loop starts
for(time in seq(1, total_time, dt)){   #dt is the defined interval between successive plots
        
        gt <- delete_edges(g,which(E(g)$time > time)) #remove absent edges
        
        layout.new <- layout_with_fr(gt,coords=layout.old,niter=10,start.temp=0.05,grid="nogrid") #jitter layout
        
        plot(gt,layout=layout.new,vertex.label="",vertex.size=1+2*log(degree(gt)),
             vertex.frame.color=V(g)$color,edge.width=1.5,asp=9/16,margin=-0.15) #make plot
        
        layout.old <- layout.new #keep layout for next image
}

g <- read.graph("graph_growth/g2500","graphml")
g1 <- induced_subgraph(g, v = degree(g)!=0)
ly.new <- as.matrix(layout.df[pmatch(V(g1)$name,V(conet)$name),1:2])
ly.new[is.na(ly.new)] <- 0
set.seed(999)
plot(g1, 
     layout=ly.new,
     vertex.label=NA,
     vertex.size=2,
     vertex.frame.color=NA,
     edge.color=grey(level = .6,.4),
     edge.width=.2)

for (i in 2900:1){
                print(i)
        
                g0 <- read.graph(paste("graph_growth/g",i+20,sep=""),"graphml") # read history graph
                g1 <- read.graph(paste("graph_growth/g",i,sep=""),"graphml") # read current graph
                #g1 <- induced_subgraph(g, v = degree(g)!=0) # remove seperated vertex
                ly.new <- as.matrix(layout.df[pmatch(V(g1)$name,V(conet)$name),1:2]) # the layout for current graph
                
                ly.new[is.na(ly.new)] <- 0 # assign posistion for vertices without layout
                # assign color to new vertices 
                V(g1)$color <- "white" 
                V(g1)$color[!V(g1)$name %in% V(g0)$name] <- getPalette(1)
                # assign size to new vertices
                V(g1)$size <- 1.5
                V(g1)$size[!V(g1)$name %in% V(g0)$name] <- 3.5
                # assign color to new edges
                E(g1)$color <- grey(level = .6,
                                    alpha = .6)
                E(g1)[from(c(1:vcount(g1))[V(g1)$color != "white"])]$color <- getPalette(1)
                # assign size to new edges
                E(g1)$size <- .5
                E(g1)[from(c(1:vcount(g1))[V(g1)$color != "white"])]$size <- 2
                
                # export PNG file
                png(filename = paste("animation/g",
                                     substr(as.character(10000+2901-i),
                                            start = 2,
                                            stop = 5),
                                     ".png",
                                     sep=""),
                    width = 500,
                    height = 500,
                    res = 128)
                
                par(bg="black",mar=rep(0,4))
                plot(x = g1, 
                     layout = ly.new,
                     vertex.label = NA,
                     vertex.size = V(g1)$size,
                     vertex.frame.color = NA,
                     edge.width = E(g1)$size
                     )
                dev.off()
                }
