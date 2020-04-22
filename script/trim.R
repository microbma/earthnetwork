# generate trimmed OTU-table
## 
otutb.name <- list.files("otutable/")

ot.size <- NULL
for(i in 1:14){
        df <- read_csv(paste("otutable/",otutb.name[i],sep=""))
        ot.size <- rbind(ot.size,dim(df))
        print(i)
}

# select 360 samples and 400 OTUs after remove corpus of animal and plants

for(i in c(2:5,7:14)){
        df <- read_csv(paste("otutable/",otutb.name[i],sep=""))[,-1]
        set.seed(i)
        df1 <- df[rowSums(df)>0,]
        df.sample <- df1[sample(1:dim(df1)[1],360,replace = FALSE),]
        otu.abu <- colSums(df.sample)
        otu.abu[is.na(otu.abu)] <- 0
        otu.abu.trim <- sort(colSums(df.sample),decreasing = TRUE)[400]
        set.seed(i*100)
        df2 <- df.sample[, otu.abu >= otu.abu.trim]
        set.seed(i*100)
        df3 <- rowSums(df.sample[,otu.abu < otu.abu.trim])
        cond <- factor(rep("A",360))
        df.deseq <- DESeqDataSetFromMatrix(t(as.matrix(cbind(df2,df3)+1)), DataFrame(cond), ~1)
        df.dss <- DESeq(df.deseq,
                        fitType = "parametric",
                        test = "Wald")
        df.norm <- counts(df.dss,normalize=T)[1:400,]
        
        write.csv(df.norm, paste("trim_otu/",otutb.name[i],sep = ""))
        print(i)
}

# permutation

bray_p_matrix <- function(otutable=paste("trim_otu/",otutb.name[2],sep = "")){
        perm <- function(x){return(sample(x))}
        ot <- read.csv(otutable,header = TRUE,row.names = 1)
        bray.mean.perm <- matrix(0, nrow=nrow(ot), ncol=nrow(ot))
        for(l in 1:999){
                ot.perm <- apply(ot, MARGIN = 1, FUN = perm)
                bray.perm <- vegdist(t(ot.perm),method = "bray")
                bray.mean.perm  <- bray.mean.perm + as.matrix(bray.perm)  
               # print(x = paste("l =",l))
        }
        bray.mean.perm <- round(bray.mean.perm/999,4)
        return(1-bray.mean.perm)
}

bray_boot_matrix <- function(otutable="otutable/rhizosphere.csv"){
        ot <- read.csv(otutable,header = TRUE,row.names = 1)
        bray.mean.boot <- matrix(0, nrow=400, ncol=400)
        bray.sd.boot <- matrix(0, nrow=400, ncol=400)
        for(j in 1:999){
                bray.boot <- vegdist(ot[,sample(ncol(ot),replace = TRUE)],method = "bray")
                bray.mean.boot  <- bray.mean.boot + as.matrix(bray.boot)  
               # print(x = paste("j =",j))
        }
        for(k in 1:999){
                bray.boot <- vegdist(ot[,sample(ncol(ot),replace = TRUE)],method = "bray")
                bray.sd.boot  <- bray.sd.boot+(as.matrix(bray.boot)- bray.mean.boot/999)^2
               # print(x = paste("k =",k))
        }
        bray.mean.boot <- round(bray.mean.boot/999,4)
        bray.sd.boot <- round(sqrt(bray.sd.boot/999),4)
        return(list(mean = bray.mean.boot,
                    sd = bray.sd.boot))
}

for(i in c(11:14)){
        ot.trim.bray.boot <- bray_boot_matrix(paste("trim_otu/",otutb.name[i],sep = ""))
        ot.trim.bray.perm <- bray_p_matrix(paste("trim_otu/",otutb.name[i],sep = ""))
        write.csv(ot.trim.bray.perm, paste("trim_otu/perm/",otutb.name[i],sep = ""))
        write.csv(ot.trim.bray.boot[[1]],paste("trim_otu/boot_mean/",otutb.name[i],sep = ""))
        write.csv(ot.trim.bray.boot[[2]],paste("trim_otu/boot_sd/",otutb.name[i],sep = ""))
        print(paste("i =",i))
        }

sp_p_matrix <- function(otutable=paste("trim_otu/",otutb.name[2],sep = "")){
        perm <- function(x){return(sample(x))}
        ot <- read.csv(otutable,header = TRUE, row.names = 1)
        bray.mean.perm <- matrix(0, nrow=400, ncol=400)
        for(l in 1:999){
                ot.perm <- apply(ot, MARGIN = 1, FUN = perm)
                bray.perm <- cor(ot.perm,method = "spearman")
                bray.mean.perm  <- bray.mean.perm + as.matrix(bray.perm)  
                # print(x = paste("l =",l))
        }
        bray.mean.perm <- round(bray.mean.perm/999,4)
        return(bray.mean.perm)
}

sp_boot_matrix <- function(otutable="otutable/rhizosphere.csv"){
        ot <- t(read.csv(otutable,header = TRUE,row.names = 1))
        bray.mean.boot <- matrix(0, nrow=400, ncol=400)
        bray.sd.boot <- matrix(0, nrow=400, ncol=400)
        
        for(j in 1:999){
                bray.boot <- cor(ot[sample(nrow(ot),replace = TRUE),],
                                 method = "spearman")
                bray.mean.boot  <- bray.mean.boot + as.matrix(bray.boot)  
                # print(x = paste("j =",j))
        }
        for(k in 1:999){
                bray.boot <- cor(ot[sample(nrow(ot),replace = TRUE),],
                                 method = "spearman")
                bray.sd.boot  <- bray.sd.boot+(as.matrix(bray.boot) - bray.mean.boot/999)^2
                # print(x = paste("k =",k))
        }
        bray.mean.boot <- round(bray.mean.boot/999,4)
        bray.sd.boot <- round(sqrt(bray.sd.boot/999),4)
        return(list(mean = bray.mean.boot,
                    sd = bray.sd.boot))
}

for(i in c(2:5,7:14)){
        ot.trim.sp.boot <- sp_boot_matrix(paste("trim_otu/",otutb.name[i],sep = ""))
        ot.trim.sp.perm <- sp_p_matrix(paste("trim_otu/",otutb.name[i],sep = ""))
        write.csv(ot.trim.sp.perm, paste("trim_otu/perm/sp_",otutb.name[i],sep = ""))
        write.csv(ot.trim.sp.boot[[1]],paste("trim_otu/boot_mean/sp_",otutb.name[i],sep = ""))
        write.csv(ot.trim.sp.boot[[2]],paste("trim_otu/boot_sd/sp_",otutb.name[i],sep = ""))
        print(paste("i =",i))
}

## estimate p values

for(k in c(2:5,7:14)){
        sp.mean.perm <- read.csv(paste("trim_otu/perm/sp_",otutb.name[k],sep=""),row.names = 1)
        sp.mean.boot <- read.csv(paste("trim_otu/boot_mean/sp_",otutb.name[k],sep=""),row.names = 1)
        sp.sd.boot <- read.csv(paste("trim_otu/boot_sd/sp_",otutb.name[k],sep=""),row.names = 1)
        sp.mean.perm[is.na(sp.mean.perm)] <- 0
        sp.sd.boot[is.na(sp.sd.boot)] <- 0
        sp.mean.boot[is.na(sp.mean.boot)] <- 0
        if (sp.mean.perm[1,1]==0){sp.mean.perm <- 1- sp.mean.perm}
        if (sp.mean.boot[1,1]==0){sp.mean.boot <- 1- sp.mean.boot}
        
        cor.p <- matrix(0, nrow=400, ncol=400)
        for(i in 1:400) {
                for (j in 1:400){        
                        p <- ks.test(x = abs(sp.mean.perm[i,j]),
                                     "pnorm",
                                     mean = abs(sp.mean.boot[i,j]),
                                     sd = sp.sd.boot[i,j],
                                     alternative = "greater")
                        cor.p[i,j] <- p$p.value
                }
        }
        write.csv(cor.p,paste("trim_otu/pmat/sp_",otutb.name[k],sep=""))
        print(paste("k = ", k))
}

for(k in c(2:5,7:14)){
        bray.mean.perm <- read.csv(paste("trim_otu/perm/",otutb.name[k],sep=""),row.names = 1)
        bray.mean.boot <- read.csv(paste("trim_otu/boot_mean/",otutb.name[k],sep=""),row.names = 1)
        bray.sd.boot <- read.csv(paste("trim_otu/boot_sd/",otutb.name[k],sep=""),row.names = 1)
        if (bray.mean.perm[1,1]==0){bray.mean.perm <- 1- bray.mean.perm}
        if (bray.mean.boot[1,1]==0){bray.mean.boot <- 1- bray.mean.boot}
        
        cor.p <- matrix(0, nrow=400, ncol=400)
        for(i in 1:400) {
                for (j in 1:400){        
                        p <- ks.test(x = abs(bray.mean.perm[i,j]),
                                     "pnorm",
                                     mean = abs(bray.mean.boot[i,j]),
                                     sd = bray.sd.boot[i,j],
                                     alternative = "greater")
                        cor.p[i,j] <- p$p.value
                }
        }
        write.csv(cor.p,paste("trim_otu/pmat/",otutb.name[k],sep=""))
        print(paste("k = ", k))
}

# Reading a P matrix

cor.threshold.trim <- c(0,.36,0.48,.41,.38,0,.42,.44,.52,.44,.38,.48,.36,.46)
bray.threshold.trim <- c(0,.89,0.94,.93,.9,0,.91,.92,.95,.94,.9,.96,.91,.92)

# Network
for(i in c(2:5,7:14)) {
        p.mat.corr = 1-read.csv(paste("trim_otu/pmat/sp_",otutb.name[i],sep=""),
                                        header = TRUE,row.names = 1)
        p.mat.bray = 1-read.csv(paste("trim_otu/pmat/",otutb.name[i],sep=""),
                              header = TRUE,row.names = 1)
        cor.mat = read.csv(paste("trim_otu/boot_mean/sp_",otutb.name[i],sep=""),
                           header = TRUE, row.names = 1)
        bray.mat = read.csv(paste("trim_otu/boot_mean/",otutb.name[i],sep=""),
                            header = TRUE, row.names = 1)
        otu.tb = read.csv(paste("trim_otu/",otutb.name[i],sep=""),
                          header = TRUE,row.names = 1)
        graph.name = gsub("csv","graphml",otutb.name[i])
        
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
        adj.mat[] <- 0
        adj.mat[p.mat < 0.05] <- 1
        diag(adj.mat) <- 0
        adj.mat[abs(cor.mat) < .3] <- 0
        adj.mat[bray.mat < .9] <- 0
        g.adj <- graph_from_adjacency_matrix(adj.mat,
                                             diag = FALSE,
                                             weighted = TRUE,
                                             mode = "undirected")
        # annotation the nodes information with OTU number
        V(g.adj)$name <- rownames(otu.tb)
        cor.mat[is.na(cor.mat)] <- 0
        #cor.dir <- cor.mat[adj.mat == 1]
        E(g.adj)$dir <- cor.mat[lower.tri(adj.mat,diag = FALSE)][adj.mat[lower.tri(adj.mat,diag = FALSE)] == 1]
        print(g.adj)
        write.graph(g.adj,
                    file = paste("trim_otu/trim_graph/", 
                                 graph.name,
                                 sep = ""),
                    format = "graphml")
        print(i)
                  
}

# Check the direction

for(i in c(2:5,7:14)){
        g.cor <- read_graph(paste("trim_otu/trim_graph/",
                                  gsub("csv","graphml",otutb.name[i]),
                                  sep = ""),format = "graphml")
        E(g.cor)$dir -> a
        print(length(a[a<0]))
}

# vertex list
vertex.name.trim <- NULL
for(i in c(2:5,7:14)){
        g.cor <- read_graph(paste("trim_otu/trim_graph/",
                                  gsub("csv","graphml",otutb.name[i]),
                                  sep = ""),format = "graphml")
        vertex.name.trim <- cbind(vertex.name.trim,V(g.cor)$name)
}

vertex.unique <- unique(c(vertex.name.trim))

vertex.trim.df <-  matrix(0, nrow = 1737, ncol = 12)

for(i in 1:12){
        vertex.trim.df[vertex.unique %in% vertex.name.trim[,i],i] <- 1
}

table(rowSums(vertex.trim.df))

# Env-network

## Jaccard distance
edgelist.trim <- list(get.edgelist(read_graph(paste("trim_otu/trim_graph/",
                                                    gsub("csv","graphml",otutb.name[2]),
                                                    sep = ""),format = "graphml")),
                      get.edgelist(read_graph(paste("trim_otu/trim_graph/",
                                                    gsub("csv","graphml",otutb.name[3]),
                                                    sep = ""),format = "graphml")),
                      get.edgelist(read_graph(paste("trim_otu/trim_graph/",
                                                    gsub("csv","graphml",otutb.name[4]),
                                                    sep = ""),format = "graphml")),
                      get.edgelist(read_graph(paste("trim_otu/trim_graph/",
                                                    gsub("csv","graphml",otutb.name[5]),
                                                    sep = ""),format = "graphml")),
                      get.edgelist(read_graph(paste("trim_otu/trim_graph/",
                                                    gsub("csv","graphml",otutb.name[7]),
                                                    sep = ""),format = "graphml")),
                      get.edgelist(read_graph(paste("trim_otu/trim_graph/",
                                                    gsub("csv","graphml",otutb.name[8]),
                                                    sep = ""),format = "graphml")),
                      get.edgelist(read_graph(paste("trim_otu/trim_graph/",
                                                    gsub("csv","graphml",otutb.name[9]),
                                                    sep = ""),format = "graphml")),
                      get.edgelist(read_graph(paste("trim_otu/trim_graph/",
                                                    gsub("csv","graphml",otutb.name[10]),
                                                    sep = ""),format = "graphml")),
                      get.edgelist(read_graph(paste("trim_otu/trim_graph/",
                                                    gsub("csv","graphml",otutb.name[11]),
                                                    sep = ""),format = "graphml")),
                      get.edgelist(read_graph(paste("trim_otu/trim_graph/",
                                                    gsub("csv","graphml",otutb.name[12]),
                                                    sep = ""),format = "graphml")),
                      get.edgelist(read_graph(paste("trim_otu/trim_graph/",
                                                    gsub("csv","graphml",otutb.name[13]),
                                                    sep = ""),format = "graphml")),
                      get.edgelist(read_graph(paste("trim_otu/trim_graph/",
                                                    gsub("csv","graphml",otutb.name[14]),
                                                    sep = ""),format = "graphml"))
                      )



edgelist.trim.p <- lapply(edgelist.trim, FUN = function(x){paste(x[,1],x[,2],sep="-")})

jac.mat <- NULL
for(i in 1:12){
        jac.i <- NULL
        for(j in 1:12){
                a <- table(edgelist.trim.p[[i]] %in% edgelist.trim.p[[j]])[2]
                b <- length(unique(c(edgelist.trim.p[[i]], edgelist.trim.p[[j]])))
                jac <- a/b
                jac.i <- c(jac.i, jac)
        }
        jac.mat <- rbind(jac.mat,jac.i) 
}
diag(jac.mat) <- 0

row.names(jac.mat) <- env.name1[c(2:5,7:14)]
colnames(jac.mat) <- env.name1[c(2:5,7:14)]

jac.mat1 <- jac.mat
jac.mat1[jac.mat<0.02] <- 0

g.inh.mod <- graph_from_adjacency_matrix(jac.mat1, 
                                         diag = F,
                                         mode = "undirected", 
                                         weighted = TRUE)

plot(g.inh.mod)
write.graph(g.inh.mod,"inh_mod_trim.graphml",format = "graphml")

# Calculate network property #######
topo_trim <- NULL

for (i in c(2:5,7:14)){
        a <- read.graph(paste("trim_otu/trim_graph/",
                                    gsub("csv","graphml",otutb.name[i]),
                                    sep = ""),format = "graphml")
        a <- subgraph(a, v = igraph::degree(a)>0)
        topo <- c(vnum = length(V(a)$name),
                  vedge = length(E(a)),
                  ave = length(E(a))/length(V(a)$name),
                  d = diameter(a),
                  c = transitivity(a),
                  s = mean_distance(a),
                  b = mean(betweenness(a)),
                  m = modularity(a,membership(cluster_fast_greedy(a)))
        )
        topo_trim <- rbind(topo_trim,topo)
}

row.names(topo_trim) <- env.name1[c(2:5,7:14)]
colnames(topo_trim) <- c(topo.name[1:2],"Average degree", topo.name[3:7])

ggplot(melt(topo_trim[order(topo_trim[,2])[c(12:1)],-c(1,3)]),aes(x=Var1,y=value,fill=Var1))+
        geom_bar(stat = "identity")+
        facet_wrap(~Var2,scales = "free_y",nrow = 3)+
        theme_bw()+labs(x="",y="Topological property")+
        scale_fill_manual(values = morandi.color[24:35])+
        theme(legend.position = "none",
              panel.grid = element_blank(),
              axis.text.x = element_text(angle = -80,
                                         hjust = 0))

ggsave("topo_subnetwork.pdf",height = 5,width = 5)

# Keystone nodes #########

hub.df <- NULL
for(i in c(2:5,7:14)){
        g.cor <- read_graph(paste("trim_otu/trim_graph/",
                                  gsub("csv","graphml",otutb.name[i]),
                                  sep = ""),format = "graphml")
        d <- degree(g.cor)
        hub.df <- cbind(hub.df,V(g.cor)$name[order(d,decreasing = TRUE)[1:10]])
}

hub.dis.df <- matrix(0,nrow = 60,ncol=12)
for(i in 1:12){
        hub.dis.df[unique(c(hub.df)) %in% hub.df[,i],i] <- 1
}

trim.otu.name <- as.numeric(gsub("otu","",vertex.unique))
trim.otu.net <- DNAStringSet(otu.id[trim.otu.name])
names(trim.otu.net) <- vertex.unique
writeXStringSet(trim.otu.net,file = "trim_otu_net.fasta")

trim.taxonomy <- read.table("trim_otu.sintax",sep="\t")

source("sintax.R")
otu.df <- taxa.df(trim.taxonomy)

vertex.trim.df[grep("Latescibacteria",otu.df$taxa.df[,6]) ,]

row.names(hub.dis.df) <- gsub("_genera_incertae_sedis","",otu.df$taxa.df[unique(c(hub.df)),6])
colnames(hub.dis.df) <- env.name1[c(2:5,7:14)]
        
hub.dis.df1 <- hub.dis.df[order(rowSums(hub.dis.df),decreasing = TRUE),]

paste(paste("italic(",row.names(hub.dis.df1),")",sep=""),collapse = ",")

pdf("hub_heatmap.pdf")
pheatmap(t(hub.dis.df1),
         labels_col =  expression(italic(Gp2),italic(Nisaea),italic(Sphingobacterium),italic(Latescibacteria),italic(Acidobacteria_gp10),italic(Treponema),italic(Micrococcus),italic(Methanobrevibacter),italic(Intestinimonas),italic(Blastopirellula),italic(Chlorobium),italic(Leptospira),italic(WPS-1),italic(Pseudobacteroides),italic(Hyphomonas),italic(Acetohalobium),italic(Acidobacteria_gp4),italic(Sandaracinobacter),italic(Zavarzinella),italic(Ralstonia),italic(Ethanoligenens),italic(Desulfovibrio),italic(Melioribacter),italic(Persicobacter),italic(Parasporobacterium),italic(Kofleria),italic(Acidobacteria_gp6),italic(Acidobacteria_gp25),italic(Anaerobranca),italic(Cronobacter),italic(Acidobacteria_gp4),italic(Desulfuromonas),italic(Effusibacillus),italic(Gluconacetobacter),italic(Microgenomates),italic(Parcubacteria),italic(Entomoplasma),italic(Legionella),italic(Aminicenantes),italic(Saccharofermentans),italic(Plantactinospora),italic(Gemmatimonas),italic(Flexibacter),italic(Woesearchaeota),italic(Armatimonadetes_gp2),italic(Eubacterium),italic(Pontimonas),italic(Omnitrophica),italic(Parcubacteria),italic(Citrobacter),italic(Polyangium),italic(Rhodobaca),italic(Acidobacteria_gp5),italic(Cystobacter),italic(Acidobacteria_gp12),italic(Thermomarinilinea),italic(Anaerostipes),italic(Omnitrophica),italic(Hymenobacter),italic(Actinobacillus)),
         cellwidth = 6,
         cellheight = 6,
         legend = FALSE,
         color = c("black","darkgreen"),
         cutree_rows = 2, 
         #cutree_cols = 2,
         treeheight_row = 3,
         treeheight_col = 3,
         cluster_rows = TRUE, 
         cluster_cols = FALSE,
         fontsize_col = 6,
         fontsize_row = 6)
dev.off()

colSum(hub.dis.df)




## OS test 
### edgelist.trim.df

otu.trim.list <- list(read.csv(paste("trim_otu/",otutb.name[-c(1,6)][1],sep=""),row.names = 1),
                      read.csv(paste("trim_otu/",otutb.name[-c(1,6)][2],sep=""),row.names = 1),
                      read.csv(paste("trim_otu/",otutb.name[-c(1,6)][3],sep=""),row.names = 1),
                      read.csv(paste("trim_otu/",otutb.name[-c(1,6)][4],sep=""),row.names = 1),
                      read.csv(paste("trim_otu/",otutb.name[-c(1,6)][5],sep=""),row.names = 1),
                      read.csv(paste("trim_otu/",otutb.name[-c(1,6)][6],sep=""),row.names = 1),
                      read.csv(paste("trim_otu/",otutb.name[-c(1,6)][7],sep=""),row.names = 1),
                      read.csv(paste("trim_otu/",otutb.name[-c(1,6)][8],sep=""),row.names = 1),
                      read.csv(paste("trim_otu/",otutb.name[-c(1,6)][9],sep=""),row.names = 1),
                      read.csv(paste("trim_otu/",otutb.name[-c(1,6)][10],sep=""),row.names = 1),
                      read.csv(paste("trim_otu/",otutb.name[-c(1,6)][11],sep=""),row.names = 1),
                      read.csv(paste("trim_otu/",otutb.name[-c(1,6)][12],sep=""),row.names = 1))

edge.present.sum <- rowSums(edgelist.trim.df)
edgelist.trim.unique <- data.frame(left,right)

for(n in 2:12){
        multi.edge <- c(1:14201)[edge.present.sum == n]
        
        p <- NULL
        os <- NULL
        for(j in 1:length(multi.edge)){
                graph.p2 <- c(1:12)[edgelist.trim.df[multi.edge[j],]==1]
                ot.os <- NULL
                ot <- NULL
                for(k in 1:n){
                        ot.os <- c(ot.os,list(otu.trim.list[[graph.p2[k]]][as.character(unlist(edgelist.trim.unique[multi.edge[j],])),]))
                        ot <- rbind(ot,ot.os[[k]])
                }
                
                cor.org <- cor(ot[1,],ot[2,],method = "spearman")
                
                cor.os <- NULL
                for(k in 1:n){
                        os1 <- lapply(ot.os, FUN = function(x){x[1,]})
                        os2 <- lapply(ot.os, FUN = function(x){x[2,]})
                        cor.os <- c(cor.os, 
                                    cor(unlist(os1), 
                                        unlist(os2),
                                        method = "spearman"))
                }
                
                p.os <- NULL
                for(k in 1:n){
                        p.os1 <- 0
                        for(i in 1:500){
                                cor.random <- cor(ot[,-sample(1:dim(ot)[2],dim(ot.os[[k]])[1])])[1,2]/cor.os[k]
                                if (cor.random < 1){ p.os1 <- p.os1+1 }
                        }
                        p.os <- c(p.os, p.os1)
                }
                
                p <- rbind(p,p.os/500)
                os <- rbind(os, c(cor.os/cor.org))
        }
        
        write.csv(p,file = paste("os_test_trim/p_",n,sep=""),row.names = TRUE,col.names = TRUE)
        write.csv(os,file = paste("os_test_trim/os_",n,sep=""),row.names = TRUE,col.names = TRUE)
        print(n)
}

local <- NULL
local1 <- NULL
for(i in c(2:6,8:10)){
        p.file <- read.csv(file = paste("os_test_trim/p_",i,sep=""))
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

# Dataset 2

edge.list.comb <- rbind(data.frame(type = "Generalist edge",generalist.edge.list),
                        data.frame(type = "Specialist edge (generalist vertex pair)",sp.edge.list),
                        data.frame(type = "Specialist edge (specialist vertex pair)",sp1.edge.list)
                        )
write.csv(edge.list.comb ,file = "manuscript/MS_V11/Dataset2.csv",row.names = FALSE)

# Dataset 3
write.csv(neg.list ,file = "manuscript/MS_V11/Dataset3.csv",row.names = FALSE)
