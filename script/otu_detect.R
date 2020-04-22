#! /usr/bin/env Rscript

# Read Biom file
require(biomformat)
require(rhdf5)
require(SparseM)


counts <- h5read("emp_deblur_90bp.qc_filtered.biom",name = "observation/matrix")
sample <- h5read("emp_deblur_90bp.qc_filtered.biom",name = "sample/ids")

# Read mappling file

mapping <- read.table("emp_qiime_mapping_release1.tsv",
                      stringsAsFactors = FALSE,
                      header=TRUE,
                      sep = "\t",
                      quote = "",
                      fill = TRUE)
soil.sample <- mapping$SampleID[mapping$empo_3 == "Soil (non-saline)"]

soil.sample.position <- pmatch(soil.sample,sample)
soil.sample.position <- soil.sample.position[!is.na(soil.sample.position)]

data.soil <- counts$data[counts$indices %in% soil.sample.position]
indices.soil <- counts$indices[counts$indices %in% soil.sample.position]
level.pos <- levels(factor(indices.soil))

indptr.soil <- NULL

for (i in 1:(length(indices.soil)-1)){
        if (indices.soil[i] >= indices.soil[i+1]){
                indptr.soil <- c(indptr.soil,i)}}


n1 <- 10
p <- 5
a <- rnorm(n1*p)
a[abs(a)<0.5] <- 0
A <- matrix(a,n1,p)
B <- t(A)%*%A
A.csr <- as.matrix.csr(A)

slot(A.csr,"ra") <- as.numeric(data.soil)
slot(A.csr,"ja") <- as.integer(c(factor(indices.soil)))
slot(A.csr,"ia") <- as.integer(indptr.soil)
slot(A.csr,"dimension") <- as.integer(c(length(indptr.soil),length(soil.sample.position)))

soil.matrix <- as.matrix(A.csr)
row.names(soil.matrix) <- paste("D",1:dim(soil.matrix)[1])
write.csv(soil.matrix,"soil_otu_table.csv")

#soil.matrix <- read.csv("soil_otu_table.csv")
soil.matrix.rare <- soil.matrix[rowSums(soil.matrix)>220,]
soil.matrix.present <- soil.matrix.rare
soil.matrix.present[soil.matrix.present!=0] <- 1
soil.matrix.rare.present <- soil.matrix.rare[rowSums(soil.matrix.present)>100,]

otu.detect <- function(otu.table = soil.matrix.rare.present[1,-1]){
        profile <- otu.table
        profile1 <- profile[profile!=0]
        col.num <- c(1:4270)[profile!=0]
        
        pos1 <- counts$indices[as.numeric(counts$data)== profile1[1]] # The first number values
        pos2 <- which(pos1 == level.pos[col.num[1]]) # values in the levels of factor(indices.soil)
        pos3 <- which(as.numeric(counts$data) == profile1[1])
        pos4 <- pos3[pos2]
        
        out1 <- sapply(pos4, FUN= function(x){table(x<counts$indptr)[1]})
        
        pos5 <- counts$indices[as.numeric(counts$data)== profile1[2]] # The first number values
        pos6 <- which(pos5==level.pos[col.num[2]]) # values in the levels of factor(indices.soil)
        pos7 <- which(as.numeric(counts$data) ==  profile1[2])
        pos8 <- pos7[pos6]
        
        out2 <- sapply(pos8, FUN= function(x){table(x<counts$indptr)[1]})
        
        pos9 <- counts$indices[as.numeric(counts$data)==  profile1[3]] # The first number values
        pos10 <- which(pos9==level.pos[col.num[3]]) # values in the levels of factor(indices.soil)
        pos11 <- which(as.numeric(counts$data) ==  profile1[3])
        pos12 <- pos11[pos10]
        
        out3 <- sapply(pos12, FUN= function(x){table(x<counts$indptr)[1]})
        
        pos13 <- counts$indices[as.numeric(counts$data)==  profile1[4]] # The first number values
        pos14 <- which(pos13==level.pos[col.num[4]]) # values in the levels of factor(indices.soil)
        pos15 <- which(as.numeric(counts$data) ==  profile1[4])
        pos16 <- pos15[pos14]
        
        out4 <- sapply(pos16, FUN= function(x){table(x<counts$indptr)[1]})
        otu.num <- out3[out3 %in% out2[out2 %in% out1]] 
        return(otu.num)
}

otu.num <- apply(soil.matrix.rare.present,MARGIN = 1, FUN = otu.detect)
write.table(otu.num,file = "otu_num.csv")
              
              