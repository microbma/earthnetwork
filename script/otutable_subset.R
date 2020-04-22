# Read Biom file

require(rhdf5)
require(SparseM)

h5ls("deblur/deblur.biom")

counts <- h5read("deblur/deblur.biom",name = "observation/matrix")
sample <- h5read("deblur/deblur.biom",name = "sample/ids")

sample.matrix <- h5read("deblur/deblur.biom",name = "sample/")
otu.id <- h5read("deblur/deblur.biom",name = "observation/ids")


taxon <- h5read("deblur/deblur.biom",name = "observation/metadata/taxonomy")
write.csv(taxon,"taxon.csv")

#taxdf <- h5read("deblur/deblur.biom",name = "observation/metadata/")
otu.number <- counts$indices
for(i in 1:(length(counts$indptr)-1)){
        otu.number[counts$indptr[i]:(counts$indptr[i+1]-1)] <- i
}
otu.number[length(otu.number)] <- 317315
require(hdf5r)
# Read mappling file

mapping <- read.table("mapping_files/emp_qiime_mapping_release1.tsv",
                      stringsAsFactors = FALSE,
                      header=TRUE,
                      sep = "\t",
                      quote = "",
                      fill = TRUE)

otutable_subset <- function(subset = "Plant surface", samples = 65){
        subset.sample <- mapping$SampleID[mapping$empo_3 == subset]

        # find the positions of animal corpus samples AND remove NA values
        subset.sample.position <- pmatch(subset.sample,sample)
        subset.sample.position <- subset.sample.position[!is.na(subset.sample.position)]

        # find the data 
        data.subset <- counts$data[counts$indices %in% subset.sample.position]
        # find the indices
        indices.subset <- counts$indices[counts$indices %in% subset.sample.position]
        level.pos <- levels(factor(indices.subset))

        otu.subset <- otu.number[counts$indices %in% subset.sample.position]

        indptr.subset <- NULL

        for (i in 1:(length(indices.subset)-1)){
                 if (otu.subset[i] != otu.subset[i+1]){
                         indptr.subset <- c(indptr.subset,i)}}

        indptr.subset <- c(indptr.subset,i+1)

        n1 <- 10
        p <- 5
        a <- rnorm(n1*p)
        a[abs(a)<0.5] <- 0
        A <- matrix(a,n1,p)
        B <- t(A)%*%A
        A.csr <- as.matrix.csr(A)

        slot(A.csr,"ra") <- as.numeric(data.subset)
        slot(A.csr,"ja") <- as.integer(c(factor(indices.subset)))
        slot(A.csr,"ia") <- as.integer(indptr.subset)
        slot(A.csr,"dimension") <- as.integer(c(length(indptr.subset),length(subset.sample.position)))

        subset.matrix <- as.matrix(A.csr)
        colnames(subset.matrix) <- paste("D",1:dim(subset.matrix)[2],sep = "")
        row.names(subset.matrix) <- paste("otu",unique(otu.subset),sep="")

        #subset.matrix <- read.csv("correlation_matrix/subset/subset_otu_table.csv",row.names = 1)
        subset.matrix.rare <- subset.matrix[rowSums(subset.matrix) > sum(subset.matrix)/100000,-1]
        subset.matrix.present <- subset.matrix.rare
        subset.matrix.present[subset.matrix.present!=0] <- 1
        subset.matrix.rare.present <- t(subset.matrix.rare[rowSums(subset.matrix.present) > samples,])
        return(subset.matrix.rare.present)
}

# generate OTU-Tables
group <- sort(table(mapping$empo_3))
group.name <- c("N","H","Mock","Aerosol","Surface","plant_corpus","black","animal_pgut","sediment_nonsaline",
                "animal_corpus","rhizosphere","sediment_saline","water_saline","animal_secretion",
                "surface_nonsaline","plant_surface","animal_surface","animal_dgut","soil","water_nonsaline")
for(i in 8:16){
        otutable <- otutable_subset(subset = names(group)[i], samples = round(group[i]/10,1))
        write.csv(otutable,paste("otutable/",group.name[i],".csv",sep=""))
        print(group.name[i])
}

for(i in 17:19){
        otutable <- otutable_subset(subset = names(group)[i], samples = 200)
        write.csv(otutable,paste("otutable/",group.name[i],".csv",sep=""))
        print(group.name[i])
}

otutable.names <- list.files("otutable/")
num <- NULL
for(i in 1:14){
        ot <- read.csv(paste("otutable/",otutable.names[i],sep=""))
        num <- rbind(num,dim(ot))
}
