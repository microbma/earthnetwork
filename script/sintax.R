sintax.info <- function(taxa.inf = taxa.inf){
        taxa.split <- strsplit(as.character(taxa.inf[,2]), ",")
        taxa.df <- data.frame(Domain = gsub(pattern = "d:",
                                            replacement = "" ,
                                            as.character(sapply(X = taxa.split, function(x){x[grep(pattern = "d:",x)]}))),
                              Phylum = gsub(pattern = "p:",
                                            replacement = "" ,
                                            as.character(sapply(X = taxa.split, function(x){x[grep(pattern = "p:",x)]}))),
                              Class = gsub(pattern = "c:",
                                           replacement = "" ,
                                           as.character(sapply(X = taxa.split, function(x){x[grep(pattern = "c:",x)]}))),
                              Order = gsub(pattern = "o:",
                                           replacement = "" ,
                                           as.character(sapply(X = taxa.split, function(x){x[grep(pattern = "o:",x)]}))),
                              Family = gsub(pattern = "f:",
                                            replacement = "" ,
                                            as.character(sapply(X = taxa.split, function(x){x[grep(pattern = "f:",x)]}))),
                              Genus = gsub(pattern = "g:",
                                           replacement = "" ,
                                           as.character(sapply(X = taxa.split, function(x){x[grep(pattern = "g:",x)]})))
        )
        
        taxa.conf <- apply(taxa.df,
                           MARGIN = 2,
                           function(x){
                                   step1 <- gsub(pattern = ".*\\(",
                                                 replacement = "" ,
                                                 x = x,
                                                 perl = T)
                                   step2 <- gsub(pattern = "\\)",
                                                 replacement = "" ,
                                                 x = step1,
                                                 perl = T)
                                   return(as.numeric(step2))
                           })
        
        taxa.name <- apply(taxa.df,
                           MARGIN = 2,
                           function(x){gsub(pattern = "\\([0-9.]*\\)",
                                            replacement = "" ,
                                            x = x,
                                            perl = T)})
        
        taxa.name.unknown <- taxa.name
        taxa.name.unknown[taxa.conf < .8] <- "Unknown"
        return(list(taxa.name.unknown,taxa.name))
}
