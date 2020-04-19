#  Earth microbial network project

The repository provides R script for the microbial network analysis of the Earth Microbiome Project (EMP).

There network file is stored as *entirenetwork.graphml*.

## Dataset preparation 

### Abundance table from the EMP

- We then extracted 14 count matrices for 14 environmental categories at level-3 of the EMP ontology (Table S1) from the 90-bp Deblur BIOM table. We filtered the OTUs with relative abundance less than 0.001% and presenting in less than 10% of samples in corresponding count matrices of environments.
- The corresponding  R script is *dataset.R*. 

### Trimmed dataset

- We kept 400 top-abundant ESVs and randomly selected 360 samples in the trimmed microbial community matrices. 
- The corresponding  R script is *trim.R*. 

## Network inference

- Microbial taxon-taxon association networks were constructed by selecting Spearman correlation and Bray-Curtis dissimilarity measures. 
- The corresponding  R script is *netinfer.R*. 

## Omission score

- The impact of environmental categories on the Spearman correlation of each edge in the network was assessed through dividing the absolute omission score (OS) (Spearman correlation without the environmental categories) by the absolute original Spearman score. 
- The corresponding  R script is *os.R*. 


## Overrepresentation analysis

- Taxon–taxon counts at high taxonomic ranks were assessed for significance using the hypergeometric distribution in the R stats::phyper. 
- Mutual exclusion versus co-presence analysis was performed using the binomial distribution implemented in the R stats::pbinom, with background probability estimated by the frequency of edges in the network.
- The corresponding  R script is *overrepresentation.R*. 

## Topological features

- Topological features were estimated with igraph package.
- The corresponding  R script is *nettopo.R*. 

## Generalist and specialist edges

- Edges present in only one subnetwork were specialist edges, which were further clustered into two groups: a specialist edge linking a specialist vertex pair or the same linking a generalist vertex pair.
- The corresponding  R script is *generalist.R*. 

## Hubs

- We identified ten hubs at the top-degree from each subnetwork inferred from the 12 trimmed datasets.
- The corresponding  R script is *hub.R*. 


## Negative edges

- We counted the number and percentage of negative edges in the subnetworks inferred from the 12 trimmed datasets
- The corresponding  R script is *negative.R*. 


