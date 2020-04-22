# Molecular dating

require(ape)
entire.tree <- read.tree("deblur/emp90.5000_1000_rxbl_placement_pruned75.tog.tre")
chr <- chronos(entire.tree)
