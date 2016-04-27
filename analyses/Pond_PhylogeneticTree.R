################################################################################
#                                                                              #
#  Pond Microbe Stability Experiment: Microbial Community Phylogenetic Tree    #
#                                                                              #
################################################################################
#                                                                              #
#  Written by: Mario Muscarella                                                #
#                                                                              #
#	Last update: 2015/05/19                                                      #
#                                                                              #
################################################################################
#                                                                              #
# Notes: This code conducts the calculations needed for microbial              #
#        phylogenetic tree as part of the Pond Experiment (2009). For general  #
#        information about this project see: Lennon et al. 2013 Plos One,      #
#        Jones and Lennon 2015 Ecology.                                        #
#                                                                              #
# Dependencies: 1) R 3.1.3 (Smooth Sidewalk)                                   #
#               2) Vegan v2.2-1                                                #
#               3) Picante v1.6-2                                              #
#               4) Ape v3.2                                                    #
#               3) DiversityFunctions.R (2014/09/04; included in Repo;         #
#                    https://github.com/mmuscarella/Diversity_Functions)       #
#               4) MothurTools.R (2015/02/22; included in Repo)                #
#                                                                              #
# Issues: None Identified (yet)                                                #
#         1.                                                                   #
#                                                                              #
# Recent Changes:                                                              #
#         1.                                                                   #
#                                                                              #
# Future Changes (To-Do List):                                                 #
#         1.                                                                   #
#                                                                              #
################################################################################

# Setup Work Environment (Directory, Packages, Source Code, Functions)
rm(list=ls())
setwd("~/GitHub/SubsidyMicrobialStability")
source("./bin/DiversityFunctions.R")
source("./bin/MothurTools.R")
se <- function(x, ...){sd(x, ...)/sqrt(length(na.omit(x)))}
require("vegan")
require("seqRFLP")
require("stringr")
require("ape")
require("phytools")
require("phylotools")

# Import Responder Tables
DNAcorr <- read.table("./data/DNA_corr.txt", sep = ",", as.is = TRUE, header=T)
RNAcorr <- read.table("./data/RNA_corr.txt", sep = ",", as.is = TRUE, header=T)

# Rename Sequences in Tree
pond.fasta <- seqRFLP::read.fasta(file="./mothur/output/Pond_all.final.FW.tx.1.rep.fasta")
old.names <- seqRFLP::gnames.fas(pond.fasta)
new.names <- gsub("\\|.*$", "", gsub("^.*?\\t", "", old.names))
new.names <- as.character(str_pad(new.names, width=3, pad="0"))
new.names <- paste("OTU", new.names, sep="")
names <- cbind(old.names, new.names)
pond.fasta <- phylotools::rename.fasta(pond.fasta, names)
write.fasta(pond.fasta, "./data/Pond.responders.fasta")

# Outside of R I merged the responders with the outgroup (Aquifex)
# The resulting file was saved as Pond.merged.fasta")
# using Mega a created a tree

# Import Raw Tree from Mega
phy <- read.tree("./data/FW.tree.consensus.nwk")

# Identify Outgroup Sequence
outgroup <- match("Aquifex", phy$tip.label)

# Root the Tree {ape}
phy <- root(phy, outgroup, resolve.root = TRUE)

# Identify responding taxa
responders <- union(DNAcorr$X, RNAcorr$X)

responders <- gsub("Otu", "OTU", responders)

# Drop Tips of Zero-Occurrence OTUs (Removes Taxa Only Found via RNA Sequencing)
phy.res <- drop.tip(phy, phy$tip.label[!phy$tip.label %in% c(responders, "Aquifex")])

# Initial Plot
plot(phy.res)
# Save tree
write.tree(phy.res, file="./data/Responders.tree.nwk")

# I used the tree above to make the tree figure in the paper.
# Not very reproducible :(
