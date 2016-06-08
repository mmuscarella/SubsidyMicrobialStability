################################################################################
#                                                                              #
#	Pond Microbe Stability Experiment:  Community Assembly                       #
#                                                                              #
################################################################################
#                                                                              #
#	Written by: Mario Muscarella                                                 #
#                                                                              #
#	Last update: 2015/05/18                                                      #
#                                                                              #
################################################################################
#                                                                              #
# Notes: This code conducts the calculations needed determining selection      #
#        using correlations as part of the Pond Experiment (2009). For general #
#        information about this project see: Lennon et al. 2013 Plos One,      #
#        Jones and Lennon 2015 Ecology.                                        #
#                                                                              #
# Dependencies: 1) R 3.1.3 (Smooth Sidewalk)                                   #
#               2) Vegan v2.2-1                                                #
#               3) DiversityFunctions.R (2014/09/04; included in Repo;         #
#                    https://github.com/mmuscarella/Diversity_Functions)       #
#               4) MothurTools.R (2015/02/22; included in Repo)                #
#                                                                              #
# Issues: None Identified (yet)                                                #
#         1.                                                                   #
#                                                                              #
# Recent Changes:                                                              #
#         1. I was calculationg selective pressures using RDA and now I'm      #
#            moving to using just spearman's rank correlations. I think the    #
#            story is easier to interpret                                      #
#         2.                                                                   #
#         3.                                                                   #
#                                                                              #
# Future Changes (To-Do List):                                                 #
#         1.                                                                   #
#         2.                                                                   #
#                                                                              #
################################################################################

# Setup Work Environment (Directory, Packages, Source Code, Functions)
rm(list=ls())
setwd("~/GitHub/SubsidyMicrobialStability")
require("vegan")
source("./bin/DiversityFunctions.R")
source("./bin/MothurTools.R")
se <- function(x, ...){sd(x, ...)/sqrt(length(na.omit(x)))}

# Save Standard Plot Settings
opar <- par(no.readonly = TRUE)  # Saves plot defaults

# Define Inputs
# Design       = general design file for experiment
# shared       = OTU table from mothur with sequence similarity clustering
# shared.phylo = OTU table from mothur using phylotypes (Silva Database)
# tax          = Taxonomy for 97% similarity OTUs
# tax.phylo    = Taxonomy for genus level phylotype OTUs
design       <- "./data/Ponds_Design.txt"
shared       <- "./mothur/output/Pond_all.final.an.shared"
shared.phylo <- "./mothur/output/Pond_all.final.FW.tx.shared"
tax          <- "./mothur/output/Pond_all.final.an.0.03.cons.taxonomy"
tax.phylo    <- "./mothur/output/Pond_all.final.FW.tx.1.cons.taxonomy"

# Import Design
design <- read.delim(design, header=T, row.names=1)

# Import Shared Files
Pond97    <- read.otu(shared = shared, cutoff = "0.03")         # 97% Similarity
PondPhylo <- read.otu(shared = shared.phylo, cutoff = "1")   # Genus Level

# Remove OTUs with less than two occurances across all sites
Pond97    <- Pond97[, which(colSums(Pond97) >= 2)]
PondPhylo <- PondPhylo[, which(colSums(PondPhylo) >= 2)]

# Pull Out Pre Data Only
pre         <- rownames(design[which(design$Time == "Pre"),])
design.p    <- design[pre,]
PondPhylo.p <- PondPhylo[pre,]

# Make Relative Abundence Matrices
PondsREL.raw <- PondPhylo.p
  for(i in 1:dim(PondPhylo.p)[1]){
    PondsREL.raw[i,]<-PondPhylo.p[i,]/sum(PondPhylo.p[i,])
	  }

# Sample Normalization
# I am doing a few things here. I am using decostand (in vegan) to
# log transform the abundance data. I am using a relative abund matrix so
# it first divides by the lowest value. I am using this transformation
# to help with the bias against low abundance (rare species).
# See Anderson et al 2006 or Legendre & Legendre 2012 p 327
PondsREL <- decostand(PondsREL.raw, method="log")

# Only DNA
DNA <- rownames(design.p[which(design.p$Molecule == "DNA"),])
PondsREL_DNA <- PondsREL[DNA,]
design_DNA   <- design[DNA,]

# Only RNA
RNA <- rownames(design.p[which(design.p$Molecule == "RNA"),])
PondsREL_RNA <- PondsREL[RNA,]
design_RNA   <- design[RNA,]

# PERMANOVA ####################################################################
adonis(PondsREL ~ design.p$Loading/design.p$Molecule, method = "bray",
       permutations = 999)
adonis(PondsREL_DNA ~ design_DNA$Loading, method = "bray", permutations = 999)
adonis(PondsREL_RNA ~ design_RNA$Loading, method = "bray", permutations = 999)


# DNA Correlations #############################################################
DNA_corr            <- as.data.frame(matrix(NA, dim(PondsREL_DNA)[2], 2))
row.names(DNA_corr) <- colnames(PondsREL_DNA)
colnames(DNA_corr)  <- c("rho", "p.value")

for (i in 1:dim(PondsREL_DNA)[2]){
  out.i <- cor.test(PondsREL_DNA[,i], design_DNA$Loading, method="spearman",
                    exact=FALSE)
  DNA_corr[i,1] <- out.i$estimate
  DNA_corr[i,2] <- out.i$p.value
}

DNA_corr     <- na.omit(DNA_corr)
DNA_corr2    <- DNA_corr[DNA_corr$p.value < 0.05, ]
DNA_corr2$BH <- p.adjust(DNA_corr2$p.value, method="BH")
DNA_corr2

write.csv(DNA_corr2, "./data/DNA_corr.txt")

sum(DNA_corr2$rho > 0)
sum(DNA_corr2$rho < 0)

# RNA Correlations #############################################################
RNA_corr            <- as.data.frame(matrix(NA, dim(PondsREL_RNA)[2], 2))
row.names(RNA_corr) <- colnames(PondsREL_RNA)
colnames(RNA_corr)  <- c("rho", "p.value")

for (i in 1:dim(PondsREL_RNA)[2]){
  out.i <- cor.test(PondsREL_RNA[,i], design_RNA$Loading, method="spearman",
                    exact=FALSE)
  RNA_corr[i,1] <- out.i$estimate
  RNA_corr[i,2] <- out.i$p.value
}

RNA_corr     <- na.omit(RNA_corr)
RNA_corr2    <- RNA_corr[RNA_corr$p.value < 0.05, ]
RNA_corr2$BH <- p.adjust(RNA_corr2$p.value, method="BH")
RNA_corr2

write.csv(RNA_corr2, "./data/RNA_corr.txt")

sum(RNA_corr2$rho > 0)
sum(RNA_corr2$rho < 0)

# Percent of Total Seqs ########################################################
DNA.names <- rownames(DNA_corr2)
RNA.names <- rownames(RNA_corr2)
responders <- union(DNA.names, RNA.names)
total.seqs <- sum(colSums(PondPhylo.p))
responder.matrix <- PondPhylo.p[, responders]
respond.seqs <- sum(colSums(responder.matrix))
respond.seqs/total.seqs
total.otu <- dim(PondPhylo.p)[2]

respond.seqsDNA <- PondPhylo.p[which(design$MolTime == "DNA_Pre"), DNA.names]
total.seqsDNA <- sum(colSums(PondPhylo.p[which(design$MolTime == "DNA_Pre"), ]))
respond.seqsDNA <- sum(colSums(respond.seqsDNA))
respond.seqsDNA / total.seqsDNA

respond.seqsRNA <- PondPhylo.p[which(design$MolTime == "RNA_Pre"), RNA.names]
total.seqsRNA <- sum(colSums(PondPhylo.p[which(design$MolTime == "RNA_Pre"), ]))
respond.seqsRNA <- sum(colSums(respond.seqsRNA))
respond.seqsRNA / total.seqsRNA

# Exmaple Graphs ###############################################################
png(filename="./figures/Pond_Increase.png",
    width = 1400, height = 1200, res = 96*2)
par(opar)

par(mar=c(6,7,1,1))
plot(log(PondsREL.raw[which(design.p$Molecule == "DNA"),"Otu008"])
     ~ design_DNA$Loading, ylim = log(c(0.0008, 0.11)),
     bg="gray", col="black", pch=21, cex=2.5, lwd=2,
     las = 1, axes = FALSE, xlab="", ylab="")
  lines(lowess(log(PondsREL.raw[which(design.p$Molecule == "DNA"),"Otu008"])
               ~ design_DNA$Loading), lwd = 4, lty = 2)
  axis(side = 1, labels=T, lwd.ticks=2, cex.axis=1.5)
  axis(side = 2, lwd.ticks=2, cex.axis=1.5, las = 1,
       labels=c(0.001, 0.01, 0.1), at=log(c(0.001, 0.01, 0.1)))
  axis(side = 1, tck=0.01, labels=F, lwd.ticks=2)
  axis(side = 3, tck=0.01, labels=F, lwd.ticks=2)
  axis(side = 2, tck=0.01, labels=F, lwd.ticks=2, at=log(c(0.001, 0.01, 0.1)))
  axis(side = 4, tck=0.01, labels=F, lwd.ticks=2, at=log(c(0.001, 0.01, 0.1)))
  axis(side = 2, tck=-0.015, labels=F, lwd.ticks=2,
       at=log(c(seq(0.002, 0.009, 0.001), seq(0.01, 0.1, 0.01))))
  mtext(expression(paste("tDOC Supply Rate (g m"^"-2",")")),
        side=1, line=4, cex=2)
  mtext("Relative Abundance", side=2, line=4, cex=2)
  legend("topleft", expression("OTU008 ("*italic("Methylomonas sp.")*")"), bty = 'n', cex = 1.5, x.intersp = -0.5)
  box(lwd=2)
dev.off() # this writes plot to folder
graphics.off() # shuts down open devices 
par(opar)

png(filename="./figures/Pond_Decrease.png",
    width = 1400, height = 1200, res = 96*2)
par(opar)

par(mar=c(6,7,1,1))
plot(log(PondsREL.raw[which(design.p$Molecule == "DNA"),"Otu011"])
     ~ design_DNA$Loading, ylim = log(c(0.0008, 0.175)),
     bg="gray", col="black", pch=21, cex=2.5, lwd=2,
     las = 1, axes = FALSE, xlab="", ylab="")
  lines(lowess(log(PondsREL.raw[which(design.p$Molecule == "DNA"),"Otu011"])
               ~ design_DNA$Loading), lwd = 4, lty = 2)
  axis(side = 1, labels=T, lwd.ticks=2, cex.axis=1.5)
  axis(side = 2, lwd.ticks=2, cex.axis=1.5, las = 1,
       labels=c(0.001, 0.01, 0.1), at=log(c(0.001, 0.01, 0.1)))
  axis(side = 1, tck=0.01, labels=F, lwd.ticks=2)
  axis(side = 3, tck=0.01, labels=F, lwd.ticks=2)
  axis(side = 2, tck=0.01, labels=F, lwd.ticks=2, at=log(c(0.001, 0.01, 0.1)))
  axis(side = 4, tck=0.01, labels=F, lwd.ticks=2, at=log(c(0.001, 0.01, 0.1)))
  axis(side = 2, tck=-0.015, labels=F, lwd.ticks=2,
       at=log(c(seq(0.002, 0.009, 0.001), seq(0.01, 0.1, 0.01), seq(0.1, 1, 0.1))))
  mtext(expression(paste("tDOC Supply Rate (g m"^"-2",")")),
        side=1, line=4, cex=2)
  mtext("Relative Abundance", side=2, line=4, cex=2)
  legend("topright", expression("OTU011 ("*italic("Rhodococcus sp.")*")"), bty = 'n', cex = 1.5)
  box(lwd=2)
dev.off() # this writes plot to folder
graphics.off() # shuts down open devices 
par(opar)
