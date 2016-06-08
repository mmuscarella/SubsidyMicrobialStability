################################################################################
#                                                                              #
#  Pond Microbe Stability Experiment: Microbial Community Stability            #
#                                     Analysis Using Euclidiean distance PCoA  #
#                                     (See Brown 2003 Ecology Letters)         #
#                                                                              #
################################################################################
#                                                                              #
#	Written by: Mario Muscarella                                                 #
#                                                                              #
#	Last update: 2015/05/18                                                      #
#                                                                              #
################################################################################
#                                                                              #
# Notes: This code conducts the calculations needed for microbial stability    #
#        and responsivenses as part of the Pond Experiment (2009). For general #
#        information about this project see: Lennon et al 2013 Plos One,       #
#        Jones and Lennon (2015).                                              #
#                                                                              #
# Dependencies: 1) Vegan v2.2-1                                                #
#               2) DiversityFunctions.R (2014/09/04)                           #
#                                                                              #
# Issues: None Identified (yet)                                                #
#         1.                                                                   #
#                                                                              #
# Recent Changes:                                                              #
#         1. Now using Phylotype instead of 97% Sim                            #
#         2. New Vegan Package                                                 #
#         3. Added Diversity Estimation                                        #
#         4.                                                                   #
#                                                                              #
# Future Changes (To-Do List):                                                 #
#         1. Check phylotype designations                                      #
#         2.                                                                   #
#                                                                              #
################################################################################

# Setup Work Environment
rm(list=ls())
setwd("~/GitHub/SubsidyMicrobialStability")
source("./bin/DiversityFunctions.r")
require("vegan")
se <- function(x, ...){sd(x, ...)/sqrt(length(na.omit(x)))}

# Save Standard Plot Settings
opar <- par(no.readonly = TRUE)  # Saves plot defaults

# Define Inputs
# Design = general design file for experiment
# shared = OTU table from mothur with sequence similarity clustering
# shared.phylo = OTU table from mothur using phylotypes (Freshwater Database)
design <- "./data/Ponds_Design.txt"
shared <- "./mothur/output/Pond_all.final.an.shared"
shared.phylo <- "./mothur/output/Pond_all.final.FW.tx.shared"

# Import Design
design <- read.delim(design, header=T, row.names=1)

# Import Shared Files
Pond97 <- read.otu(shared = shared, cutoff = "0.03")         # 97% Similarity
PondPhylo <- read.otu(shared = shared.phylo, cutoff = "1")   # Genus Level

# Remove OTUs with less than two occurances across all sites
Pond97 <- Pond97[, which(colSums(Pond97) >= 2)]
PondPhylo <- PondPhylo[, which(colSums(PondPhylo) >= 2)]

# Coverage Stats (the two should be the same, notice 9&3 PrecDNA are <1000)
count.groups(Pond97)
count.groups(PondPhylo)

# Make Presence Absence Matrices
PondsPA <- (PondPhylo > 0) * 1

# Make Relative Abundence Matrices
PondsREL <- PondPhylo
  for(i in 1:44){
	  PondsREL[i,]<-PondPhylo[i,]/sum(PondPhylo[i,])
	  }

# I am doing a few things here.
# I am using decostand (in vegan) to log transform the abundance data.
# I am using a rel abund matrix so it first divides by the lowest value.
# I am using this transformation to help with the bias against low abundance
# (rare species). See Anderson et al 2006 or Legendre & Legendre 2012 p 327
PondsREL.dist <- vegdist(decostand(PondsREL, method="log"),method="bray")

# I am using cmdscale (vegan) to do PcoA.
PondsREL.pcoa <- cmdscale(PondsREL.dist, eig=T)
per.exp <- PondsREL.pcoa$eig/sum(PondsREL.pcoa$eig)
sum(per.exp[1:6])

# Here I am just using 3 axes to explain the data (explains 70% of variation)
PondsREL.pcoa <- cmdscale(PondsREL.dist, k=6)
Pond_Pcoa <- as.data.frame(PondsREL.pcoa)
Pond_Pcoa <- merge(Pond_Pcoa, design, by="row.names")

# I am using metaMDS (vegan) to do NonMetric Multi Demensional Scalling.
# Here I am just using 5 axes to explain the data
# PondsREL.nmds <- metaMDS(PondsREL.dist, k=5)
# Pond_nmds <- as.data.frame(PondsREL.pcoa)
# Pond_nmds <- merge(Pond_nmds, design, by="row.names")

# PCoA Plot if you are interested
# par(mar=c(5,6,2,2))
# fill.color <- ifelse(Pond_Pcoa$Molecule=="DNA", "gray", "black")
# symbol <- ifelse(Pond_Pcoa$Time=="Pre", 24, 21)
#
# plot(PondsREL.pcoa,  bg=fill.color, col="black", pch=symbol, cex=1.5,
#      cex.lab=1.5, cex.axis=1.2, xlab="PCoA 1", ylab="PCoA 2")
# text(PondsREL.pcoa, labels=Pond_Pcoa$Pond, pos=3)
# legend(-0.45, -0.10, c("", ""), pch=c(24), pt.bg=c("gray", "black"),
#        col="black", bty='n', cex=1.5, horiz=T)
# legend(-0.45, -0.12, c("", ""), pch=c(21), pt.bg=c("gray", "black"),
#        col="black", bty='n', cex=1.5, horiz=T)
# text(-0.46, -0.12, "Pre")
# text(-0.46, -0.14, "Post")
# text(-0.42, -0.10, "Total")
# text(-0.34, -0.10, "Active")
# box(lwd=2)


# Here I am calculating the Euclidean distance
# and I am making the output a symetrical matrix
matrixPcoA <- as.matrix(PondsREL.pcoa)
distPcoA <- as.matrix(vegdist(matrixPcoA, method="euclidean", diag=T, upper=T))

# I now need to extract the distances to represent pre and post perturbation
# for each pond/molecule type
Pond10DNA <- distPcoA[3,1]; Pond10cDNA <- distPcoA[4,2]
Pond13DNA <- distPcoA[7,5]; Pond13cDNA <- distPcoA[8,6]
Pond14DNA <- distPcoA[11,9]; Pond14cDNA <- distPcoA[12,10]
Pond15DNA <- distPcoA[15,13]; Pond15cDNA <- distPcoA[16,14]
Pond16DNA <- distPcoA[19,17]; Pond16cDNA <- distPcoA[20,18]
Pond3DNA <- distPcoA[23,21]; Pond3cDNA <- distPcoA[24,22]
Pond4DNA <- distPcoA[27,25]; Pond4cDNA <- distPcoA[28,26]
Pond5DNA <- distPcoA[31,29]; Pond5cDNA <- distPcoA[32,30]
Pond7DNA <- distPcoA[35,33]; Pond7cDNA <- distPcoA[36,34]
Pond8DNA <- distPcoA[39,37]; Pond8cDNA <- distPcoA[40,38]
Pond9DNA <- distPcoA[43,41]; Pond9cDNA <- distPcoA[44,42]

# Now I am going to organize the data
StabDNA <- c(Pond3DNA, Pond4DNA, Pond5DNA, Pond7DNA, Pond8DNA, Pond9DNA,
             Pond10DNA, Pond13DNA, Pond14DNA, Pond15DNA, Pond16DNA)
StabRNA <- c(Pond3cDNA, Pond4cDNA, Pond5cDNA, Pond7cDNA, Pond8cDNA, Pond9cDNA,
             Pond10cDNA, Pond13cDNA, Pond14cDNA, Pond15cDNA, Pond16cDNA)
Stab <- c(StabDNA, StabRNA)
PondLoad <- c(25.72, 115.2, 0, 151.63, 84.12, 55.72, 0, 184.31,
              201.46, 133.41, 212.71)
PondLoad2 <- c(PondLoad, PondLoad)
Molecule <- c(rep("DNA", 11), rep("RNA", 11))

# Linear Models
Stab_Model <- lm(Stab ~ Molecule/PondLoad2)
summary(Stab_Model)

# New Matrix
stab_data <- as.data.frame(cbind(stab = round(Stab, 2), load = PondLoad2))
stab_data$molecule <- Molecule

stab_data_DNA <- subset(stab_data, molecule=="DNA")
stab_data_RNA <- subset(stab_data, molecule=="RNA")
model.dna <- lm(stab ~ load, data=stab_data_DNA)
model.rna <- lm(stab ~ load, data=stab_data_RNA)

################################################################################
#### Stability Plot ############################################################
################################################################################

fill.color <- ifelse(stab_data$molecule=="DNA", "gray", "black")
symbol <- ifelse(stab_data$molecule=="DNA", 21, 17)
windows.options(width=6, height=6)
par(mar=c(5,5,1,1))

plot(stab_data$load, (1 - stab_data$stab), type='n',
     ylim=c(0.55, 1.0), xlim=c(0, 225), axes = FALSE, xlab="", ylab="")

pred.frame <- data.frame(load = seq(0,225,5))
CI.U_dna <- predict(model.dna, interval = "c", newdata=pred.frame)[, "upr"]
CI.L_dna <- predict(model.dna, interval = "c", newdata=pred.frame)[, "lwr"]
CI.U_rna <- predict(model.rna, interval = "c", newdata=pred.frame)[, "upr"]
CI.L_rna <- predict(model.rna, interval = "c", newdata=pred.frame)[, "lwr"]
pred.frame2 <- unlist(pred.frame)
X.Vec_dna <- c(pred.frame2, tail(pred.frame2, 1), rev(pred.frame2),
               head(pred.frame2, 1))
Y.Vec_dna <- c(CI.U_dna, tail(CI.L_dna, 1), rev(CI.L_dna), head(CI.U_dna,1))
X.Vec_rna <- c(pred.frame2, tail(pred.frame2, 1), rev(pred.frame2),
               head(pred.frame2, 1))
Y.Vec_rna <- c(CI.U_rna, tail(CI.L_rna, 1), rev(CI.L_rna), head(CI.U_rna,1))
polygon(X.Vec_dna, 1 - Y.Vec_dna, col = "gray90", border = NA)
polygon(X.Vec_rna, 1 - Y.Vec_rna, col = "gray90", border = NA)

matlines(pred.frame, 1 - predict(model.rna, interval = "c", newdata=pred.frame),
         lty=c(2,3,3), lwd=c(4,3,3), col="black")
matlines(pred.frame, 1 - predict(model.dna, interval = "c", newdata=pred.frame),
         lty=c(2,3,3), lwd=c(4,3,3), col="black")

points(stab_data$load, (1 - stab_data$stab), bg=fill.color, col="black",
       pch=symbol, cex=2.5, lwd=2)

axis(side = 1, labels=T, lwd.ticks=2, cex=1.5)
axis(side = 2, labels=T, lwd.ticks=2, at=seq(0.5,1.0,0.1), las=1, cex=1.5)
axis(side = 1, tck=0.01, labels=F, lwd.ticks=2)
axis(side = 2, tck=0.01, labels=F, lwd.ticks=2, at=seq(0.5,1.0,0.1))
axis(side = 3, tck=0.01, labels=F, lwd.ticks=2)
axis(side = 4, tck=0.01, labels=F, lwd.ticks=2, at=seq(0.5,1.0,0.1))
mtext(expression(paste("Carbon Loading (g" %.% "m"^"-2",")")),
      side=1, line=3.5, cex=2)
mtext("Stability (1 - EDCoA)", side=2, line=3, cex=2)
legend("bottomright", c("Total", "Active"), col="black",
       pt.bg=c("gray", "black"), pch=c(21, 24), bty='n', cex=1.5)
box(lwd=2)


################################################################################
#### Responsiveness Plot #######################################################
################################################################################

fill.color <- ifelse(stab_data$molecule=="DNA", "gray", "black")
symbol <- ifelse(stab_data$molecule=="DNA", 21, 17)

png(filename="./figures/Pond_Responsiveness.png",
    width = 1400, height = 1200, res = 96*2)
par(opar)
par(mar=c(6,6,1,1))

plot(stab_data$load, (stab_data$stab), type='n',
     ylim=c(0, 0.45), xlim=c(0, 225), axes = FALSE, xlab="", ylab="")

pred.frame <- data.frame(load = seq(0,225,5))
CI.U_dna <- predict(model.dna, interval = "c", newdata=pred.frame)[, "upr"]
CI.L_dna <- predict(model.dna, interval = "c", newdata=pred.frame)[, "lwr"]
CI.U_rna <- predict(model.rna, interval = "c", newdata=pred.frame)[, "upr"]
CI.L_rna <- predict(model.rna, interval = "c", newdata=pred.frame)[, "lwr"]
pred.frame2 <- unlist(pred.frame)
X.Vec_dna <- c(pred.frame2, tail(pred.frame2, 1), rev(pred.frame2),
               head(pred.frame2, 1))
Y.Vec_dna <- c(CI.U_dna, tail(CI.L_dna, 1), rev(CI.L_dna), head(CI.U_dna,1))
X.Vec_rna <- c(pred.frame2, tail(pred.frame2, 1), rev(pred.frame2),
               head(pred.frame2, 1))
Y.Vec_rna <- c(CI.U_rna, tail(CI.L_rna, 1), rev(CI.L_rna), head(CI.U_rna,1))
polygon(X.Vec_dna, Y.Vec_dna, col = "gray90", border = NA)
polygon(X.Vec_rna, Y.Vec_rna, col = "gray90", border = NA)

matlines(pred.frame, predict(model.rna, interval = "c", newdata=pred.frame),
         lty=c(2,3,3), lwd=c(4,3,3), col="black")
matlines(pred.frame, predict(model.dna, interval = "c", newdata=pred.frame),
         lty=c(2,3,3), lwd=c(4,3,3), col="black")

points(stab_data$load, stab_data$stab, bg=fill.color, col="black",
       pch=symbol, cex=2.5, lwd=2)

axis(side = 1, labels=T, lwd.ticks=2, cex.axis=1.5)
axis(side = 2, labels=T, lwd.ticks=2, at=seq(0,0.5,0.1), las=1, cex.axis=1.5)
axis(side = 1, tck=0.01, labels=F, lwd.ticks=2)
axis(side = 2, tck=0.01, labels=F, lwd.ticks=2, at=seq(0,0.5,0.1))
axis(side = 3, tck=0.01, labels=F, lwd.ticks=2)
axis(side = 4, tck=0.01, labels=F, lwd.ticks=2, at=seq(0,0.5,0.1))
mtext(expression(paste("tDOC Supply Rate (g m"^"-2",")")),
      side=1, line=4, cex=2)
mtext("Responsiveness", side=2, line=3, cex=2)
legend("topright", c("Total", "Active"), col="black",
       pt.bg=c("gray", "black"), pch=c(21, 24), bty='n', cex=1.5)
box(lwd=2)

dev.off() # this writes plot to folder
graphics.off() # shuts down open devices 
par(opar)
