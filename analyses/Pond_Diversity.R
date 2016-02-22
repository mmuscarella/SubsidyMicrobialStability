################################################################################
#                                                                              #
#	Pond Microbe Stability Experiment: Microbial Community Diversity Graphs      #
#                                                                              #
################################################################################
#                                                                              #
#	Written by: Mario Muscarella                                                 #
#                                                                              #
#	Last update: 2015/06/24                                                      #
#                                                                              #
################################################################################
#                                                                              #
# Notes: This code conducts the calculations needed for microbial diversity    #
#        estimation as part of the Pond Experiment (2009). For general         #
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
#         1. Evenness linear model seems inapprppriate.                        #
#                                                                              #
# Recent Changes:                                                              #
#         1. Change from Peilou's Evenness to Simpson's Evenness               #
#         2. New Vegan Package                                                 #
#         3. Finalizing Figures                                                #
#                                                                              #
# Future Changes (To-Do List):                                                 #
#         1. Add Phylogenetic Diversity                                        #
#         2. Add Phylogenetic Signals                                          #
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
Pond97 <- read.otu(shared = shared, cutoff = "0.03")         # 97% Similarity
PondPhylo <- read.otu(shared = shared.phylo, cutoff = "1")   # Genus Level

# Remove OTUs with less than two occurances across all sites
Pond97 <- Pond97[, which(colSums(Pond97) >= 2)]
PondPhylo <- PondPhylo[, which(colSums(PondPhylo) >= 2)]

# Coverage Stats (the two should be the same, notice 9&3 PrecDNA are <1000)
count.groups(Pond97)
count.groups(PondPhylo)

# Pull Out Pre Data Only
pre <- rownames(design[which(design$Time == "Pre"),])
design.p <- design[pre,]
Pond97.p <- Pond97[pre,]
PondPhylo.p <- PondPhylo[pre,]

# Alpha Diversity with Resampling
rich <- round(richness.iter(input = PondPhylo.p, size = 2000,
              iters = 1000, shared="FALSE"), 3)
even <- round(evenness.iter(input = PondPhylo.p, size = 2000,
              iters = 1000, shared="FALSE", method="simp_even"), 3)
# shan <- round(diversity.iter(input = PondPhylo.p, size = 1000,
#               iters = 1000, shared="FALSE", method="shan"), 3)

# Richness Summary
rich_data        <- merge(design.p, rich, by="row.names")
rich_data$mean   <- round(apply(rich_data[7:1006], 1, mean, na.rm = TRUE),3)
rich_data$se     <- round(apply(rich_data[7:1006], 1, se, na.rm = TRUE), 3)

# Evenness Summary
even_data        <- merge(design.p, even, by="row.names")
even_data$mean   <- round(apply(even_data[7:1006], 1, mean, na.rm = TRUE),3)
even_data$se     <- round(apply(even_data[7:1006], 1, se, na.rm = TRUE),3)

# # Shannon Diversity Summary
# shan_data        <- merge(design.p, shan, by="row.names")
# shan_data$mean   <- round(apply(shan_data[7:105], 1, mean, na.rm = TRUE),3)
# shan_data$se     <- round(apply(shan_data[7:105], 1, se, na.rm = TRUE),3)

# Linear Models
Rich_Model <- lm(mean ~ Molecule/Loading, data = rich_data)
Even_Model <- lm(mean ~ Molecule/Loading, data = even_data)

summary(Rich_Model)
summary(Even_Model)


################################################################################
#### Alpha Diversity Plots #####################################################
################################################################################

################################################################################
# Pond Richness Plot
# DNA and RNA ##################################################################

rich_data_DNA <- subset(rich_data, Molecule=="DNA")
rich_data_RNA <- subset(rich_data, Molecule=="RNA")
model.dna <- lm(mean ~ Loading, data=rich_data_DNA)
model.rna <- lm(mean ~ Loading, data=rich_data_RNA)

fill.color <- ifelse(rich_data$Molecule=="DNA", "gray", "black")
symbol <- ifelse(rich_data$Molecule=="DNA", 21, 17)

png(filename="./figures/Pond_Richness.png",
    width = 1400, height = 1200, res = 96*2)
par(opar)

par(mar=c(6,6,1,1))
plot(rich_data$Loading, rich_data$mean, type='n',
     ylim=c(20, 95), xlim=c(0, 225), axes = FALSE, xlab="", ylab="")
pred.frame <- data.frame(Loading = seq(0,225,5))
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
         lty=c(2,3,3), lwd=c(4,2,2), col="black")
matlines(pred.frame, predict(model.dna, interval = "c", newdata=pred.frame),
         lty=c(2,3,3), lwd=c(4,2,2), col="black")
points(rich_data$Loading, rich_data$mean,
       bg=fill.color, col="black", pch=symbol, cex=2.5, lwd=2)
axis(side = 1, labels=T, lwd.ticks=2, cex=1.75)
axis(side = 2, labels=T, lwd.ticks=2, at=seq(0,100,20), las=1, cex=1.75)
axis(side = 1, tck=0.01, labels=F, lwd.ticks=2)
axis(side = 2, tck=0.01, labels=F, lwd.ticks=2, at=seq(0,100,20))
axis(side = 3, tck=0.01, labels=F, lwd.ticks=2)
axis(side = 4, tck=0.01, labels=F, lwd.ticks=2, at=seq(0,100,20))
mtext(expression(paste("tDOC Supply Rate (g m"^"-2",")")),
      side=1, line=4, cex=2)
mtext("Taxa Richness", side=2, line=3, cex=2)
legend("bottomleft", c("Total", "Active"), col="black",
       pt.bg=c("gray", "black"), pch=c(21, 24), bty='n', cex=1.5)
box(lwd=2)

dev.off() # this writes plot to folder
graphics.off() # shuts down open devices 
par(opar)

################################################################################
# Pond Evenness Plot
# DNA and RNA ##################################################################

even_data_DNA <- subset(even_data, Molecule=="DNA")
even_data_RNA <- subset(even_data, Molecule=="RNA")
even_data_RNA <- subset(even_data_RNA, Loading!="201")
model.dna <- lm(mean ~ Loading, data=even_data_DNA)
model.rna <- lm(mean ~ Loading, data=even_data_RNA)

fill.color <- ifelse(even_data$Molecule=="DNA", "gray", "black")
symbol <- ifelse(even_data$Molecule=="DNA", 21, 17)

png(filename="./figures/Pond_Evenness.png",
    width = 1400, height = 1200, res = 96*2)
par(opar)

par(mar=c(6,6,1,1))
plot(even_data$Loading, even_data$mean, type='n',
     ylim=c(0, 0.38), xlim=c(0, 225), axes = FALSE, xlab="", ylab="")
pred.frame <- data.frame(Loading = seq(0,225,5))
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
polygon(X.Vec_rna, Y.Vec_rna, col = "gray90", border = NA)
polygon(X.Vec_dna, Y.Vec_dna, col = "gray90", border = NA)
matlines(pred.frame, predict(model.rna, interval = "c", newdata=pred.frame),
         lty=c(2,3,3), lwd=c(4,2,2), col="black")
matlines(pred.frame, predict(model.dna, interval = "c", newdata=pred.frame),
         lty=c(2,3,3), lwd=c(4,2,2), col="black")
points(even_data$Loading, even_data$mean,
       bg=fill.color, col="black", pch=symbol, cex=2.5, lwd=2)
axis(side = 1, labels=T, lwd.ticks=2, cex=1.75)
axis(side = 2, labels=T, lwd.ticks=2, at=seq(0,0.5,0.1), las=1, cex=1.75)
axis(side = 1, tck=0.01, labels=F, lwd.ticks=2)
axis(side = 2, tck=0.01, labels=F, lwd.ticks=2, at=seq(0,0.5,0.1))
axis(side = 3, tck=0.01, labels=F, lwd.ticks=2)
axis(side = 4, tck=0.01, labels=F, lwd.ticks=2, at=seq(0,0.5,0.1))
mtext(expression(paste("tDOC Supply Rate (g m"^"-2",")")),
      side=1, line=4, cex=2)
mtext("Simpson's Evenness", side=2, line=3, cex=2)
legend("topleft", c("Total", "Active"), col="black",
       pt.bg=c("gray", "black"), pch=c(21, 24), bty='n', cex=1.5)
box(lwd=2)

dev.off() # this writes plot to folder
graphics.off() # shuts down open devices 
par(opar)

################################################################################
# Pond Diversity Plot
# DNA and RNA ##################################################################

# shan_data_DNA <- subset(shan_data, Molecule=="DNA")
# shan_data_RNA <- subset(shan_data, Molecule=="RNA")
# model.dna <- lm(mean ~ Loading, data=shan_data_DNA)
# model.rna <- lm(mean ~ Loading, data=shan_data_RNA)
#
# fill.color <- ifelse(shan_data$Molecule=="DNA", "gray", "black")
# symbol <- ifelse(shan_data$Molecule=="DNA", 21, 17)
# windows.options(width=6, height=6)
# par(mar=c(5,5,1,1))
# plot(rich_data$Loading, rich_data$mean, type='n', ylim=c(0, 6),
#      xlim=c(0, 225), axes = FALSE, xlab="", ylab="")
# pred.frame <- data.frame(Loading = seq(0,225,5))
# CI.U_dna <- predict(model.dna, interval = "c", newdata=pred.frame)[, "upr"]
# CI.L_dna <- predict(model.dna, interval = "c", newdata=pred.frame)[, "lwr"]
# CI.U_rna <- predict(model.rna, interval = "c", newdata=pred.frame)[, "upr"]
# CI.L_rna <- predict(model.rna, interval = "c", newdata=pred.frame)[, "lwr"]
# pred.frame2 <- unlist(pred.frame)
# X.Vec_dna <- c(pred.frame2, tail(pred.frame2, 1), rev(pred.frame2),
#                head(pred.frame2, 1))
# Y.Vec_dna <- c(CI.U_dna, tail(CI.L_dna, 1), rev(CI.L_dna), head(CI.U_dna,1))
# X.Vec_rna <- c(pred.frame2, tail(pred.frame2, 1), rev(pred.frame2),
#                head(pred.frame2, 1))
# Y.Vec_rna <- c(CI.U_rna, tail(CI.L_rna, 1), rev(CI.L_rna), head(CI.U_rna,1))
# polygon(X.Vec_dna, Y.Vec_dna, col = "gray80", border = NA)
# polygon(X.Vec_rna, Y.Vec_rna, col = "gray90", border = NA)
# matlines(pred.frame, predict(model.rna, interval = "c", newdata=pred.frame),
#          lty=c(2,3,3), lwd=c(4,3,3), col="black")
# matlines(pred.frame, predict(model.dna, interval = "c", newdata=pred.frame),
#          lty=c(2,3,3), lwd=c(4,3,3), col="black")
# points(shan_data$Loading, shan_data$mean, bg=fill.color, col="black",
#        pch=symbol, cex=2.5, lwd=2)
# axis(side = 1, labels=T, lwd.ticks=2, cex=1.5)
# axis(side = 2, labels=T, lwd.ticks=2, at=seq(0,8,2), las=1, cex=1.5)
# axis(side = 1, tck=0.01, labels=F, lwd.ticks=2)
# axis(side = 2, tck=0.01, labels=F, lwd.ticks=2, at=seq(0,8,2))
# axis(side = 3, tck=0.01, labels=F, lwd.ticks=2)
# axis(side = 4, tck=0.01, labels=F, lwd.ticks=2, at=seq(0,8,2))
# mtext(expression(paste("Carbon Loading (g" %.% "m"^"-2",")")),
#       side=1, line=3.5, cex=2)
# mtext("Shannon Diversity", side=2, line=3, cex=2)
# legend("bottomleft", c("Total", "Active"), col="black",
#        pt.bg=c("gray", "black"), pch=c(21, 24), bty='n', cex=1.5)
# box(lwd=2)

################################################################################
# Summary Table
# All Plots ####################################################################

# Linear Models
Rich_Model <- lm(mean ~ Molecule/Loading + -1, data = rich_data)
Even_Model <- lm(mean ~ Molecule/Loading + -1, data = even_data)

summary(Rich_Model)
summary(Even_Model)


# Linear Models
Rich_Model2 <- lm(mean ~ Molecule/Loading, data = rich_data)
Even_Model2 <- lm(mean ~ Molecule/Loading, data = even_data)

summary(Rich_Model2)
summary(Even_Model2)
