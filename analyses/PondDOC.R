################################################################################
#                                                                              #
#	Pond Microbial Stability: DOC Analysis and Graph                             #
#                                                                              #
################################################################################
#                                                                              #
#	Written by: Mario Muscarella                                                 #
#                                                                              #
#	Last update: 2014/08/26                                                      #
#                                                                              #
################################################################################

# Setup Work Environment
rm(list=ls())
setwd("~./GitHub/pond-microbes/")

# Import Data
doc.data <- read.delim("./data/PondOrganicCarbon_10-23-09.txt", header=T)
design   <- read.delim("./data/Pond_DOC_Loading.txt")

doc.data2 <- subset(doc.data, doc.data$Total.Dissolved == "Dissolved")
doc.data2$Date <- as.Date(strptime(doc.data2$Date, format="%m/%d/%y"))

day60s <- min(which(doc.data2$Date == "2009-08-05"))
day100e <- max(which(doc.data2$Date == "2009-09-14"))

doc.data3 <- doc.data2[day60s:day100e,c(2,1,5)]
colnames(doc.data3) <- c("Pond", "Date", "ppm")
doc.data3$Pond <- as.character(doc.data3$Pond)

design <- design[,2:3]
colnames(design) <- c("Pond", "loading")

DOC <- merge(doc.data3, design, by="Pond")

# Plot
windows.options(width=8, height=4)
par(mar=c(3,5,0.5,0.5))
#par(mfrow=c(3,2), mar=c(0.25,4.5,0.25,1), oma=c(3,1,1,1)+0.1 )
boxplot(DOC$ppm ~ DOC$loading, axes = FALSE)
box(lwd=2)
axis(side=1, labels=F, tick=F)
axis(side=2, labels=T, las=2)
mtext("Pond (ordered by tDOC supply rate)", side=1, line=1, cex=1.5)
mtext(expression(paste("DOC (mg C L"^" -1",")", sep="")), side=2, line=2.5, cex=1.5)


png(file="./figures/DOC_Loading.png", width=800, height=400, antialias = "cleartype")
windows.options(width=8, height=4)
par(mar=c(3,5,0.5,0.5))
#par(mfrow=c(3,2), mar=c(0.25,4.5,0.25,1), oma=c(3,1,1,1)+0.1 )
boxplot(DOC$ppm ~ DOC$loading, axes = FALSE)
box(lwd=2)
axis(side=1, labels=F, tick=F)
axis(side=2, labels=T, las=2)
mtext("Pond (by tDOC supply rate)", side=1, line=1, cex=1.5)
mtext(expression(paste("DOC (mg C L"^" -1",")", sep="")), side=2, line=2.5, cex=1.5)
dev.off()


