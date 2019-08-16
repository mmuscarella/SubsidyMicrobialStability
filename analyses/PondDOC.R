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


################################################################################
#### BGE - from Plos Paper ##############################################
################################################################################
DOC[which(DOC$Date == "2009-08-19"), ]
bge <- c(0.560, 0.467, 0.534, 0.499, 0.486, 0.503, 0.497, 0.462, 0.381, 0.424)
mean.supply <- aggregate(DOC$ppm, by = list(DOC$Pond), FUN = mean)
doc <- c(692, 781, 823, 950, 959, 1114, 1362, 1540, 1460, 1756)
pond <- c(10, 9, 3, 4, 8, 15, 7, 13, 16, 14)
loading <- c(0, 55.72, 25.72, 115.2, 84.12, 133.41, 151.63, 184.31, 201.46, 212.71)

model.bge <- lm(bge ~ loading)

png(filename="./figures/Pond_BGE.png",
    width = 1400, height = 1200, res = 96*2)
par(opar)
par(mar=c(6,6,1,1))

plot(loading, bge, type='n',
     ylim=c(0.35, 0.6), xlim=c(0, 225), axes = FALSE, xlab="", ylab="")

pred.frame <- data.frame(loading = seq(0,225,5))

CI.U <- predict(model.bge, interval = "c", newdata=pred.frame)[, "upr"]
CI.L <- predict(model.bge, interval = "c", newdata=pred.frame)[, "lwr"]
pred.frame2 <- unlist(pred.frame)

X.Vec <- c(pred.frame2, tail(pred.frame2, 1), rev(pred.frame2),
               head(pred.frame2, 1))
Y.Vec <- c(CI.U, tail(CI.L, 1), rev(CI.L), head(CI.U,1))
polygon(X.Vec, Y.Vec, col = "gray90", border = NA)

matlines(pred.frame, predict(model.bge, interval = "c", newdata=pred.frame),
         lty=c(2,3,3), lwd=c(4,3,3), col="black")

points(loading, bge, bg="gray", col="black",
       pch=21, cex=2.5, lwd=2)

axis(side = 1, labels=T, lwd.ticks=2, cex.axis=1.5)
axis(side = 2, labels=T, lwd.ticks=2, at=seq(0,0.5,0.1), las=1, cex.axis=1.5)
axis(side = 1, tck=0.01, labels=F, lwd.ticks=2)
axis(side = 2, tck=0.01, labels=F, lwd.ticks=2, at=seq(0,0.5,0.1))
axis(side = 3, tck=0.01, labels=F, lwd.ticks=2)
axis(side = 4, tck=0.01, labels=F, lwd.ticks=2, at=seq(0,0.5,0.1))
mtext(expression(paste("tDOC Supply Rate (g m"^"-2",")")),
      side=1, line=4, cex=2)
mtext("Growth Efficiency", side=2, line=3, cex=2)

box(lwd=2)

dev.off() # this writes plot to folder
graphics.off() # shuts down open devices 
par(opar)
