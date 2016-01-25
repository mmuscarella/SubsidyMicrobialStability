################################################################################
#                                                                              #
#	Pond Microbial Stability: Nutrient Time Series Graph                         #
#                                                                              #
################################################################################
#                                                                              #
#	Written by: Mario Muscarella                                                 #
#                                                                              #
#	Last update: 2014/09/17                                                      #
#                                                                              #
################################################################################

# Setup Work Environment
rm(list=ls())
setwd("~/GitHub/SubsidyMicrobialStability/")
require(reshape)

# Import Data
TN.data <- read.delim("./data/PondNitrogen_10-23-09.txt", header=T)
TP.data <- read.delim("./data/PondPhosphorus_10-23-09.txt", header=T)
SRP.data <- read.delim("./data/PondPhosphorus_10-23-09.txt", header=T)
design  <- read.delim("./data/Pond_DOC_Loading.txt", header=T)

# Subset/Format Data
TN.data2 <- subset(TN.data, TN.data$Total.Dissolved == "Total")
TN.data2$Date <- as.Date(strptime(TN.data2$Date, format="%m/%d/%y"))

TP.data2 <- subset(TP.data, TP.data$Total.Dissolved == "Total")
TP.data2$Date <- as.Date(strptime(TP.data2$Date, format="%m/%d/%Y"))

SRP.data2 <- subset(SRP.data, SRP.data$Analyte == "SRP")
SRP.data2$Date <- as.Date(strptime(SRP.data2$Date, format="%m/%d/%Y"))

# Organize Data
TN <- TN.data2[,c(2,1,5)]
TP <- TP.data2[,c(2,1,5)]
SRP <- SRP.data2[,c(2,1,5)]
colnames(TN) <- c("Pond", "Date", "ppm")
colnames(TP) <- c("Pond", "Date", "ppm")
colnames(SRP) <- c("Pond", "Date", "ppm")
TN <- TN[order(TN[,2]),]
TP <- TP[order(TP[,2]),]
SRP <- SRP[order(SRP[,2]),]

# SRP Calculations
SRP.m <- melt(SRP, id.vars=c("Pond","Date"), measure.vars="ppm")
SRP.w <- na.omit(as.data.frame(cast(SRP.m, ... ~ Pond)))
SRP.w[,3:13][SRP.w[,3:13] < 0] <- 0
SRP.w <- SRP.w[rowMeans(SRP.w[,3:13]) > 0,]

SRP.w$mn <- apply(SRP.w[,3:13], 1, mean, na.rm=T)
SRP.w$mid <- apply(SRP.w[,3:13], 1, quantile, na.rm=T)[3,]
SRP.w$low <- apply(SRP.w[,3:13], 1, quantile, na.rm=T, type=5)[2,]
SRP.w$high <- apply(SRP.w[,3:13], 1, quantile, na.rm=T, type=5)[4,]

# SRP Timeseries Plot
windows.options(width=8, height=2.5)
windows()
par(mfrow=c(1,1), mar=c(3,5,0.5,0.5))
plot(SRP.w$mn ~ SRP.w$Date,
     type='n', lwd=2, xlab="Time", xaxt='n',
     yaxt='n', ylab="", ylim=c(0,60))
X.Vec <- c(SRP.w$Date, tail(SRP.w$Date, 1), rev(SRP.w$Date), head(SRP.w$Date, 1))
Y.Vec <- c(SRP.w$high, tail(SRP.w$low, 1), rev(SRP.w$low), head(SRP.w$high,1))
polygon(X.Vec, Y.Vec, col = "gray85", border = NA)
points(SRP.w$mid ~ SRP.w$Date, type='l', lwd=2, lty=1)
points(SRP.w$low ~ SRP.w$Date, type='l', lwd=1, lty=2)
points(SRP.w$high ~ SRP.w$Date, type='l', lwd=1, lty=2)
axis(side=2, labels=T, las=2, at=seq(0,60,15))
axis(side=1, labels=c("Day 55", "Day 70", "Day 85", "Day 100"), las=1,
     at = c(as.Date("2009-08-01"), as.Date("2009-08-15"),
            as.Date("2009-09-01"), as.Date("2009-09-15")))
mtext(expression(paste("SRP (µg P L"^"-1",")", sep="")),
      side=2, line=2.5, cex=1.5)
abline(v=as.Date("2009-08-28"), lwd=2, lty=3)
box(lwd=2)

# SRP Timeseries Peak
peak <- as.Date("2009-08-31")
peak.srp <- subset(SRP, SRP$Date == peak)
peak.srp <- merge(peak.srp, design, by="Pond")
colnames(peak.srp) <- c("Pond", "Date", "ppm", "plank", "supply")

# SRP Peak Plot
windows.options(width=8, height=4)
par(mar=c(5,5,0.5,0.5))
plot(peak.srp$ppm ~ peak.srp$supply,
     type='p', xlab='', yaxt='n',
     ylab='', ylim=c(0,75), pch=17, cex=1.5)
axis(side=2, labels=T, las=2, at=seq(0,75,25))
mtext(expression(paste("SRP (µg P L"^" -1",")", sep="")),
      side=2, line=2.5, cex=1.5)
mtext(expression(paste("DOC Supply (mg C L"^" -1",")", sep="")),
      side=1, line=3, cex=1.5)
box(lwd=2)

# TP vs Carbon
TP.pre <- TP[TP$Date > "2009-08-05" & TP$Date < "2009-08-25", ]
TP.pre2 <- merge(TP.pre, design, by="Pond")
colnames(TP.pre2) <- c("Pond", "Date", "ppm", "plank", "supply")
TP.model <- lm(TP.pre2$ppm ~ TP.pre2$supply)
summary(TP.model)

# TN vs Carbon
TN.pre <- TN[TN$Date > "2009-08-05" & TN$Date < "2009-08-25", ]
TN.pre2 <- merge(TN.pre, design, by="Pond")
colnames(TN.pre2) <- c("Pond", "Date", "ppm", "plank", "supply")
TN.model <- lm(TN.pre2$ppm ~ TN.pre2$supply)
summary(TN.model)

# Useful Dates
day60s <- min(which(TN.data2$Date == "2009-08-05"))
day100e <- max(which(TN.data2$Date == "2009-09-14"))
day80 <- as.Date("2009-08-25")
day82 <- as.Date("2009-08-27")
