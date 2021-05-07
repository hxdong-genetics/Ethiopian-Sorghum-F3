setwd("~/OneDrive - University of Georgia/Postdoc Projects/Ethiopia_Mapping_Pops/1_Separate")
library(qtl)

# Plot the mapA
mapA <- read.csv("PopA map.csv")
mapA.pos <- split(mapA[,3], mapA[,2])
mapA.mar <- split(mapA[,1], mapA[,2])
for(i in seq(along=mapA.pos))
  names(mapA.pos[[i]]) <- mapA.mar[[i]]
mapA.map <- mapA.pos
class(mapA.map) <- "map"

plot(mapA.map, main="")


# Plot the mapB
mapB <- read.csv("PopB map.csv")
mapB.pos <- split(mapB[,3], mapB[,2])
mapB.mar <- split(mapB[,1], mapB[,2])
for(i in seq(along=mapB.pos))
  names(mapB.pos[[i]]) <- mapB.mar[[i]]
mapB.map <- mapB.pos
class(mapB.map) <- "map"

plot(mapB.map, main="")


# Plot the mapC
mapC <- read.csv("PopC map.csv")
mapC.pos <- split(mapC[,3], mapC[,2])
mapC.mar <- split(mapC[,1], mapC[,2])
for(i in seq(along=mapC.pos))
  names(mapC.pos[[i]]) <- mapC.mar[[i]]
mapC.map <- mapC.pos
class(mapC.map) <- "map"

plot(mapC.map, main="")


## Composite map
Composite <- read.csv("Pop CompositeABC.csv", header = TRUE)
CMplot(Composite,plot.type="d",bin.size=1e6,col=c("darkgreen", "yellow", "red"),file="jpg",memo="",dpi=300)

length(unique(Composite[,1]))

table(mapA[,1]%in%mapB[,1])
table(mapA[,1]%in%mapC[,1])
table(mapB[,1]%in%mapC[,1])
table(mapA[,1]%in%mapB[,1])
## Order comparison between genetic map and reference genome


LG1=mapB[mapB[,2]==1,]
LG2=mapB[mapB[,2]==2,]
LG3=mapB[mapB[,2]==3,]
LG4=mapB[mapB[,2]==4,]
LG5=mapB[mapB[,2]==5,]
LG6=mapB[mapB[,2]==6,]
LG7=mapB[mapB[,2]==7,]
LG8=mapB[mapB[,2]==8,]
LG9=mapB[mapB[,2]==9,]
LG10=mapB[mapB[,2]==10,]


tiff("PopB collinearity between genetic map and physical reference genome.tiff",
     units = "in", width = 7, height = 7, res = 300)
par(mfrow=c(2,5))
### LG1 ###
plot(c(rep(1,nrow(LG1))),LG1[,3]*(1000/max(LG1[,3])),type = "l",col="blue",main="Chr 1",tck=0,axes = FALSE,xlab = "",ylab = "",lwd=10,xlim = c(0.5,2.5),ylim = c(0,1100))
lines(c(rep(2,nrow(LG1))),LG1[,5]*(1000/max(LG1[,5])),lwd=10)

for (i in 1:length(LG1[,2])){
  if (LG1[,2][i]==LG1[,4][i]) {
    segments(x0=1,y0=as.numeric(LG1[i,3]*(1000/max(LG1[,3]))),x1=2,y1=as.numeric(LG1[i,5]*(1000/max(LG1[,5]))),col = "orange",lty = 1,lwd=1)}
  }


### LG2 ###
plot(c(rep(1,nrow(LG2))),LG2[,3]*(1000/max(LG2[,3])),type = "l",col="blue",main="Chr 2",tck=0,axes = FALSE,xlab = "",ylab = "",lwd=10,xlim = c(0.5,2.5),ylim = c(0,1100))
lines(c(rep(2,nrow(LG2))),LG2[,5]*(1000/max(LG2[,5])),lwd=10)

for (i in 1:length(LG2[,2])){
  if (LG2[,2][i]==LG2[,4][i]) {
    segments(x0=1,y0=as.numeric(LG2[i,3]*(1000/max(LG2[,3]))),x1=2,y1=as.numeric(LG2[i,5]*(1000/max(LG2[,5]))),col = "orange",lty = 1,lwd=1)}
}



### LG3 ###
plot(c(rep(1,nrow(LG3))),LG3[,3]*(1000/max(LG3[,3])),type = "l",col="blue",main="Chr 3",tck=0,axes = FALSE,xlab = "",ylab = "",lwd=10,xlim = c(0.5,2.5),ylim = c(0,1100))
lines(c(rep(2,nrow(LG3))),LG3[,5]*(1000/max(LG3[,5])),lwd=10)

for (i in 1:length(LG3[,2])){
  if (LG3[,2][i]==LG3[,4][i]) {
    segments(x0=1,y0=as.numeric(LG3[i,3]*(1000/max(LG3[,3]))),x1=2,y1=as.numeric(LG3[i,5]*(1000/max(LG3[,5]))),col = "orange",lty = 1,lwd=1)}
}


### LG4 ###
plot(c(rep(1,nrow(LG4))),LG4[,3]*(1000/max(LG4[,3])),type = "l",col="blue",main="Chr 4",tck=0,axes = FALSE,xlab = "",ylab = "",lwd=10,xlim = c(0.5,2.5),ylim = c(0,1100))
lines(c(rep(2,nrow(LG4))),LG4[,5]*(1000/max(LG4[,5])),lwd=10)

for (i in 1:length(LG4[,2])){
  if (LG4[,2][i]==LG4[,4][i]) {
    segments(x0=1,y0=as.numeric(LG4[i,3]*(1000/max(LG4[,3]))),x1=2,y1=as.numeric(LG4[i,5]*(1000/max(LG4[,5]))),col = "orange",lty = 1,lwd=1)}
}

### LG5 ###
plot(c(rep(1,nrow(LG5))),LG5[,3]*(1000/max(LG5[,3])),type = "l",col="blue",main="Chr 5",tck=0,axes = FALSE,xlab = "",ylab = "",lwd=10,xlim = c(0.5,2.5),ylim = c(0,1100))
lines(c(rep(2,nrow(LG5))),LG5[,5]*(1000/max(LG5[,5])),lwd=10)

for (i in 1:length(LG5[,2])){
  if (LG5[,2][i]==LG5[,4][i]) {
    segments(x0=1,y0=as.numeric(LG5[i,3]*(1000/max(LG5[,3]))),x1=2,y1=as.numeric(LG5[i,5]*(1000/max(LG5[,5]))),col = "orange",lty = 1,lwd=1)}
}


### LG6 ###
plot(c(rep(1,nrow(LG6))),LG6[,3]*(1000/max(LG6[,3])),type = "l",col="blue",main="Chr 6",tck=0,axes = FALSE,xlab = "",ylab = "",lwd=10,xlim = c(0.5,2.5),ylim = c(0,1100))
lines(c(rep(2,nrow(LG6))),LG6[,5]*(1000/max(LG6[,5])),lwd=10)

for (i in 1:length(LG6[,2])){
  if (LG6[,2][i]==LG6[,4][i]) {
    segments(x0=1,y0=as.numeric(LG6[i,3]*(1000/max(LG6[,3]))),x1=2,y1=as.numeric(LG6[i,5]*(1000/max(LG6[,5]))),col = "orange",lty = 1,lwd=1)}
}

### LG7 ###
plot(c(rep(1,nrow(LG7))),LG7[,3]*(1000/max(LG7[,3])),type = "l",col="blue",main="Chr 7",tck=0,axes = FALSE,xlab = "",ylab = "",lwd=10,xlim = c(0.5,2.5),ylim = c(0,1100))
lines(c(rep(2,nrow(LG7))),LG7[,5]*(1000/max(LG7[,5])),lwd=10)

for (i in 1:length(LG7[,2])){
  if (LG7[,2][i]==LG7[,4][i]) {
    segments(x0=1,y0=as.numeric(LG7[i,3]*(1000/max(LG7[,3]))),x1=2,y1=as.numeric(LG7[i,5]*(1000/max(LG7[,5]))),col = "orange",lty = 1,lwd=1)}
}

### LG8 ###
plot(c(rep(1,nrow(LG8))),LG8[,3]*(1000/max(LG8[,3])),type = "l",col="blue",main="Chr 8",tck=0,axes = FALSE,xlab = "",ylab = "",lwd=10,xlim = c(0.5,2.5),ylim = c(0,1100))
lines(c(rep(2,nrow(LG8))),LG8[,5]*(1000/max(LG8[,5])),lwd=10)

for (i in 1:length(LG8[,2])){
  if (LG8[,2][i]==LG8[,4][i]) {
    segments(x0=1,y0=as.numeric(LG8[i,3]*(1000/max(LG8[,3]))),x1=2,y1=as.numeric(LG8[i,5]*(1000/max(LG8[,5]))),col = "orange",lty = 1,lwd=1)}
}


### LG9 ###
plot(c(rep(1,nrow(LG9))),LG9[,3]*(1000/max(LG9[,3])),type = "l",col="blue",main="Chr 9",tck=0,axes = FALSE,xlab = "",ylab = "",lwd=10,xlim = c(0.5,2.5),ylim = c(0,1100))
lines(c(rep(2,nrow(LG9))),LG9[,5]*(1000/max(LG9[,5])),lwd=10)

for (i in 1:length(LG9[,2])){
  if (LG9[,2][i]==LG9[,4][i]) {
    segments(x0=1,y0=as.numeric(LG9[i,3]*(1000/max(LG9[,3]))),x1=2,y1=as.numeric(LG9[i,5]*(1000/max(LG9[,5]))),col = "orange",lty = 1,lwd=1)}
}


### LG10 ###
plot(c(rep(1,nrow(LG10))),LG10[,3]*(1000/max(LG10[,3])),type = "l",col="blue",main="Chr 10",tck=0,axes = FALSE,xlab = "",ylab = "",lwd=10,xlim = c(0.5,2.5),ylim = c(0,1100))
lines(c(rep(2,nrow(LG10))),LG10[,5]*(1000/max(LG10[,5])),lwd=10)

for (i in 1:length(LG10[,2])){
  if (LG10[,2][i]==LG10[,4][i]) {
    segments(x0=1,y0=as.numeric(LG10[i,3]*(1000/max(LG10[,3]))),x1=2,y1=as.numeric(LG10[i,5]*(1000/max(LG10[,5]))),col = "orange",lty = 1,lwd=1)}
}

dev.off()