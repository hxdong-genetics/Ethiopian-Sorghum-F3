######################################################################################
setwd("C:/Users/hdong/Desktop/ABHGenotypeR")

library(qtl)

######################################################################################
# Plot the map
PopA <- read.csv("genetic map of mapA.csv")
pos <- split(PopA[,3], PopA[,2])
mar <- split(PopA[,1], PopA[,2])
for(i in seq(along=pos))
  names(pos[[i]]) <- mar[[i]]
PopA <- pos
class(PopA) <- "map"
plot.map(PopA, show.marker.names = FALSE, main="Genetic map of population A (Baji x 76T1-23)")

######################################################################################
# Plot the map
PopB <- read.csv("genetic map of mapB.csv")
pos <- split(PopB[,3], PopB[,2])
mar <- split(PopB[,1], PopB[,2])
for(i in seq(along=pos))
  names(pos[[i]]) <- mar[[i]]
PopB <- pos
class(PopB) <- "map"
plot.map(PopB, show.marker.names = FALSE, main="Genetic map of population B (Birmash x Meco-1)")

######################################################################################
# Plot the map
PopC <- read.csv("genetic map of mapC.csv")
pos <- split(PopC[,3], PopC[,2])
mar <- split(PopC[,1], PopC[,2])
for(i in seq(along=pos))
  names(pos[[i]]) <- mar[[i]]
PopC <- pos
class(PopC) <- "map"
plot.map(PopC, show.marker.names = FALSE, main="Genetic map of population C (Birmash x 76T1-23)")

