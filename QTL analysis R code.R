setwd("C:/Users/hdong/Desktop/QTL analysis/mapA")

library(qtl)

mapa=read.cross("csv", , "Rqtl mapA.csv")
summary(mapa)

mapa=calc.genoprob(mapa,step=1,error.prob=0.001)
mapa=sim.geno(mapa,step=1,n.draws=100,error.prob=0.001)
####################################################
##    QTL analysis for Kobo data   ##
###################################################

## No QTL
Kobo.DF=cim(mapa, pheno.col=1, window = 10)
summary(Kobo.DF)
plot(Kobo.DF, ylab="LOD score")
add.cim.covar(Kobo.DF, col="green")


## No QTL
Kobo.DM=cim(mapa, pheno.col=2, window = 10)
summary(Kobo.DM)
plot(Kobo.DM, ylab="LOD score")
add.cim.covar(Kobo.DM, col="green")


## 
Kobo.GYP=cim(mapa, pheno.col=3, window = 10)
summary(Kobo.GYP)
plot(Kobo.GYP, ylab="LOD score")
add.cim.covar(Kobo.GYP, col="green")
lodint(Kobo.GYP, 2, 1.5)
find.marker(mapa, 2, 125.8)
effectplot(mapa, pheno.col = 3, mname1 = "S2_56650097")$Means
qtl <- makeqtl(mapa, chr = 2, pos = 125.8)
summary(fitqtl(mapa, qtl = qtl, pheno.col=3, get.ests=TRUE))
  
## 
Kobo.HE=cim(mapa, pheno.col=4, window = 10)
summary(Kobo.HE)
plot(Kobo.HE, ylab="LOD score")
add.cim.covar(Kobo.HE, col="green")
lodint(Kobo.HE, 1, 1)
find.marker(mapa, 1, 59)
lodint(Kobo.HE, 4, 1)
find.marker(mapa, 4, 43)
effectplot(mapa, pheno.col = 4, mname1 = "S1_55771897")$Means
effectplot(mapa, pheno.col = 4, mname1 = "S4_58856579")$Means
qtl <- makeqtl(mapa, chr = c(1,4), pos = c(59,43))
summary(fitqtl(mapa, qtl = qtl, pheno.col=4, get.ests=TRUE))


## 
Kobo.LS=cim(mapa, pheno.col=5, window = 10)
summary(Kobo.LS,chr=1)
plot(Kobo.LS, ylab="LOD score")
add.cim.covar(Kobo.LS, col="green")
lodint(Kobo.LS, 1, 1)
find.marker(mapa, 1, 48)
effectplot(mapa, pheno.col = 5, mname1 = "S1_54707166")$Means
qtl <- makeqtl(mapa, chr = 1, pos = 48)
summary(fitqtl(mapa, qtl = qtl, pheno.col=5, get.ests=TRUE))



## 
Kobo.NL=cim(mapa, pheno.col=6, window = 10)
summary(Kobo.NL)
plot(Kobo.NL, ylab="LOD score")
add.cim.covar(Kobo.NL, col="green")


## No QTL
Kobo.NPP=cim(mapa, pheno.col=7, window = 10)
summary(Kobo.NPP)
plot(Kobo.NPP, ylab="LOD score")
add.cim.covar(Kobo.NPP, col="green")
lodint(Kobo.NPP, 2, 1.5)
find.marker(mapa, 2, 86)
effectplot(mapa, pheno.col = 7, mname1 = "S2_3154794")


## No QTL
Kobo.NT=cim(mapa, pheno.col=8, window = 10)
summary(Kobo.NT)
plot(Kobo.NT, ylab="LOD score")
add.cim.covar(Kobo.NT, col="green")


## 
Kobo.PH=cim(mapa, pheno.col=9, window = 10)
summary(Kobo.PH)
plot(Kobo.PH, ylab="LOD score")
add.cim.covar(Kobo.PH, col="green")
lodint(Kobo.PH, 2, 1.5)
find.marker(mapa, 2, 160.41)
effectplot(mapa, pheno.col = 9, mname1 = "S2_60462641")$Means

lodint(Kobo.PH, 6, 1.5)
find.marker(mapa, 6, 4.94)
effectplot(mapa, pheno.col = 9, mname1 = "S6_42495445")$Means

lodint(Kobo.PH, 7, 1.5)
find.marker(mapa, 7, 85)
effectplot(mapa, pheno.col = 9, mname1 = "S7_53119949")$Means

qtl <- makeqtl(mapa, chr = c(2,6,7), pos = c(160.41,4.94,85))
summary(fitqtl(mapa, qtl = qtl, pheno.col=9, get.ests=TRUE))



## No QTL
Kobo.PL=cim(mapa, pheno.col=10, window = 10)
summary(Kobo.PL)
plot(Kobo.PL, ylab="LOD score")
add.cim.covar(Kobo.PL, col="green")


## No QTL
Kobo.PW=cim(mapa, pheno.col=11, window = 10)
summary(Kobo.PW)
plot(Kobo.PW, ylab="LOD score")
add.cim.covar(Kobo.PW, col="green")


## No QTL
Kobo.PWD=cim(mapa, pheno.col=12, window = 10)
summary(Kobo.PWD)
plot(Kobo.PWD, ylab="LOD score")
add.cim.covar(Kobo.PWD, col="green")


## 
Kobo.TSW=cim(mapa, pheno.col=13, window = 10)
summary(Kobo.TSW)
plot(Kobo.TSW, ylab="LOD score")
add.cim.covar(Kobo.TSW, col="green")

lodint(Kobo.TSW, 3, 1.5)
find.marker(mapa, 3, 67)
effectplot(mapa, pheno.col = 13, mname1 = "S3_13049407")$Means

qtl <- makeqtl(mapa, chr = 3, pos = 67)
summary(fitqtl(mapa, qtl = qtl, pheno.col=13, get.ests=TRUE))







####################################################
##    QTL analysis for Meiso data   ##
###################################################

## No QTL
Meiso.DF=cim(mapa, pheno.col=14, window = 10)
summary(Meiso.DF)
plot(Meiso.DF, ylab="LOD score")
add.cim.covar(Meiso.DF, col="green")

## No QTL
Meiso.DM=cim(mapa, pheno.col=15, window = 10)
summary(Meiso.DM)
plot(Meiso.DM, ylab="LOD score")
add.cim.covar(Meiso.DM, col="green")


## No QTL
Meiso.GYP=cim(mapa, pheno.col=16, window = 10)
summary(Meiso.GYP)
plot(Meiso.GYP, ylab="LOD score")
add.cim.covar(Meiso.GYP, col="green")


## No QTL
Meiso.HE=cim(mapa, pheno.col=17, window = 10)
summary(Meiso.HE)
plot(Meiso.HE, ylab="LOD score")
add.cim.covar(Meiso.HE, col="green")


## No QTL
Meiso.LS=cim(mapa, pheno.col=18, window = 10)
summary(Meiso.LS)
plot(Meiso.LS, ylab="LOD score")
add.cim.covar(Meiso.LS, col="green")


## No QTL
Meiso.NL=cim(mapa, pheno.col=19, window = 10)
summary(Meiso.NL)
plot(Meiso.NL, ylab="LOD score")
add.cim.covar(Meiso.NL, col="green")


## No QTL
Meiso.NPP=cim(mapa, pheno.col=20, window = 10)
summary(Meiso.NPP)
plot(Meiso.NPP, ylab="LOD score")
add.cim.covar(Meiso.NPP, col="green")


## No QTL
Meiso.NT=cim(mapa, pheno.col=21, window = 10)
summary(Meiso.NT)
plot(Meiso.NT, ylab="LOD score")
add.cim.covar(Meiso.NT, col="green")


## Near Dw1?
Meiso.PH=cim(mapa, pheno.col=22, window = 10)
summary(Meiso.PH)
plot(Meiso.PH, ylab="LOD score")
add.cim.covar(Meiso.PH, col="green")

lodint(Meiso.PH, 3, 1.5)
find.marker(mapa, 3, 139.02)
effectplot(mapa, pheno.col = 22, mname1 = "S3_58535398")$Means

qtl <- makeqtl(mapa, chr = 3, pos = 139.02)
summary(fitqtl(mapa, qtl = qtl, pheno.col=22, get.ests=TRUE))

## 
Meiso.PL=cim(mapa, pheno.col=23, window = 10)
summary(Meiso.PL)
plot(Meiso.PL, ylab="LOD score")
add.cim.covar(Meiso.PL, col="green")


## No QTL
Meiso.PW=cim(mapa, pheno.col=24, window = 10)
summary(Meiso.PW)
plot(Meiso.PW, ylab="LOD score")
add.cim.covar(Meiso.PW, col="green")


## No QTL
Meiso.PWD=cim(mapa, pheno.col=25, window = 10)
summary(Meiso.PWD)
plot(Meiso.PWD, ylab="LOD score")
add.cim.covar(Meiso.PWD, col="green")


## 
Meiso.TSW=cim(mapa, pheno.col=26, window = 10)
summary(Meiso.TSW)
plot(Meiso.TSW, ylab="LOD score")
add.cim.covar(Meiso.TSW, col="green")

lodint(Meiso.TSW, 8, 1.5)
lodint(Meiso.TSW, 10, 1.5)
find.marker(mapa, 8, 117.41)
find.marker(mapa, 10, 96.02)
effectplot(mapa, pheno.col = 26, mname1 = "S8_51144389")$Means
effectplot(mapa, pheno.col = 26, mname1 = "S10_46553796")$Means

qtl <- makeqtl(mapa, chr = c(8, 10), pos = c(117.41, 96.02))
summary(fitqtl(mapa, qtl = qtl, pheno.col=26, get.ests=TRUE))


