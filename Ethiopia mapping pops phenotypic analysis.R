##### Ethiopia mapping populations (3 pops), Phenotypic data analysis
setwd("C:/Users/hdong/OneDrive - University of Georgia/Postdoc Projects/Ethiopia_Mapping_Pops/Phenotypic Data")

## PopA: 76T1-23 x Baji
## PopB: 76T1-23 x Birmash
## PopC: Meco-1  x Birmash

## 13 traits, 2 locations (Kobo, Meiso), Kobo had more drought than Meiso

## DF: Days to 50% flowering
## DM: Days to maturity
## PH: Plant height (cm)
## NT: Number of tillers per plant
## NL: Number of leaves per plant
## HE: Head exsertion (cm)
## NPP: Number of panicle per plant
## LS: Leaf sensence (score)
## PW: Panicle weight (g)
## PL: Panicle length (cm)
## PWD: Panicle width (cm)
## GYP: Grain yield per plant (g)
## TSW: Thousand seed weight (g)


# Load raw data
Pheno <- read.csv("Ethiopia_Mapping_Pops_Pheno.csv", header = TRUE)
Pheno[,7] <- as.factor(Pheno[,7])
Pheno[,6] <- as.factor(Pheno[,6])
library(lsmeans)
library(lattice)



################################################################################
## Calculate LSMeans with lsmeans package
################################################################################
Pheno.A <- Pheno[c(Pheno[,1]=="A"),]
Pheno.B <- Pheno[c(Pheno[,1]=="B"),]
Pheno.C <- Pheno[c(Pheno[,1]=="C"),]

Pheno.A.Kobo <- Pheno.A[c(Pheno.A[,3]=="Kobo"),]
Pheno.B.Kobo <- Pheno.B[c(Pheno.B[,3]=="Kobo"),]
Pheno.C.Kobo <- Pheno.C[c(Pheno.C[,3]=="Kobo"),]

Pheno.A.Meiso <- Pheno.A[c(Pheno.A[,3]=="Meiso"),]
Pheno.B.Meiso <- Pheno.B[c(Pheno.B[,3]=="Meiso"),]
Pheno.C.Meiso <- Pheno.C[c(Pheno.C[,3]=="Meiso"),]


# Calculate LSMeans for each population, at each location (Kobo and Meiso).
## Population C, kobo
C.Kobo.DF.lm <- lm(DF ~ Genotypes, data = Pheno.C.Kobo)
C.Kobo.DF.lsmean <- as.data.frame(lsmeans(C.Kobo.DF.lm, ~ Genotypes))
hist(C.Kobo.DF.lsmean[,2])

C.Kobo.DM.lm <- lm(DM ~ Genotypes, data = Pheno.C.Kobo)
C.Kobo.DM.lsmean <- as.data.frame(lsmeans(C.Kobo.DM.lm, ~ Genotypes))
hist(C.Kobo.DM.lsmean[,2])

C.Kobo.PH.lm <- lm(PH ~ Genotypes, data = Pheno.C.Kobo)
C.Kobo.PH.lsmean <- as.data.frame(lsmeans(C.Kobo.PH.lm, ~ Genotypes))
hist(C.Kobo.PH.lsmean[,2])

C.Kobo.NT.lm <- lm(NT ~ Genotypes, data = Pheno.C.Kobo)
C.Kobo.NT.lsmean <- as.data.frame(lsmeans(C.Kobo.NT.lm, ~ Genotypes))
hist(C.Kobo.NT.lsmean[,2])

C.Kobo.NL.lm <- lm(NL ~ Genotypes, data = Pheno.C.Kobo)
C.Kobo.NL.lsmean <- as.data.frame(lsmeans(C.Kobo.NL.lm, ~ Genotypes))
hist(C.Kobo.NL.lsmean[,2])

C.Kobo.HE.lm <- lm(HE ~ Genotypes, data = Pheno.C.Kobo)
C.Kobo.HE.lsmean <- as.data.frame(lsmeans(C.Kobo.HE.lm, ~ Genotypes))
hist(C.Kobo.HE.lsmean[,2]) # All genotypes have value of 0, NO variation for head exsertion.

C.Kobo.NPP.lm <- lm(NPP ~ Genotypes, data = Pheno.C.Kobo)
C.Kobo.NPP.lsmean <- as.data.frame(lsmeans(C.Kobo.NPP.lm, ~ Genotypes))
hist(C.Kobo.NPP.lsmean[,2])

C.Kobo.LS.lm <- lm(LS ~ Genotypes, data = Pheno.C.Kobo)
C.Kobo.LS.lsmean <- as.data.frame(lsmeans(C.Kobo.LS.lm, ~ Genotypes))
hist(C.Kobo.LS.lsmean[,2])

C.Kobo.PW.lm <- lm(PW ~ Genotypes, data = Pheno.C.Kobo)
C.Kobo.PW.lsmean <- as.data.frame(lsmeans(C.Kobo.PW.lm, ~ Genotypes))
hist(C.Kobo.PW.lsmean[,2])

C.Kobo.PL.lm <- lm(PL ~ Genotypes, data = Pheno.C.Kobo)
C.Kobo.PL.lsmean <- as.data.frame(lsmeans(C.Kobo.PL.lm, ~ Genotypes))
hist(C.Kobo.PL.lsmean[,2])

C.Kobo.PWD.lm <- lm(PWD ~ Genotypes, data = Pheno.C.Kobo)
C.Kobo.PWD.lsmean <- as.data.frame(lsmeans(C.Kobo.PWD.lm, ~ Genotypes))
hist(C.Kobo.PWD.lsmean[,2])

C.Kobo.GYP.lm <- lm(GYP ~ Genotypes, data = Pheno.C.Kobo)
C.Kobo.GYP.lsmean <- as.data.frame(lsmeans(C.Kobo.GYP.lm, ~ Genotypes))
hist(C.Kobo.GYP.lsmean[,2])

C.Kobo.TSW.lm <- lm(TSW ~ Genotypes, data = Pheno.C.Kobo)
C.Kobo.TSW.lsmean <- as.data.frame(lsmeans(C.Kobo.TSW.lm, ~ Genotypes))
hist(C.Kobo.TSW.lsmean[,2])

C.Kobo.Trait.LSMeans <- cbind(C.Kobo.DF.lsmean[,1:2], C.Kobo.DM.lsmean[,2],
                              C.Kobo.GYP.lsmean[,2], C.Kobo.HE.lsmean[,2],
                              C.Kobo.LS.lsmean[,2], C.Kobo.NL.lsmean[,2],
                              C.Kobo.NPP.lsmean[,2], C.Kobo.NT.lsmean[,2],
                              C.Kobo.PH.lsmean[,2], C.Kobo.PL.lsmean[,2],
                              C.Kobo.PW.lsmean[,2], C.Kobo.PWD.lsmean[,2],
                              C.Kobo.TSW.lsmean[,2])
colnames(C.Kobo.Trait.LSMeans) <- c("Genotypes", "Kobo.DF", "Kobo.DM", "Kobo.GYP", "Kobo.HE",
                                    "Kobo.LS", "Kobo.NL", "Kobo.NPP", "Kobo.NT", "Kobo.PH",
                                    "Kobo.PL", "Kobo.PW", "Kobo.PWD", "Kobo.TSW")
write.csv(C.Kobo.Trait.LSMeans, "mapC Kobo trait lsmeans.csv")





## Population C, Meiso
C.Meiso.DF.lm <- lm(DF ~ Genotypes, data = Pheno.C.Meiso)
C.Meiso.DF.lsmean <- as.data.frame(lsmeans(C.Meiso.DF.lm, ~ Genotypes))
hist(C.Meiso.DF.lsmean[,2])

C.Meiso.DM.lm <- lm(DM ~ Genotypes, data = Pheno.C.Meiso)
C.Meiso.DM.lsmean <- as.data.frame(lsmeans(C.Meiso.DM.lm, ~ Genotypes))
hist(C.Meiso.DM.lsmean[,2])

C.Meiso.PH.lm <- lm(PH ~ Genotypes, data = Pheno.C.Meiso)
C.Meiso.PH.lsmean <- as.data.frame(lsmeans(C.Meiso.PH.lm, ~ Genotypes))
hist(C.Meiso.PH.lsmean[,2])

C.Meiso.NT.lm <- lm(NT ~ Genotypes, data = Pheno.C.Meiso)
C.Meiso.NT.lsmean <- as.data.frame(lsmeans(C.Meiso.NT.lm, ~ Genotypes))
hist(C.Meiso.NT.lsmean[,2])

C.Meiso.NL.lm <- lm(NL ~ Genotypes, data = Pheno.C.Meiso)
C.Meiso.NL.lsmean <- as.data.frame(lsmeans(C.Meiso.NL.lm, ~ Genotypes))
hist(C.Meiso.NL.lsmean[,2])

C.Meiso.HE.lm <- lm(HE ~ Genotypes, data = Pheno.C.Meiso)
C.Meiso.HE.lsmean <- as.data.frame(lsmeans(C.Meiso.HE.lm, ~ Genotypes))
hist(C.Meiso.HE.lsmean[,2]) # All genotypes have value of 0, NO variation for head exsertion.

C.Meiso.NPP.lm <- lm(NPP ~ Genotypes, data = Pheno.C.Meiso)
C.Meiso.NPP.lsmean <- as.data.frame(lsmeans(C.Meiso.NPP.lm, ~ Genotypes))
hist(C.Meiso.NPP.lsmean[,2])

C.Meiso.LS.lm <- lm(LS ~ Genotypes, data = Pheno.C.Meiso)
C.Meiso.LS.lsmean <- as.data.frame(lsmeans(C.Meiso.LS.lm, ~ Genotypes))
hist(C.Meiso.LS.lsmean[,2])

C.Meiso.PW.lm <- lm(PW ~ Genotypes, data = Pheno.C.Meiso)
C.Meiso.PW.lsmean <- as.data.frame(lsmeans(C.Meiso.PW.lm, ~ Genotypes))
hist(C.Meiso.PW.lsmean[,2])

C.Meiso.PL.lm <- lm(PL ~ Genotypes, data = Pheno.C.Meiso)
C.Meiso.PL.lsmean <- as.data.frame(lsmeans(C.Meiso.PL.lm, ~ Genotypes))
hist(C.Meiso.PL.lsmean[,2])

C.Meiso.PWD.lm <- lm(PWD ~ Genotypes, data = Pheno.C.Meiso)
C.Meiso.PWD.lsmean <- as.data.frame(lsmeans(C.Meiso.PWD.lm, ~ Genotypes))
hist(C.Meiso.PWD.lsmean[,2])

C.Meiso.GYP.lm <- lm(GYP ~ Genotypes, data = Pheno.C.Meiso)
C.Meiso.GYP.lsmean <- as.data.frame(lsmeans(C.Meiso.GYP.lm, ~ Genotypes))
hist(C.Meiso.GYP.lsmean[,2])

C.Meiso.TSW.lm <- lm(TSW ~ Genotypes, data = Pheno.C.Meiso)
C.Meiso.TSW.lsmean <- as.data.frame(lsmeans(C.Meiso.TSW.lm, ~ Genotypes))
hist(C.Meiso.TSW.lsmean[,2])

C.Meiso.Trait.LSMeans <- cbind(C.Meiso.DF.lsmean[,1:2], C.Meiso.DM.lsmean[,2],
                              C.Meiso.GYP.lsmean[,2], C.Meiso.HE.lsmean[,2],
                              C.Meiso.LS.lsmean[,2], C.Meiso.NL.lsmean[,2],
                              C.Meiso.NPP.lsmean[,2], C.Meiso.NT.lsmean[,2],
                              C.Meiso.PH.lsmean[,2], C.Meiso.PL.lsmean[,2],
                              C.Meiso.PW.lsmean[,2], C.Meiso.PWD.lsmean[,2],
                              C.Meiso.TSW.lsmean[,2])
colnames(C.Meiso.Trait.LSMeans) <- c("Genotypes", "Meiso.DF", "Meiso.DM", "Meiso.GYP", "Meiso.HE",
                                    "Meiso.LS", "Meiso.NL", "Meiso.NPP", "Meiso.NT", "Meiso.PH",
                                    "Meiso.PL", "Meiso.PW", "Meiso.PWD", "Meiso.TSW")
write.csv(C.Meiso.Trait.LSMeans, "mapC Meiso trait lsmeans.csv")
















# Calculate LSMeans for each population, at each location (Kobo and Meiso).
## Population B, kobo
B.Kobo.DF.lm <- lm(DF ~ Genotypes, data = Pheno.B.Kobo)
B.Kobo.DF.lsmean <- as.data.frame(lsmeans(B.Kobo.DF.lm, ~ Genotypes))
hist(B.Kobo.DF.lsmean[,2])

B.Kobo.DM.lm <- lm(DM ~ Genotypes, data = Pheno.B.Kobo)
B.Kobo.DM.lsmean <- as.data.frame(lsmeans(B.Kobo.DM.lm, ~ Genotypes))
hist(B.Kobo.DM.lsmean[,2])

B.Kobo.PH.lm <- lm(PH ~ Genotypes, data = Pheno.B.Kobo)
B.Kobo.PH.lsmean <- as.data.frame(lsmeans(B.Kobo.PH.lm, ~ Genotypes))
hist(B.Kobo.PH.lsmean[,2])

B.Kobo.NT.lm <- lm(NT ~ Genotypes, data = Pheno.B.Kobo)
B.Kobo.NT.lsmean <- as.data.frame(lsmeans(B.Kobo.NT.lm, ~ Genotypes))
hist(B.Kobo.NT.lsmean[,2])

B.Kobo.NL.lm <- lm(NL ~ Genotypes, data = Pheno.B.Kobo)
B.Kobo.NL.lsmean <- as.data.frame(lsmeans(B.Kobo.NL.lm, ~ Genotypes))
hist(B.Kobo.NL.lsmean[,2])

B.Kobo.HE.lm <- lm(HE ~ Genotypes, data = Pheno.B.Kobo)
B.Kobo.HE.lsmean <- as.data.frame(lsmeans(B.Kobo.HE.lm, ~ Genotypes))
hist(B.Kobo.HE.lsmean[,2]) # All genotypes have value of 0, NO variation for head exsertion.

B.Kobo.NPP.lm <- lm(NPP ~ Genotypes, data = Pheno.B.Kobo)
B.Kobo.NPP.lsmean <- as.data.frame(lsmeans(B.Kobo.NPP.lm, ~ Genotypes))
hist(B.Kobo.NPP.lsmean[,2])

B.Kobo.LS.lm <- lm(LS ~ Genotypes, data = Pheno.B.Kobo)
B.Kobo.LS.lsmean <- as.data.frame(lsmeans(B.Kobo.LS.lm, ~ Genotypes))
hist(B.Kobo.LS.lsmean[,2])

B.Kobo.PW.lm <- lm(PW ~ Genotypes, data = Pheno.B.Kobo)
B.Kobo.PW.lsmean <- as.data.frame(lsmeans(B.Kobo.PW.lm, ~ Genotypes))
hist(B.Kobo.PW.lsmean[,2])

B.Kobo.PL.lm <- lm(PL ~ Genotypes, data = Pheno.B.Kobo)
B.Kobo.PL.lsmean <- as.data.frame(lsmeans(B.Kobo.PL.lm, ~ Genotypes))
hist(B.Kobo.PL.lsmean[,2])

B.Kobo.PWD.lm <- lm(PWD ~ Genotypes, data = Pheno.B.Kobo)
B.Kobo.PWD.lsmean <- as.data.frame(lsmeans(B.Kobo.PWD.lm, ~ Genotypes))
hist(B.Kobo.PWD.lsmean[,2])

B.Kobo.GYP.lm <- lm(GYP ~ Genotypes, data = Pheno.B.Kobo)
B.Kobo.GYP.lsmean <- as.data.frame(lsmeans(B.Kobo.GYP.lm, ~ Genotypes))
hist(B.Kobo.GYP.lsmean[,2])

B.Kobo.TSW.lm <- lm(TSW ~ Genotypes, data = Pheno.B.Kobo)
B.Kobo.TSW.lsmean <- as.data.frame(lsmeans(B.Kobo.TSW.lm, ~ Genotypes))
hist(B.Kobo.TSW.lsmean[,2])

B.Kobo.Trait.LSMeans <- cbind(B.Kobo.DF.lsmean[,1:2], B.Kobo.DM.lsmean[,2],
                              B.Kobo.GYP.lsmean[,2], B.Kobo.HE.lsmean[,2],
                              B.Kobo.LS.lsmean[,2], B.Kobo.NL.lsmean[,2],
                              B.Kobo.NPP.lsmean[,2], B.Kobo.NT.lsmean[,2],
                              B.Kobo.PH.lsmean[,2], B.Kobo.PL.lsmean[,2],
                              B.Kobo.PW.lsmean[,2], B.Kobo.PWD.lsmean[,2],
                              B.Kobo.TSW.lsmean[,2])
colnames(B.Kobo.Trait.LSMeans) <- c("Genotypes", "Kobo.DF", "Kobo.DM", "Kobo.GYP", "Kobo.HE",
                                    "Kobo.LS", "Kobo.NL", "Kobo.NPP", "Kobo.NT", "Kobo.PH",
                                    "Kobo.PL", "Kobo.PW", "Kobo.PWD", "Kobo.TSW")
write.csv(B.Kobo.Trait.LSMeans, "mapB Kobo trait lsmeans.csv")





## Population B, Meiso
B.Meiso.DF.lm <- lm(DF ~ Genotypes, data = Pheno.B.Meiso)
B.Meiso.DF.lsmean <- as.data.frame(lsmeans(B.Meiso.DF.lm, ~ Genotypes))
hist(B.Meiso.DF.lsmean[,2])

B.Meiso.DM.lm <- lm(DM ~ Genotypes, data = Pheno.B.Meiso)
B.Meiso.DM.lsmean <- as.data.frame(lsmeans(B.Meiso.DM.lm, ~ Genotypes))
hist(B.Meiso.DM.lsmean[,2])

B.Meiso.PH.lm <- lm(PH ~ Genotypes, data = Pheno.B.Meiso)
B.Meiso.PH.lsmean <- as.data.frame(lsmeans(B.Meiso.PH.lm, ~ Genotypes))
hist(B.Meiso.PH.lsmean[,2])

B.Meiso.NT.lm <- lm(NT ~ Genotypes, data = Pheno.B.Meiso)
B.Meiso.NT.lsmean <- as.data.frame(lsmeans(B.Meiso.NT.lm, ~ Genotypes))
hist(B.Meiso.NT.lsmean[,2])

B.Meiso.NL.lm <- lm(NL ~ Genotypes, data = Pheno.B.Meiso)
B.Meiso.NL.lsmean <- as.data.frame(lsmeans(B.Meiso.NL.lm, ~ Genotypes))
hist(B.Meiso.NL.lsmean[,2])

B.Meiso.HE.lm <- lm(HE ~ Genotypes, data = Pheno.B.Meiso)
B.Meiso.HE.lsmean <- as.data.frame(lsmeans(B.Meiso.HE.lm, ~ Genotypes))
hist(B.Meiso.HE.lsmean[,2]) # All genotypes have value of 0, NO variation for head exsertion.

B.Meiso.NPP.lm <- lm(NPP ~ Genotypes, data = Pheno.B.Meiso)
B.Meiso.NPP.lsmean <- as.data.frame(lsmeans(B.Meiso.NPP.lm, ~ Genotypes))
hist(B.Meiso.NPP.lsmean[,2])

B.Meiso.LS.lm <- lm(LS ~ Genotypes, data = Pheno.B.Meiso)
B.Meiso.LS.lsmean <- as.data.frame(lsmeans(B.Meiso.LS.lm, ~ Genotypes))
hist(B.Meiso.LS.lsmean[,2])

B.Meiso.PW.lm <- lm(PW ~ Genotypes, data = Pheno.B.Meiso)
B.Meiso.PW.lsmean <- as.data.frame(lsmeans(B.Meiso.PW.lm, ~ Genotypes))
hist(B.Meiso.PW.lsmean[,2])

B.Meiso.PL.lm <- lm(PL ~ Genotypes, data = Pheno.B.Meiso)
B.Meiso.PL.lsmean <- as.data.frame(lsmeans(B.Meiso.PL.lm, ~ Genotypes))
hist(B.Meiso.PL.lsmean[,2])

B.Meiso.PWD.lm <- lm(PWD ~ Genotypes, data = Pheno.B.Meiso)
B.Meiso.PWD.lsmean <- as.data.frame(lsmeans(B.Meiso.PWD.lm, ~ Genotypes))
hist(B.Meiso.PWD.lsmean[,2])

B.Meiso.GYP.lm <- lm(GYP ~ Genotypes, data = Pheno.B.Meiso)
B.Meiso.GYP.lsmean <- as.data.frame(lsmeans(B.Meiso.GYP.lm, ~ Genotypes))
hist(B.Meiso.GYP.lsmean[,2])

B.Meiso.TSW.lm <- lm(TSW ~ Genotypes, data = Pheno.B.Meiso)
B.Meiso.TSW.lsmean <- as.data.frame(lsmeans(B.Meiso.TSW.lm, ~ Genotypes))
hist(B.Meiso.TSW.lsmean[,2])

B.Meiso.Trait.LSMeans <- cbind(B.Meiso.DF.lsmean[,1:2], B.Meiso.DM.lsmean[,2],
                               B.Meiso.GYP.lsmean[,2], B.Meiso.HE.lsmean[,2],
                               B.Meiso.LS.lsmean[,2], B.Meiso.NL.lsmean[,2],
                               B.Meiso.NPP.lsmean[,2], B.Meiso.NT.lsmean[,2],
                               B.Meiso.PH.lsmean[,2], B.Meiso.PL.lsmean[,2],
                               B.Meiso.PW.lsmean[,2], B.Meiso.PWD.lsmean[,2],
                               B.Meiso.TSW.lsmean[,2])
colnames(B.Meiso.Trait.LSMeans) <- c("Genotypes", "Meiso.DF", "Meiso.DM", "Meiso.GYP", "Meiso.HE",
                                     "Meiso.LS", "Meiso.NL", "Meiso.NPP", "Meiso.NT", "Meiso.PH",
                                     "Meiso.PL", "Meiso.PW", "Meiso.PWD", "Meiso.TSW")
write.csv(B.Meiso.Trait.LSMeans, "mapB Meiso trait lsmeans.csv")






# Calculate LSMeans for each population, at each location (Kobo and Meiso).
## Population B, kobo
A.Kobo.DF.lm <- lm(DF ~ Genotypes, data = Pheno.A.Kobo)
A.Kobo.DF.lsmean <- as.data.frame(lsmeans(A.Kobo.DF.lm, ~ Genotypes))
hist(A.Kobo.DF.lsmean[,2])

A.Kobo.DM.lm <- lm(DM ~ Genotypes, data = Pheno.A.Kobo)
A.Kobo.DM.lsmean <- as.data.frame(lsmeans(A.Kobo.DM.lm, ~ Genotypes))
hist(A.Kobo.DM.lsmean[,2])

A.Kobo.PH.lm <- lm(PH ~ Genotypes, data = Pheno.A.Kobo)
A.Kobo.PH.lsmean <- as.data.frame(lsmeans(A.Kobo.PH.lm, ~ Genotypes))
hist(A.Kobo.PH.lsmean[,2])

A.Kobo.NT.lm <- lm(NT ~ Genotypes, data = Pheno.A.Kobo)
A.Kobo.NT.lsmean <- as.data.frame(lsmeans(A.Kobo.NT.lm, ~ Genotypes))
hist(A.Kobo.NT.lsmean[,2])

A.Kobo.NL.lm <- lm(NL ~ Genotypes, data = Pheno.A.Kobo)
A.Kobo.NL.lsmean <- as.data.frame(lsmeans(A.Kobo.NL.lm, ~ Genotypes))
hist(A.Kobo.NL.lsmean[,2])

A.Kobo.HE.lm <- lm(HE ~ Genotypes, data = Pheno.A.Kobo)
A.Kobo.HE.lsmean <- as.data.frame(lsmeans(A.Kobo.HE.lm, ~ Genotypes))
hist(A.Kobo.HE.lsmean[,2]) # All genotypes have value of 0, NO variation for head exsertion.

A.Kobo.NPP.lm <- lm(NPP ~ Genotypes, data = Pheno.A.Kobo)
A.Kobo.NPP.lsmean <- as.data.frame(lsmeans(A.Kobo.NPP.lm, ~ Genotypes))
hist(A.Kobo.NPP.lsmean[,2])

A.Kobo.LS.lm <- lm(LS ~ Genotypes, data = Pheno.A.Kobo)
A.Kobo.LS.lsmean <- as.data.frame(lsmeans(A.Kobo.LS.lm, ~ Genotypes))
hist(A.Kobo.LS.lsmean[,2])

A.Kobo.PW.lm <- lm(PW ~ Genotypes, data = Pheno.A.Kobo)
A.Kobo.PW.lsmean <- as.data.frame(lsmeans(A.Kobo.PW.lm, ~ Genotypes))
hist(A.Kobo.PW.lsmean[,2])

A.Kobo.PL.lm <- lm(PL ~ Genotypes, data = Pheno.A.Kobo)
A.Kobo.PL.lsmean <- as.data.frame(lsmeans(A.Kobo.PL.lm, ~ Genotypes))
hist(A.Kobo.PL.lsmean[,2])

A.Kobo.PWD.lm <- lm(PWD ~ Genotypes, data = Pheno.A.Kobo)
A.Kobo.PWD.lsmean <- as.data.frame(lsmeans(A.Kobo.PWD.lm, ~ Genotypes))
hist(A.Kobo.PWD.lsmean[,2])

A.Kobo.GYP.lm <- lm(GYP ~ Genotypes, data = Pheno.A.Kobo)
A.Kobo.GYP.lsmean <- as.data.frame(lsmeans(A.Kobo.GYP.lm, ~ Genotypes))
hist(A.Kobo.GYP.lsmean[,2])

A.Kobo.TSW.lm <- lm(TSW ~ Genotypes, data = Pheno.A.Kobo)
A.Kobo.TSW.lsmean <- as.data.frame(lsmeans(A.Kobo.TSW.lm, ~ Genotypes))
hist(A.Kobo.TSW.lsmean[,2])

A.Kobo.Trait.LSMeans <- cbind(A.Kobo.DF.lsmean[,1:2], A.Kobo.DM.lsmean[,2],
                              A.Kobo.GYP.lsmean[,2], A.Kobo.HE.lsmean[,2],
                              A.Kobo.LS.lsmean[,2], A.Kobo.NL.lsmean[,2],
                              A.Kobo.NPP.lsmean[,2], A.Kobo.NT.lsmean[,2],
                              A.Kobo.PH.lsmean[,2], A.Kobo.PL.lsmean[,2],
                              A.Kobo.PW.lsmean[,2], A.Kobo.PWD.lsmean[,2],
                              A.Kobo.TSW.lsmean[,2])
colnames(A.Kobo.Trait.LSMeans) <- c("Genotypes", "Kobo.DF", "Kobo.DM", "Kobo.GYP", "Kobo.HE",
                                    "Kobo.LS", "Kobo.NL", "Kobo.NPP", "Kobo.NT", "Kobo.PH",
                                    "Kobo.PL", "Kobo.PW", "Kobo.PWD", "Kobo.TSW")
write.csv(A.Kobo.Trait.LSMeans, "mapA Kobo trait lsmeans.csv")





## Population A, Meiso
A.Meiso.DF.lm <- lm(DF ~ Genotypes, data = Pheno.A.Meiso)
A.Meiso.DF.lsmean <- as.data.frame(lsmeans(A.Meiso.DF.lm, ~ Genotypes))
hist(A.Meiso.DF.lsmean[,2])

A.Meiso.DM.lm <- lm(DM ~ Genotypes, data = Pheno.A.Meiso)
A.Meiso.DM.lsmean <- as.data.frame(lsmeans(A.Meiso.DM.lm, ~ Genotypes))
hist(A.Meiso.DM.lsmean[,2])

A.Meiso.PH.lm <- lm(PH ~ Genotypes, data = Pheno.A.Meiso)
A.Meiso.PH.lsmean <- as.data.frame(lsmeans(A.Meiso.PH.lm, ~ Genotypes))
hist(A.Meiso.PH.lsmean[,2])

A.Meiso.NT.lm <- lm(NT ~ Genotypes, data = Pheno.A.Meiso)
A.Meiso.NT.lsmean <- as.data.frame(lsmeans(A.Meiso.NT.lm, ~ Genotypes))
hist(A.Meiso.NT.lsmean[,2])

A.Meiso.NL.lm <- lm(NL ~ Genotypes, data = Pheno.A.Meiso)
A.Meiso.NL.lsmean <- as.data.frame(lsmeans(A.Meiso.NL.lm, ~ Genotypes))
hist(A.Meiso.NL.lsmean[,2])

A.Meiso.HE.lm <- lm(HE ~ Genotypes, data = Pheno.A.Meiso)
A.Meiso.HE.lsmean <- as.data.frame(lsmeans(A.Meiso.HE.lm, ~ Genotypes))
hist(A.Meiso.HE.lsmean[,2]) # All genotypes have value of 0, NO variation for head exsertion.

A.Meiso.NPP.lm <- lm(NPP ~ Genotypes, data = Pheno.A.Meiso)
A.Meiso.NPP.lsmean <- as.data.frame(lsmeans(A.Meiso.NPP.lm, ~ Genotypes))
hist(A.Meiso.NPP.lsmean[,2])

A.Meiso.LS.lm <- lm(LS ~ Genotypes, data = Pheno.A.Meiso)
A.Meiso.LS.lsmean <- as.data.frame(lsmeans(A.Meiso.LS.lm, ~ Genotypes))
hist(A.Meiso.LS.lsmean[,2])

A.Meiso.PW.lm <- lm(PW ~ Genotypes, data = Pheno.A.Meiso)
A.Meiso.PW.lsmean <- as.data.frame(lsmeans(A.Meiso.PW.lm, ~ Genotypes))
hist(A.Meiso.PW.lsmean[,2])

A.Meiso.PL.lm <- lm(PL ~ Genotypes, data = Pheno.A.Meiso)
A.Meiso.PL.lsmean <- as.data.frame(lsmeans(A.Meiso.PL.lm, ~ Genotypes))
hist(A.Meiso.PL.lsmean[,2])

A.Meiso.PWD.lm <- lm(PWD ~ Genotypes, data = Pheno.A.Meiso)
A.Meiso.PWD.lsmean <- as.data.frame(lsmeans(A.Meiso.PWD.lm, ~ Genotypes))
hist(A.Meiso.PWD.lsmean[,2])

A.Meiso.GYP.lm <- lm(GYP ~ Genotypes, data = Pheno.A.Meiso)
A.Meiso.GYP.lsmean <- as.data.frame(lsmeans(A.Meiso.GYP.lm, ~ Genotypes))
hist(A.Meiso.GYP.lsmean[,2])

A.Meiso.TSW.lm <- lm(TSW ~ Genotypes, data = Pheno.A.Meiso)
A.Meiso.TSW.lsmean <- as.data.frame(lsmeans(A.Meiso.TSW.lm, ~ Genotypes))
hist(A.Meiso.TSW.lsmean[,2])

A.Meiso.Trait.LSMeans <- cbind(A.Meiso.DF.lsmean[,1:2], A.Meiso.DM.lsmean[,2],
                               A.Meiso.GYP.lsmean[,2], A.Meiso.HE.lsmean[,2],
                               A.Meiso.LS.lsmean[,2], A.Meiso.NL.lsmean[,2],
                               A.Meiso.NPP.lsmean[,2], A.Meiso.NT.lsmean[,2],
                               A.Meiso.PH.lsmean[,2], A.Meiso.PL.lsmean[,2],
                               A.Meiso.PW.lsmean[,2], A.Meiso.TSW.lsmean[,2])
colnames(A.Meiso.Trait.LSMeans) <- c("Genotypes", "Meiso.DF", "Meiso.DM", "Meiso.GYP", "Meiso.HE",
                                     "Meiso.LS", "Meiso.NL", "Meiso.NPP", "Meiso.NT", "Meiso.PH",
                                     "Meiso.PL", "Meiso.PW", "Meiso.TSW")
write.csv(A.Meiso.Trait.LSMeans, "mapA Meiso trait lsmeans.csv")
write.csv(A.Meiso.PWD.lsmean, "mapA Meiso PWD lsmeans.csv")
