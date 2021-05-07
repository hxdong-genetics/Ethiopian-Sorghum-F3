setwd("C:/Users/hdong/Desktop/ABHGenotypeR")

library(ABHgenotypeR)

# Raw data
mapa <- readABHgenotypes("ABH mapA snp data.csv", nameA = "Baji", nameB = "76T1-23", readPos = TRUE)
mapb <- readABHgenotypes("ABH mapB snp data.csv", nameA = "Birmash", nameB = "Meco-1", readPos = TRUE)
mapc <- readABHgenotypes("ABH mapC snp data.csv", nameA = "Birmash", nameB = "76T1-23", readPos = TRUE)

# Correct short miscalled stretches based on flanking alleles
mapa.correct <- correctStretches(mapa, maxHapLength = 1)
mapb.correct <- correctStretches(mapb, maxHapLength = 1)
mapc.correct <- correctStretches(mapc, maxHapLength = 1)

# Correct undercalled heterozygous sites based on flanking alleles
mapa.correct.het <- correctUndercalledHets(mapa.correct, maxHapLength = 1)
mapb.correct.het <- correctUndercalledHets(mapb.correct, maxHapLength = 1)
mapc.correct.het <- correctUndercalledHets(mapc.correct, maxHapLength = 1)

# Impute missing genotypes based on flanking alleles
mapa.imp <- imputeByFlanks(mapa.correct.het)
mapb.imp <- imputeByFlanks(mapb.correct.het)
mapc.imp <- imputeByFlanks(mapc.correct.het)

#plotGenos(mapa.imp)
#plotGenos(mapb.imp)
#plotGenos(mapc.imp)

#plotMarkerDensity(mapa.imp)
#plotMarkerDensity(mapb.imp)
#plotMarkerDensity(mapc.imp)

# Plot genotype data overview across all markers and all individuals
tiff("Population A genotype data overview.tiff", height = 3, width = 6, units = 'in', res=300)
plotGenos(mapa.imp)
dev.off()
par(mfrow = c(1,1))

tiff("Population B genotype data overview.tiff", height = 3, width = 6, units = 'in', res=300)
plotGenos(mapb.imp)
dev.off()
par(mfrow = c(1,1))

tiff("Population C genotype data overview.tiff", height = 3, width = 6, units = 'in', res=300)
plotGenos(mapc.imp)
dev.off()
par(mfrow = c(1,1))


# Plot marker density
tiff("Population A marker density.tiff", height = 9, width = 5, units = 'in', res=300)
plotMarkerDensity(mapa.imp)
dev.off()
par(mfrow = c(1,1))

tiff("Population B marker density.tiff", height = 9, width = 5, units = 'in', res=300)
plotMarkerDensity(mapb.imp)
dev.off()
par(mfrow = c(1,1))

tiff("Population C marker density.tiff", height = 9, width = 5, units = 'in', res=300)
plotMarkerDensity(mapc.imp)
dev.off()
par(mfrow = c(1,1))



# Plot allele frequency
tiff("Population A allele frequency.tiff", height = 9, width = 5, units = 'in', res=300)
plotAlleleFreq(mapa.imp)
dev.off()
par(mfrow = c(1,1))

tiff("Population B allele frequency.tiff", height = 9, width = 5, units = 'in', res=300)
plotAlleleFreq(mapb.imp)
dev.off()
par(mfrow = c(1,1))

tiff("Population C allele frequency.tiff", height = 9, width = 5, units = 'in', res=300)
plotAlleleFreq(mapc.imp)
dev.off()
par(mfrow = c(1,1))

