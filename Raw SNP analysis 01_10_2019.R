setwd("C:/Users/hdong/OneDrive - University of Georgia/Postdoc Projects/Ethiopia_Mapping_Pops/2_Combined")
# Load RData
load(".RData")

####################################################
# import hmp (genotype) file
SNP = read.table("PopABC_SNP.hmp.txt", header=T,sep="\t",comment.char="",stringsAsFactors=F )
# Remove mono-morphic markers
SNP2 <- SNP[c(SNP[,2]=='A' | SNP[,2]=='G' | SNP[,2]=='C' | SNP[,2]=='T'),]
SNP3 <- SNP[c(!SNP[,2]%in%SNP2[,2]),] # Retain 11327 SNPs

SNP_Numeric=read.table("PopABC_SNP_Numeric.txt", header = TRUE)
SNP_Numeric=t(SNP_Numeric)
SNP_Numeric=SNP_Numeric[c(rownames(SNP_Numeric)%in%SNP3[,1]),]
dim(SNP_Numeric)


SNP3 <- SNP3[,-c(5:11)] # Remove useless columns

# Remove one in pairs of adjacent SNPs (< 150 bp)
dist_bp = c()
for (i in 2:nrow(SNP3)){
  dist_bp[i] = SNP3[i,4] - SNP3[i-1,4]
  }
dist_bp[1] <- 80989

SNP3 <- SNP3[c(dist_bp > 150 | dist_bp < 0),] # Retain 6819 SNPs
SNP_Numeric <- SNP_Numeric[c(rownames(SNP_Numeric)%in%SNP3[,1]),]
dim(SNP_Numeric)

# Remove individuals with high missing data rate (> 20%)
miss_by_ind <- apply(SNP3[,5:ncol(SNP3)], 
                     2, function(x) length(grep("N", x)))
hist(miss_by_ind/nrow(SNP3)) # plot missing 
table(miss_by_ind/nrow(SNP3) < 0.2)
SNP4 <- cbind(SNP3[,1:4], SNP3[,5:ncol(SNP3)][,miss_by_ind/nrow(SNP3)<0.2])

SNP_Numeric <- SNP_Numeric[,c(colnames(SNP_Numeric)%in%colnames(SNP4))]


# Remove SNPs with high missing data rate (> 20%)
miss_by_snp <- apply(SNP4[,5:ncol(SNP4)], 
                     1, function(x) length(grep("N", x)))
hist(miss_by_snp/575) # plot missing 
table(miss_by_snp/575 < 0.1)
SNP5 <- SNP4[miss_by_snp/575 < 0.1,]

SNP_Numeric <- SNP_Numeric[rownames(SNP_Numeric)%in%SNP5[,1],]
dim(SNP_Numeric)


T.Individuals <- t(SNP_Numeric)

# replaces all NA values with last non-NA values
na.lomf <- function(x) {
  
  na.lomf.0 <- function(x) {
    non.na.idx <- which(!is.na(x))
    if (is.na(x[1L])) {
      non.na.idx <- c(1L, non.na.idx)
    }
    rep.int(x[non.na.idx], diff(c(non.na.idx, length(x) + 1L)))
  }
  
  dim.len <- length(dim(x))
  
  if (dim.len == 0L) {
    na.lomf.0(x)
  } else {
    apply(x, dim.len, na.lomf.0)
  }
}

T.Individuals.noNA <- na.lomf(T.Individuals)
rownames(T.Individuals.noNA) <- rownames(T.Individuals)
  
pca_matrix = scale(T.Individuals)
eig_value = prcomp(T.Individuals.noNA, center = T) 
scores_pca = eig_value$x 
names(eig_value)

Pop.name <- c(rep("A",162),rep("B",205),rep("C",208))
PC14 <- cbind(scores_pca[,1:4], as.matrix(Pop.name))

plot(PC14[,1], PC14[,2], pch=19, cex=1, xlab = "PC1", ylab = "PC2",
     col=ifelse(PC14[,5]=="A", "orange", ifelse(PC14[,5]=="B", "blue", "green"))) 



# Split SNP3 into three populations
PopA.SNP <- SNP3[,c(1:216)]
# Remove SNPs that are homozygous for the same allele in both parents
PopA.SNP <- PopA.SNP[!c(PopA.SNP[,5]=='A' & PopA.SNP[,6]=='A' | PopA.SNP[,5]=='G' & PopA.SNP[,6]=='G' |
                          PopA.SNP[,5]=='C' & PopA.SNP[,6]=='C' | PopA.SNP[,5]=='T' & PopA.SNP[,6]=='T' |
                          PopA.SNP[,5]=='M' & PopA.SNP[,6]=='M' | PopA.SNP[,5]=='R' & PopA.SNP[,6]=='R' |
                          PopA.SNP[,5]=='W' & PopA.SNP[,6]=='W' | PopA.SNP[,5]=='S' & PopA.SNP[,6]=='S' |
                          PopA.SNP[,5]=='Y' & PopA.SNP[,6]=='Y' | PopA.SNP[,5]=='K' & PopA.SNP[,6]=='K' |
                          PopA.SNP[,5]=='N' & PopA.SNP[,6]=='N'),]
# Remove individuals with high missing data rate
miss_by_ind <- apply(PopA.SNP[,5:ncol(PopA.SNP)], 
                     2, function(x) length(grep("N", x)))
hist(miss_by_ind/nrow(PopA.SNP)) # plot missing 
table(miss_by_ind/nrow(PopA.SNP) < 0.4)
PopA.SNP2 <- cbind(PopA.SNP[,1:4], PopA.SNP[,5:ncol(PopA.SNP)][,miss_by_ind/nrow(PopA.SNP)<0.4])


PopB.SNP <- SNP3[,c(1:4,217:426)]
# Remove SNPs that are homozygous for the same allele in both parents
PopB.SNP <- PopB.SNP[!c(PopB.SNP[,5]=='A' & PopB.SNP[,6]=='A' | PopB.SNP[,5]=='G' & PopB.SNP[,6]=='G' |
                          PopB.SNP[,5]=='C' & PopB.SNP[,6]=='C' | PopB.SNP[,5]=='T' & PopB.SNP[,6]=='T' |
                          PopB.SNP[,5]=='M' & PopB.SNP[,6]=='M' | PopB.SNP[,5]=='R' & PopB.SNP[,6]=='R' |
                          PopB.SNP[,5]=='W' & PopB.SNP[,6]=='W' | PopB.SNP[,5]=='S' & PopB.SNP[,6]=='S' |
                          PopB.SNP[,5]=='Y' & PopB.SNP[,6]=='Y' | PopB.SNP[,5]=='K' & PopB.SNP[,6]=='K' |
                          PopB.SNP[,5]=='N' & PopB.SNP[,6]=='N'),]
# Remove individuals with high missing data rate
miss_by_ind <- apply(PopB.SNP[,5:ncol(PopB.SNP)], 
                     2, function(x) length(grep("N", x)))
hist(miss_by_ind/nrow(PopB.SNP)) # plot missing 
table(miss_by_ind/nrow(PopB.SNP) < 0.3)
PopB.SNP2 <- cbind(PopB.SNP[,1:4], PopB.SNP[,5:ncol(PopB.SNP)][,miss_by_ind/nrow(PopB.SNP)<0.3])


PopC.SNP <- SNP3[,c(1:4,427:638)]
# Remove SNPs that are homozygous for the same allele in both parents
PopC.SNP <- PopC.SNP[!c(PopC.SNP[,5]=='A' & PopC.SNP[,6]=='A' | PopC.SNP[,5]=='G' & PopC.SNP[,6]=='G' |
                          PopC.SNP[,5]=='C' & PopC.SNP[,6]=='C' | PopC.SNP[,5]=='T' & PopC.SNP[,6]=='T' |
                          PopC.SNP[,5]=='M' & PopC.SNP[,6]=='M' | PopC.SNP[,5]=='R' & PopC.SNP[,6]=='R' |
                          PopC.SNP[,5]=='W' & PopC.SNP[,6]=='W' | PopC.SNP[,5]=='S' & PopC.SNP[,6]=='S' |
                          PopC.SNP[,5]=='Y' & PopC.SNP[,6]=='Y' | PopC.SNP[,5]=='K' & PopC.SNP[,6]=='K' |
                          PopC.SNP[,5]=='N' & PopC.SNP[,6]=='N'),]
# Remove individuals with high missing data rate
miss_by_ind <- apply(PopC.SNP[,5:ncol(PopC.SNP)], 
                     2, function(x) length(grep("N", x)))
hist(miss_by_ind/nrow(PopC.SNP)) # plot missing 
table(miss_by_ind/nrow(PopC.SNP) < 0.3)
PopC.SNP2 <- cbind(PopC.SNP[,1:4], PopC.SNP[,5:ncol(PopC.SNP)][,miss_by_ind/nrow(PopC.SNP)<0.3])

