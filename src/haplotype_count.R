setwd("~/Documents/test/files/bindingfiles")

haplotypes <- read.csv("globalhaps.csv")
genotypes <- read.csv("individual_data.csv")

hla_cleavednet <- read.csv("covid_netmhcpan_scores_all_alleles.csv")
hla_cleavednet_filt <- hla_cleavednet[as.integer(hla_cleavednet$Binding_affinity) < 500,]

haplotypes$A <- paste("HLA-", gsub("\\*","",haplotypes$A), sep="")
haplotypes$B <- paste("HLA-", gsub("\\*","",haplotypes$B), sep="")
haplotypes$C <- paste("HLA-", gsub("\\*","",haplotypes$C), sep="")

library(dplyr)
library(ggplot2)
#haplotypes$count <- length(unique(hla_cleavednet_filt$Peptide[haplotypes$A %in% hla_cleavednet_filt$Allele, ]))

haplotypes[500, 1:3]
haplotypes$count <- rep(1,nrow(haplotypes))
for (i in 1:dim(haplotypes)[1])
{
  haplotypes$count[i] <- length(unique(union(union(hla_cleavednet_filt[hla_cleavednet_filt$Allele == haplotypes[i,1], ]$Peptide, 
                              hla_cleavednet_filt[hla_cleavednet_filt$Allele == haplotypes[i,2], ]$Peptide), 
                              hla_cleavednet_filt[hla_cleavednet_filt$Allele == haplotypes[i,3], ]$Peptide)))/8409
}

#haplotypes$count <- length(unique(union(union(hla_cleavednet_filt[hla_cleavednet_filt$Allele == haplotypes[500,1], ]$Peptide, 
#      hla_cleavednet_filt[hla_cleavednet_filt$Allele == haplotypes[500,2], ]$Peptide), 
#      hla_cleavednet_filt[hla_cleavednet_filt$Allele == haplotypes[500,3], ]$Peptide)))
haplotypes_4 <- haplotypes
haplotypes_4 <- subset(haplotypes_4, nchar(as.character(A)) == 10)
haplotypes_4 <- subset(haplotypes_4, nchar(as.character(B)) == 10)
haplotypes_4 <- subset(haplotypes_4, nchar(as.character(C)) == 10)

haplotypes_4 <- haplotypes_4[]

haplotypes_4 <- haplotypes_4[haplotypes_4$A %in% hla_cleavednet_filt$Allele,]
haplotypes_4 <- haplotypes_4[haplotypes_4$B %in% hla_cleavednet_filt$Allele,]
haplotypes_4 <- haplotypes_4[haplotypes_4$C %in% hla_cleavednet_filt$Allele,]

haplotypes_4 <- 

#intersect(haplotypes_4$A, hla_cleavednet_filt$Allele)

#qplot(haplotypes$A, haplotypes$count)

ggplot(haplotypes_4, aes(x=A, y=count)) + geom_boxplot()  + coord_flip()
ggplot(haplotypes_4, aes(x=B, y=count)) + geom_boxplot() + coord_flip()
ggplot(haplotypes_4, aes(x=C, y=count)) + geom_boxplot() + coord_flip()



ggplot(haplotypes_4, aes(x=count, y=freq)) + geom_jitter()




genotypes$A1 <- paste("HLA-", gsub("\\*","",genotypes$A1), sep="")
genotypes$B1 <- paste("HLA-", gsub("\\*","",genotypes$B1), sep="")
genotypes$C1 <- paste("HLA-", gsub("\\*","",genotypes$C1), sep="")
genotypes$A2 <- paste("HLA-", gsub("\\*","",genotypes$A2), sep="")
genotypes$B2 <- paste("HLA-", gsub("\\*","",genotypes$B2), sep="")
genotypes$C2 <- paste("HLA-", gsub("\\*","",genotypes$C2), sep="")

genotypes <- subset(genotypes, nchar(as.character(A1)) == 10)
genotypes <- subset(genotypes, nchar(as.character(B1)) == 10)
genotypes <- subset(genotypes, nchar(as.character(C1)) == 10)
genotypes <- subset(genotypes, nchar(as.character(A2)) == 10)
genotypes <- subset(genotypes, nchar(as.character(B2)) == 10)
genotypes <- subset(genotypes, nchar(as.character(C2)) == 10)

