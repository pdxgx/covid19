## TCGA_filtered$Allele is a data frame column with the allele associated with an epitope, each row in the data frame represents an epitope
## allele_freqs is a data frame column with one row per HLA allele, where the 'Allele' column contains the allele name and the 'Frequency' column contains the frequency
TCGA_filtered <- read.csv("allele_freqs.tsv", sep='\t')
TCGA_alleles_filtered <- as.data.frame(rev(sort(table(TCGA_filtered$allele))))
colors <- NULL
A_counts <- NULL
B_counts <- NULL
C_counts <- NULL
for (i in 1:length(TCGA_alleles_filtered$Var1)){
  if (grepl("A", TCGA_alleles_filtered$Var1[i]) == TRUE){
    colors <- append(colors, "green")
    A_counts <- append(A_counts, TCGA_alleles_filtered$Freq[i])
  } else if (grepl("B", TCGA_alleles_filtered$Var1[i]) == TRUE) {
    colors <- append(colors, "blue")
    B_counts <- append(B_counts, TCGA_alleles_filtered$Freq[i])
  } else if (grepl("C", TCGA_alleles_filtered$Var1[i]) == TRUE){
    colors <- append(colors, "orange")
    C_counts <- append(C_counts, TCGA_alleles_filtered$Freq[i])
  }
}
TCGA_alleles_filtered$Allele_colors <- colors

TCGA_alleles_filtered <- TCGA_alleles_filtered[order(TCGA_alleles_filtered$Var1),]





frequencies <- NULL
allele_freqs <- read.table("allele_freqs.tsv", header=FALSE, sep="\t")
colnames(allele_freqs) <- c("Allele", "Frequency")
for (i in 1:length(TCGA_alleles_filtered$Var1)){
  allele_f <- allele_freqs$Frequency[as.character(allele_freqs$Allele) == TCGA_alleles_filtered$Var1[i]]
  frequencies <- append(frequencies, allele_f)
}
TCGA_alleles_filtered$Freq2 <- frequencies

layout(matrix(c(1,2), nrow = 2, byrow=TRUE))
par(mar=c(0,6,2,0.5))
par(lty = 0)
barplot(TCGA_alleles_filtered$Freq, col=TCGA_alleles_filtered$Allele_colors, axes = F, ylim=c(0,70000), ylab="", cex.lab=1.2)
par(lty = "solid")
mtext("Novel binding neoepitopes (#)", side=2, line=4, cex=1.5)
#legend("topright", c("HLA-A", "HLA-B", "HLA-C"), fill=c("green", "blue", "orange"), cex=0.9)
rect(xleft=125, xright=175, ytop=69000, ybottom=54000, col="white", border="black")
rect(xleft=130, xright=140, ytop=68000, ybottom=64500, col="green", border="black")
rect(xleft=130, xright=140, ytop=63250, ybottom=59750, col="blue", border="black")
rect(xleft=130, xright=140, ytop=58500, ybottom=55000, col="orange", border="black")
text(157, 66200, labels="HLA-A",cex=1.5, col="black")
text(157, 61450, labels="HLA-B",cex=1.5, col="black")
text(157, 56700, labels="HLA-C",cex=1.5, col="black")
abline(h=0)
axis(2, lwd = 0, tick = FALSE, at = seq(0,70000,10000), labels = seq(0,70000,10000), las = 2, cex.axis=1)
box(lty="solid", col="black")
par(mar=c(4,6,0.9,0.5))
par(lty = 0)
barplot(-100*TCGA_alleles_filtered$Freq2, col=TCGA_alleles_filtered$Allele_colors, names.arg="", xlab="", axes = FALSE, ylim=c(-20, 0), ylab="", cex.lab=1.2)
par(lty = "solid")
axis(2, lwd = 0, tick = FALSE, at = seq(-20,0,5), labels = rev(c(0, 5, 10, 15, 20)), las = 2, cex.axis=1)
mtext("Allele frequency (%)", side=2, line=4, cex=1.5)
#mtext("HLA-A*02:01", side=1, line=-1.5, at=15, cex=0.5)
text(18, -17, labels="HLA-A*02:01",cex=1, col="black")
text(53, -12, labels="HLA-C*07:02",cex=1, col="black")
text(91, -12, labels="HLA-A*24:02",cex=1, col="black")
text(162, -12, labels="HLA-C*04:01",cex=1, col="black")
mtext("Allele", side=1, line=1, at=90, cex=1.5)
box(lty="solid", col="black")
par(mar=c(5.1,4.1,4.1,2.1))
