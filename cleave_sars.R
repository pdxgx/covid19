setwd("~/Documents/Coding/covid/covid19")
library(ggplot2)
library(dplyr)
#nugs
sars_hla_cleavednugs <- read.csv("cleaved_sars_nuggets.csv")
typeof(sars_hla_cleavednugs$file_name)
gsub('.{4}$', '', sars_hla_cleavednugs$file_name) -> sars_hla_cleavednugs$file_name
sub('.............', '', sars_hla_cleavednugs$file_name) -> sars_hla_cleavednugs$file_name
colnames(sars_hla_cleavednugs)[1] <- 'allele'
sars_hla_cleavednugs_filt <- sars_hla_cleavednugs[is.na(as.numeric(as.character(sars_hla_cleavednugs$allele))),]
sars_hla_cleavednugs_filt <- sars_hla_cleavednugs_filt[as.integer(sars_hla_cleavednugs_filt$ic50) < 500,]
sars_hla_cleavednugs_filt$less_than_50 <- ifelse(sars_hla_cleavednugs_filt$ic50<50,'50','500')
sars_hla_cleavednugs_filt$occ <- rep(1,nrow(sars_hla_cleavednugs_filt))
sars_hla_cleavednugs_filt$allele <- reorder(sars_hla_cleavednugs_filt$allele,sars_hla_cleavednugs_filt$allele,FUN=length)

table(sars_hla_cleavednugs_filt$allele)
sars_hla_cleavednugs_filt$class <- rep(1,nrow(sars_hla_cleavednugs_filt))
sars_hla_cleavednugs_filt$class <- ifelse(grepl("LA-A", sars_hla_cleavednugs_filt$allele), 'A', ifelse(grepl("LA-B", sars_hla_cleavednugs_filt$allele), 'B', 'C'))
ggplot(sars_hla_cleavednugs_filt, aes(fill=class, y=occ, x=reorder(desc(allele), -occ))) + 
  geom_bar(position="stack", stat="identity") + theme_classic() + 
  theme(axis.text.x = element_blank()) + 
  scale_fill_manual(values=c("#003366","#bc5090", "#ffa600"))

sars_hla_cleavednugs_filt$newvar <- paste(sars_hla_cleavednugs_filt$less_than_50, sars_hla_cleavednugs_filt$class)
sars_hla_cleavednugs_filt$newvar[which(sars_hla_cleavednugs_filt$less_than_50 == '50')] <- '50'
ggplot(sars_hla_cleavednugs_filt, aes(fill=newvar, y=occ, x=reorder(desc(allele), -occ))) + 
  geom_bar(position="stack", stat="identity")


ggplot(sars_hla_cleavednugs_filt, aes(fill=less_than_50, y=occ, x=reorder(desc(allele), -occ))) + 
  geom_bar(position="stack", stat="identity")

barplot(table(sars_hla_cleavednugs_filt$allele))
#flurry

sars_hla_cleavedflurry <- read.csv("cleaved_sars_flurry.csv")
sars_hla_cleavedflurry_filt <- sars_hla_cleavedflurry[as.integer(sars_hla_cleavedflurry$mhcflurry_prediction) < 500,]
sars_hla_cleavedflurry_filt$less_than_50 <- ifelse(sars_hla_cleavedflurry_filt$mhcflurry_prediction<50,'50','500')
sars_hla_cleavedflurry_filt$occ <- rep(1,nrow(sars_hla_cleavedflurry_filt))
sars_hla_cleavedflurry_filt$allele <- reorder(sars_hla_cleavedflurry_filt$allele,sars_hla_cleavedflurry_filt$allele,FUN=length)
ggplot(sars_hla_cleavedflurry_filt, aes(fill=less_than_50, y=occ, x=reorder(desc(allele), -occ))) + 
  geom_bar(position="stack", stat="identity")
table(sars_hla_cleavedflurry_filt$allele)

sars_hla_cleavedflurry_filt$class <- rep(1,nrow(sars_hla_cleavedflurry_filt))
sars_hla_cleavedflurry_filt$class <- ifelse(grepl("LA-A", sars_hla_cleavedflurry_filt$allele), 'A', ifelse(grepl("LA-B", sars_hla_cleavedflurry_filt$allele), 'B', 'C'))
ggplot(sars_hla_cleavedflurry_filt, aes(fill=less_than_50, y=occ, x=reorder(desc(allele), -occ))) + 
  geom_bar(position="stack", stat="identity")
ggplot(sars_hla_cleavedflurry_filt, aes(fill=class, y=occ, x=reorder(desc(allele), -occ))) + 
  geom_bar(position="stack", stat="identity") + theme_classic() + 
  theme(axis.text.x = element_blank()) + 
  scale_fill_manual(values=c("#003366","#bc5090", "#ffa600")) + 
  scale_y_continuous(expand=c(0,0))

write.csv(sars_hla_cleavedflurry_filt, "sars_hla_cleavedflurry_filt.csv", row.names = FALSE)
write.csv(sars_hla_cleavednugs_filt, "sars_hla_cleavednugs_filt.csv", row.names = FALSE)

#50 cutoff
sars_hla_cleavednugs_filt_50 <- sars_hla_cleavednugs_filt[as.integer(sars_hla_cleavednugs_filt$ic50) < 50,]
sars_hla_cleavedflurry_filt_50 <- sars_hla_cleavedflurry[as.integer(sars_hla_cleavedflurry$mhcflurry_prediction) < 50,]
write.csv(sars_hla_cleavednugs_filt_50, "sars_hla_cleavednugs_filt_50.csv", row.names = FALSE)
write.csv(sars_hla_cleavedflurry_filt_50, "sars_hla_cleavedflurry_filt_50.csv", row.names = FALSE)
