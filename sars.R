setwd("~/Documents/Coding/covid/covid19")
library(ggplot2)
#nugs
hla_nugs <- read.csv("nugget_sars.csv")
typeof(hla_nugs$file_name)
gsub('.{4}$', '', hla_nugs$file_name) -> hla_nugs$file_name
sub('.........', '', hla_nugs$file_name) -> hla_nugs$file_name
colnames(hla_nugs)[1] <- 'allele'
hla_nugs_filt <- hla_nugs[is.na(as.numeric(as.character(hla_nugs$allele))),]
hla_nugs_filt <- hla_nugs_filt[as.integer(hla_nugs_filt$ic50) < 500,]
hla_nugs_filt$less_than_50 <- ifelse(hla_nugs_filt$ic50<50,'50','500')
hla_nugs_filt$occ <- rep(1,nrow(hla_nugs_filt))
hla_nugs_filt$allele <- reorder(hla_nugs_filt$allele,hla_nugs_filt$allele,FUN=length)
ggplot(hla_nugs_filt, aes(fill=less_than_50, y=occ, x=reorder(desc(allele), -occ))) + 
  geom_bar(position="stack", stat="identity")
table(hla_nugs_filt$allele)
#flurry
hla_flurry <- read.csv("flurry_sars.csv")
hla_flurry_filt <- hla_flurry[as.integer(hla_flurry$mhcflurry_prediction) < 500,]
hla_flurry_filt$less_than_50 <- ifelse(hla_flurry_filt$mhcflurry_prediction<50,'50','500')
hla_flurry_filt$occ <- rep(1,nrow(hla_flurry_filt))
hla_flurry_filt$allele <- reorder(hla_flurry_filt$allele,hla_flurry_filt$allele,FUN=length)
ggplot(hla_flurry_filt, aes(fill=less_than_50, y=occ, x=reorder(desc(allele), -occ))) + 
  geom_bar(position="stack", stat="identity") 
table(hla_flurry_filt$allele)

