setwd("~/Documents/Coding/covid/covid19")
library(ggplot2)
#nugs
hla_nugs <- read.csv("mhc_nug_conserved.csv")
typeof(hla_nugs$file_name)
gsub('.{4}$', '', hla_nugs$file_name) -> hla_nugs$file_name
sub('..........', '', hla_nugs$file_name) -> hla_nugs$file_name
colnames(hla_nugs)[1] <- 'allele'
hla_nugs_filt <- hla_nugs[is.na(as.numeric(as.character(hla_nugs$allele))),]
hla_nugs_filt <- hla_nugs_filt[as.integer(hla_nugs_filt$ic50) < 500,]
hla_nugs_filt$less_than_50 <- ifelse(hla_nugs_filt$ic50<50,'50','500')
hla_nugs_filt$occ <- rep(1,nrow(hla_nugs_filt))
hla_nugs_filt$allele <- reorder(hla_nugs_filt$allele,hla_nugs_filt$allele,FUN=length)
ggplot(hla_nugs_filt, aes(fill=less_than_50, y=occ, x=reorder(allele, -occ))) + 
  geom_bar(position="stack", stat="identity")

write.csv(table(hla_nugs_filt$allele), "conserved_table_nuggets_covid.csv", row.names = FALSE)
#flurry
hla_flurry_conserved <- read.csv("mhc_flurry_conserved.csv")
hla_flurry_conserved_filt <- hla_flurry_conserved[as.integer(hla_flurry_conserved$mhcflurry_prediction) < 500,]
hla_flurry_conserved_filt$less_than_50 <- ifelse(hla_flurry_conserved_filt$mhcflurry_prediction<50,'50','500')
hla_flurry_conserved_filt$occ <- rep(1,nrow(hla_flurry_conserved_filt))
hla_flurry_conserved_filt$allele <- reorder(hla_flurry_conserved_filt$allele,hla_flurry_conserved_filt$allele,FUN=length)
ggplot(hla_flurry_conserved_filt, aes(fill=less_than_50, y=occ, x=reorder(desc(allele), -occ))) + 
  geom_bar(position="stack", stat="identity")

write.csv(table(hla_flurry_conserved_filt$allele), "conserved_table_flurry_covid.csv", row.names = FALSE)

write.csv(hla_flurry_conserved_filt, "covid_flurry_conserved.csv", row.names=FALSE)


