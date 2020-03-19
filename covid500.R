setwd("~/Documents/Coding/covid/covid19")
hla <- read.csv("covid_hla_filt_nuggetswflurryalleles.csv")
typeof(hla$file_name)
gsub('.{4}$', '', hla$file_name) -> hla$file_name
sub('..........', '', hla$file_name) -> hla$file_name
colnames(hla)[1] <- 'allele'
hla_filt <- hla[is.na(as.numeric(as.character(hla$allele))),]
hla_filt <- hla_filt[as.integer(hla_filt$ic50) < 500,]
hla_filt$less_than_50 <- ifelse(hla_filt$ic50<50,'50','500')
hla_filt$occ <- rep(1,nrow(hla_filt))
#hla_50 <- hla_filt[as.integer(hla_filt$ic50) < 50,]
#hla_450 <- hla_filt[as.integer(hla_filt$ic50) > 50,]

typeof(as.character(hla_filt$allele))
library(ggplot2)

table(hla_filt$allele)

#hla_50_tab <- data.frame(table(hla_50$allele))
#hla_450_tab <- data.frame(table(hla_450$allele))
#merged <- merge(hla_50_tab, hla_450_tab, by='Var1')


library(forcats)
hla_filt$allele <- reorder(hla_filt$allele,hla_filt$allele,FUN=length)
# Grouped
#ggplot(data, aes(fill=condition, y=value, x=specie)) + 
#  geom_bar(position="stack", stat="identity")

ggplot(hla_filt, aes(fill=less_than_50, y=occ, x=reorder(allele, -occ))) + 
  geom_bar(position="stack", stat="identity")

table(hla_filt$allele) -> hla_list2

write.csv(rownames(hla_list2), "mhc_nugs_order.csv")
write.csv(hla_filt, file="peps_nugs_filt.csv", row.names=FALSE)
