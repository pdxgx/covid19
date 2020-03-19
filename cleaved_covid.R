setwd("~/Documents/Coding/covid/covid19")
library(ggplot2)
#nugs
hla_cleavednugs <- read.csv("cleaved_covid_nuggets.csv")
typeof(hla_cleavednugs$file_name)
gsub('.{4}$', '', hla_cleavednugs$file_name) -> hla_cleavednugs$file_name
sub('.............', '', hla_cleavednugs$file_name) -> hla_cleavednugs$file_name
colnames(hla_cleavednugs)[1] <- 'allele'
hla_cleavednugs_filt <- hla_cleavednugs[is.na(as.numeric(as.character(hla_cleavednugs$allele))),]
hla_cleavednugs_filt <- hla_cleavednugs_filt[as.integer(hla_cleavednugs_filt$ic50) < 500,]
hla_cleavednugs_filt$less_than_50 <- ifelse(hla_cleavednugs_filt$ic50<50,'50','500')
hla_cleavednugs_filt$occ <- rep(1,nrow(hla_cleavednugs_filt))
hla_cleavednugs_filt$allele <- reorder(hla_cleavednugs_filt$allele,hla_cleavednugs_filt$allele,FUN=length)
ggplot(hla_cleavednugs_filt, aes(fill=less_than_50, y=occ, x=reorder(desc(allele), -occ))) + 
  geom_bar(position="stack", stat="identity")
table(hla_cleavednugs_filt$allele)
hla_cleavednugs_filt$class <- rep(1,nrow(hla_cleavednugs_filt))
hla_cleavednugs_filt$class <- ifelse(grepl("LA-A", hla_cleavednugs_filt$allele), 'A', ifelse(grepl("LA-B", hla_cleavednugs_filt$allele), 'B', 'C'))
ggplot(hla_cleavednugs_filt, aes(fill=class, y=occ, x=reorder(desc(allele), -occ))) + 
  geom_bar(position="stack", stat="identity") + theme_classic() + 
  theme(axis.text.x = element_blank()) + 
  scale_fill_manual(values=c("#003366","#bc5090", "#ffa600"))
ggsave("covid19_netchopfilt_035_nuggets.tiff", device="tiff")
ggsave("covid19_netchopfilt_035_nuggets.png", device="png")





#flurry
hla_cleavedflurry <- read.csv("cleaved_covid_flurry.csv")
hla_cleavedflurry_filt <- hla_cleavedflurry[as.integer(hla_cleavedflurry$mhcflurry_prediction) < 500,]
hla_cleavedflurry_filt$less_than_50 <- ifelse(hla_cleavedflurry_filt$mhcflurry_prediction<50,'50','500')
hla_cleavedflurry_filt$occ <- rep(1,nrow(hla_cleavedflurry_filt))
hla_cleavedflurry_filt$allele <- reorder(hla_cleavedflurry_filt$allele,hla_cleavedflurry_filt$allele,FUN=length)
ggplot(hla_cleavedflurry_filt, aes(fill=less_than_50, y=occ, x=reorder(desc(allele), -occ))) + 
  geom_bar(position="stack", stat="identity")
table(hla_cleavedflurry_filt$allele)
hla_cleavedflurry_filt$class <- rep(1,nrow(hla_cleavedflurry_filt))
hla_cleavedflurry_filt$class <- ifelse(grepl("LA-A", hla_cleavedflurry_filt$allele), 'A', ifelse(grepl("LA-B", hla_cleavedflurry_filt$allele), 'B', 'C'))
ggplot(hla_cleavedflurry_filt, aes(fill=class, y=occ, x=reorder(desc(allele), -occ))) + 
  geom_bar(position="stack", stat="identity") + theme_classic() + 
  theme(axis.text.x = element_blank()) + 
  scale_fill_manual(values=c("#003366","#bc5090", "#ffa600"))
ggsave("covid19_netchopfilt_035_flurry.tiff", device="tiff")
ggsave("covid19_netchopfilt_035_flurry.png", device="png")



#netmhcpan
hla_cleavednet <- read.csv("covid_scores.csv")
hla_cleavednet_filt <- hla_cleavednet[as.integer(hla_cleavednet$Binding_affinity) < 500,]
hla_cleavednet_filt$less_than_50 <- ifelse(hla_cleavednet_filt$Binding_affinity<50,'50','500')
hla_cleavednet_filt$occ <- rep(1,nrow(hla_cleavednet_filt))
colnames(hla_cleavednet_filt)[2] <- "allele"
hla_cleavednet_filt$allele <- reorder(hla_cleavednet_filt$allele,hla_cleavednet_filt$allele,FUN=length)

table(hla_cleavednet_filt$allele)
hla_cleavednet_filt$class <- rep(1,nrow(hla_cleavednet_filt))
hla_cleavednet_filt$class <- ifelse(grepl("LA-A", hla_cleavednet_filt$allele), 'A', ifelse(grepl("LA-B", hla_cleavednet_filt$allele), 'B', 'C'))
ggplot(hla_cleavednet_filt, aes(fill=class, y=occ, x=reorder(desc(allele), -occ))) + 
  geom_bar(position="stack", stat="identity") + theme_classic() + 
  theme(axis.text.x = element_blank()) + 
  scale_fill_manual(values=c("#003366","#bc5090", "#ffa600"))
ggsave("covid19_netchopfilt_035_netmhcpan.tiff", device="tiff")
ggsave("covid19_netchopfilt_035_netmhcpan.png", device="png")




#writing
write.csv(hla_cleavednugs_filt, "hla_cleavednugs_filt.csv", row.names = FALSE)
write.csv(hla_cleavedflurry_filt, "hla_cleavedflurry_filt.csv", row.names = FALSE)
write.csv(hla_cleavednet_filt, "hla_cleavednet_filt.csv", row.names = FALSE)

#50 cutoff
hla_cleavednugs_filt_50 <- hla_cleavednugs_filt[as.integer(hla_cleavednugs_filt$ic50) < 50,]
hla_cleavedflurry_filt_50 <- hla_cleavedflurry[as.integer(hla_cleavedflurry$mhcflurry_prediction) < 50,]
write.csv(hla_cleavednugs_filt_50, "hla_cleavednugs_filt_50.csv", row.names = FALSE)
write.csv(hla_cleavedflurry_filt_50, "hla_cleavedflurry_filt_50.csv", row.names = FALSE)

