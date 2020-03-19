setwd("~/Documents/Coding/covid/covid19")
library(ggplot2)
library(dplyr)
#nugs
hla_cleavednugs <- read.csv("covid_sort.csv")
typeof(hla_cleavednugs$file_name)
gsub('.{4}$', '', hla_cleavednugs$file_name) -> hla_cleavednugs$file_name
sub('.............', '', hla_cleavednugs$file_name) -> hla_cleavednugs$file_name
colnames(hla_cleavednugs)[1] <- 'allele'
hla_cleavednugs_filt <- hla_cleavednugs[is.na(as.numeric(as.character(hla_cleavednugs$allele))),]
hla_cleavednugs_filt <- hla_cleavednugs_filt[as.integer(hla_cleavednugs_filt$ic50) < 500,]
hla_cleavednugs_filt$less_than_50 <- ifelse(hla_cleavednugs_filt$ic50<50,'50','500')
hla_cleavednugs_filt$occ <- rep(1,nrow(hla_cleavednugs_filt))
hla_cleavednugs_filt$allele <- reorder(hla_cleavednugs_filt$allele,hla_cleavednugs_filt$allele,FUN=length)
table(hla_cleavednugs_filt$allele)
hla_cleavednugs_filt$class <- rep(1,nrow(hla_cleavednugs_filt))
hla_cleavednugs_filt$class <- ifelse(grepl("LA-A", hla_cleavednugs_filt$allele), 'A', ifelse(grepl("LA-B", hla_cleavednugs_filt$allele), 'B', 'C'))
ggplot(hla_cleavednugs_filt, aes(fill=class, y=occ, x=reorder(desc(allele), -occ))) + 
  geom_bar(position="stack", stat="identity") + theme_classic() + 
  theme(axis.text.x = element_blank()) + 
  scale_fill_manual(values=c("#003366","#bc5090", "#ffa600")) + coord_flip()
ggsave("covid19_netchopfilt_010_nuggets.tiff", device="tiff")
ggsave("covid19_netchopfilt_010_nuggets.png", device="png")





#flurry
hla_cleavedflurry <- read.csv("covid_0.10.csv")
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
  scale_fill_manual(values=c("#003366","#bc5090", "#ffa600")) + coord_flip()
ggsave("covid19_netchopfilt_010_flurry.tiff", device="tiff")
ggsave("covid19_netchopfilt_010_flurry.png", device="png")



#netmhcpan
hla_cleavednet <- read.csv("covid_netmhc_0.1.csv")
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
  scale_fill_manual(values=c("#003366","#bc5090", "#ffa600")) + coord_flip()
ggsave("covid19_netchopfilt_010_netmhcpan.tiff", device="tiff")
ggsave("covid19_netchopfilt_010_netmhcpan.png", device="png")




#writing
write.csv(hla_cleavednugs_filt, "covid0.10_cleavednugs_filt.csv", row.names = FALSE)
write.csv(hla_cleavedflurry_filt, "covid0.10_cleavedflurry_filt.csv", row.names = FALSE)
write.csv(hla_cleavednet_filt, "hla_cleavednet_filt.csv", row.names = FALSE)

write.csv(gsub(".*HLA-", "", rownames(table(hla_cleavednugs_filt$allele))), "allele_nug_cleave.csv", row.names = FALSE, col.names = FALSE)
write.csv(gsub(".*HLA-", "", rownames(table(hla_cleavedflurry_filt$allele))), "allele_flurry_cleave.csv", row.names = FALSE, col.names = FALSE)
write.csv(gsub(".*HLA-", "", rownames(table(hla_cleavednet_filt$allele))), "allele_net_cleave.csv", row.names = FALSE, col.names = FALSE)


#50 cutoff
hla_cleavednugs_filt_50 <- hla_cleavednugs_filt[as.integer(hla_cleavednugs_filt$ic50) < 50,]
hla_cleavedflurry_filt_50 <- hla_cleavedflurry[as.integer(hla_cleavedflurry$mhcflurry_prediction) < 50,]
hla_cleavednet_filt_50 <- hla_cleavednet[as.integer(hla_cleavednet$Binding_affinity) < 50,]

hla_cleavednugsfilt_50_table <- as.data.frame(table(hla_cleavednugs_filt_50$peptide))
hla_cleavednugsfilt_50_table_filt <- hla_cleavednugsfilt_50_table[hla_cleavednugsfilt_50_table$Freq > 5,]
hla_cleavednugsfilt_50_table_filt_order <- hla_cleavednugsfilt_50_table_filt[order(hla_cleavednugsfilt_50_table_filt$Freq, decreasing=TRUE),]

hla_cleavedflurryfilt_50_table <- as.data.frame(table(hla_cleavedflurry_filt_50$peptide))
hla_cleavedflurryfilt_50_table_filt <- hla_cleavedflurryfilt_50_table[hla_cleavedflurryfilt_50_table$Freq > 5,]
hla_cleavedflurryfilt_50_table_filt_order <- hla_cleavedflurryfilt_50_table_filt[order(hla_cleavedflurryfilt_50_table_filt$Freq, decreasing=TRUE),]


hla_cleavednetfilt_50_table <- as.data.frame(table(hla_cleavednet_filt_50$Peptide))

hla_cleavednetfilt_50_table_filt <- hla_cleavednetfilt_50_table[hla_cleavednetfilt_50_table$Freq > 5,]
hla_cleavednetfilt_50_table_filt_order <- hla_cleavednetfilt_50_table_filt[order(hla_cleavednetfilt_50_table_filt$Freq, decreasing=TRUE),]

write.csv(hla_cleavednugsfilt_50_table_filt_order, "covid_nug_peps.csv", row.names = FALSE)
write.csv(hla_cleavedflurryfilt_50_table_filt_order, "covid_flurry_peps.csv", row.names = FALSE)
write.csv(hla_cleavednetfilt_50_table_filt_order, "covid_net_peps.csv", row.names = FALSE)

merged_peps_covid <- merge(hla_cleavednugsfilt_50_table_filt_order, merge(hla_cleavedflurryfilt_50_table_filt_order, hla_cleavednetfilt_50_table_filt_order, by = "Var1"))

for (i in 1:dim(merged_peps_covid)[1])
{
  merged_peps_covid$average[i] <- mean(as.numeric(merged_peps_covid[i,2:4]))
}

merged_peps_covid <- merged_peps_covid[order(merged_peps_covid$average, decreasing=FALSE),]
write.csv(merged_peps_covid, "covid_peps_common.csv", row.names=FALSE)




#sars


setwd("~/Documents/Coding/covid/covid19")
library(ggplot2)
#nugs
sars_hla_cleavednugs <- read.csv("sars_sort.csv")
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
  scale_fill_manual(values=c("#003366","#bc5090", "#ffa600")) + coord_flip()
ggsave("sars_netchopfilt_010_nuggets.png", device="png")

#flurry

sars_hla_cleavedflurry <- read.csv("sars_0.10.csv")
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
  scale_fill_manual(values=c("#003366","#bc5090", "#ffa600")) +coord_flip()
ggsave("sars_netchopfilt_010_flurry.png", device="png")


#net
sars_hla_cleavednet <- read.csv("sars_netmhc_0.1.csv")
sars_hla_cleavednet_filt <- sars_hla_cleavednet[as.integer(sars_hla_cleavednet$Binding_affinity) < 500,]
sars_hla_cleavednet_filt$less_than_50 <- ifelse(sars_hla_cleavednet_filt$Binding_affinity<50,'50','500')
sars_hla_cleavednet_filt$occ <- rep(1,nrow(sars_hla_cleavednet_filt))
colnames(sars_hla_cleavednet_filt)[2] <- "allele"
sars_hla_cleavednet_filt$allele <- reorder(sars_hla_cleavednet_filt$allele,sars_hla_cleavednet_filt$allele,FUN=length)


sars_hla_cleavednet_filt$occ_50 <- ifelse(sars_hla_cleavednet_filt$less_than_50 == '50', 1, 0)
table(sars_hla_cleavednet_filt$allele)
sars_hla_cleavednet_filt$class <- rep(1,nrow(sars_hla_cleavednet_filt))
sars_hla_cleavednet_filt$class <- ifelse(grepl("LA-A", sars_hla_cleavednet_filt$allele), 'A', ifelse(grepl("LA-B", sars_hla_cleavednet_filt$allele), 'B', 'C'))
ggplot(sars_hla_cleavednet_filt, aes(fill=class, y=occ, x=allele)) + 
  geom_bar(position="stack", stat="identity") + theme_classic() + 
  #theme(axis.text.x = element_blank()) + 
  scale_fill_manual(values=c("#003366","#bc5090", "#ffa600")) + coord_flip()

ggplot(sars_hla_cleavednet_filt, aes(fill=class, y=occ_50, x=allele)) + 
  geom_bar(position="stack", stat="identity") + theme_classic() + 
  #theme(axis.text.x = element_blank()) + 
  scale_fill_manual(values=c("#003366","#bc5090", "#ffa600")) + coord_flip()
ggsave("sars19_netchopfilt_010_netmhcpan.tiff", device="tiff")
ggsave("sars19_netchopfilt_010_netmhcpan.png", device="png")



write.csv(sars_hla_cleavedflurry_filt, "sars_hla_cleavedflurry_filt.csv", row.names = FALSE)
write.csv(sars_hla_cleavednugs_filt, "sars_hla_cleavednugs_filt.csv", row.names = FALSE)




#50 cutoff
sars_hla_cleavednugs_filt_50 <- sars_hla_cleavednugs_filt[as.integer(sars_hla_cleavednugs_filt$ic50) < 50,]
sars_hla_cleavedflurry_filt_50 <- sars_hla_cleavedflurry[as.integer(sars_hla_cleavedflurry$mhcflurry_prediction) < 50,]
sars_hla_cleavednetfilt_50 <- sars_hla_cleavednet[as.integer(sars_hla_cleavednet$Binding_affinity) < 50,]

sars_hla_cleavednetfilt_50_table <- as.data.frame(table(sars_hla_cleavednetfilt_50$Peptide))
sars_hla_cleavednetfilt_50_table_filt <- sars_hla_cleavednetfilt_50_table[sars_hla_cleavednetfilt_50_table$Freq > 5,]
sars_hla_cleavednetfilt_50_table_filt_order <- sars_hla_cleavednetfilt_50_table_filt[order(sars_hla_cleavednetfilt_50_table_filt$Freq, decreasing=TRUE),]



write.csv(gsub(".*HLA-", "", rownames(table(sars_hla_cleavednugs_filt$allele))), "sars_allele_nug_cleave.csv", row.names = FALSE, col.names = FALSE)
write.csv(gsub(".*HLA-", "", rownames(table(sars_hla_cleavedflurry_filt$allele))), "sars_allele_flurry_cleave.csv", row.names = FALSE, col.names = FALSE)
write.csv(gsub(".*HLA-", "", rownames(table(sars_hla_cleavednet_filt$allele))), "sars_allele_net_cleave.csv", row.names = FALSE, col.names = FALSE)










###plotting


library(gridExtra)

a <- ggplot(hla_cleavednugs_filt, aes(fill=class, y=occ, x=allele)) + 
  geom_bar(position="stack", stat="identity") + theme_classic() + 
#  theme(axis.text.x = element_blank()) + 
  scale_fill_manual(values=c("#003366","#bc5090", "#ffa600")) + 
  coord_flip() + ggtitle("Covid-19 peptides binding <500nm vs allele, mhcnuggets") + scale_x_discrete(expand=c(0,0))


b <- ggplot(hla_cleavedflurry_filt, aes(fill=class, y=occ, x=allele)) + 
  geom_bar(position="stack", stat="identity") + theme_classic() + 
#  theme(axis.text.x = element_blank()) + 
  scale_fill_manual(values=c("#003366","#bc5090", "#ffa600")) + 
  coord_flip() + ggtitle("Covid-19 peptides binding <500nm vs allele, mhcflurry") + scale_x_discrete(expand=c(0,0))

c <- ggplot(hla_cleavednet_filt, aes(fill=class, y=occ, x=allele)) + 
  geom_bar(position="stack", stat="identity") + theme_classic() + 
#  theme(axis.text.x = element_blank()) + 
  scale_fill_manual(values=c("#003366","#bc5090", "#ffa600")) + 
  coord_flip() + ggtitle("Covid-19 peptides binding <500nm vs allele, netmhcpan") + scale_x_discrete(expand=c(0,0))

d <- ggplot(sars_hla_cleavednugs_filt, aes(fill=class, y=occ, x=allele)) + 
  geom_bar(position="stack", stat="identity") + theme_classic() + 
#  theme(axis.text.x = element_blank()) + 
  scale_fill_manual(values=c("#003366","#bc5090", "#ffa600")) + 
  coord_flip() + ggtitle("SARS peptides binding <500nm vs allele, mhcnuggets") + scale_x_discrete(expand=c(0,0))

e <- ggplot(sars_hla_cleavedflurry_filt, aes(fill=class, y=occ, x=allele)) + 
  geom_bar(position="stack", stat="identity") + theme_classic() + 
#  theme(axis.text.x = element_blank()) + 
  scale_fill_manual(values=c("#003366","#bc5090", "#ffa600")) + 
  coord_flip() + ggtitle("SARS peptides binding <500nm vs allele, mhcflurry") + scale_x_discrete(expand=c(0,0))

f <- ggplot(sars_hla_cleavednet_filt, aes(fill=class, y=occ, x=allele)) + 
  geom_bar(position="stack", stat="identity") + theme_classic() + 
#  theme(axis.text.x = element_blank()) + 
  scale_fill_manual(values=c("#003366","#bc5090", "#ffa600")) + 
  coord_flip() + ggtitle("SARS peptides binding <500nm vs allele, netmhcpan") + scale_x_discrete(expand=c(0,0))

grid.arrange(a,b,c,d,e,f, ncol = 2)


#plotting
sars_hla_cleavednet_filt$combine <- paste(sars_hla_cleavednet_filt$class,"<",sars_hla_cleavednet_filt$less_than_50, sep="")

ggplot(sars_hla_cleavednet_filt) + 
  geom_bar(aes(fill=combine, y=occ, x=allele), position="stack", stat="identity") + 
  theme_classic() + 
  theme(axis.title.y=element_blank()) + 
  scale_fill_manual(values=c("#0061c2", "#003366", "#d48cb7", "#bc5090", "#ffdd9e", "#ffa600")) + 
  coord_flip() + ggtitle("SARS <500nm- netmhcpan") + scale_y_discrete(expand=c(0,0))



sars_recode <- sars_hla_cleavednet_filt %>%
  arrange(desc(less_than_50))


allele_freqs <- as.data.frame(read.csv("allele_freqs.csv", header=FALSE))
colnames(allele_freqs) <- c("allele", "frequency")
allele_freqs$class <- ifelse(grepl("A", allele_freqs$allele), 'A', ifelse(grepl("B", allele_freqs$allele), 'B', 'C'))
allele_freqs$allele <- paste("HLA-", allele_freqs$allele, sep="")
allele_freqs$allele <- gsub("(\\d{2})(?=\\d{2})", "\\1:", allele_freqs$allele, perl = TRUE)
allele_freqs_filt <- subset(allele_freqs, allele_freqs$allele %in% sars_hla_cleavednet_filt$allele)

new_allele_df <- as.data.frame(table(sars_recode$allele))
colnames(new_allele_df)[1] <- "allele"
merged_allele_df <- merge(new_allele_df, allele_freqs_filt, by="allele", all=TRUE)

write.csv(merged_allele_df, "nick2.csv", row.names=FALSE)
#merged_allele_df2 <- read.csv("nick2.csv")

a <- ggplot(sars_recode, aes(x = allele, y = occ, fill = class, alpha = as.integer(less_than_50)))+ 
  geom_bar(position="stack", stat="identity")+
  theme_classic()+
  theme(axis.title.y=element_blank()) + 
  scale_fill_manual(values=c('#1B9E77', '#D95F02', '#7570B3'))+
  coord_flip()+
  ggtitle("SARS <500nm- netmhcpan")+
  scale_y_discrete(expand=c(0,0))+
  scale_alpha(range = c(0.6, 1), guide = FALSE) + 
  theme(axis_ticks = element_line(size = 2))
b <- ggplot(merged_allele_df, aes(x = allele, y = as.numeric(as.character(frequency)), fill = class))+
  geom_bar(position="stack", stat="identity")+
  theme_classic() + scale_x_discrete(position='top')+
  theme(axis.title.y=element_blank(), axis.text.y=element_blank())+
  scale_fill_manual(values=c('#1B9E77', '#D95F02', '#7570B3'))+
  coord_flip()+
  ggtitle("Allele frequencies") + scale_y_reverse() + theme(legend.position = "none")+
  theme(axis.ticks = element_line(size = 2))

grid.arrange(b,a, ncol=2)


#ggplot(sars_hla_cleavednet_filt) + 
#  geom_bar(aes(fill=less_than_50, y=occ, x=allele), position="stack", stat="identity", alpha=0.8) + 
#  theme_classic() + scale_y_reverse() + scale_x_discrete(position='top')+
#  theme(axis.title.y=element_blank(),axis.text.y=element_blank()) +
#  scale_fill_manual(values=c("#003366","#bc5090", "#ffa600")) + 
#  coord_flip() + ggtitle("SARS <50nm- netmhcpan") + 
#  theme(legend.position = "none")
  
#grid.arrange(b,a, ncol=2)











