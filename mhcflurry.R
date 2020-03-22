require(dplyr)

#SARS-CoV-2
hla_cleavedflurry <- read.csv("covid_mhcflurry.csv")
hla_cleavedflurry_filt <- hla_cleavedflurry[as.integer(hla_cleavedflurry$mhcflurry_prediction) < 500,]
hla_cleavedflurry_filt$less_than_50 <- ifelse(hla_cleavedflurry_filt$mhcflurry_prediction<50,'50','500')
hla_cleavedflurry_filt$occ <- rep(1,nrow(hla_cleavedflurry_filt))
colnames(hla_cleavedflurry_filt)[2] <- "allele"
hla_cleavedflurry_filt$allele <- reorder(hla_cleavedflurry_filt$allele,hla_cleavedflurry_filt$allele,FUN=length)

table(hla_cleavedflurry_filt$allele)
hla_cleavedflurry_filt$class <- rep(1,nrow(hla_cleavedflurry_filt))
hla_cleavedflurry_filt$class <- ifelse(grepl("LA-A", hla_cleavedflurry_filt$allele), 'A', ifelse(grepl("LA-B", hla_cleavedflurry_filt$allele), 'B', 'C'))


#flurrysars
sars_hla_cleavedflurry <- read.csv("sars_mhcflurry.csv")
sars_hla_cleavedflurry_filt <- sars_hla_cleavedflurry[as.integer(sars_hla_cleavedflurry$mhcflurry_prediction) < 500,]
sars_hla_cleavedflurry_filt$less_than_50 <- ifelse(sars_hla_cleavedflurry_filt$mhcflurry_prediction<50,'50','500')
sars_hla_cleavedflurry_filt$occ <- rep(1,nrow(sars_hla_cleavedflurry_filt))
colnames(sars_hla_cleavedflurry_filt)[2] <- "allele"
sars_hla_cleavedflurry_filt$allele <- reorder(sars_hla_cleavedflurry_filt$allele,sars_hla_cleavedflurry_filt$allele,FUN=length)
sars_hla_cleavedflurry_filt$occ_50 <- ifelse(sars_hla_cleavedflurry_filt$less_than_50 == '50', 1, 0)
table(sars_hla_cleavedflurry_filt$allele)
sars_hla_cleavedflurry_filt$class <- rep(1,nrow(sars_hla_cleavedflurry_filt))
sars_hla_cleavedflurry_filt$class <- ifelse(grepl("LA-A", sars_hla_cleavedflurry_filt$allele), 'A', ifelse(grepl("LA-B", sars_hla_cleavedflurry_filt$allele), 'B', 'C'))


#plotting
covid_recode_flurry <- hla_cleavedflurry_filt %>%
  arrange(desc(less_than_50))
covid_recode_flurry$allele <- reorder(covid_recode_flurry$allele,covid_recode_flurry$allele,FUN=length)

allele_freqs <- as.data.frame(read.csv("allele_freqs.csv", header=TRUE))
colnames(allele_freqs) <- c("allele", "frequency")
allele_freqs$class <- ifelse(grepl("A", allele_freqs$allele), 'A', ifelse(grepl("B", allele_freqs$allele), 'B', 'C'))
allele_freqs$allele <- paste("HLA-", allele_freqs$allele, sep="")
allele_freqs$allele <- gsub("(\\d{2})(?=\\d{2})", "\\1:", allele_freqs$allele, perl = TRUE)
allele_freqs$allele <- gsub("HLA-A", "HLA-A*", allele_freqs$allele)
allele_freqs$allele <- gsub("HLA-B", "HLA-B*", allele_freqs$allele)
allele_freqs$allele <- gsub("HLA-C", "HLA-C*", allele_freqs$allele)


allele_freqs_filt <- subset(allele_freqs, allele_freqs$allele %in% hla_cleavedflurry_filt$allele)

new_allele_df <- as.data.frame(table(covid_recode_flurry$allele))
colnames(new_allele_df)[1] <- "allele"
merged_allele_df <- merge(new_allele_df, allele_freqs_filt, by="allele", all=TRUE)
merged_allele_df <- merged_allele_df[merged_allele_df$Freq>0,]

a <- ggplot(covid_recode_flurry, aes(x = allele, y = occ, fill = class, alpha = as.integer(less_than_50)))+ 
  geom_bar(position="stack", stat="identity")+
  theme_classic()+
  ylab("")+
  xlab("HLA alleles")+
  scale_fill_manual(values=c('#1B9E77', '#D95F02', '#7570B3'),name = "Gene", labels = c("HLA-A", "HLA-B", "HLA-C"))+
  coord_flip()+
  ggtitle("SARS-CoV-2 Presentation-flurry")+
  scale_y_discrete(expand=c(0,0))+
  scale_alpha(range = c(0.5, 1), guide = FALSE)+
  theme(axis.ticks = element_line(size = .2), text = element_text(size=20), axis.text.y=element_blank(), legend.key.size=unit(1,"cm"), panel.border=element_rect(colour="black", fill=NA, size=1))

b <- ggplot(merged_allele_df, aes(x = reorder(allele,Freq), y = as.numeric(as.character(frequency)), fill = class))+
  geom_bar(position="stack", stat="identity")+
  theme_classic() + scale_x_discrete(position='top')+
  theme(axis.title.y=element_blank(), axis.text.y=element_blank())+
  ylab("Global allele frequency")+
  scale_fill_manual(values=c('#1B9E77', '#D95F02', '#7570B3'))+
  coord_flip()+
  #ggtitle("Global allele frequency")+
  scale_y_reverse(expand=(c(0,0)))+
  scale_alpha(range = c(0.5, 1), guide = FALSE)+
  theme(legend.position = "none")+
  theme(axis.ticks = element_line(size = .2), text = element_text(size=20),panel.border=element_rect(colour="black", fill=NA, size=1))

#grid.arrange(b,a, ncol=2)
gA <- ggplotGrob(a)
gB <- ggplotGrob(b)
pdf(file="HLA_distrib_flurry.pdf",height=7,width=9)
grid::grid.newpage()
grid::grid.draw(cbind(gB, gA))
dev.off()






sars_recode_flurry <- sars_hla_cleavedflurry_filt %>%
  arrange(desc(less_than_50))
sars_recode_flurry$allele <- reorder(sars_recode_flurry$allele,sars_recode_flurry$allele,FUN=length)

new_allele_df <- as.data.frame(table(sars_recode_flurry$allele))
colnames(new_allele_df)[1] <- "allele"
merged_allele_df <- merge(new_allele_df, allele_freqs_filt, by="allele", all=TRUE)
merged_allele_df <- merged_allele_df[merged_allele_df$Freq>0,]


a <- ggplot(sars_recode_flurry, aes(x = allele, y = occ, fill = class, alpha = as.integer(less_than_50)))+ 
  geom_bar(position="stack", stat="identity")+
  theme_classic()+
  ylab("")+
  xlab("HLA alleles")+
  scale_fill_manual(values=c('#1B9E77', '#D95F02', '#7570B3'),name = "Gene", labels = c("HLA-A", "HLA-B", "HLA-C"))+
  coord_flip()+
  ggtitle("SARS-CoV Presentation-flurry")+
  scale_y_discrete(expand=c(0,0))+
  scale_alpha(range = c(0.5, 1), guide = FALSE)+
  theme(axis.ticks = element_line(size = .2), text = element_text(size=20), axis.text.y=element_blank(), legend.key.size=unit(1,"cm"), panel.border=element_rect(colour="black", fill=NA, size=1))

b <- ggplot(merged_allele_df, aes(x = reorder(allele,Freq), y = as.numeric(as.character(frequency)), fill = class))+
  geom_bar(position="stack", stat="identity")+
  theme_classic() + scale_x_discrete(position='top')+
  theme(axis.title.y=element_blank(), axis.text.y=element_blank())+
  ylab("Global allele frequency")+
  scale_fill_manual(values=c('#1B9E77', '#D95F02', '#7570B3'))+
  coord_flip()+
  #ggtitle("Global allele frequency")+
  scale_y_reverse(expand=(c(0,0)))+
  scale_alpha(range = c(0.5, 1), guide = FALSE)+
  theme(legend.position = "none")+
  theme(axis.ticks = element_line(size = .2), text = element_text(size=20),panel.border=element_rect(colour="black", fill=NA, size=1))

#grid.arrange(b,a, ncol=2)
gA <- ggplotGrob(a)
gB <- ggplotGrob(b)
pdf(file="HLA_distrib_SARS_flurry.pdf",height=7,width=9)
grid::grid.newpage()
grid::grid.draw(cbind(gB, gA))
dev.off()



#conserved
conserved_peps <- read.csv("conserved-peptides.csv")
covid_recode_flurry_conserved <- covid_recode_flurry[covid_recode_flurry$Peptide %in% conserved_peps$peptide,]
covid_recode_flurry_conserved$allele <- reorder(covid_recode_flurry_conserved$allele,covid_recode_flurry_conserved$allele,FUN=length)
write.csv(covid_recode_flurry_conserved, "covid_recode_flurry_conserved.csv", row.names = FALSE)


a <- ggplot(covid_recode_flurry_conserved, aes(x = allele, y = occ, fill = class, alpha = as.integer(less_than_50)))+ 
  geom_bar(position="stack", stat="identity")+
  theme_classic()+
  theme(axis.title.y=element_blank()) + 
  scale_fill_manual(values=c('#1B9E77', '#D95F02', '#7570B3'))+
  coord_flip()+
  ggtitle("covid <500nm- flurry")+
  scale_y_discrete(expand=c(0,0))+
  scale_alpha(range = c(0.6, 1), guide = FALSE)+
  theme(axis.ticks = element_line(size = .2), text = element_text(size=20))
new_allele_df_con <- as.data.frame(table(covid_recode_flurry_conserved$allele))
colnames(new_allele_df_con)[1] <- "allele"
new_allele_df_con <- merge(new_allele_df_con, allele_freqs_filt, by="allele", all=TRUE)
new_allele_df_con <- new_allele_df_con[new_allele_df_con$Freq>0,]

b <- ggplot(new_allele_df_con, aes(x = reorder(allele,Freq), y = as.numeric(as.character(frequency)), fill = class))+
  geom_bar(position="stack", stat="identity")+
  theme_classic() + scale_x_discrete(position='top')+
  #theme(axis.title.y=element_blank(), axis.text.y=element_blank())+
  theme(axis.title.y=element_blank())+
  scale_fill_manual(values=c('#1B9E77', '#D95F02', '#7570B3'))+
  coord_flip()+
  ggtitle("Allele frequencies") + scale_y_reverse()+
  theme(legend.position = "none")+
  theme(axis.ticks = element_line(size = .2))

#grid.arrange(b,a, ncol=2)
gA <- ggplotGrob(a)
gB <- ggplotGrob(b)

grid::grid.newpage()
grid::grid.draw(cbind(gB, gA))


#####Plot new conserved
a <- ggplot(covid_recode_flurry_conserved, aes(x = allele, y = occ, fill = class, alpha = as.integer(less_than_50)))+ 
  geom_bar(position="stack", stat="identity")+
  theme_classic()+
  ylab("")+
  xlab("HLA alleles")+
  scale_fill_manual(values=c('#1B9E77', '#D95F02', '#7570B3'),name = "Gene", labels = c("HLA-A", "HLA-B", "HLA-C"))+
  coord_flip()+
  ggtitle("SARS-CoV-2 Conserved")+
  scale_y_discrete(expand=c(0,0))+
  scale_alpha(range = c(0.5, 1), guide = FALSE)+
  theme(axis.ticks = element_line(size = .2), text = element_text(size=20), axis.text.y=element_blank(), legend.key.size=unit(1,"cm"), panel.border=element_rect(colour="black", fill=NA, size=1))

b <- ggplot(new_allele_df_con, aes(x = reorder(allele,Freq), y = as.numeric(as.character(frequency)), fill = class))+
  geom_bar(position="stack", stat="identity")+
  theme_classic() + scale_x_discrete(position='top')+
  theme(axis.title.y=element_blank(), axis.text.y=element_blank())+
  ylab("Global allele frequency")+
  scale_fill_manual(values=c('#1B9E77', '#D95F02', '#7570B3'))+
  coord_flip()+
  #ggtitle("Global allele frequency")+
  scale_y_reverse(expand=(c(0,0)))+
  scale_alpha(range = c(0.5, 1), guide = FALSE)+
  theme(legend.position = "none")+
  theme(axis.ticks = element_line(size = .2), text = element_text(size=20),panel.border=element_rect(colour="black", fill=NA, size=1))

#grid.arrange(b,a, ncol=2)
gA <- ggplotGrob(a)
gB <- ggplotGrob(b)
pdf(file="Conserved_HLA_distrib.pdf",height=7,width=9)
grid::grid.newpage()
grid::grid.draw(cbind(gB, gA))
dev.off()




###
#write.csv(table(covid_recode_flurry$allele), "covid_full_order.csv", row.names=FALSE)
#write.csv(table(covid_recode_flurry_conserved$allele), "covid_conserved_order.csv", row.names=FALSE)