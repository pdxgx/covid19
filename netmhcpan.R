require(dplyr)

#SARS-CoV-2
hla_cleavednet <- read.csv("covid_netmhcpan_scores_all_alleles.csv")
hla_cleavednet_filt <- hla_cleavednet[as.integer(hla_cleavednet$Binding_affinity) < 500,]
hla_cleavednet_filt$less_than_50 <- ifelse(hla_cleavednet_filt$Binding_affinity<50,'50','500')
hla_cleavednet_filt$occ <- rep(1,nrow(hla_cleavednet_filt))
colnames(hla_cleavednet_filt)[2] <- "allele"
hla_cleavednet_filt$allele <- reorder(hla_cleavednet_filt$allele,hla_cleavednet_filt$allele,FUN=length)

table(hla_cleavednet_filt$allele)
hla_cleavednet_filt$class <- rep(1,nrow(hla_cleavednet_filt))
hla_cleavednet_filt$class <- ifelse(grepl("LA-A", hla_cleavednet_filt$allele), 'A', ifelse(grepl("LA-B", hla_cleavednet_filt$allele), 'B', 'C'))


#netsars
sars_hla_cleavednet <- read.csv("sars_netmhcpan_scores_all_alleles.csv")
sars_hla_cleavednet_filt <- sars_hla_cleavednet[as.integer(sars_hla_cleavednet$Binding_affinity) < 500,]
sars_hla_cleavednet_filt$less_than_50 <- ifelse(sars_hla_cleavednet_filt$Binding_affinity<50,'50','500')
sars_hla_cleavednet_filt$occ <- rep(1,nrow(sars_hla_cleavednet_filt))
colnames(sars_hla_cleavednet_filt)[2] <- "allele"
sars_hla_cleavednet_filt$allele <- reorder(sars_hla_cleavednet_filt$allele,sars_hla_cleavednet_filt$allele,FUN=length)
sars_hla_cleavednet_filt$occ_50 <- ifelse(sars_hla_cleavednet_filt$less_than_50 == '50', 1, 0)
table(sars_hla_cleavednet_filt$allele)
sars_hla_cleavednet_filt$class <- rep(1,nrow(sars_hla_cleavednet_filt))
sars_hla_cleavednet_filt$class <- ifelse(grepl("LA-A", sars_hla_cleavednet_filt$allele), 'A', ifelse(grepl("LA-B", sars_hla_cleavednet_filt$allele), 'B', 'C'))


#plotting
covid_recode <- hla_cleavednet_filt %>%
  arrange(desc(less_than_50))
covid_recode$allele <- reorder(covid_recode$allele,covid_recode$allele,FUN=length)

allele_freqs <- as.data.frame(read.csv("allele_freqs.csv", header=FALSE))
colnames(allele_freqs) <- c("allele", "frequency")
allele_freqs$class <- ifelse(grepl("A", allele_freqs$allele), 'A', ifelse(grepl("B", allele_freqs$allele), 'B', 'C'))
allele_freqs$allele <- paste("HLA-", allele_freqs$allele, sep="")
allele_freqs$allele <- gsub("(\\d{2})(?=\\d{2})", "\\1:", allele_freqs$allele, perl = TRUE)
allele_freqs_filt <- subset(allele_freqs, allele_freqs$allele %in% hla_cleavednet_filt$allele)

new_allele_df <- as.data.frame(table(covid_recode$allele))
colnames(new_allele_df)[1] <- "allele"
merged_allele_df <- merge(new_allele_df, allele_freqs_filt, by="allele", all=TRUE)
merged_allele_df <- merged_allele_df[merged_allele_df$Freq>0,]

a <- ggplot(covid_recode, aes(x = allele, y = occ, fill = class, alpha = as.integer(less_than_50)))+ 
  geom_bar(position="stack", stat="identity")+
  theme_classic()+
  theme(axis.title.y=element_blank()) + 
  scale_fill_manual(values=c('#1B9E77', '#D95F02', '#7570B3'))+
  coord_flip()+
  ggtitle("covid <500nm- netmhcpan")+
  scale_y_discrete(expand=c(0,0))+
  scale_alpha(range = c(0.6, 1), guide = FALSE)+
  theme(axis.ticks = element_line(size = .2), text = element_text(size=20))

b <- ggplot(merged_allele_df, aes(x = reorder(allele,Freq), y = as.numeric(as.character(frequency)), fill = class))+
  geom_bar(position="stack", stat="identity")+
  theme_classic() + scale_x_discrete(position='top')+
  theme(axis.title.y=element_blank(), axis.text.y=element_blank())+
  #theme(axis.title.y=element_blank())+
  scale_fill_manual(values=c('#1B9E77', '#D95F02', '#7570B3'))+
  coord_flip()+
  ggtitle("Allele frequencies") + scale_y_reverse(expand=(c(0,0)))+
  theme(legend.position = "none")+
  theme(axis.ticks = element_line(size = .2))

#grid.arrange(b,a, ncol=2)
gA <- ggplotGrob(a)
gB <- ggplotGrob(b)

grid::grid.newpage()
grid::grid.draw(cbind(gB, gA))







sars_recode <- sars_hla_cleavednet_filt %>%
  arrange(desc(less_than_50))
sars_recode$allele <- reorder(sars_recode$allele,sars_recode$allele,FUN=length)

new_allele_df <- as.data.frame(table(sars_recode$allele))
colnames(new_allele_df)[1] <- "allele"
merged_allele_df <- merge(new_allele_df, allele_freqs_filt, by="allele", all=TRUE)
merged_allele_df <- merged_allele_df[merged_allele_df$Freq>0,]


#merged_allele_df <- read.csv("nick2.csv")
a <- ggplot(sars_recode, aes(x = allele, y = occ, fill = class, alpha = as.integer(less_than_50)))+ 
  geom_bar(position="stack", stat="identity")+
  theme_classic()+
  theme(axis.title.y=element_blank()) + 
  scale_fill_manual(values=c('#1B9E77', '#D95F02', '#7570B3'))+
  coord_flip()+
  ggtitle("SARS <500nm- netmhcpan")+
  scale_y_discrete(expand=c(0,0))+
  scale_alpha(range = c(0.6, 1), guide = FALSE)+
  theme(axis.ticks = element_line(size = .2), text = element_text(size=20))

b <- ggplot(merged_allele_df, aes(x = reorder(allele,Freq), y = as.numeric(as.character(frequency)), fill = class))+
  geom_bar(position="stack", stat="identity")+
  theme_classic() + scale_x_discrete(position='top')+
  theme(axis.title.y=element_blank(), axis.text.y=element_blank())+
  #theme(axis.title.y=element_blank())+
  scale_fill_manual(values=c('#1B9E77', '#D95F02', '#7570B3'))+
  coord_flip()+
  ggtitle("Allele frequencies") + scale_y_reverse(expand=c(0,0))+
  theme(legend.position = "none")+
  theme(axis.ticks = element_line(size = .2))

#grid.arrange(b,a, ncol=2)
gA <- ggplotGrob(a)
gB <- ggplotGrob(b)

grid::grid.newpage()
grid::grid.draw(cbind(gB, gA))



#conserved
conserved_peps <- read.csv("conserved-peptides.csv")
covid_recode_conserved <- covid_recode[covid_recode$Peptide %in% conserved_peps$peptide,]
covid_recode_conserved$allele <- reorder(covid_recode_conserved$allele,covid_recode_conserved$allele,FUN=length)
write.csv(covid_recode_conserved, "covid_recode_conserved.csv", row.names = FALSE)


a <- ggplot(covid_recode_conserved, aes(x = allele, y = occ, fill = class, alpha = as.integer(less_than_50)))+ 
  geom_bar(position="stack", stat="identity")+
  theme_classic()+
  theme(axis.title.y=element_blank()) + 
  scale_fill_manual(values=c('#1B9E77', '#D95F02', '#7570B3'))+
  coord_flip()+
  ggtitle("covid <500nm- netmhcpan")+
  scale_y_discrete(expand=c(0,0))+
  scale_alpha(range = c(0.6, 1), guide = FALSE)+
  theme(axis.ticks = element_line(size = .2), text = element_text(size=20))

new_allele_df_con <- as.data.frame(table(covid_recode_conserved$allele))
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



###
write.csv(table(covid_recode$allele), "covid_full_order.csv", row.names=FALSE)
write.csv(table(covid_recode_conserved$allele), "covid_conserved_order.csv", row.names=FALSE)