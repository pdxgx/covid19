#sars_hla_cleavednet_filt <- read.csv("nickfilekekw.csv")

setwd("~/Documents/Coding/covid/covid19/csv")

covid_recode <- hla_cleavednet_filt %>%
  arrange(desc(less_than_50))
covid_recode$allele <- reorder(covid_recode$allele,covid_recode$allele,FUN=length)

new_allele_df <- as.data.frame(table(covid$allele))
colnames(new_allele_df)[1] <- "allele"
merged_allele_df <- merge(new_allele_df, allele_freqs_filt, by="allele", all=TRUE)


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
  ggtitle("Allele frequencies") + scale_y_reverse()+
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
  ggtitle("Allele frequencies") + scale_y_reverse()+
  theme(legend.position = "none")+
  theme(axis.ticks = element_line(size = .2))

#grid.arrange(b,a, ncol=2)
gA <- ggplotGrob(a)
gB <- ggplotGrob(b)

grid::grid.newpage()
grid::grid.draw(cbind(gB, gA))



#####Conserved
conserved_peps <- read.csv("conserved-peptides.csv")
covid_recode_conserved <- covid_recode[covid_recode$Peptide %in% conserved_peps$peptide,]
covid_recode_conserved$allele <- reorder(covid_recode_conserved$allele,covid_recode_conserved$allele,FUN=length)

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
