individual <- read.csv("individual_data.csv")
haplotypes <- read.csv("globalhaps.csv")

b46 <- haplotypes[haplotypes$B == "B*46:01",]
b46 <- b46[complete.cases(b46),]


#A*02:07 B*46:01    C*01:02

iedb <- read.csv("iedb.csv")
iedbprots <- iedb$Description

iedb_netmhc_pan <- hla_cleavednet_filt[hla_cleavednet_filt$Peptide %in% iedbprots,]
covid_recode <- iedb_netmhc_pan %>%
  arrange(desc(less_than_50))
covid_recode$allele <- reorder(covid_recode$allele,covid_recode$allele,FUN=length)

allele_freqs <- as.data.frame(read.csv("allele_freqs.csv", header=TRUE))
colnames(allele_freqs) <- c("allele", "frequency")
allele_freqs$class <- ifelse(grepl("A", allele_freqs$allele), 'A', ifelse(grepl("B", allele_freqs$allele), 'B', 'C'))
allele_freqs$allele <- paste("HLA-", allele_freqs$allele, sep="")
allele_freqs$allele <- gsub("(\\d{2})(?=\\d{2})", "\\1:", allele_freqs$allele, perl = TRUE)
allele_freqs_filt <- subset(allele_freqs, allele_freqs$allele %in% hla_cleavednet_filt$allele)

new_allele_df <- as.data.frame(table(covid_recode$allele))
colnames(new_allele_df)[1] <- "allele"
merged_allele_df <- merge(new_allele_df, allele_freqs_filt, by="allele", all=TRUE)
merged_allele_df <- merged_allele_df[merged_allele_df$Freq>0,]
library(ggplot2)
a <- ggplot(covid_recode, aes(x = allele, y = occ, fill = class, alpha = as.integer(less_than_50)))+ 
  geom_bar(position="stack", stat="identity")+
  theme_classic()+
  ylab("")+
  xlab("HLA alleles")+
  scale_fill_manual(values=c('#1B9E77', '#D95F02', '#7570B3'),name = "Gene", labels = c("HLA-A", "HLA-B", "HLA-C"))+
  coord_flip()+
  ggtitle("iedb SARS-CoV-2 Presentation")+
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
pdf(file="iedb_HLA_distrib_netmhcpan.pdf",height=7,width=9)
grid::grid.newpage()
grid::grid.draw(cbind(gB, gA))
dev.off()

alleles <- table(covid_recode$allele)

