sars_hla_cleavednet_filt <- read.csv("nickfilekekw.csv")

sars_recode <- sars_hla_cleavednet_filt %>%
  arrange(desc(less_than_50))
sars_recode$allele <- reorder(sars_recode$allele,sars_recode$allele,FUN=length)

merged_allele_df <- read.csv("nick2.csv")
a <- ggplot(sars_recode, aes(x = allele, y = occ, fill = class, alpha = as.integer(less_than_50)))+ 
  geom_bar(position="stack", stat="identity")+
  theme_classic()+
  theme(axis.title.y=element_blank()) + 
  scale_fill_manual(values=c('#1B9E77', '#D95F02', '#7570B3'))+
  coord_flip()+
  ggtitle("SARS <500nm- netmhcpan")+
  scale_y_discrete(expand=c(0,0))+
  scale_alpha(range = c(0.6, 1), guide = FALSE)+
  theme(axis.ticks = element_line(size = .2))

b <- ggplot(merged_allele_df, aes(x = allele, y = as.numeric(as.character(frequency)), fill = class))+
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
