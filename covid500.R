hla <- read.csv("covid_hla_500.csv")
typeof(hla$file_name)
gsub('.{4}$', '', hla$file_name) -> hla$file_name
sub('..........', '', hla$file_name) -> hla$file_name
colnames(hla)[1] <- 'allele'
hla_filt <- hla[is.na(as.numeric(as.character(hla$allele))),]
hla_filt$less_than_50 <- ifelse(hla_filt$ic50<50,'50','500')
hla_filt$occ <- rep(1,nrow(hla_filt))
#hla_50 <- hla_filt[as.integer(hla_filt$ic50) < 50,]
#hla_450 <- hla_filt[as.integer(hla_filt$ic50) > 50,]

typeof(as.character(hla_filt$allele))
library(ggplot2)

ggplot(hla_filt, aes(x=allele)) + geom_bar()
ggplot(hla_50, aes(x=allele)) + geom_bar()
table(hla_filt$allele)
#hla_50_tab <- data.frame(table(hla_50$allele))
#hla_450_tab <- data.frame(table(hla_450$allele))
#merged <- merge(hla_50_tab, hla_450_tab, by='Var1')

# create a dataset
specie <- c(rep("sorgho" , 3) , rep("poacee" , 3) , rep("banana" , 3) , rep("triticum" , 3) )
condition <- rep(c("normal" , "stress" , "Nitrogen") , 4)
value <- abs(rnorm(12 , 0 , 15))
data <- data.frame(specie,condition,value)

# Grouped
ggplot(data, aes(fill=condition, y=value, x=specie)) + 
  geom_bar(position="stack", stat="identity")

ggplot(hla_filt, aes(fill=less_than_50, y=occ, x=reorder(allele, -occ))) + 
  geom_bar(position="stack", stat="identity")

library(forcats)
hla_filt$allele <- reorder(hla_filt$allele,hla_filt$allele,FUN=length)

