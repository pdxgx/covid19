hla_flurry <- read.csv("flurry_preds.csv")
hla_flurry_filt <- hla_flurry[as.integer(hla_flurry$mhcflurry_prediction) < 500,]
hla_flurry_filt$less_than_50 <- ifelse(hla_flurry_filt$mhcflurry_prediction<50,'50','500')
hla_flurry_filt$occ <- rep(1,nrow(hla_flurry_filt))
#hla_50 <- hla_flurry_filt[as.integer(hla_flurry_filt$ic50) < 50,]
#hla_450 <- hla_flurry_filt[as.integer(hla_flurry_filt$ic50) > 50,]

typeof(as.character(hla_flurry_filt$allele))
library(ggplot2)

library(forcats)
hla_flurry_filt$allele <- reorder(hla_flurry_filt$allele,hla_flurry_filt$allele,FUN=length)


ggplot(hla_flurry_filt, aes(fill=less_than_50, y=occ, x=reorder(allele, -occ))) + 
  geom_bar(position="stack", stat="identity")

write.csv(rownames(table(hla_flurry_filt$allele)), "mhc_flurry_order.csv")
write.csv(hla_flurry_filt, file="peps_flurry_filt.csv", row.names=FALSE)

library(genefilter)
