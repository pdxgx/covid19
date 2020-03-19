list <- read.csv("sars_rank.csv", stringsAsFactors = FALSE)


nugs <- as.data.frame(list$nugs)
nugs$order_nugs <- rownames(nugs)
colnames(nugs)[1] <- "allele"
net <- as.data.frame(list$net)
net$order_flurry <- rownames(net)
colnames(net)[1] <- "allele"
flurry <- as.data.frame(list$flurry)
flurry$order_flurry <- rownames(flurry)
colnames(flurry)[1] <- "allele"

combined_sars <- merge(merge(nugs, net, by="allele"), flurry, by="allele")
combined_sars$average[1] <- mean(as.numeric(combined_sars[1,2:4]))

for (i in 1:dim(combined_sars)[1])
{
  combined_sars$average[i] <- mean(as.numeric(combined_sars[i,2:4]))
}

write.csv(combined_sars, "sars_rank_alleles.csv", row.names=FALSE)
