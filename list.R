covid_list <- read.csv("combined.csv", stringsAsFactors = FALSE)
#order(covid_list$Nuggets)

nugs <- as.data.frame(covid_list$Nuggets)
nugs$order_nugs <- rownames(nugs)
colnames(nugs)[1] <- "allele"
net <- as.data.frame(covid_list$Net)
net$order_flurry <- rownames(net)
colnames(net)[1] <- "allele"
flurry <- as.data.frame(covid_list$Flurry)
flurry$order_flurry <- rownames(flurry)
colnames(flurry)[1] <- "allele"

combined_covid <- merge(merge(nugs, net, by="allele"), flurry, by="allele")
combined_covid$average[1] <- mean(as.numeric(combined_covid[1,2:4]))

for (i in 1:dim(combined_covid)[1])
{
  combined_covid$average[i] <- mean(as.numeric(combined_covid[i,2:4]))
}

write.csv(combined_covid, "rank_alleles.csv", row.names=FALSE)

combined_both <- merge(combined_covid, combined_sars, by="allele")
write.csv(combined_both, "bothrank_alleles.csv", row.names=FALSE)
