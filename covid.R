setwd("~/Documents/Coding/covid/covid19")
hla_binding <- read.csv("covid_hla_filt.csv")
gsub('.{4}$', '', hla_binding$file_name) -> hla_binding$file_name
sub('..........', '', hla_binding$file_name) -> hla_binding$file_name
hla_binding_500 <- hla_binding[as.integer(hla_binding$ic50) < 500,]
hla_binding_500 <- hla_binding_500[!is.na(as.integer(hla_binding_500$ic50)),]
colnames(hla_binding_500)[1] <- "allele"
write.csv(hla_binding_500, file="covid_500.csv", row.names = FALSE)
library(ggplot2)

count <- hla_binding_500$allele
t


