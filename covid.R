setwd("~/Documents/Coding/covid/covid19")
hla_binding <- read.csv("covid_hla.csv", header=FALSE)
hla_binding_500 <- hla_binding[as.integer(hla_binding$V3) < 500,]
