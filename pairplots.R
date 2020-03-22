panel.cor <- function(x, y, digits=2, prefix="") 
{
  usr <- par("usr"); on.exit(par(usr)) 
  par(usr = c(0, 1, 0, 1)) 
  r <- abs(cor(x, y)) 
  txt <- format(c(r, 0.123456789), digits=digits)[1] 
  txt <- paste(prefix, txt, sep="") 
  
  test <- cor.test(x,y) 
  # borrowed from printCoefmat
  Signif <- symnum(test$p.value, corr = FALSE, na = FALSE, 
                   cutpoints = c(0, 0.001, 0.01, 0.05, 0.1, 1),
                   symbols = c("***", "**", "*", ".", " ")) 
  
  text(0.5, 0.5, txt, cex=1.25) 
  text(.8, .8, Signif, cex=1.25, col=2) 
}

hla_a <- read.csv("covid.HLA-A.full.csv")
hla_b <- read.csv("covid.HLA-B.full.csv")
hla_c <- read.csv("covid.HLA-C.full.csv")
tempdf <- data.frame(matrix(NA, nrow = 9755, ncol=0))
tempdf$hla_a <- hla_a
tempdf$hla_b <- hla_b
tempdf$hla_c <- hla_c
tempdf$sum <- tempdf$hla_b+  tempdf$hla_a+ tempdf$hla_c
# Corplot for mutational burden across callers
pairs(cbind(tempdf$hla_a, tempdf$hla_b, tempdf$hla_c), 
      lower.panel=panel.smooth, upper.panel=panel.cor,
      labels = c('HLA-A', 'HLA-B', 'HLA-C'),
      cex.labels = 1.25)
