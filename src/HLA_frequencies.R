#!/usr/bin/env R

# Geographic maps of HLA data 
# Helpful resources
# https://www.r-spatial.org/r/2018/10/25/ggplot2-sf.html
# https://www.r-bloggers.com/map-visualization-of-covid-19-across-the-world-with-r/
# https://cran.r-project.org/web/packages/rworldmap/index.html

require(rvest)
require(rworldmap)

alleles.outfile <- "HLA-freqs.rda"
iso3.outfile <- "HLA-ISO3-freqs.rda"
global.outfile <- "HLA-freqs-global.rda"
global.haplotype.outfile <- "haplotype-freqs-global.rda"
individuals.outfile <- "HLA-individuals.rda"
haplotypes.outfile <- "HLA-haplotypes.rda"
haplotype.iso3.outfile <- "haplotype-ISO3-freqs.rda"
MHC.alleles.infile <- "supporting_data/netmhcpan_allele_list.txt"
plotdir <- "plots"

#--------------------
# get individual allele frequency data
#--------------------
alleles <- data.frame(HLA=character(), pop.ID=numeric(), pop.name=character(), freq=numeric(), sample.size=numeric())
for (HLA in LETTERS[1:3]) {
	page <- 1
	while (page > 0) {
		print(paste0("HLA-",HLA," (page ",page,")"))
		qurl <- paste0("http://www.allelefrequencies.net/hla6006a.asp?", "hla_locus=", HLA, "&page=", page)
 		data <- qurl %>%
   			read_html()
   		recs <- data %>%
   			html_nodes(xpath = '//*[@id="divGenNavig2"]/table') %>%
   			html_table()
   		recs <- recs[[1]][1]
   		if (sub(".*to ([0-9,]+)[^0-9]*[(]from (.*)[)].*", "\\1", recs) ==
   			sub(".*to ([0-9,]+)[^0-9]*[(]from (.*)[)].*", "\\2", recs)) {
   			page <- 0
   		} else {
   			page <- page + 1
   		}
		data <- data %>%
  			html_nodes(xpath = '//*[@id="divGenDetail"]/table') 
    	popIDs <- data %>% 
    		html_nodes("a") %>% 
    		html_attr("href") %>% 
    		grep(pattern="pop6001c", value=TRUE) %>%
    		sub(pattern=".*[=]",replacement="") %>%
    		as.numeric()
    	data <- data %>% 
    		html_table() %>%
    		unlist(recursive=FALSE)
    	alleles <- rbind(alleles, data.frame(HLA=data[[2]], pop.ID=popIDs, pop.name=data[[4]], freq=sub("[^0-9,.]+","",data[[6]]), sample.size=gsub(",","",data[[8]])))
    }
	save(alleles,file=alleles.outfile)
}

#--------------------
# get population data
#--------------------
popIDs <- unique(alleles[,"pop.ID"])
N <- nrow(alleles)
alleles <- cbind(alleles, data.frame(lat=rep(NA, N), long=rep(NA, N), region=rep(NA, N), ethnicity=rep(NA, N), country=rep(NA, N)))
for (pop in popIDs) {
	print(paste0("Population ID: ",pop))
	data <- paste0("http://www.allelefrequencies.net/pop6001c.asp?pop_id=",pop) %>%
		read_html() %>%
		html_nodes(xpath='//*[@style="width: 700px; border-top: 1px solid #cccccc; border-left: 1px solid #cccccc; border-bottom: 3px solid #aaaaaa; border-right: 3px solid #aaaaaa;"]/table')
	country <- data[[1]] %>% 
		html_nodes("img") %>% 
		html_attr("src") %>%
		sub(pattern=".*[/]([A-Z]+)[.]bmp", replacement="\\1")
	data <- data[[2]] %>%
		html_table()
	lat <- data[which(data[,2]=="Latitude:"),3]
	long <- data[which(data[,2]=="Longitude:"),3]
	region <- data[which(data[,2]=="Geographic Region:"),3]
	ethnicity <- data[which(data[,2]=="Ethnic origin:"),3]
	alleles[which(alleles[,"pop.ID"]==pop), c("lat", "long", "region", "ethnicity", "country")] <- list(lat=lat, long=long, region=region, country=country)
}
# fix ISO3 code match issues...
# England -> Great Britain
alleles[which(alleles[,"country"]=="ENG"), "country"] <- "GBR"
# Scotland -> Great Britain
alleles[which(alleles[,"country"]=="SCO"), "country"] <- "GBR"
# Wales -> Great Britain
alleles[which(alleles[,"country"]=="WAL"), "country"] <- "GBR"
# Northern Ireland -> Great Britain
alleles[which(alleles[,"country"]=="NIR"), "country"] <- "GBR"
# Gaza -> Israel
alleles[which(alleles[,"country"]=="GAZ"), "country"] <- "ISR"
# Romania -> Romania
alleles[which(alleles[,"country"]=="ROM"), "country"] <- "ROU"
# Iraqi Kurdistan -> Iraq
alleles[which(alleles[,"country"]=="KUR"), "country"] <- "IRQ"
save(alleles,file=alleles.outfile)

#--------------------
# get per-individual HLA combinations
#--------------------
individuals <- data.frame(A1=character(), A2=character(), B1=character(), B2=character(), C1=character(), C2=character(), pop.ID=numeric())
for (pop in popIDs) {
	data <- (paste0("http://www.allelefrequencies.net/viewrawdata.asp?pop_id=",pop) %>%
		read_html() %>%
  		html_nodes(xpath = '//*[@style="width: 700px; border-top: 1px solid #cccccc; border-left: 1px solid #cccccc; border-bottom: 3px solid #aaaaaa; border-right: 3px solid #aaaaaa;"]/table'))[[2]]
  	if (data %>% html_text() %>% grepl(pattern="no public genotype")) {
  		print(paste0("No genotype data available for population ID: ",pop)) 
  	} else {
		print(paste0("Population ID: ",pop))
	  	data <- data %>%
	  		html_table()
	  	A <- which(apply(data,2,function(x) { any(grepl(pattern="A*",x, fixed=TRUE))}))
	  	if (length(A) < 2) {
	  		A1 <- A2 <- rep(NA, nrow(data))
	  	} else {
	  		A1 <- data[,A[1]]
	  		if (length(C) < 1) {
	  			A2 <- rep(NA, nrow(data))
	  		} else {
		  		A2 <- data[,A[2]]	  			
	  		}
	  	}
	  	B <- which(apply(data,2,function(x) { any(grepl(pattern="B*",x, fixed=TRUE))}))
	  	if (length(B) < 2) {
	  		B1 <- B2 <- rep(NA, nrow(data))
	  	} else {
	  		B1 <- data[,B[1]]
	  		if (length(B) < 1) {
	  			B2 <- rep(NA, nrow(data))
	  		} else {
		  		B2 <- data[,B[2]]	  			
	  		}
	  	}
	  	C <- which(apply(data,2,function(x) { any(grepl(pattern="C*",x, fixed=TRUE))}))
	  	if (length(C) < 2) {
	  		C1 <- C2 <- rep(NA, nrow(data))
	  	} else {
	  		C1 <- data[,C[1]]
	  		if (length(C) < 1) {
	  			C2 <- rep(NA, nrow(data))
	  		} else {
		  		C2 <- data[,C[2]]	  			
	  		}
	  	}
	  	individuals <- rbind(individuals, data.frame(A1=A1, A2=A2, B1=B1, B2=B2, C1=C1, C2=C2, pop.ID=rep(pop, nrow(data))))
  	}
}
individuals[individuals=="untyped"] <- NA
individuals <- individuals[which(unlist(apply(individuals,1,function(x) {!all(is.na(x))}))),]
save(individuals, file=individuals.outfile)

#--------------------
# get HLA haplotype frequency data
#--------------------
haplotypes <- data.frame(A=character(), B=character(), C=character(), pop.ID=numeric(), freq=numeric())
for (pop in popIDs) {
	page <- 1
	while (page > 0) {
		print(paste0("Population ID: ",pop," (page ",page,")"))
		qurl <- paste0("http://www.allelefrequencies.net/hla6003a.asp?hla_population=", pop, "&page=", page)
		data <- qurl %>%
			read_html()
   		recs <- data %>%
   			html_nodes(xpath = '//*[@id="divGenNavig"]/table') %>%
   			html_table()
   		if (length(recs) < 1) {
   			page <- 0
   			next
   		}
   		recs <- recs[[1]][1]
   		if (sub(".*to ([0-9,]+)[^0-9]*[(]from (.*)[)].*", "\\1", recs) ==
   			sub(".*to ([0-9,]+)[^0-9]*[(]from (.*)[)].*", "\\2", recs)) {
   			page <- 0
   		} else {
   			page <- page + 1
   		}
		data <- data %>%
  			html_nodes(xpath = '//*[@id="divGenDetail"]/table') %>% 
    		html_table() %>%
    		unlist(recursive=FALSE)
    	A <- sub(".*(A[*][0-9:]+).*", "\\1", data$Haplotype)
    	A[!grepl("A*", fixed=TRUE, A)] <- NA
    	B <- sub(".*(B[*][0-9:]+).*", "\\1", data$Haplotype)
    	B[!grepl("B*", fixed=TRUE, B)] <- NA
    	C <- sub(".*(C[*][0-9:]+).*", "\\1", data$Haplotype)
    	C[!grepl("C*", fixed=TRUE, C)] <- NA
    	haplotypes <- rbind(haplotypes, data.frame(A=A, B=B, C=C, pop.ID=rep(pop, length(A)), freq=as.numeric(sub("[^0-9,.]+","",data[[5]]))))
	}
}
save(haplotypes,file=haplotypes.outfile)

#-----------------------
# aggregate individual allele data by country (ISO3)
#-----------------------
# FUTURE DEVELOPMENT: first summarize allele frequency by ethnic population, THEN collapse to country level based on ethnic frequency in population
# --> will achieve more accurate population-level estimates because some small ethnic groups are otherwise over-represented in individual countries, e.g. Aboriginal data is >>20% of HLA info for Australia, but they are only ~3% of country's population
# UN-level data with relative size of ethnic populations in each country: http://data.un.org/Data.aspx?d=POP&f=tableCode:26
countries <- get(data(countryExData))
hla.iso3 <- data.frame(ISO3=character(), HLA=character(), freq=numeric(), popsize=numeric())
for (ISO3 in unique(alleles[,"country"])) {
	data <- alleles[which(alleles[,"country"] == ISO3),]
	popsize <- switch(ISO3,
		# 2005 population statistics from World Bank, http://api.worldbank.org/v2/en/indicator/SP.POP.TOTL?downloadformat=csv & https://www.worldometers.info/world-population/
		MTQ = 397.19,
		SRB = 7441,
		HKG = 6813.2,
		STP = 155.63,
		SGP = 4265.762,
		CPV = 474.567,
		LBY = 5798.614,
		ASM = 59.118,
		NCL = 232.25,
		GNQ = 757.317,
		TWN = 22705.713,
		# Canada HLA data only available for BC First Nations (Athabaskan and Penutian)
		# Canada First Nations are ~4% of national population 
		CAN = 32243.753 * 0.04,
		countries[which(countries[,1]==ISO3),"Population2005"]
	)
	for (HLA in unique(as.character(data[,"HLA"]))) {
		# aggregate codes (A*01:01 will include contributions from e.g. A*01:01:01 and A*01:01:02 in addition to A*01:01 itself)
		hla.data <- data[grepl(HLA,data[,"HLA"], fixed=TRUE),c("freq","sample.size")]
		# alternative approach is to only lump by exact allele match
		# hla.data <- data[which(data[,"HLA"]==HLA),c("freq","sample.size")]
		hla.iso3 <- rbind(hla.iso3, data.frame(ISO3=ISO3, HLA=HLA, freq=weighted.mean(as.numeric(hla.data[,1]), w=as.numeric(hla.data[,2])), popsize=popsize))
	}
}
save(hla.iso3,file=iso3.outfile)

#-----------------------
# aggregate haplotype data by country (ISO3)
#-----------------------
# FUTURE DEVELOPMENT: first summarize haplotype frequency by ethnic population, THEN collapse to country level based on ethnic frequency in population
# --> will achieve more accurate population-level estimates because some small ethnic groups are otherwise over-represented in individual countries, e.g. Aboriginal data is >>20% of HLA info for Australia, but they are only ~3% of country's population
# UN-level data with relative size of ethnic populations in each country: http://data.un.org/Data.aspx?d=POP&f=tableCode:26
haplotype.iso3 <- data.frame(ISO3=character(), A=character(), B=character(), C=character(), freq=numeric(), popsize=numeric())
for (ISO3 in unique(alleles[,"country"])) {
	data <- unique(alleles[which(alleles[,"country"] == ISO3),c("pop.ID", "sample.size")])
	popsize <- switch(ISO3,
		# 2005 population statistics from World Bank, http://api.worldbank.org/v2/en/indicator/SP.POP.TOTL?downloadformat=csv & https://www.worldometers.info/world-population/
		MTQ = 397.19,
		SRB = 7441,
		HKG = 6813.2,
		STP = 155.63,
		SGP = 4265.762,
		CPV = 474.567,
		LBY = 5798.614,
		ASM = 59.118,
		NCL = 232.25,
		GNQ = 757.317,
		TWN = 22705.713,
		# Canada HLA data only available for BC First Nations (Athabaskan and Penutian)
		# Canada First Nations are ~4% of national population 
		CAN = 32243.753 * 0.04,
		countries[which(countries[,1]==ISO3),"Population2005"]
	)
	haps <- haplotypes[which(haplotypes[,"pop.ID"] %in% data[,"pop.ID"]),]
	if (nrow(haps) < 1) { next }
	apply(unique(haps[,1:3]),1,function(hap) {
			if (!is.na(hap[1])) {
				which.hap <- which(haps[,"A"]==hap[1])
			}
			else {
				which.hap <- which(is.na(haps[,"A"]))
			}
			if (!is.na(hap[2])) {
				which.hap <- intersect(which.hap, which(haps[,"B"]==hap[2]))
			}
			else {
				which.hap <- intersect(which.hap, which(is.na(haps[,"B"])))
			}
			if (!is.na(hap[3])) {
				which.hap <- intersect(which.hap, which(haps[,"C"]==hap[3]))
			}
			else {
				which.hap <- intersect(which.hap, which(is.na(haps[,"C"])))
			}
			sample.sizes <- as.numeric(data[match(haps[which.hap,"pop.ID"], data[,"pop.ID"]), "sample.size"])
			haplotype.iso3 <<- rbind(haplotype.iso3, data.frame(ISO3=ISO3, A=hap[1], B=hap[2], C=hap[3], freq=weighted.mean(as.numeric(haps[which.hap,"freq"]), w=sample.sizes), popsize=popsize))
		})
}
save(haplotype.iso3,file=haplotype.iso3.outfile)

#-----------------------
# estimate global allele and haplotype frequency
#-----------------------
global_allele_freqs <- data.frame(HLA=character(), freq=numeric())
for (hla in as.character(unique(hla.iso3[,"HLA"]))) {
	hla.data <- hla.iso3[grepl(hla,hla.iso3[,"HLA"], fixed=TRUE),c("freq","popsize")]
	global_allele_freqs <- rbind(global_allele_freqs, data.frame(HLA=hla, freq=weighted.mean(as.numeric(hla.data[,1]), w=as.numeric(hla.data[,2]))))
}
save(global_allele_freqs, file=global.outfile)

global_haplotype_freqs <- data.frame(A=character(), B=character(), C=character(), freq=numeric())
apply(unique(haplotype.iso3[,2:4]),1,function(hap) {
	if (!is.na(hap[1])) {
		which.hap <- which(haplotype.iso3[,"A"]==hap[1])
	}
	else {
		which.hap <- which(is.na(haplotype.iso3[,"A"]))
	}
	if (!is.na(hap[2])) {
		which.hap <- intersect(which.hap, which(haplotype.iso3[,"B"]==hap[2]))
	}
	else {
		which.hap <- intersect(which.hap, which(is.na(haplotype.iso3[,"B"])))
	}
	if (!is.na(hap[3])) {
		which.hap <- intersect(which.hap, which(haplotype.iso3[,"C"]==hap[3]))
	}
	else {
		which.hap <- intersect(which.hap, which(is.na(haplotype.iso3[,"C"])))
	}
	global_haplotype_freqs <<- rbind(global_haplotype_freqs, data.frame(A=hap[1], B=hap[2], C=hap[3], freq=weighted.mean(as.numeric(haplotype.iso3[which.hap,"freq"]), w=as.numeric(haplotype.iso3[which.hap,"popsize"]))))
})
save(global_haplotype_freqs, file=global.haplotype.outfile)

#----------
# plot/map data function
#----------
HLAgeoPlot <- function(HLA, outdir=NULL,unified.freq=FALSE) {
	data <- hla.iso3[which(hla.iso3[,"HLA"] %in% HLA),]
	maxfreq <- ceiling(max(data[,"freq"])*100)
	for (hla in HLA) {
		sPDF <- joinCountryData2Map(data[which(data[,"HLA"]==hla),], joinCode="ISO3", nameJoinColumn="ISO3", verbose=TRUE)
		if (!is.null(outdir)) {
			dir.create(outdir, showWarnings = FALSE)
		  	pdf(paste0(outdir, "/", 
	             gsub("\\:", "-", gsub("\\*", "", hla)), ".pdf"), 
	    	   width=8, height=6)
		}
		if (unified.freq) {
			freq <- maxfreq
		}
		else {
			freq <- ceiling(max(data[which(data[,"HLA"]==hla),"freq"])*100)
		}
		mapParams <- mapCountryData(sPDF, nameColumnToPlot="freq", missingCountryCol="lightgray", 
	  			borderCol="black", oceanCol="lightblue", catMethod=0:freq/100, addLegend=FALSE,
	  			mapTitle=paste0(hla,"\n(~", format(global_allele_freqs[which(global_allele_freqs[,1]==hla),2]*100, digits=2),"% globally)"))
		do.call(addMapLegend, c(mapParams, legendWidth=0.5, legendMar=6.5))
		if (!is.null(outdir)) {
			dev.off()
		}
	}
}
#----------
# generate figures for paper
#----------
MHC_alleles <- as.character(unlist(read.csv(MHC.alleles.infile, header=TRUE))) %>%
	sub(pattern="HLA[-]([A-C])",replacement="\\1*")

# Figure 5
pdf(file=paste0(plotdir, "/best+worst_HLA_all_peptides.pdf"), width=9,height=8.5)
layout(matrix(1:6,nr=3,nc=2))
par(mar=c(3.1,0.1,2.1,0.1))
HLAgeoPlot(c("A*02:06", "A*02:07", "B*15:03", "B*46:01", "C*12:03", "C*01:02"))
dev.off()

# Supplementary Figure S1
pdf(file=paste0(plotdir, "/best+worst_HLA_conserved_peptides.pdf"), width=9,height=8.5)
layout(matrix(1:6,nr=3,nc=2))
par(mar=c(3.1,0.1,2.1,0.1))
HLAgeoPlot(c("A*02:06", "A*02:07", "B*08:01", "B*46:01", "C*12:03", "C*01:02"))
dev.off()

# Appendix 2
for (HLA in MHC_alleles) {
	HLAgeoPlot(HLA, outdir=plotdir)
}
