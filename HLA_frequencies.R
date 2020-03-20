#!/usr/bin/env R

# Geographic maps of HLA data 
# Helpful resources
# https://www.r-spatial.org/r/2018/10/25/ggplot2-sf.html
# https://www.r-bloggers.com/map-visualization-of-covid-19-across-the-world-with-r/
# https://cran.r-project.org/web/packages/rworldmap/index.html

require(rvest)
require(rworldmap)

#--------------------
# get allele frequency data
#--------------------
alleles.outfile <- "HLA-freqs.rda"
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

#-----------------------
# aggregate allele data by country (ISO3)
#-----------------------
# fix ISO3 code match issues (England -> Great Britain, Northern Ireland -> Great Britain, Gaza -> Israel, Romania -> Romania)
alleles[which(alleles[,"country"]=="ENG"), "country"] <- "GBR"
alleles[which(alleles[,"country"]=="NIR"), "country"] <- "GBR"
alleles[which(alleles[,"country"]=="GAZ"), "country"] <- "ISR"
alleles[which(alleles[,"country"]=="ROM"), "country"] <- "ROU"
save(alleles,file=alleles.outfile)
countries <- get(data(countryExData))
hla.iso3 <- data.frame(ISO3=character(), HLA=character(), freq=numeric(), popsize=numeric())
#ISO3 <- unique(alleles[,"country"])
for (ISO3 in unique(alleles[,"country"])) {
	data <- alleles[which(alleles[,"country"] == ISO3),]
	switch(ISO3,
		# 2005 population statistics from World Bank, http://api.worldbank.org/v2/en/indicator/SP.POP.TOTL?downloadformat=csv & https://www.worldometers.info/world-population/
		MTQ = popsize <- 397.19,
		SRB = popsize <- 7441,
		HKG = popsize <- 6813.2,
		STP = popsize <- 155.63,
		SGP = popsize <- 4265.762,
		CPV = popsize <- 474.567,
		LBY = popsize <- 5798.614,
		ASM = popsize <- 59.118,
		NCL = popsize <- 232.25,
		GNQ = popsize <- 757.317,
		TWN = popsize <- 22705.713,
		# Canada HLA data only available for BC First Nations (Athabaskan and Penutian)
		# Canada First Nations are ~4% of national population 
		CAN = popsize <- 32243.753 * 0.04,
		popsize <- countries[which(countries[,1]==ISO3),"Population2005"]
	)
	for (HLA in unique(as.character(data[,"HLA"]))) {
		# aggregate codes (A*01:01 will include contributions from e.g. A*01:01:01 and A*01:01:02 in addition to A*01:01 itself)
		hla.data <- data[grepl(HLA,data[,"HLA"], fixed=TRUE),c("freq","sample.size")]
		# alternative approach is to only lump by exact allele match
		# hla.data <- data[which(data[,"HLA"]==HLA),c("freq","sample.size")]
		hla.iso3 <- rbind(hla.iso3, data.frame(ISO3=ISO3, HLA=HLA, freq=weighted.mean(as.numeric(hla.data[,1]), w=as.numeric(hla.data[,2])), popsize=popsize))
	}
}
save(hla.iso3,file="hla.iso3.rda")

#-----------------------
# estimate global allele frequency
#-----------------------
global_allele_freqs <- data.frame(HLA=character(), freq=numeric())
for (hla in as.character(unique(hla.iso3[,"HLA"]))) {
	hla.data <- hla.iso3[grepl(hla,hla.iso3[,"HLA"], fixed=TRUE),c("freq","popsize")]
	global_allele_freqs <- rbind(global_allele_freqs, data.frame(HLA=hla, freq=weighted.mean(as.numeric(hla.data[,1]), w=as.numeric(hla.data[,2]))))
}
save(global_allele_freqs, file="global.hla.freqs.rda")

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
	  			mapTitle=paste0(HLA," (~", format(global_allele_freqs[which(global_allele_freqs[,1]==hla),2]*100, digits=3),"%)"))
		do.call(addMapLegend, c(mapParams, legendWidth=0.5, legendMar=6.5))
		if (!is.null(outdir)) {
			dev.off()
		}
	}
	# avoid plotting messages (dev.off results) upon function return
	rm(data)
}
#----------
# generate figures for paper
#----------
pdf(file="plots/best+worst_HLA_all_peptides.pdf", width=9,height=8)
layout(matrix(1:6,nr=3,nc=2))
par(mar=c(1.1,0.1,1.1,0.1))
HLAgeoPlot(c("A*02:06", "A*02:07", "B*15:03", "B*46:01", "C*12:03", "C*01:02"))
dev.off()

pdf(file="plots/best+worst_HLA_conserved_peptides.pdf", width=9,height=8)
layout(matrix(1:6,nr=3,nc=2))
par(mar=c(1.1,0.1,1.1,0.1))
HLAgeoPlot(c("A*02:06", "A*02:07", "B*08:01", "B*46:01", "C*12:03", "C*01:02"))
dev.off()