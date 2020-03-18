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
countries <- get(data(countryExData))
popIDs <- unique(alleles[,"pop.ID"])
alleles <- cbind(alleles, data.frame(lat=rep(NA, dim(alleles)[1]), long=rep(NA, dim(alleles)[1]), region=rep(NA, dim(alleles)[1]), country=rep(NA, dim(alleles)[1])))
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
	alleles[which(alleles[,"pop.ID"]==pop), c("lat", "long", "region", "country")] <- list(lat=lat, long=long, region=region, country=country)
}
save(alleles,file=alleles.outfile)

#-----------------------
# aggregate allele data by country (ISO3)
#-----------------------
# fix ISO3 code errors (England -> Great Britain)
alleles[which(alleles[,"country"]=="ENG"), "country"] <- "GBR"
hla.iso3 <- data.frame(ISO3=character(), HLA=character(), freq=numeric())
#ISO3 <- unique(alleles[,"country"])
for (ISO3 in unique(alleles[,"country"])) {
	data <- alleles[which(alleles[,"country"] == ISO3),]
	for (HLA in unique(as.character(data[,"HLA"]))) {
		# aggregate codes (A*01:01 will include contributions from e.g. A*01:01:01 and A*01:01:02 in addition to A*01:01 itself)
		hla.data <- data[grepl(HLA,data[,"HLA"]),c("freq","sample.size")]
		# alternative approach is to only lump by exact allele match
		# hla.data <- data[which(data[,"HLA"]==HLA),c("freq","sample.size")]
		hla.iso3 <- rbind(hla.iso3, data.frame(ISO3=ISO3, HLA=HLA, freq=weighted.mean(as.numeric(hla.data[,1]), w=as.numeric(hla.data[,2]))))
	}
}
save(hla.iso3,file="hla.iso3.rda")

#----------
# plot data
#----------
# note this will generate >4,000 plots!!
#for (HLA in unique(hla.iso3[,"HLA"])) {
#	sPDF <- joinCountryData2Map(hla.iso3[which(hla.iso3[,"HLA"]==HLA),], joinCode="ISO3", nameJoinColumn="ISO3")
#  	jpeg(paste0(pdn, "/", 
#             gsub("\\:", "-", gsub("\\*", "", HLA)), ".jpg"), 
#       4, 3, units = "in", res = 400)
#  	mapParams <- mapCountryData(sPDF, nameColumnToPlot="freq", missingCountryCol="gray", borderCol="black", oceanCol
#="lightblue", catMethod="pretty",addLegend=FALSE, mapTitle=HLA)
#	do.call( addMapLegend, c(mapParams, legendWidth=0.5, legendMar = 2))
#  	  dev.off()
#}

HLAgeoPlot <- function(HLA, outdir=NULL) {
	sPDF <- joinCountryData2Map(hla.iso3[which(hla.iso3[,"HLA"]==HLA),], joinCode="ISO3", nameJoinColumn="ISO3")
	if (!is.null(outdir)) {
		dir.create(outdir, showWarnings = FALSE)
	  	jpeg(paste0(outdir, "/", 
             gsub("\\:", "-", gsub("\\*", "", HLA)), ".jpg"), 
    	   4, 3, units = "in", res = 400)
	}
  	mapParams <- mapCountryData(sPDF, nameColumnToPlot="freq", missingCountryCol="gray", 
  			borderCol="black", oceanCol="lightblue", catMethod="pretty",
  			addLegend=FALSE, mapTitle=HLA)
	do.call(addMapLegend, c(mapParams, legendWidth=0.5, legendMar=2))
	if (!is.null(outdir)) {
		dev.off()
	}
}
