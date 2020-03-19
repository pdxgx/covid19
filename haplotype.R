get_hapsearch <- function(aid = "", popminsize = "", orderby = "hapfreq",
                          xpathval = '//*[@id="divGenDetail"]/table'){
  require(rvest) # lib for web scraping
  urlstr <- c("http://www.allelefrequencies.net/hla6003a_scr.asp?",
              "&hla_population=&hla_country=&hla_dataset=&hla_region=&hla_ethnic=&hla_study=",
              "&hla_sample_size_pattern=equal",
              "&hla_sample_year_pattern=equal&hla_sample_year=&hla_loci=")
  l1 <- ifelse(grepl("A", aid), aid, "A*") # first locus
  l2 <- ifelse(grepl("B", aid), aid, "B*")
  l3 <- ifelse(grepl("C", aid), aid, "C*")
  l4 <- ifelse(grepl("DRB1", aid), aid, "")
  l5 <- ifelse(grepl("DPA1", aid), aid, "")
  l6 <- ifelse(grepl("DPB1", aid), aid, "")
  l7 <- ifelse(grepl("DOA1", aid), aid, "")
  l8 <- ifelse(grepl("DOB1", aid), aid, "") # final locus
  # set the order variable
  ordervar <- ifelse(orderby == "hapfreq", "order_3", "order_1")
  # make query url
  qurl <- paste0(urlstr[1],
                 "hla_locus1=", l1,
                 "&hla_locus2=", l2,
                 "&hla_locus3=", l3,
                 "&hla_locus4=", l4,
                 "&hla_locus5=", l5,
                 "&hla_locus6=", l6,
                 "&hla_locus7=", l7,
                 "&hla_locus8=", l8,
                 urlstr[2],
                 "&hla_sample_size=", popminsize,
                 urlstr[3],
                 "&hla_order=", ordervar,
                 urlstr[4])
  # scrape data and make table
  hd <- qurl %>%
    html() %>%
    html_nodes(xpath = xpathval) %>%
    html_table()
  ht <- hd[[1]] # extract the table data
  return(ht)
}

test <- get_hapsearch("B*46:01")
