kmerConservation <- function(kmer.file, peptide.files) {
	kmer.data <- read.csv(kmer.file, header=TRUE, as.is=TRUE)
	kmers <- unique(kmer.data[,"peptide"])
	kmer.scores <- list()
	for (peptide.file in peptide.files) {
		data <- read.csv(peptide.file,header=FALSE,as.is=TRUE,row.names=1)
		sequences <- unlist(lapply(data[1,],"as.character"))[1:length(data[1,])]
		peptide <- gsub("[-]","",paste(sequences, collapse=""))
		ungapped <- which(sequences != "-")
		matches <- mapply(gregexpr, kmers, peptide)
		for (kmer in kmers) {
			kmer.matches <- matches[[kmer]] 
			if (any(kmer.matches < 0)) {
				next
			}
			N <- length(kmer.matches)		
			Qscore <- Bcons <- Hcons <- Cons <- 0
			for (j in 1:N) {
				Q <- B <- H <- C <- 0
				indices <- ungapped[kmer.matches[j]:(kmer.matches[j]+attr(kmer.matches, "match.length")[j]-1)]
				gaps <- setdiff(min(indices):max(indices), indices)
				#add scores for each successive amino acid in linear epitope
				for (k in indices) {
					if (grepl(data["SARS-CoV-2",k], gsub("[[]?([A-Z]*)[]]? [0-9]*[%]", "\\1", data["Consensus_seq", k]))) {
						Q <- Q + as.numeric(data["Quality", k])
						B <- B + as.numeric(data["Beta-CoV cons", k])
						H <- H + as.numeric(data["Human-CoV cons", k])
						C <- C + as.numeric(data["Conservation", k])
					}
				}
				#penalize for any gaps in sequence that may disrupt linear epitope among a subset of aligned sequences
				Qgap <- Bgap <- Hgap <- Cgap <- 0
				for (k in gaps) {
					Qgap <- max(Qgap, as.numeric(data["Quality", k]))
					Bgap <- max(Bgap, as.numeric(data["Beta-CoV cons", k]))
					Hgap <- max(Hgap, as.numeric(data["Human-CoV cons", k]))
					Cgap <- max(Cgap, as.numeric(data["Conservation", k]))
				}
				Qscore <- max(Q-Qgap, Qscore, na.rm=TRUE)
				Bcons <- max(B-Bgap, Bcons, na.rm=TRUE)
				Hcons <- max(H-Hgap, Hcons, na.rm=TRUE)
				Cons <- max(C-Cgap, Cons, na.rm=TRUE)
			}
			# if peptide occurs in multiple proteins it's already scored, addend with another score
			if (hasName(kmer.scores, kmer)) {
				kmer.scores[[kmer]][["quality"]] <- c(kmer.scores[[kmer]][["quality"]], Qscore/nchar(kmer))
				kmer.scores[[kmer]][["conservation.beta"]] <- c(kmer.scores[[kmer]][["conservation.beta"]], Bcons/nchar(kmer))
				kmer.scores[[kmer]][["conservation.human"]] <- c(kmer.scores[[kmer]][["conservation.human"]], Hcons/nchar(kmer))
				kmer.scores[[kmer]][["conservation.combined"]] <- c(kmer.scores[[kmer]][["conservation.combined"]], Cons/nchar(kmer))
			}
			else {
				kmer.scores[[kmer]] <- list(quality=Qscore/nchar(kmer), conservation.beta=Bcons/nchar(kmer), conservation.human=Hcons/nchar(kmer), conservation.combined=Cons/nchar(kmer))
			}
		}
	}
	kmer.data <- cbind(kmer.data, matrix(nr=dim(kmer.data)[1], ncol=4, dimnames=list(c(), c("quality", "conservation.beta", "conservation.human", "conservation.combined"))))
	for (k in 1:(dim(kmer.data)[1])) {
		scores <- kmer.scores[[kmer.data[k, "peptide"]]]
		if (is.null(scores)) { next }
		kmer.data[k, "quality"] <- max(scores[["quality"]])
		kmer.data[k, "conservation.beta"] <- max(scores[["conservation.beta"]])
		kmer.data[k, "conservation.human"] <- max(scores[["conservation.human"]])
		kmer.data[k, "conservation.combined"] <- max(scores[["conservation.combined"]])
	}
	return(kmer.data)
}