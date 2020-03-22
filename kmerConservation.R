kmerConservation <- function(kmer.file, peptide.files, seqname="SARS-CoV-2") {
	kmer.data <- read.csv(kmer.file, header=TRUE, as.is=TRUE, comment.char="#")
	kmers <- unique(kmer.data[,"Peptide"])
	kmer.scores <- new.env(hash=TRUE)
	for (peptide.file in peptide.files) {
		data <- read.csv(peptide.file,header=FALSE,as.is=TRUE,row.names=1)
		sequences <- unlist(lapply(data[seqname,],"as.character"))[1:length(data[seqname,])]
		peptide <- gsub("[-]","",paste(sequences, collapse=""))
		ungapped <- which(sequences != "-")
		matches <- mapply(gregexpr, kmers, peptide)
		names(matches) <- kmers
		for (kmer in kmers) {
			kmer.matches <- matches[[kmer]] 
			if (any(kmer.matches < 0)) { next }
			N <- length(kmer.matches)	
			Qscore <- Bcons <- Hcons <- Cons <- 0
			for (j in 1:N) {
				Q <- B <- H <- C <- 0
				indices <- ungapped[kmer.matches[j]:(kmer.matches[j]+attr(kmer.matches, "match.length")[j]-1)]
				gaps <- setdiff(min(indices):max(indices), indices)
				#add scores for each successive amino acid in linear epitope
				for (k in indices) {
					if (grepl(data[seqname,k], gsub("[[]?([A-Z]*)[]]? [0-9]*[%]", "\\1", data["Consensus_seq", k]))) {
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
			if (!is.null(kmer.scores[[kmer]])) {
				scores <- kmer.scores[[kmer]]
				scores[["quality"]] <- c(scores[["quality"]], Qscore/nchar(kmer))
				scores[["conservation.beta"]] <- c(scores[["conservation.beta"]], Bcons/nchar(kmer))
				scores[["conservation.human"]] <- c(scores[["conservation.human"]], Hcons/nchar(kmer))
				scores[["conservation.combined"]] <- c(scores[["conservation.combined"]], Cons/nchar(kmer))
				kmer.scores[[kmer]] <- scores
			}
			else {
				kmer.scores[[kmer]] <- list(quality=Qscore/nchar(kmer), conservation.beta=Bcons/nchar(kmer), conservation.human=Hcons/nchar(kmer), conservation.combined=Cons/nchar(kmer))
			}
		}
	}
	kmer.data <- t(apply(kmer.data, 1, function (x) {
			scores <- kmer.scores[[x["Peptide"]]]
			if (is.null(scores)) { return(c(x, rep(NA, 4))) }
			return(c(x, max(scores[["quality"]]),
					max(scores[["conservation.beta"]]),
					max(scores[["conservation.human"]]),
					max(scores[["conservation.combined"]])))	
		}))
	colnames(kmer.data)[ncol(kmer.data)-3:0] <- c("quality", "conservation.beta", "conservation.human", "conservation.combined")
	kmer.data <- as.data.frame(kmer.data)
	kmer.data$Binding_affinity <- as.numeric(as.character(kmer.data$Binding_affinity))
	kmer.data$quality <- as.numeric(as.character(kmer.data$quality))
	kmer.data$conservation.beta <- as.numeric(as.character(kmer.data$conservation.beta))
	kmer.data$conservation.human <- as.numeric(as.character(kmer.data$conservation.human))
	kmer.data$conservation.combined <- as.numeric(as.character(kmer.data$conservation.combined))
	return(kmer.data)
}


#----------
# generate figures for paper
#----------
# Supplementary Figure XXX
peptide.files <- dir("alignments", pattern=".csv", full.names=TRUE)
data <- kmerConservation("Appendix_4.csv", peptide.files)
cor.test(data$Binding_affinity, data$quality) # -0.01748769 
cor.test(data$Binding_affinity, data$conservation.combined) # -0.02822426 
# group by individual peptide (summarize as best binding affinity across all HLA types)
data.min <- aggregate(data[,c("Binding_affinity","quality","conservation.combined")],by=list(data[,"Peptide"]),FUN=min)
plot(data.min$Binding_affinity, data.min$quality,pch=20, cex=0.25,xlab="Best peptide-HLA binding affinity (nM)",ylab="Peptide conservation (alignment quality)")
plot(data.min$Binding_affinity, data.min$conservation.combined,pch=20, cex=0.25,xlab="Best peptide-HLA binding affinity (nM)",ylab="Peptide conservation (BLOSUM62-based)")
# group by individual peptide (summarize as median binding affinity across all HLA types)
data.med <- aggregate(data[,c("Binding_affinity","quality","conservation.combined")],by=list(data[,"Peptide"]),FUN=median)
plot(data.med$Binding_affinity, data.med$quality,pch=20, cex=0.25,xlab="Median peptide-HLA binding affinity (nM)",ylab="Peptide conservation (alignment quality)")
plot(data.med$Binding_affinity, data.med$conservation.combined,pch=20, cex=0.25,xlab="Median peptide-HLA binding affinity (nM)",ylab="Peptide conservation (BLOSUM62-based)")
