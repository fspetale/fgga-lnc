###############################################################################
# FGGA-lnc                                          
# Dicember 2020
# Spetale, F. E. and Murillo J.
# Method for characterizing lncRNA sequences based on their secondary structure
###############################################################################
library(jsonlite)
library(seqinr)

bulge_loop <- function (ctFile) {
	RNAstructure <- as.matrix(ctFile[, 5])
	bulge_loops <- list()
	n <- 1
	arr0 <- RNAstructure[, 1]
	for (i in 1:(length(arr0) - 1)) {
  		loop_length <- abs(arr0[i] - arr0[i + 1])
		num_min <- min(arr0[i], arr0[i + 1]) + 1
		num_max <- max(arr0[i], arr0[i + 1]) - 1
		
		if (arr0[i] != 0 && arr0[i + 1] != 0 && loop_length != 1 && length(which(arr0[num_min:num_max] != 0)) == 0) {
			bulge_loops[[n]] <- num_min:num_max
			n <- n + 1
		}
	}
	if (length(bulge_loops) == 0) {
		return(list())
  	}  else {
		bulge_number<- length(bulge_loops)
		bulge_max 	<- max(sapply(bulge_loops,length)) 
		bulge_min 	<- min(sapply(bulge_loops,length))
    	attr(bulge_loops, "number of bases in bulge loops") <- length(unlist(bulge_loops)) 
    	attr(bulge_loops, "number of bulge loops") 			<- bulge_number
    	attr(bulge_loops, "Maximum length of bulge loops") 	<- bulge_max
    	attr(bulge_loops, "Minimum length of bulge loops") 	<- bulge_min
    	attr(bulge_loops, "Average length of bulge loops") 	<- length(unlist(bulge_loops))/bulge_number
    	return(bulge_loops)
  	}
}

hairpin_loop <- function(ctFile){
	RNAstructure <- matrix(c(as.numeric(ctFile[, 5]), as.numeric(ctFile[, 6])), ncol = 2, byrow = F)
	arr_line <- c()
	arr0 <- RNAstructure[, 1]
	hairpin_loops <- list()
	n <- 1
	for (i in 1:dim(RNAstructure)[1]) {
		num_min <- min(RNAstructure[i, 1], RNAstructure[i, 2])
		num_max <- max(RNAstructure[i, 1], RNAstructure[i, 2])
		arr <- arr0[(num_min + 1):(num_max - 1)]
		if (num_min != 0 && num_max != 0 && length(which(arr != 0)) == 0) {
			line 	<- paste0(RNAstructure[i, 1], "_", RNAstructure[i, 2])
			line2 	<- paste0(RNAstructure[i, 2], "_", RNAstructure[i, 1])
			arr_line<- c(arr_line, line)
			if (length(which(arr_line == line2)) == 0) {
				hairpin_loops[[n]] <- (num_min + 1):(num_max - 1)
				n <- n + 1
			}
		}
	}
  	if (length(hairpin_loops) == 0) {
		return(list())
	} else {
		hairpin_number 	<- length(hairpin_loops)
		hairpin_max 	<-max(sapply(hairpin_loops,length))
		hairpin_min 	<-min(sapply(hairpin_loops,length))
		attr(hairpin_loops, "number of bases in hairpin loops") <- length(unlist(hairpin_loops))
		attr(hairpin_loops, "number of hairpin loops") 			<- hairpin_number
		attr(hairpin_loops, "Maximum length of hairpin loops") 	<- hairpin_max
		attr(hairpin_loops, "Minimum length of hairpin loops") 	<- hairpin_min
		attr(hairpin_loops, "Average length of hairpin loops") 	<- length(unlist(hairpin_loops)) / hairpin_number
		return(hairpin_loops)
	}
}

internal_loop <- function (ctFile){
	RNAstructure <- matrix(c(as.numeric(ctFile[, 5]), as.numeric(ctFile[, 6])), ncol = 2, byrow = F)
	internal_loops <- list()
	boundStart <- 0
	boundEnd <- 0
	boundStart2 <- 0
	boundEnd2 <- 0
	n <- 1
	loops <- c()
	for (i in 1:(dim(RNAstructure)[1] - 1)) {
    	if (RNAstructure[i, 1] != 0 && RNAstructure[i + 1, 1] == 0) boundStart <- i 
    	if (RNAstructure[i, 1] == 0 && RNAstructure[i + 1, 1] != 0) boundEnd <- i + 1
		if (boundEnd != 0) {
			if (boundStart == 0) boundEnd <- 0 
			else {
				boundStart2 <- RNAstructure[boundStart, 1]
				boundEnd2 	<- RNAstructure[boundEnd, 1]
				num_min 	<- min(boundStart2, boundEnd2) + 1
				num_max 	<- max(boundStart2, boundEnd2) - 1
				loop 		<- RNAstructure[num_min:num_max, 1]
				if (length(which(loop != 0)) == 0 && (boundStart + 1) != num_min && (boundEnd - 1) != num_max) {
					if (length(which(loops == (min(boundStart, boundEnd) + 1))) == 0) {
						loops <- c(loops, boundStart, boundEnd, num_min, num_max)
						internal_loops[[n]] <- sort(c(num_min:num_max, (min(boundStart, boundEnd) + 1):(max(boundStart, boundEnd) - 1)))
						n <- n + 1
					}
				}
				boundStart 	<- 0
				boundEnd 	<- 0
        		boundStart2 <- 0
        		boundEnd2 	<- 0
      		}
		}
  	}
	if (length(internal_loops) == 0) {
		return(internal_loops)
	} else {
		internal_number <- length(internal_loops)
		internal_max <- max(sapply(internal_loops,length)) 
 		internal_min <- min(sapply(internal_loops,length)) 
		attr(internal_loops, "number of bases in internal loops") <- length(unlist(internal_loops))
		attr(internal_loops, "number of internal loops") <- internal_number
		attr(internal_loops, "Maximum length of internal loops") <- internal_max
		attr(internal_loops, "Minimum length of internal loops") <- internal_min
		attr(internal_loops, "Average length of internal loops") <- length(unlist(internal_loops))/internal_number
		return(internal_loops)
	}
}

multi_branch_loop <-function(ctFile, hairpin_bases, bulge_bases, internal_bases){
	hairpin_bases 	<- unlist(hairpin_bases)
	bulge_bases 	<- unlist(bulge_bases)
	internal_bases 	<- unlist(internal_bases)
	ctFile 			<- matrix(c(as.numeric(ctFile[, 5]), as.numeric(ctFile[, 6])), ncol = 2, byrow = F)
	pair_bases		<- which(ctFile[, 1] != 0)
	known_bases 	<- c(hairpin_bases, bulge_bases, internal_bases, pair_bases)
	multi_bases 	<- setdiff(ctFile[, 2], known_bases)
	if (length(multi_bases) != 0) {
		if (length(which(multi_bases == 1)) != 0) {
			num <- 1
			while (length(which(multi_bases == num)) != 0) {
				multi_bases <- multi_bases[-which(multi_bases == num)]
				num <- num + 1
			}
		}
		if (length(which(multi_bases == dim(ctFile)[1])) != 0) {
			num <- dim(ctFile)[1]
			while (length(which(multi_bases == num)) != 0) {
				multi_bases <- multi_bases[-which(multi_bases == num)]
				num <- num - 1
			}
		}
		if (length(multi_bases) != 0) {
			links_index 	<- which(ctFile[, 1] != 0)
			links 			<- lapply(1:length(links_index), function(x,y,z) c(y[x[z], 1], y[x[z], 2]), y=ctFile, x=links_index)
			links_index		<- links_index2<-length_links<-length(links_index)
			min_dist 		<- sapply(1:links_index2, function(x,y) min(y[[x]]), y=links)
			max_dist 		<- sapply(1:links_index2, function(x,y) max(y[[x]]), y=links)
			iter.conditional<- T
			while (iter.conditional) {
 				if (links_index != links_index2){
					min_dist[links_index2:links_index]<-sapply(links_index2:links_index, function(x,y) min(y[[x]]), y=links)
					max_dist[links_index2:links_index]<-sapply(links_index2:links_index, function(x,y) max(y[[x]]), y=links)
					length_links <-length(links)
				}
				for (j in 1:length_links) {
					init.index<-1
					if (links_index2!=length_links && j < links_index2) init.index<- links_index2
					for (k in init.index:length_links) {
						if ((abs(min_dist[j] - min_dist[k]) == 1 && abs(max_dist[j] - max_dist[k]) != 1) || 
							(abs(min_dist[j] - max_dist[k]) == 1 && abs(max_dist[j] - min_dist[k]) != 1) ||
							(abs(max_dist[j] - min_dist[k]) == 1 && abs(min_dist[j] - max_dist[k]) != 1) ||
							(abs(max_dist[j] - max_dist[k]) == 1 && abs(min_dist[j] - min_dist[k]) != 1)) {
							links_index 			<- links_index + 1
							links[[links_index]] 	<- c(min_dist[j], max_dist[j], min_dist[k], max_dist[k])
							links 					<- unique(links)
							if (links_index-length(links)==1) links_index <- links_index - 1
						}
					}
				}
				links_index2 <- length_links
				if (length(links) == links_index2) iter.conditional<-F
			}
			links_2 <- list()
			base1 <- multi_bases[1]
			n <- 1
			if (length(multi_bases)>1){
				for (i in 2:length(multi_bases)){
					if (multi_bases[i] == (max(base1) + 1)) {
						base1 <- c(base1, multi_bases[i])
					} else {
						links_2[[n]] <- base1
						n <- n + 1
 						base1 <- multi_bases[i]
					}
				}
				links_2[[n]] <- base1
				links_4 <- list(c(0))
        
				while (!identical(links_4, links_2) && length(links_2) > 1) {
					min_dist2 <- sapply(1:length(links_2), function(x,y) min(y[[x]]), y=links_2)
					max_dist2 <- sapply(1:length(links_2), function(x,y) max(y[[x]]), y=links_2)
					links_4 <- links_2
					for (i in 1:(length(links_2) - 1)) {
						for (j in (i + 1):length(links_2)) {
							for (k in 1:length(links)) {
								if ((max_dist2[i] == (min_dist[k] - 1) && min_dist2[j] == (max_dist[k] + 1)) || 
									(min_dist2[i] == (min_dist[k] + 1) && max_dist2[j] == (max_dist[k] - 1))) {
									links_2[[j]] <- c(links_2[[i]], links_2[[j]])
									if (length(which(links_2[[j]]==0)>0)) links_2[[j]] <-links_2[[j]][which(links_2[[j]]>0)]
									links_2[[i]] <- c(0)
								}
							}
						}
					}
				}
				links_nonzero <-sapply(1:length(links_2), function(x,y) length(which(y[[x]]==0)>0), y=links_2)
				links_nonzero <- which(links_nonzero==0)
				if (length(links_nonzero) != 0) {
					links_2 <- links_2[links_nonzero]
				}
				links_min <- min(unlist(links))
				links_max <- max(unlist(links))
				min_dist2 <- sapply(1:length(links_2), function(x,y) min(y[[x]]), y=links_2)
				max_dist2 <- sapply(1:length(links_2), function(x,y) max(y[[x]]), y=links_2)
				for (i in 1:length(links_2)) {
					for (j in 1:length(links)) {
						if ((max_dist2[i] == min_dist[j] - 1) || (min_dist2[i] == (max_dist[j] + 1))) {
							if ((links_max %in% links[[j]]) || (links_min %in% links[[j]])) {
								links_2[[i]] <- NA
							}
						}
					}
				}
			}
			links_2 <- links_2[!is.na(links_2)]
			multi_bases <- unique(links_2)
			if (length(links_2) == 0) {
				return(links_2)
			} else {
				multi_branch_number <- length(links_2)
				multi_branch_max   <- max(sapply(links_2,length))		
				multi_branch_min   <- min(sapply(links_2,length))		
				attr(links_2, "number of bases in multi_branch loops") <- length(unlist(links_2))
				attr(links_2, "number of multi_branch loops") <- multi_branch_number
				attr(links_2, "Maximum length of multi_branch loops") <- multi_branch_max
				attr(links_2, "Minimum length of multi_branch loops") <- multi_branch_min
				attr(links_2, "Average length of multi_branch loops") <- length(unlist(links_2)) / multi_branch_number
				return(links_2)
			}
		} else {
			return(list())
		}
	} else {
		return(list())
	}
}

stem <-function (ctFile) {
	RNAstructure <- matrix(c(as.numeric(ctFile[, 5]), as.numeric(ctFile[, 6])), ncol = 2, byrow = F)
	stem_list <- list()
	n <- 0
	stem_arr <- c()
	indexPair <- which(RNAstructure[, 1] != 0)
	numPair <- length(indexPair)
	if (numPair == 0) {
		print("There is no stems")
		return(stem_list)
	} else {
		if (numPair == 1) {
			n <- n + 1
			stem_arr <- c(stem_arr, RNAstructure[indexPair, 1], RNAstructure[indexPair, 2])
			stem_list[[n]] <- sort(stem_arr)
		} else {
			for (index in 2:numPair) {
				if (abs(RNAstructure[indexPair[index], 1] - RNAstructure[indexPair[index - 1], 1]) == 1 &&
					abs(RNAstructure[indexPair[index], 2] - RNAstructure[indexPair[index - 1], 2]) == 1) {
            		stem_arr <- c(stem_arr, RNAstructure[indexPair[index - 1], 1], RNAstructure[indexPair[index - 1], 2])
				} else {
					listL <- length(stem_list)
					if (listL == 0) {
						n <- n + 1
						stem_arr <- sort(unique(stem_arr))
						stem_list[[n]] <- stem_arr
						stem_arr <- c(RNAstructure[indexPair[index], 1], RNAstructure[indexPair[index], 2])
					} else {
						count <- 0
						for (i in 1:listL) {
							if (length(setdiff(stem_arr, stem_list[[i]])) == 0 && length(setdiff(stem_list[[i]], stem_arr)) == 0) {
								count <- count + 1
								break
							}
						}
						if (count == 0) {
							n <- n + 1
							stem_arr <- sort(unique(stem_arr))
							stem_list[[n]] <- stem_arr
							stem_arr <- c(RNAstructure[indexPair[index], 1], RNAstructure[indexPair[index], 2])
						} else {
							stem_arr <- c(RNAstructure[indexPair[index], 1], RNAstructure[indexPair[index], 2])
						}
					}
				}
			}
		}
		if (listL>2) stem_nonzero <-sapply(1:listL, function(x,y) length(y[[x]]), y=stem_list) 
		else stem_nonzero<-length(stem_list[[1]])
		stem_nonzero <- which(stem_nonzero>0)
		if (length(stem_nonzero) != 0) {
			stem_list <- stem_list[stem_nonzero]
		}
		stem_number <- length(stem_list)
	    stem_max	<- max(sapply(stem_list,length))		
	    stem_min	<- min(sapply(stem_list,length))		
		attr(stem_list, "number of bases in stems") <- length(unlist(stem_list))
		attr(stem_list, "number of stems") 			<- stem_number
		attr(stem_list, "Maximum length of stems") 	<- stem_max
		attr(stem_list, "Minimum length of stems") 	<- stem_min
		attr(stem_list, "Average length of stems") 	<- length(unlist(stem_list))/stem_number
		return(stem_list)
	}
}

external_loop<-function(ctFile, Hairpin_loop, Internal_loop, Bulge_loop, Multi_branch_loop, Stem){
	Hairpin_loop		<- unlist(Hairpin_loop)
	Internal_loop		<- unlist(Internal_loop)
	Bulge_loop			<- unlist(Bulge_loop)
	Multi_branch_loop 	<- unlist(Multi_branch_loop)
	Stem 				<- unlist(Stem)
	loopArr 			<- c(Hairpin_loop, Internal_loop, Bulge_loop, Multi_branch_loop, Stem)
	externalArr 		<- setdiff(1:dim(ctFile)[1], loopArr)
	if (length(externalArr) == 0) {
		print("There is no external loops")
		externalList <- list()
	} else {
		externalList <- list()
		tempArr <- externalArr[1]
		n <- 0
		if (length(externalArr) >= 2) {
			for (i in 2:length(externalArr)) {
				if (externalArr[i - 1] == (externalArr[i] - 1)) {
					tempArr <- c(tempArr, externalArr[i])
				} else {
					n <- n + 1
					externalList[[n]] <- tempArr
					tempArr <- externalArr[i]}
				}
		}
		n <- n + 1
		externalList[[n]] <- tempArr
		external_Num <- length(externalList)
		external_Max <-   max(sapply(externalList,length))
		external_Min <-   min(sapply(externalList,length))
		attr(externalList, "number of bases in external loops") <- length(unlist(externalList))
		attr(externalList, "number of external loops") 			<- external_Num
		attr(externalList, "Maximum length of external loops") 	<- external_Max
		attr(externalList, "Minimum length of external loops") 	<- external_Min
		attr(externalList, "Average length of external loops") 	<- length(unlist(externalList))/external_Num
    }
	return(externalList)}
	ct2dot <- function (ctFile) {
	ctFile <- as.matrix(ctFile)
	base <- ctFile[, 2]
	dot <- list()
	for (i in 1:dim(ctFile)[1]) {
		if (as.numeric(ctFile[i, 5]) == 0) {
			dot[[i]] <- "." 
		} else if (as.numeric(ctFile[i, 6]) < as.numeric(ctFile[i, 5])) {
      		dot[[i]] <- "("
      		dot[[as.numeric(ctFile[i, 5])]] <- ")"
    	}
	}
	dot <- paste(unlist(dot), collapse = "", sep = "")
	base <- paste(base, collapse = "", sep = "")
	re_list <- list()
	re_list[[1]] <- dot
	re_list[[2]] <- base
	return(re_list)
}

convertRNAKmers<-function(seq.RNA, sec.struc.RNA, type_SS = NULL){
	Matrix.15symbols			<-matrix(0,3,5)
	colnames(Matrix.15symbols)	<-c("a", "c", "g", "u", "n")
	rownames(Matrix.15symbols)	<-c("(",")", ".")
	Matrix.15symbols[1,]		<-c("A", "D", "E", "I", "F")
	Matrix.15symbols[2,]		<-c("R", "C", "G", "L", "S")
	Matrix.15symbols[3,]		<-c("N", "Q", "H", "K", "P")
	seq.RNA <- s2c(seq.RNA)
	NucleotideT<-which(seq.RNA=="t")
	if (length(NucleotideT)>0) seq.RNA[NucleotideT]<-"u"
	sec.struc.RNA <- s2c(sec.struc.RNA)
	seq.new <- array("", length(seq.RNA))
	for (i in 1:length(seq.new)) seq.new[i] <- Matrix.15symbols[sec.struc.RNA[i],seq.RNA[i]]
	if (length(type_SS)>0) {
		Matrix.15symbols<-as.vector(Matrix.15symbols)
		Matrix.SS.symbols<-matrix(seq(1, 6*15),6,15,byrow = T)
		colnames(Matrix.SS.symbols)<-Matrix.12symbols
		rownames(Matrix.SS.symbols)<-c("h","s", "b", "i", "e", "m")
		seq.new2 <- array(0, length(seq.RNA))
		for (i in 1:length(seq.new2)) seq.new2[i] <- Matrix.SS.symbols[type_SS[i],seq.new[i]]
		RNA.Dipeptide<-2*count(seq.new2, 2, alphabet =seq(1, 6*15))/length(seq.new2)
		RNA.single<-count(seq.new2, 1, alphabet = seq(1, 6*15))/length(seq.new2)
	} else {
		RNA.Dipeptide<-2*count(seq.new, 2, alphabet = s2c("ADEIRCGLNQHKFSP"))/length(seq.new)
		RNA.single<-count(seq.new, 1, alphabet = s2c("ADEIRCGLNQHKFSP"))/length(seq.new)
	}
	return(c(RNA.single,RNA.Dipeptide))
}

lncRNA_SS2 <- function(path, namelncRNA){
	Mypath<-paste0(path, "/", namelncRNA, sep="", collapse="" )
	rnaInfo <- fromJSON(paste0("https://rnacentral.org/api/v1/rna/", namelncRNA, ".json", sep=""), simplifyVector = T)
	RNA.new<-array(0,258)
	if(rnaInfo$length<6000){
		system(paste0("mkdir ", Mypath,  sep="", collapse="" ))
		setwd(Mypath)
		fasta.file<- paste0(Mypath,"/Seq_", namelncRNA , ".fasta", sep="", collapse = "")
		write.fasta(sequences = rnaInfo$sequence, names=namelncRNA, as.string = TRUE, nbchar=60, file.out = fasta.file)
		system(paste0("mfold SEQ='",fasta.file, "' NA=RNA LC='circular' ANN='p-num' MAX=2", sep="", collapse = ""), intern = T)
		fasta.file	<- paste0(Mypath, "/Seq_",namelncRNA, "_1.ct", sep="", collapse = "")
		fasta.file2	<- paste0(Mypath,"/Seq_",namelncRNA, "_3.ct", sep="", collapse = "")
		system(paste0("sed -E  \'s/\\s+/ /g\' ", fasta.file," > ", fasta.file2, sep="", collapse = ""), intern=F)
		ctFile<-read.table(fasta.file2, sep=" ", dec=".", skip=1, stringsAsFactors =F)
  		#ctFile<-ctFile[,-1]
		ctFile[,2]		<- system("tr [:lower:] [:upper:]",input=as.character(ctFile[,2]),intern=TRUE)
		bulge 			<- bulge_loop(ctFile)
		hairpin 		<- hairpin_loop(ctFile)
		internal 		<- internal_loop(ctFile)
		multi_branch	<- multi_branch_loop(ctFile, hairpin, bulge, internal)
		stems 			<- stem(ctFile)
		external		<- external_loop(ctFile,hairpin, bulge, internal, multi_branch, stems)
		dot				<- ct2dot(ctFile)
		dot[[2]]		<- system("tr [:upper:] [:lower:]",input=dot[[2]],intern=TRUE)
		RNA.new[1:240]	<- convertRNAKmers(dot[[2]],dot[[1]])
		if (length(bulge)>0)   RNA.new[241:243]	<-	c( attr(bulge, "number of bases in bulge loops"), attr(bulge, "number of bulge loops"), attr(bulge, "Average length of bulge loops") ) / dim(ctFile)[1]
  		if (length(hairpin)>0) RNA.new[244:246] <- c( attr(hairpin, "number of bases in hairpin loops"), attr(hairpin, "number of hairpin loops"), attr(hairpin, "Average length of hairpin loops") ) / dim(ctFile)[1]
		if (length(internal)>0) RNA.new[247:249] <- c(attr(internal, "number of bases in internal loops"), attr(internal, "number of internal loops"), attr(internal, "Average length of internal loops") ) / dim(ctFile)[1]
  		if (length(multi_branch)>0) RNA.new[250:252] <- c( attr(multi_branch, "number of bases in multi_branch loops"), attr(multi_branch, "number of multi_branch loops"), attr(multi_branch, "Average length of multi_branch loops") ) / dim(ctFile)[1]
  		RNA.new[253:255] <- c(attr(stems, "number of bases in stems"), attr(stems, "number of stems"), attr(stems, "Average length of stems"))/dim(ctFile)[1]
  		if (length(external)>0) RNA.new[256:258]<- c( attr(external, "number of bases in external loops"), attr(external, "number of external loops"), attr(external, "Average length of external loops") ) / dim(ctFile)[1]
	}
	return(RNA.new)
}
