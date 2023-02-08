set.seed(120)
options(dplyr.summarise.inform = FALSE)
options(tidyverse.quiet = TRUE)



library(igraph)
library(Biostrings)
library(matrixStats)

library(data.table)
library(purrr)
library(dplyr)
library(tidyr)
library(data.table)
library(stringr)
library(cluster)
library(igraph)
library(Matrix)
library(dbscan)
library(ggplot2)
library(NbClust)
library(factoextra)
library(Rtsne)
library(factoextra)
library(furrr)
library(matrixStats)
library(Biostrings)
library(future)
library(purrr)

library(data.table)
#library(tidyverse)
library(stringdist)
library(kmer)
library(magrittr)

#------------------------------------------------------
get_scores_multi<-function(qfasta,qlen){
HMM_append<-list()
for (k in (1:qlen)){
	fasta_score=numeric(length(lPHMM))
	sc <- numeric(length(lPHMM)) # scores (log probabilities)
	threshold_score<- numeric(length(lPHMM)) 

	for (count in (1:length(lPHMM ))){
		print(paste0("count:",count,"      ","k:",k))
		fasta_score[count]<-aphid::forward(lPHMM[[count]],(qfasta[k]),residues="DNA",odds=TRUE)$score
		sc[count] <- aphid::forward(lPHMM[[count]],qfasta[k],residues="DNA",odds=F)$score
		threshold_score[count]<-aphid::forward(AllGPHMM,qfasta[k],residues="DNA",odds=TRUE)$score
	}
	total_score <- aphid::logsum(sc)
	#fasta_score<-lapply(fasta_score, function(x) replace(x, !is.finite(x), 0))
	#threshold_score<-lapply(threshold_score, function(x) replace(x, !is.finite(x), 0))
	
	lo<-which.max(fasta_score)
	selectedlphmm<-fasta_score[lo]
	ss<-softmax(fasta_score)
	best_ss<-which.max(ss)
	selectedsoftmax<-ss[best_ss]
	akwgts <- exp(sc- total_score)
  	best_model <- which.max(akwgts)
  	newakw <- akwgts[best_model]
  	logoddsbyakw <- fasta_score[best_model]

	df_score<-cbind(lo,selectedlphmm,best_ss,selectedsoftmax,newakw,best_model,logoddsbyakw )
	HMM_append<-c(HMM_append,list(df_score))
}
HMM_append1<-purrr::map(HMM_append,~as.data.frame(.x))%>% rbindlist(.,idcol="source")

return(HMM_append1)
}
#------------------------------------------------------

#softmax function
softmax <- function(par){
  n.par <- length(par)
  par1 <- sort(par, decreasing = TRUE)
  Lk <- par1[1]
  for (k in 1:(n.par-1)) {
    Lk <- max(par1[k+1], Lk) + log1p(exp(-abs(par1[k+1] - Lk))) 
  }
  val <- exp(par - Lk)
  return(val)
}
#------------------------------------------------------
maxn <- function(x, n) {
  partial <- length(x) - n + 1
  x[x >= sort(x, partial = partial)[partial]]
}
#------------------------------------------------------



## convert binary DNAbin object to character easily
dna2char <- function(x){
  cbytes <- as.raw(c(65, 84, 71, 67, 83, 87, 82, 89, 75, 77, 66, 86, 72, 68, 78, 45, 63))
  indices <- c(136, 24, 72, 40, 96, 144, 192, 48, 80 ,160, 112, 224, 176, 208, 240, 4, 2)
  vec <- raw(240)
  vec[indices] <- cbytes
  if(is.list(x)){
    res <- sapply(x, function(s) rawToChar(vec[as.integer(s)]))
    if(!is.null(attr(x[[1]], "quality"))){
      attr(res, "quality") <- unname(sapply(x, function(s) .qual2char(attr(s, "quality"))))
    }
  }else{
    res <- rawToChar(vec[as.integer(x)])
    if(!is.null(attr(x, "quality"))){
      attr(res, "quality") <- unname(.qual2char(attr(x, "quality")))
    }
  }
  attr(res, "rerep.names") <- attr(x, "rerep.names")
  attr(res, "rerep.pointers") <- attr(x, "rerep.pointers")
  return(res)
}

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ 

multiple_scores<-function(sequences_list){
  ss<-future_map(sequences_list,~
      as.DNAbin.character(str_split(.x,"")%>%
                      unlist()))
result<-future_map(ss,~get_scores(.x,qlen=1)%>% as.data.frame)%>% rbindlist(idcol="source")
return(result)
}




#---------------------
#scores<-get_scores(qfasta,qlen=1)


get_scores<-function(qfasta,qlen){

#cutoff_calc=data.frame()
score_withself=data.frame()
softmax_score<-list()
HMM_score<-list()
HMM_index<-list()
HMM_append<-list()
j=0
	fasta_score=numeric(length(lPHMM))
	sc <- numeric(length(lPHMM)) # scores (log probabilities)
	threshold_score<- numeric(length(lPHMM)) 
	for (count in (1:length(lPHMM))){
		#print(paste0("count:",count,"      ","k:",k))
		fasta_score[count]<-forward(lPHMM[[count]],(qfasta),residues="DNA",odds=TRUE)$score
		sc[count] <- forward(lPHMM[[count]],qfasta,residues="DNA",odds=F)$score
		#threshold_score[count]<-forward(AllGPHMM,qfasta,residues="DNA",odds=TRUE)$score
		}
	
	total_score <- aphid::logsum(sc)
	#fasta_score<-lapply(fasta_score, function(x) replace(x, !is.finite(x), 0))
	#threshold_score<-lapply(threshold_score, function(x) replace(x, !is.finite(x), 0))
	
	lo<-which.max(fasta_score)
	selectedlphmm<-fasta_score[lo]
	ss<-softmax(fasta_score)
	best_ss<-which.max(ss)
	selectedsoftmax<-ss[best_ss]
	akwgts <- exp(sc- total_score)
  	best_model <- which.max(akwgts)
  	newakw <- akwgts[best_model]
  	logoddsbyakw <- fasta_score[best_model]

	df_score<-cbind(lo,selectedlphmm,best_ss,selectedsoftmax,newakw,best_model,logoddsbyakw )

return(df_score)
}


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

checkg4<- function(a){
 	a4<-"(?=(G{3,7}?[A|C|G|T|U|N]{1,10}?G{3,7}?[A|C|G|T|U|N]{1,10}?G{3,7}?[A|C|G|T|U|N]{1,10}?G{3,7}))"
	#a4 <- "G{2,7}?[A|C|G|T|U|N]{1,10}?G{3,7}?[A|C|G|T|U|N]{1,10}?G{3,7}?[A|C|G|T|U|N]{1,10}?G{2,7}?"
	#a5 <- gregexpr(a4,  a, ignore.case = TRUE,perl=TRUE)
	a5<-stringi::stri_match_all_regex(a,a4)
	a5<-lapply(a5,function(x) x[,2]%>% as.data.frame())%>% rbindlist(idcol = "source")%>%
		 as.data.frame()%>% filter(!is.na(.))
	colnames(a5)<-c("seqnum","sequence")

	sequence<-a5$sequence%>% as.character()
	sequence_length <- nchar(sequence)

	sequence_position <- a5$seqnum%>% as.character()

	a10 <- cbind(sequence_position, sequence, sequence_length)
  	a11 <- a10[as.integer(as.character(a10[,3]))<50,]
  
 	 if(length(a11) > 3){
    		a12 <- as.data.frame(a11, stringsAsFactors=FALSE)
    		a12 [,1] <- as.numeric(as.character(a12 [,1]))
    		a12 [,3] <- as.integer(as.character(a12 [,3]))
    		a13 <- a12[order(a12[,1]),]
    		return(a13)
  	}  
  
 	 if(length(a11) == 3){
    		if (as.integer(as.character(a11[[3]])) < 3 ){
      			a14 <- data.frame("sequence_position" = "-", "sequence" = "-", "sequence_length" = "-")
      		return(a14)
    	}	else{
      			a14 <- data.frame("sequence_position" = a11[[1]], "sequence" = a11[[2]], "sequence_length" = a11[[3]])
      		return(a14)
    		}
  	}
  if(length(a11)  < 3 ){
    		a14 <- data.frame("sequence_position" = "-", "sequence" = "-", "sequence_length" = "-")
    		return(a14)
  		}
	}
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
gquad_main  <- function(b){

  if(length(b) == 1){
    b1 <- checkg4(b)
    return(b1)
  }else{
    input_pos = 0
    q <- data.frame("input_ID" = integer(0), "sequence_position" = character(0), "sequence" = character(0), "sequence_length" = character(0))
    for(i in b){
      b1 <- checkg4(i)
      input_pos = input_pos + 1
      b2 <- cbind(input_ID = input_pos, b1)
      b2[,c(2,4)] <- sapply(b2[,c(2,4)],as.character)
      q <- rbind(q, b2)
	q<-q%>% distinct(input_ID,sequence, .keep_all=TRUE)

    }
    return(q)
  }
}
#-------------------------------------------------------------------------------------------------

msaTrim <- function(msa, gap.end = 0.95, gap.mid = 0.99){
  nc <- unique(str_length(msa$Sequence))
  if(length(nc) != 1) stop("This is not a multiple alignment, sequences have different lengths!")
  cmat <- str_split(msa$Sequence, pattern = "", simplify = T)
  gap.frac <- colSums(cmat == "-")/nrow(cmat)
  if(min(gap.frac) > min(gap.end, gap.mid)){
    warning("All positions have more gaps than arguments allow, please increase gap.end value")
    msa$Sequence <- ""
    return(msa)
  }
  idx.keep <- (gap.frac <= gap.mid)
  j <- 1
  while(gap.frac[j] > gap.end){
    idx.keep[j] <- F
    j <- j + 1
  }
  j <- nc
  while(gap.frac[j] > gap.end){
    idx.keep[j] <- F
    j <- j - 1
  }
  msa$Sequence <- apply(cmat[,idx.keep], 1, paste, collapse="")
  return(msa)
}

readFasta <- function(in.file){
  if(!is.character(in.file) | length(in.file) > 1) stop("The argument in.file must be a single text (a filename)")
  if(file.exists(in.file)){
    in.file <- normalizePath(in.file)
    fread(in.file, header = F, sep = "\t", data.table = F, na.strings = "") %>% 
      filter(!is.na(.data$V1)) -> tbl
    hh <- grep("^>", tbl[,1])
    hhh <- c(hh, nrow(tbl) + 1)
    tibble(Header = str_remove(tbl[hh,1], "^>"),
           Sequence = sapply(1:length(hh), function(i){str_c(tbl[(hhh[i]+1):(hhh[i+1]-1),1], collapse = "")})) -> fdta
    return(fdta)
  } else {
    stop("Cannot find ", in.file, ", please correct path and/or file name")
  }
}
writeFastausedhere <- function(fdta, out.file, width = 0){
  if(!("Header" %in% colnames(fdta))){
    stop("This is not a Fasta object, column Header is lacking\n")
  }
  if(!("Sequence" %in% colnames(fdta))){
    stop("This is not a Fasta object, column Sequence is lacking\n")
  }
  out.file <- file.path(normalizePath(dirname(out.file)),
                        basename(out.file))
  lst <- vector("list", 2 * nrow(fdta))
  lst[seq(1, 2 * nrow(fdta), 2)] <- as.list(str_c(">", fdta$Header))
  if(width > 0){
    lst[seq(2, 2 * nrow(fdta), 2)] <- lapply(fdta$Sequence, chop, width)
  } else {
    lst[seq(2, 2 * nrow(fdta), 2)] <- as.list(fdta$Sequence)
  }
  fwrite(tibble(unlist(lst)), file = out.file, compress = "auto", quote = F, col.names = F)
}

chop <- function(sequence, width){
  n <- str_length(sequence)
  s1 <- seq(1, n, width)
  s2 <- pmin(s1 + width - 1, n)
  return(str_sub(sequence, s1, s2))
}

MSA_files<-function(filenames_fasta,msa_Dir){
	#input as DNAstringset
	lfasta_set<- lapply(filenames_fasta, function(x) genWS(input=x))
	lfasta_length<- lapply(lfasta_set,  length)
	
	#lfasta_set<-lfasta_set[lengths(lfasta_set)>4]
	#lfasta_keepgrt<-which(lengths(lfasta_set)>4)
	#filenames_fasta<- filenames_fasta[lfasta_keepgrt]
	
	lfasta_set1<-lapply(lfasta_set, function(x) x %>% as.character()%>% gquad_main(.) %>% 
	dplyr::select(sequence)%>% unlist(.,use.names=F)%>% as.character())
	
	lfasta_set_unique<- lapply(lfasta_set1,function(x) length(unique(x)))%>% unlist()
	lfasta_keepgrt<-which(lfasta_set_unique>1)

	lfasta_set2<-lapply(lfasta_set1, DNAStringSet)
	
	lfasta_tosttager_index <- which(lfasta_set_unique>1)
	lfasta_asitis_index<- which(lfasta_set_unique<2)
	#filenames_fasta<- filenames_fasta#[lfasta_tosttager_index]
	nlen<-length(lfasta_set)

alignsequences<-function(lfasta_tosttager_index,lfasta_asitis_index,index){
	if (index %in% lfasta_tosttager_index ){
	alignedDNA <- AlignSeqs(lfasta_set2[[index]],iterations=20,refinements=15,gapOpening=-8,gapExtension=-10, verbose = FALSE)
	alignedDNA_stgger<-StaggerAlignment(alignedDNA, tree = NULL, threshold = 4, fullLength = TRUE, processors = 16, verbose = FALSE)
	alignedDNA_adjust<-AdjustAlignment(alignedDNA_stgger,substitutionMatrix=sub_mat,gapOpening=-6,gapExtension=-8)
	alignedDNA_gapsremoved <- del.colgapsonly(as.matrix( alignedDNA_adjust),threshold=0.98)
	alignedDNA_gapsremoved <- del.rowgapsonly(as.matrix(alignedDNA_gapsremoved),threshold=0.98)
	alignedDNA<-alignedDNA_gapsremoved  %>%as.list()%>% as.character() %>% lapply(.,paste0,collapse="") %>% unlist %>% DNAStringSet
	}
	else if (index %in% lfasta_asitis_index){
	alignedDNA <- AlignSeqs(lfasta_set[[index]],iterations=20,refinements=15,gapOpening=-8,gapExtension=-10, verbose = FALSE)
	}
return(alignedDNA)
}

	
	alignedDNA<-lapply(seq_along(1:nlen),function(x) alignsequences(lfasta_tosttager_index,lfasta_asitis_index,index=x))

	#alignedDNA <- lapply(lfasta_set, function(x) AlignSeqs(x,iterations=20,refinements=15,gapOpening=-8,gapExtension=-10, verbose = FALSE))
	#alignedDNA_stgger<-lapply(alignedDNA, function(x) StaggerAlignment(x, tree = NULL, threshold = 4, fullLength = TRUE, processors = 16, verbose = FALSE) )
	#alignedDNA_adjust<-purrr::map(alignedDNA_stgger,~AdjustAlignment(.x,substitutionMatrix=sub_mat,gapOpening=-6,gapExtension=-8))
	###alignedDNA_adjust<-purrr::map(alignedDNA_stgger,~AdjustAlignment(.x,substitutionMatrix=sub_mat))

	#alignedDNA_gapsremoved <- purrr::map(alignedDNA_adjust, ~ del.colgapsonly(as.matrix(.x),threshold=0.98))
	#alignedDNA_gapsremoved <- purrr::map(alignedDNA_gapsremoved , ~ del.rowgapsonly(as.matrix(.x),threshold=0.98))
	
	#toupper(unlist(as.character(alignedDNA_gapsremoved[[1]])))
	#
	dir.create(file.path(msa_Dir), showWarnings = FALSE,recursive=TRUE)
	file.path(msa_Dir)
	msa_Dir_setwd<-paste0(msa_Dir,"/","allplots/")
	dir.create(file.path(msa_Dir_setwd), showWarnings = FALSE,recursive=TRUE)
	
	setwd(msa_Dir_setwd)
	files<-basename(filenames_fasta)
	#outfiles<-purrr::map(files,~ paste0(msa_Dir,gsub(".fasta","",.x),"_al_t",times,".fasta"))
	outfiles<-purrr::map(files,~ paste0(msa_Dir,.x))

	# write the alignment to a new FASTA file
	purrr::map2(alignedDNA, outfiles,~writeXStringSet(.x, file=.y))
		
	alignedDNA1<-purrr::map(outfiles,~readFasta(.x))
	alignedDNA_gapsremoved1 <- purrr::map(alignedDNA1, ~msaTrim(.x,gap.end=0.95,gap.mid=1))
	
	purrr::map2(alignedDNA_gapsremoved1, outfiles,~writeFastausedhere(.x, .y))
	msa_stat_p<-"/bio/home/goonerrn/g_quadruplex/annotation_database/G_quad_clustering/starcode/codes/MSA_stat.sh"
	system(paste0(msa_stat_p,"  ",msa_Dir," ", times),intern=TRUE)
	msa_statplot_p<-"/bio/home/goonerrn/g_quadruplex/annotation_database/G_quad_clustering/starcode/codes/MSTAT_plot.R"
	source(msa_statplot_p)
	
}


genWS<-function(input) {
x<-readDNAStringSet(file=input, format="fasta")
return(x)
}
check=function(x) tryCatch(if(class(x) == 'logical') 1 else 1, error=function(e) 0) 



# Function
ReadFasta_profile<-function(file) {
   # Read the file line by line
   fasta<-readLines(file)
   # Identify header lines
   ind<-grep(">", fasta)
   # Identify the sequence lines
   s<-data.frame(ind=ind, from=ind+1, to=c((ind-1)[-1], length(fasta)))
   # Process sequence lines
   seqs<-rep(NA, length(ind))
   for(i in 1:length(ind)) {
      seqs[i]<-paste(fasta[s$from[i]:s$to[i]], collapse="")
   }
	filename=tools::file_path_sans_ext(basename(file))
   # Create a data frame 
   DF<-data.frame(id=gsub(">", "", as.character(fasta[ind])), sequence=as.character(seqs),cluster=filename)
   # Return the data frame as a result object from the function
   return(DF)
}


writeFasta1<-function(data,cluster_fasta_Dir){
  fastaLines = c()
  for (rowNum in 1:nrow(data)){
    fastaLines = c(fastaLines, as.character(paste(">", data[rowNum,"id"], sep = "")))
    fastaLines = c(fastaLines,as.character(data[rowNum,"sequence"]))
  }
	filename<-paste0(cluster_fasta_Dir,"/",as.character(data[1,"cluster"]),".fasta")

  fileConn<-file(filename)
  writeLines(fastaLines, fileConn)
  close(fileConn)
}



first_MSA<-function(filenames_fasta,msa_Dir){
	
	lfasta_set<- purrr::map(filenames_fasta, ~ genWS(input=.x))
	lfasta_length<- purrr::map(lfasta_set, ~ length(.x))
	
	lfasta_set<-lfasta_set[lengths(lfasta_set)>4]
	nlen<-length(lfasta_set)	
	alignedDNA <- purrr::map(lfasta_set, ~ AlignSeqs(.x,iterations=20,refinements=15,gapOpening=-3,gapExtension=-5, verbose = TRUE))
	alignedDNA_stgger<-purrr::map(alignedDNA,~StaggerAlignment(.x, tree = NULL, threshold = 3, fullLength = TRUE, processors = 16, verbose = TRUE))
	alignedDNA_adjust<-purrr::map(alignedDNA_stgger,~AdjustAlignment(.x,substitutionMatrix=sub_mat,gapOpening=-2,gapExtension=-2))
	#alignedDNA_adjust<-purrr::map(alignedDNA_stgger,~AdjustAlignment(.x,substitutionMatrix=sub_mat))
	
	alignedDNA_gapsremoved <- purrr::map(alignedDNA_adjust, ~ del.colgapsonly(as.matrix(.x),threshold=9))
	
        alignedDNA_gapsremoved <- purrr::map(alignedDNA_gapsremoved , ~ del.rowgapsonly(as.matrix(.x),threshold=1))
	
	library(magrittr)
	
	alignedDNA<-purrr::map(alignedDNA_gapsremoved ,~.x %>%as.list()%>% as.character() %>%lapply(.,paste0,collapse="") %>% unlist %>% DNAStringSet)
	dir.create(file.path(msa_Dir), showWarnings = FALSE,recursive=TRUE)
	file.path(msa_Dir)
	msa_Dir_setwd<-paste0(msa_Dir,"/","allplots/")
	dir.create(file.path(msa_Dir_setwd), showWarnings = FALSE,recursive=TRUE)
	
	setwd(msa_Dir_setwd)
	files<-basename(filenames_fasta)
	#outfiles<-purrr::map(files,~ paste0(msa_Dir,gsub(".fasta","",.x),"_al_t",times,".fasta"))
	outfiles<-purrr::map(files,~ paste0(msa_Dir,.x))

	# write the alignment to a new FASTA file
	purrr::map2(alignedDNA, outfiles,~writeXStringSet(.x, file=.y))
	
	msa_stat_p<-"/bio/home/goonerrn/g_quadruplex/annotation_database/G_quad_clustering/starcode/codes/MSA_stat.sh"
	system(paste0(msa_stat_p," ",msa_Dir," ", times),intern=FALSE)
	msa_statplot_p<-"/bio/home/goonerrn/g_quadruplex/annotation_database/G_quad_clustering/starcode/codes/MSTAT_plot.R"
	source(msa_statplot_p)
	return(alignedDNA)
}


	sth=1
 	while ((length(filelist_fasta)>0) && (length(files))>0 && (sth<100)){
	
	if (times==1){
	fasta_Dir<-paste0(output_loc,"/","fasta_",path_fasta,"_",times,"/")
	dir.create(file.path(fasta_Dir), showWarnings = FALSE,recursive=TRUE)

	
	final_add<<-fasta_Dir
	file.copy(from=filelist_fasta, to=fasta_Dir, 
				overwrite = TRUE, recursive = FALSE, 
				copy.mode = TRUE)
	
	filenames_fasta<-(list.files(path=fasta_Dir, pattern=NULL, all.files=FALSE,
	    full.names=TRUE))
	
	}


	msa_Dir<-paste0(output_loc,"/MSA_fasta_",times,"/")
	print(paste0("------------------------------------------CHECKPOINT 2 ---",times,"(TIMES) PASSED------------------------------------------------"))
	rate <- rate_backoff(pause_base = 0.1, pause_min = 0.005, max_times = 10)

#-------------------------------------------------
	insistent_MSA_files<-insistently(MSA_files,rate,quiet=F)
	invisible(insistent_MSA_files(filenames_fasta,msa_Dir))
	print(paste0("------------------------------------------CHECKPOINT 3 ---",times,"(TIMES) PASSED------------------------------------------------"))
	
#--------------------------------------------------------------------------
	filteredMSA<-fread(paste0(edit6_path,"6_editMSTAT_final.tsv"))	

	if (times==1){
		files_msa<-filteredMSA$basenames
		filenames_msa<-filteredMSA$file
		laterstore<-filteredMSA %>% filter(score< 0)%>% select(basenames)%>% unlist()
		laterstore_filenames<-filteredMSA %>% filter(score< 0)%>% select(file)%>% unlist()
		later2<-NULL
		AllG<-rDNAbin(seq(15, 30, by=1) , base.freq =c(0,0,1,0),prefix="ALL_G")
		AllGalign<-align(AllG,progressive=TRUE,lambda = 1,residues="DNA",qe=NULL,refine="BaumWelch",maxiter=20)
		AllGPHMM <- derivePHMM(AllGalign,inserts="inherited",  refine = "BaumWelch",progressive=TRUE)
		randomG<-rDNAbin(seq(15, 30, by=1) , base.freq =c(.25,.25,.25,.25),prefix="random_G")
		randomGalign<-align(randomG,progressive=TRUE,lambda = 1,residues="DNA",qe=NULL,refine="BaumWelch",maxiter=20)
		randomGPHMM <- derivePHMM(randomGalign,inserts="inherited",  refine = "BaumWelch",progressive=TRUE)

			filesmsa<-paste0(files_msa,".fasta")
			files_2<-paste0(final_add,"/",filesmsa)	
			filenames_fasta_filt<-paste0(final_add,"/",files_msa,".fasta")		
	}else if (times>1){
		laterstore<-filteredMSA %>% filter(score< 0)%>% select(basenames)%>% unlist()
		laterstore_filenames<-filteredMSA %>% filter(score < 0)%>% select(file)%>% unlist()

		files_msa<-filteredMSA$basenames
		filenames_msa<-filteredMSA$file
		filesmsa<-paste0(files_msa,".fasta")
		filenames_fasta_filt<-paste0(clustermergedfastaDir,"/",filesmsa) 
		later2<-NULL	
		 }
	
	

	print(files_msa[1])
	print(filenames_msa[1])
	print(filenames_fasta_filt[1])
	print(length(files_msa))
	print(length(filenames_msa))
	
	
print(paste0("------------------------------------------CHECKPOINT 4 ---",times,"(TIMES) PASSED------------------------------------------------"))

	#lmsass <- purrr::map(filenames_msa, ~ genWS(input=.x))
	#lmsa <- purrr::map(filenames_msa, ~ genWS(input=.x)%>% as.character()%>% strsplit(.,split=""))
	#lmsa1<-purrr::map(lmsa,~matrix(unlist(.x),nrow=length(.x)))
	print("lmsa to be  done")
	if (times<3){
	lmsa1 <- purrr::map(filenames_msa, ~ ape::read.dna(.x,format="fasta",as.matrix=TRUE))
	lPHMM <- purrr::map(lmsa1, ~aphid::derivePHMM(.x,seqweights=NULL,pseudocounts="Laplace",residues="DNA",refine="BaumWelch",progressive=TRUE,omit.endgaps=TRUE,maxiter=100,deltaLL=0.1))
	} else if (times>2){
	lmsa1 <- purrr::map(filenames_msa, ~ ape::read.dna(.x,format="fasta",as.matrix=TRUE))	
	lPHMM <- purrr::map(lmsa1, ~aphid::derivePHMM(.x,seqweights=NULL,pseudocounts="Laplace",residues="DNA",inserts="inherited",refine="BaumWelch",progressive=FALSE,omit.endgaps=TRUE,maxiter=100,deltaLL=0.1))
	}
	print("lphmm done")

	lfasta <- purrr::map(filenames_fasta_filt, ~ape::read.dna(.x,format="fasta",as.matrix=FALSE))
	nlen<-length(filenames_fasta_filt)

	namefiles<-filteredMSA$basenames
	#cutoff_calc=data.frame()
	score_withself=data.frame()
	#toupper(unlist(as.character(lfasta[[k]][value]),use.names=F))

	j=0
	score_here<-function(k){
		HMM_outerscore=list()
		fasta_score=numeric(length(lPHMM))
		 ll<-dna2char(lfasta[[k]])%>%unlist(.,use.names = F)%>% as.character()%>% gquad_main(.)
		sequences_quadparser<-split(ll$sequence,ll$input_ID)
		c<-data.frame()
		c<-future_map(seq_along(1:length(sequences_quadparser)),~ multiple_scores(sequences_quadparser[[.x]]))%>%
		rbindlist(idcol = "col")%>% group_by(col)%>% filter(logoddsbyakw == max(logoddsbyakw))%>% ungroup()%>%
 		as.data.frame()
		c[, "conf"] <-findInterval(c$newakw, c(0.3, 0.6))
		c[, "scoreconf"] <-findInterval(c$logoddsbyakw, c(0, 3,5))
		HMM_outerscore<-c
		return(HMM_outerscore)
}

family_number<-seq_along(1:length(lfasta))		
#HMM_outerscore1<-future_map(seq_along(1:length(lfasta)),~ dna2char(lfasta[[.x]])%>% unlist(.,use.names = F)%>% as.character()%>% gquad_main(.))

HMM_outerscore<-future_map(seq_along(1:length(lfasta)),~score_here(.x))
HMM_outerscore<- mapply(cbind, HMM_outerscore, "family_num"=family_number, SIMPLIFY=F)
names(HMM_outerscore)<-family_number

length_cls<-map(lfasta,~length(.x))

#tf<-future_map(HMM_outerscore,~as.data.frame.table(prop.table(table(.x$best_model==.x$family_num))))%>% rbindlist(idcol="fam")%>%
#as.data.frame()%>% reshape2::dcast(fam~Var1,value.val="Freq",drop=FALSE,fill=0)
#these are ones we protect but only merge with ones through their values
#tf[is.na(tf)]<-0

#tf1<- tf%>%   group_by(fam)%>% mutate(cond=`TRUE`> `FALSE`)  %>% as.data.frame()%>% filter(if(any(cond==TRUE)) (cond==TRUE &  `TRUE` > 0.7))

#fam<-tf1$fam%>% unlist()


#-------------------------------------------------

great_groups<-future_map(seq_along(1:length(HMM_outerscore)),~HMM_outerscore[[.x]]%>%
 group_by(best_model,family_num)%>%summarise(dbm=median(logoddsbyakw),n=n(),confscore=median(newakw))%>% mutate(per=prop.table(n)*100)%>%
filter((n>5) & (per>95) & (best_model==family_num) & (dbm>6) & (confscore > 0.8)))%>% rbindlist()

print("HIGHLIGHTED GROUPS")
print(great_groups)

if (dim(great_groups)[1] != 0) {
gg<-great_groups%>% select(best_model)%>% unlist(use.names=F)
} else {
gg<-NA
}

##----------------
if (times>0){
		laterstore<-unique(files_msa[gg])
		laterstore_filenames<-unique(filenames_msa[gg])
		#laterstore<-filteredMSA %>% filter(score< 0.10)%>% select(basenames)%>% unlist()
		#laterstore_filenames<-filteredMSA %>% filter(score< 0.10)%>% select(file)%>% unlist()
		
		files_msa<-filteredMSA$basenames[!filteredMSA$basenames %in% laterstore]
		filenames_msa<-filteredMSA$file[!filteredMSA$file %in% laterstore_filenames]
		
		print(paste0("FAMILIES SAVED:",length(!is.na(laterstore))))
		#print(length(!is.na(laterstore)))
		print("SELECTED CLUSTERS")
		print(head(laterstore))
		if(!exists("clustermergedfastaDir"))
		{
 		clustermergedfastaDir<<-paste0(edit6_path,"fasta_",path_fasta,"_",times,"/")
		}
		if (length(laterstore)>0){	
			later_filt<-laterstore
			#later2<-paste0(paste0(final_add,"/",laterstore,".fasta"))
			filesmsa<-paste0(files_msa,".fasta")
			files_1<-filesmsa
			files_2<-paste0(clustermergedfastaDir,"/",files_1) 
			later2<-paste0(clustermergedfastaDir,"/",laterstore,".fasta")
			filenames_fasta_filt<-files_2[!files_2  %in% later2]			
		} else if (length(laterstore)==0){
			filesmsa<-paste0(files_msa,".fasta")
			files_1<-filesmsa
			files_2<-paste0(clustermergedfastaDir,"/",files_1) 
			later2<-NA
			filenames_fasta_filt<-files_2[!files_2  %in% later2]
		 }
	}


#--------------------------------------------------------------------------	


get_groups<-function(HMM_outerscore,gg){

tomerge<-future_map(HMM_outerscore,~.x%>%group_by(best_model,family_num)%>% summarise(n=n(),dbm=mean(logoddsbyakw))%>% ungroup()%>%
filter(dbm<5 & n>1)%>% dplyr::select(family_num,best_model))%>% rbindlist()
#these are the good ones
dontmerge<-map(HMM_outerscore,~.x%>% group_by(best_model)%>% mutate(n=n(),dbm=mean(logoddsbyakw),tf=best_model==family_num)%>% ungroup%>%
filter((tf==TRUE) & (dbm < -15 ))%>%  nrow())%>% unlist()
names(dontmerge)<-family_number
dm<-dontmerge[dontmerge<1]
print("DONT MERGE ELEMENTS:")
print(dm)
 tomerge1<-tomerge[!(tomerge$family_num %in% as.numeric(names(dm)))]
 tomerge1<-tomerge1[!(tomerge1$best_model %in% as.numeric(names(dm)))]

 #tomerge1<-tomerge1[!(tomerge1$best_model %in% as.numeric(fam))]
 tomerge1<-tomerge1[!(tomerge1$best_model %in% as.numeric(gg))]
 tomerge1<-tomerge1[!(tomerge1$family_num %in% as.numeric(gg))]

groups<-tomerge1%>% filter(family_num!=best_model)
filtertemp3<- lapply(as.numeric(names(dm)),function(x) filenames_fasta_filt[x])%>% unlist()
return(list(as.data.frame(groups),filtertemp3))

}
#--------------------------------------------------------------------------
#--------------------------------------------------------------------------
	
possget_groups<-possibly(.f=get_groups,otherwise=list(data.frame(family_num=numeric(0),best_model=numeric(0)),NA))
#--------------------------------------------------------------------------

groups_sth<-possget_groups(HMM_outerscore,gg)
groups<-groups_sth[[1]]
filtertemp3<-groups_sth[[2]]
#here we collect filenames of those whose ttest is same
#these are the ones not to be merged because of their low scores on themselves


print(paste0("------------------------------------------CHECKPOINT 5 ---",times,"(TIMES) PASSED------------------------------------------------"))

	
	nrow_group<-nrow(as.data.frame(groups))
	#cutoff_calc: need to get this, forgot now how i was thinking to get it, but easy
	
	groups<-as.data.frame(groups)
 	print(groups)	
	groups_file=data.frame()
	
	print(length(filenames_fasta_filt))
	for( i in rownames(groups) ){
		file1_m<-as.numeric(groups[i, 'family_num'])
		print(file1_m)
		file2_m<-as.numeric(groups[i,"best_model"])
		print(file2_m)
		print(paste0(filenames_fasta_filt[file1_m],"    ",filenames_fasta_filt[file2_m]))
		groups_file<-rbind(groups_file,data.frame(cbind(V1=as.character(filenames_fasta_filt[file1_m]),V2=as.character(filenames_fasta_filt[file2_m]))))
	}

#----------#----------#----------#----------#----------#----------
	#here we collect filenames of those whose ttest is same
	
	
	if (nrow(groups_file)>0){
		filtertemp1<-groups_file[,c("V1","V2")]
		filtertemp1$V1<-basename(as.character(groups_file$V1))
		filtertemp1$V2<-basename(as.character(groups_file$V2))
		filtertemp1$V1[!(filtertemp1$V1 %in% filtertemp3)]
		filtertemp1$V2[!(filtertemp1$V2 %in% filtertemp3)]
		groups_file<-filtertemp1[!((filtertemp1$V1 %in% filtertemp3) | (filtertemp1$V2 %in% filtertemp3)) , ]
		mergedfastaDir<-paste0(path_to_use,"merged_fasta_",times,"/")
		dir.create(file.path(mergedfastaDir ), showWarnings = FALSE,recursive=TRUE)
		if (nrow(groups_file)>0){

			groups_filep<-groups_file
			groups_filep$V1<-paste0((dirname(filenames_fasta_filt))[1],"/",groups_file$V1)
			groups_filep$V2<-paste0(dirname(filenames_fasta_filt)[1],"/",groups_file$V2)
			temp_unlist<-groups_file %>%       mutate(across(where(is.factor), as.character))%>%
			unlist()%>% unique()
			net=graph_from_data_frame(groups_filep,directed=FALSE)
			net <- simplify(net, remove.multiple = T, remove.loops = T) 
			cl<-components(net)
						
			temp0<-lapply(seq_along(cl$csize)[cl$csize > 1], function(x) 
                                         V(net)$name[cl$membership %in% x])
			temp_unlist<-lapply(temp0,basename)%>%unlist() 
			temp<-map(temp0,~paste(.x, collapse= ' ')%>% as.data.frame())%>%bind_rows()%>%as.data.frame()
			colnames(temp)<-"V1"

			temp$V1<-as.character(temp$V1)
			#temp<-groups%>% group_by(V1)%>% mutate(V2=paste0(unique(V2),collapse=" "))%>% distinct()%>%as.data.frame()
			#
			
			#remaining fasta which are unmerged
			tempsth<-basename(filenames_fasta_filt)
			temp2<-basename(temp_unlist)
			fasta_unmerged<-filenames_fasta_filt[!(tempsth %in% temp2)]
			new_fasta_merged<-"merged_"
			
			
			#print(temp)
			for (i in 1:length(temp$V1)){
				out<-paste0(mergedfastaDir,new_fasta_merged,i,".fasta")
				filenames<-system(paste0('cat   ',  temp$V1[i] , '  >  ', out))
				print(head(filenames))
			}

	 } else {
		temp_unlist<-NA
		fasta_unmerged<-filenames_fasta_filt[!filenames_fasta_filt %in% temp_unlist]	
		}
		

	} else if (nrow(groups_file)==0){
		mergedfastaDir<-paste0(path_to_use,"merged_fasta_",times,"/")
		dir.create(file.path(mergedfastaDir ), showWarnings = FALSE,recursive=TRUE)
		temp_unlist<-NA
		fasta_unmerged<-filenames_fasta_filt[!filenames_fasta_filt %in% temp_unlist]
			
	}

		
	#if (check(later2) ==1 ){
	if (length(!is.na(later2))>0){
		final_fasta<-later2 
		#fasta_unmerged.1 <- fasta_unmerged[!fasta_unmerged %in% later2]

		finalfastaDIR<-paste0(output_loc,"/final_RESULT_fasta/")	
		if (!(file.exists(finalfastaDIR))){
			dir.create(finalfastaDIR, showWarnings = FALSE,recursive=TRUE)	
			}
	
			file.copy(from=final_fasta, to=finalfastaDIR, 
				overwrite = TRUE, recursive = FALSE, 
				copy.mode = TRUE)
			
			
			file.remove(final_fasta)


	} 

	file.copy(from=fasta_unmerged, to=mergedfastaDir ,
		overwrite = TRUE,	 recursive = FALSE, copy.mode = TRUE)

    	file.remove(fasta_unmerged)


print(paste0("------------------------------------------CHECKPOINT 6 ---",times,"(TIMES) PASSED------------------------------------------------"))


		
	filenames_fasta_temp<-list.files(mergedfastaDir,full.names=TRUE)
	files<-list.files(path=mergedfastaDir, pattern="*.fasta", all.files=FALSE,
		full.names=FALSE)
		
	
	source("/bio/home/goonerrn/g_quadruplex/annotation_database/G_quad_clustering/starcode/codes/recovery_graphcluster.R")
	pairwise_dist_folder<-paste0(path_to_use,"pairwise_dist_second_",times,"/")
	dir.create(file.path(pairwise_dist_folder), showWarnings = FALSE,recursive=TRUE)
	

	fasta_cluster<-purrr::map(filenames_fasta_temp,~ ReadFasta_profile(.x))
	fasta_cluster_Dnabin<-purrr::map(filenames_fasta_temp,~ ape::read.dna(.x,format="fasta",as.matrix=FALSE))

	n=length(fasta_cluster)
	
	invisible(final_result<-lapply(1:n,function(i){

	dat_frame<-as.data.frame(fasta_cluster[[i]])

	if (length(unique(unlist(as.character(dat_frame$sequence))))<2){
	res<-dat_frame
	#res$cluster<-as.character(dat_frame$cluster)
	#res$sequence<-as.character(dat_frame$sequence)
	res$agnesclusters2<-as.character(dat_frame$cluster)
	res$dianaclusters2<-res$agnesclusters2
	res$hclusters2<-res$agnesclusters2
	res[] <- lapply(res, as.character)
	result<-res%>% select(id,sequence,cluster, agnesclusters2,dianaclusters2,hclusters2)%>% as.data.frame()

	} else {
		if (nrow(dat_frame)<1000){
			result<-second_clustering(dat_frame,pairwise_dist_folder)
			result[]<-lapply(result,as.character)
			result<-as.data.frame(result)
		}
		if (nrow(dat_frame)>999){
			dat_frame1 <-(fasta_cluster_Dnabin[[i]])
			kmerasotu<-otu(dat_frame1 , k = 6, threshold = 0.6, method = "centroid", nstart = 20)
			res<-kmerasotu%>% as.data.frame()
			#we make these same as above for uniformity, otherwise only one required. to lazy to go back and write that function again. again, three types of clustering is provided, so added info.
			colnames(res)<-"agnesclusters2"
			res$dianaclusters2<-res$agnesclusters2
			res$hclusters2<-res$agnesclusters2
			res$id<-names(dat_frame1)
			res$cluster<-as.character(fasta_cluster[[i]]$cluster)
			res$sequence<-fasta_cluster[[i]]$sequence
			result<-res%>% select(id,sequence,cluster, agnesclusters2,dianaclusters2,hclusters2)%>% as.data.frame()
			}
	}
	return(result)  
	}))
#------------------------------
	times<<-times+1
	clustermergedfastaDir<<-paste0(edit6_path,"fasta_",path_fasta,"_",times,"/")
	#clustermergedfastaDir<-paste0(path_to_use,"clusteredaftermerged_fasta_",times,"/")
	dir.create(file.path(clustermergedfastaDir), showWarnings = FALSE,recursive=TRUE)
	
	library(purrr)
	


	result3<-lapply(1:length(final_result),function(y){
	
		results3b<-final_result[[y]]%>% as.data.frame() %>%  #group_by(hclusters2)%>%
		mutate(n=n())%>%mutate(cluster2=ifelse(n >= 1000,hclusters2,"s"))%>% mutate(cluster3=paste0(cluster,cluster2,"b"))%>%
		select(id,cluster3,sequence)
	
		result4<-results3b %>% bind_rows(.id="source")
		names(result4)[grepl("source", names(result4))] <- "osource"
	
		#print(head(result4))
		#print(paste0("NROW==",nrow(result4)))
		return(results3b)
	})



	#-----------------------------------------


	
	result4<-bind_rows(result3,.id="source")
	result4$cluster<-paste0(result4$cluster3)
	data1<-split(result4,result4$cluster)
	
	data2<-purrr::map(data1,~writeFasta1(.x,clustermergedfastaDir))
	

	print('completed')


	#print(MSA_Dir)
	print("hurray")
	
print(paste0("------------------------------------------CHECKPOINT 7 ---",times,"(TIMES) PASSED------------------------------------------------"))
	
	
	end <- Sys.time()
	difftime(end, start, units="secs")
	#if (times<6) {
	##first we filter files with less than 4 sequences in cluster
	stringa<-paste0("for f in ", clustermergedfastaDir,"*.fasta; do if [ $(grep -c '>'  $f | cut -d' ' -f1) -le 4 ]; then rm -f $f; fi; done")
	#} else {
	#stringa<-paste0("for f in ", clustermergedfastaDir,"*.fasta; do if [ $(grep -c '>'  $f | cut -d' ' -f1) -le 4 ]; then rm -f $f; fi; done")
	#}

	system(stringa)
	print(clustermergedfastaDir)
	
	filenames_fasta<<-(list.files(path=clustermergedfastaDir, pattern="*.fasta", all.files=FALSE,recursive=T,
		full.names=TRUE))
	
	files<<-list.files(path=clustermergedfastaDir, pattern="*.fasta", all.files=FALSE,recursive=T,
		full.names=FALSE)
	print(length(filenames_fasta))
	final_add<<-clustermergedfastaDir
	print(paste0("------------------------------------------CHECKPOINT 8 ---",times,"(TIMES) PASSED------------------------------------------------"))
	sth=sth+1
	print(paste0("LOOP STH",sth))
	print(paste0("number of files in fasta_files,times, ", times))
	rm(laterstore)
	rm(later2,fasta_cluster,dontmerge,gg,dm,groups,groups_sth,file2_m,file1_m)
	rm(laterstore_filenames,laterstore,lmsa1,lPHMM,lfasta,filenames_fasta_filt,filteredMSA,data2,data1,result4,mergedfastaDir,groups_file,filtertemp3)
	}






















