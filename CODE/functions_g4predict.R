
if.is.empty <- function(x){
  is.null(need(x, message = FALSE))}
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  checkfasta <- function(fasta){
    
    # empty input
    if(length(fasta) == 0){
      return("NO INPUT")
    }
    
    # starts with > but not a sequence file (not containing A C G T N)
    if(startsWith(fasta[1], ">")){
      
      comp = unique(unlist(strsplit(paste(fasta[2:length(fasta)], collapse = ""), "")))
      all_symbols = c("A", "T", "C", "G", "N")
      if(sum(comp %in% all_symbols) != length(comp)){
        return(FALSE)
      } else {
        return(NULL)
      }
    } 
    
    # not starts with >, but it is a sequence file (contains only ACGTN)
    if(startsWith(fasta[1], ">") == F){
      comp = unique(unlist(strsplit(paste(fasta, collapse = ""), "")))
      all_symbols = c("A", "T", "C", "G", "N")
      if(sum(comp %in% all_symbols) != length(comp)){
        return(FALSE)
      } else {
        return(NULL)
      }
    }
    
  }
  
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~


checkg4<- function(a){
 	a4<-"(?=(G{3,7}?[A|C|G|T|U|N]{1,10}?G{3,7}?[A|C|G|T|U|N]{1,10}?G{3,7}?[A|C|G|T|U|N]{1,10}?G{2,7}))"
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


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ 
 
outputEval <- function(fa){
    
    fasta = readLines(fa)
#     print(paste("check results:", checkFasta(fasta)))
    
    # if it is null, pass the test
    if(length(checkFasta(fasta)) == 0){
      return(NULL)
    }
    
    # if it is "NOINPUT", FALSE will cause validate() fail silently
    if(checkFasta(fasta) == "NOINPUT"){
      return(FALSE)
    }
    
    # if it is not a right format, return character strings
    if(checkFasta(fasta) == FALSE){
      return("Please input valid fasta format")
    }
    
  }
  
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ 

multiple_scores<-function(sequences_list){
  ss<-future_map(sequences_list,~
      as.DNAbin.character(str_split(.x,"")%>%
                      unlist()))
result<-map(ss,~get_scores(.x,qlen=1)%>% as.data.frame)%>% rbindlist(idcol="source")
return(result)
}





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


#------------------------------------------------------

#------------------------------------------------------#------------------------------------------------------


classify_g4<-function(sequence,lPHMMtemp){

qfasta <- (sequence)
#print(sequence)
qlen<-length(qfasta)
qfasta<-as.DNAbin(unlist(strsplit(as.character(qfasta,use.names=F),"")))

#cutoff_calc=data.frame()
score_withself=data.frame()
HMM_score<-list()
HMM_index<-list()
HMM_append<-list()
j=0

#toupper((unlist(as.character(qfasta[k])),collapse=""))
# toupper(unlist(as.character(qfasta[1]),use.names=F))


for (k in (1:qlen)){
	fasta_score=list()
	threshold_score<-list()
	if.is.empty(lPHMM)
	for (count in seq_along(1:length(lPHMM))){
		print(paste0("count:",count,"      ","k:",k))
#	  print(lPHMMtemp[[count]])
#	  print(as.matrix(qfasta[k]))
		cscore<-aphid::forward(lPHMM[[count]],strsplit(sequence,""),odds=TRUE,windowspace="all",refine="BaumWelch")
		fasta_score<-c(fasta_score,cscore$score)
		}
		#print(fasta_score)
	#print(paste(round(unlist(fasta_score,3)),collapse=","))	
	fasta_score<-lapply(fasta_score, function(x) replace(x, !is.finite(x), 0))
	max_value2=max(unlist(fasta_score))	
	print(unlist(fasta_score))
	#setDT(fasta_score)
	#setorder(fasta_score,-x)
	#n_top_indice<-d$indice[1:3]
	index <- which((fasta_score)==max(unlist(fasta_score)))
	#print(max_value2)
	HMM_score<-c(HMM_score,list(max_value2))
	HMM_index<-c(HMM_index,list(index))
	HMM_append<-c(HMM_append,paste(max_value2,index,k,collapse=";"))
}


HMM_index1<-map(HMM_index,~as.data.frame(.x))
tempindex<-HMM_index1 %>% bind_rows(.id="source")



HMM_score1<-map(HMM_score,~as.data.frame(.x))%>% bind_rows(.id="source")
temp<-cbind(tempindex,HMM_score1)
return(temp)

}




#------------------------------------------------------
get_scores_multi<-function(qfasta,qlen){
HMM_append<-list()
for (k in (1:qlen)){
	fasta_score=numeric(length(lPHMM))
	sc <- numeric(length(lPHMM)) # scores (log probabilities)
	threshold_score<- numeric(length(lPHMM_filt)) 

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