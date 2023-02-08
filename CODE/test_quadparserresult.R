#!/bin/bash
start <- Sys.time()

args <- commandArgs(trailingOnly = TRUE)
path_input<-as.character(args[1])

source("/bio/home/goonerrn/g_quadruplex/annotation_database/G_quad_clustering/starcode/codes/functions_g4predict.R")

library(aphid)
library(dplyr)
library(ape)
library(purrr)
library(data.table)
require(Biostrings) 
library(furrr)
library(stringr)
library(stringi)
library(DECIPHER)
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#path_HMM<-"/bio/home/goonerrn/g_quadruplex/annotation_database/G_quad_clustering/starcode/gquad_all_70g4score/HMM_out/"
path_fasta<-"/bio/home/goonerrn/g_quadruplex/annotation_database/G_quad_clustering/starcode/allfinalfa/6_edit/final_RESULT_fasta/"
path_msa<-"/bio/home/goonerrn/g_quadruplex/annotation_database/G_quad_clustering/starcode/allfinalfa/6_edit/FINAL_MSA/"


filenames_fasta<-list.files(path=path_msa, pattern=".fasta", all.files=FALSE,   full.names=TRUE)
names<-list.files(path=path_msa, pattern=".fasta", all.files=FALSE,   full.names=F)

#filenames_fasta<-list.files(path=path_fasta, pattern=NULL, all.files=FALSE,    full.names=TRUE)

#names<-paste0("cluster_",seq_along(1:length(filenames_fasta)))


#lPHMM <- lapply(filenames, readPHMM)
#lPHMM<-lPHMM

family_lfasta <- lapply(filenames_fasta, function(x) ape::read.dna(x,format="fasta",as.matrix=TRUE))
lPHMM <- lapply(family_lfasta, function(x) aphid::derivePHMM(x, residues = "DNA",inserts="inherited",refine="BaumWelch",progressive=TRUE))	
len<-length(lPHMM )



if ((path_input)=="load"){
load("/bio/home/goonerrn/g_quadruplex/annotation_database/epimap/additional_bed/G4_expr_all_g4.RData")
} else {
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#possible inputs
#this is based on experimental G4
#gquadinput<-fread("/bio/home/goonerrn/g_quadruplex/annotation_database/epimap/additional_bed/output_header",header=FALSE)
gquadinput<-fread("/bio/home/goonerrn/g_quadruplex/annotation_database/epimap/additional_bed/GSM3003539_G4_quadparser.tsv",header=FALSE)

gquadinput<-gquadinput	%>% rowwise() %>%mutate(Ccontent=lengths(regmatches(V4, gregexpr("c|C", V4))))
gquadinput<-gquadinput	%>% rowwise()%>% mutate(Ggontent=lengths(regmatches(V4, gregexpr("g|G", V4))))


gquadinput<-setDT(gquadinput)
gquadinput[, V2rc := ifelse((Ccontent> Ggontent) , sapply( V4, function(x)
 as.character(reverseComplement(DNAString(x)))), sapply(V4,function(x) as.character(x)))]

gquadinput[, strand := ifelse((Ccontent> Ggontent) , sapply( V4, function(x) as.character(
-1)), sapply(V4,function(x) as.character(1)))]

ll<-gquadinput$V2rc%>%unlist(.,use.names = F)%>% as.character()%>% 
  gquad_main(.)
gquadinput$input_ID<-seq_along(1:nrow(gquadinput))

}
sequences<-split(as.character(ll$sequence),ll$input_ID)

time1<-Sys.time()
a<-map(seq_along(1:length(sequences)),
              ~multiple_scores(sequences[[.x]]))%>% rbindlist(idcol = "id")
time2<-Sys.time()
print(paste0("::::",time2-time1))
paste0("total time for scores of sqeuences length", length(sequences),":           ",
round(as.numeric(difftime(time1 = time2, time2 = time1, units = "secs")), 3), " Seconds")
a1<-cbind(ll,a)
a2<-a1%>%  filter(logoddsbyakw>3)%>% as.data.frame()
a3<-left_join(gquadinput,a2,by="input_ID")#%>% select(-c(V4,sequence_position)) %>%group_by(input_ID)%>% top_n(1,logoddsbyakw)
a3%>% fwrite(.,"Experimental_family_predictions.tsv",sep = "\t",quote = F,row.names = F)

 a4<-a3%>% filter(!is.na(best_model))%>% filter(logoddsbyakw > 5)

 a4$family<-names[a$best_model]
a5<-a4%>% group_by(input_ID)%>% top_n(1,logoddsbyakw)%>% as.data.frame()
a5%>% fwrite(.,"PREDICTED_TO_USE.tsv",sep="\t",sep="\t",quote=F,row.names=F)
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
multiple_scores<-function(sequences_list){
  ss<-future_map(sequences_list,~
      as.DNAbin.character(str_split(.x,"")%>%
                      unlist()))
result<-map(ss,~get_scores(.x,qlen=1)%>% as.data.frame)%>% rbindlist(idcol="source")
return(result)
}

#----------------------------------------------------------------------------------------

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




#filenames3<-"/bio/home/goonerrn/g_quadruplex/annotation_database/G_quad_clustering/starcode/pu27.fa"

#filenames3<-"/bio/home/goonerrn/g_quadruplex/annotation_database/query/g4_10added.fasta"
#filenames3<-"/bio/home/goonerrn/g_quadruplex/random/fasta files generation/FinalBlast/fullfastafile_known_duplicates.fasta"
#filenames3<-"/bio/home/goonerrn/g_quadruplex/annotation_database/G_quad_clustering/starcode/final_mixed/GGGTGGGTGGGTGGG.fa"

#qfasta <- ape::read.dna(filenames3,format="fasta")


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
