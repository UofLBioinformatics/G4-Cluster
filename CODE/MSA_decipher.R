#!/bin/bash

library(aphid)
library(dplyr)
library(ape)
library(purrr)
library(data.table)
require(Biostrings) 
library(DECIPHER)
library(aphid)
library(dplyr)
library(ape)
library(purrr)
library(data.table)
library(parallel)
#dump <- get(load("last.dump.rda"))


standard_path<-"/bio/home/goonerrn/g_quadruplex/annotation_database/G_quad_clustering/starcode/"

 setwd(standard_path)

path_fasta<<-"allfinalfa"
#path_fasta<<-"allfinalf"
#path_fasta<<-"total_again"
#path_fasta<<-"ggg25rc"
#path_fasta<<-"upd_60c"
#path_fasta<<-"upd_cdhit"
#path_fasta<<-"gggblat60rc_same"
#path_fasta<<-"gggblat60rc_copied"
#path_fasta<<-"DNAclust_cls"
#path_fasta<<-"upd_dnaclust"
#path_fasta<<-"gggmeshclust"
#path_fasta<<-"ggg_blat_c_redo_highergap"


 options(error = quote({dump.frames(); save.image(file = paste0(standard_path,path_fasta,"_","last.dump.rda"))}))


final_add<-paste0(standard_path,path_fasta,"/6_edit/fasta_",path_fasta,"/")
path_output<-paste0(standard_path,path_fasta,"/")
edit6_path<<-paste0(standard_path,path_fasta,"/6_edit/")



##first we filter files with less than 4 sequences in cluster
stringa<-paste0("for f in ", final_add,"*.fasta; do if [ $(grep -c '>'  $f | cut -d' ' -f1) -le 4 ]; then rm -f $f; fi; done")

system(stringa)

filelist_fasta<-(list.files(path=final_add, pattern=NULL, all.files=FALSE,
    full.names=TRUE))
#msa_path<-'/bio/home/goonerrn/g_quadruplex/annotation_database/G_quad_clustering/starcode/gquad_all_70g4score/6_edit/msa_out/'

files<-list.files(path=final_add, pattern=NULL, all.files=FALSE,
    full.names=FALSE)

filenames<-filelist_fasta

filelist_fasta<-filelist_fasta
files<-files
#lfasta <- purrr::map(filelist_fasta, ~ape::read.dna(.x,format="fasta"))

genWS<-function(input) {
x<-readDNAStringSet(file=input, format="fasta")
return(x)
}


sub_mat<-nucleotideSubstitutionMatrix(match = 15, mismatch = -2, baseOnly = TRUE, type = "DNA")
sub_mat['A','A']<-7
sub_mat['G','G']<-10
sub_mat['C','C']<-7
sub_mat['T','T']<-7
sub_mat['A','G']<- -15
sub_mat['C','G']<- -15
sub_mat['G','T']<- -15
sub_mat['A','C']<- 0
sub_mat['A','T']<-5
sub_mat['C','T']<- 0

sub_mat<-Matrix::forceSymmetric(sub_mat)
sub_mat<-as.matrix(sub_mat)


sub_mat<-nucleotideSubstitutionMatrix(match = 3, mismatch = -2, baseOnly = TRUE, type = "DNA")
sub_mat['A','A']<-7
sub_mat['G','G']<-10
sub_mat['C','C']<-7
sub_mat['T','T']<-7
sub_mat['A','G']<- -8
sub_mat['C','G']<- -6
sub_mat['G','T']<- -8
sub_mat['A','C']<- -8
sub_mat['A','T']<-4
sub_mat['C','T']<- -8

sub_mat<-Matrix::forceSymmetric(sub_mat)
sub_mat<-as.matrix(sub_mat)

output_loc<-dirname(final_add)
path_to_use<-edit6_path
#times=1

#source("/bio/home/goonerrn/g_quadruplex/annotation_database/G_quad_clustering/starcode/codes/merge_PHMM_ttest_try3.R")
#HMM_ttest_merge(filelist_fasta,sub_mat,files,path_to_use,output_loc)

nrow_group<-1
times<-1
start <- Sys.time()
print("------------------------------------------CHECKPOINT 1 PASSED------------------------------------------------")

source("/bio/home/goonerrn/g_quadruplex/annotation_database/G_quad_clustering/starcode/codes/while_HMM_new.R")
print("completed")
print("Thank you for your patience!")
print("Have a great Day!")
print("-------------------------------------------------------")