#!/usr/bin/env Rscript

#SBATCH --mem=20g

##there should be two arguments one character name of the folder to be created and partition of the train test ratio.

#./cluster_R.R "file_name" 70


library(stringr)
require(tidyr)
library(dplyr)
library(Biostrings)
#library(tidyverse)
library(stringdist)
library(igraph)


library(BiocParallel)


library(fastcluster)

library(cluster)
library(furrr)
library(future)
library(stringr)
library(purrr)

library(data.table)


##this is for output_fwd
source("/bio/home/goonerrn/g_quadruplex/annotation_database/G_quad_clustering/clustering_output/cluster_functions.R")
##this and SubDir should be input

args <- commandArgs(trailingOnly = TRUE)

RNAfold<-fread("/bio/home/goonerrn/g_quadruplex/annotation_database/redoing_reproduce_results/RNAFOLD_allG4fasta_var_features.txt",header=TRUE)
chromosome=paste0("chr",c(1:22,"X","Y"))

sapply(RNAfold, function(x) sum(is.na(x)))

RNAfold.1<-do.call(data.frame,lapply(RNAfold, function(x) replace(x, is.infinite(x),NA)))

RNAfold.1<-RNAfold.1[!duplicated(RNAfold.1$name), ]

numeric_var <- names(RNAfold.1)[sapply(RNAfold.1, is.numeric)]
setDT(RNAfold.1)[, summary(.SD), .SDcols = numeric_var]
RNAfold.2<-RNAfold.1[,c("name","GC","GGr","AGr","GAr","GCr","GTr","TGr","CGr","AT.GCr","MFEadj.GC","DiffMFE.EFE","EFEadj","ED","MEA","MFE")]

total.start.time <- Sys.time()

added_name_for_output_file<-args[1]
partition_per<-args[2]
partition_per<-as.numeric(partition_per)/100

mainDir<-"/bio/home/goonerrn/g_quadruplex/annotation_database/G_quad_clustering/clustering_output/"
SubDir <- paste0(mainDir,added_name_for_output_file,"/")

model_Dir <- paste0(SubDir,"models_",added_name_for_output_file,"/")

mat_Dir<-paste0(SubDir,"G4dist_",added_name_for_output_file,"/")
path_out_mat<-mat_Dir
dir.create(file.path( SubDir ), showWarnings = FALSE,recursive=TRUE)
dir.create(file.path( model_Dir ), showWarnings = FALSE,recursive=TRUE)
dir.create(file.path( mat_Dir ), showWarnings = FALSE,recursive=TRUE)

setwd(file.path( SubDir ))

setwd(SubDir )
message("Filtering Sequences...")

output_fwd<-fread("/bio/home/goonerrn/g_quadruplex/annotation_database/query/output_fwd_withCHR.txt",header=FALSE)

chromosome=paste0("chr",c(1:22,"X","Y"))

output_fwd<-output_fwd[output_fwd$V1 %in% chromosome]
output_fwd$name<-paste0(output_fwd$V1,":",output_fwd$V2)
output_fwd<-output_fwd%>%separate(V3,c("nggg","nquad","fquad"),":")
output_fwd$nggg<-as.numeric(output_fwd$nggg)
output_fwd$grp1<-ifelse(as.numeric(output_fwd$nggg>7),1,as.numeric(output_fwd$nggg))
output_fwd$grp1<-as.factor(output_fwd$grp1)

output_fwd.1<-output_fwd%>% filter(nggg==4 | nggg==5 | nggg==6)
colnames(output_fwd.1)[6]<-"sequence"

sequence.1<-gsub('(^[G])\\1+','\\1',output_fwd.1$sequence)
output_fwd.1$seq<-gsub('([G])\\1+$','\\1',sequence.1)
output_fwd.1$length<-nchar(output_fwd.1$seq)

output_fwd.1<-merge(output_fwd.1,RNAfold.2,by="name")

output_fwd.1<-output_fwd.1 %>% 
    group_by(sequence) %>% 
    mutate(ID = paste(unique(name), collapse = '|'),nggg=mean(nggg),length=mean(length))%>%distinct()

output_fwd.1$name <- seq.int(nrow(output_fwd.1))
output_fwd.1$name<-paste0("g4_",output_fwd.1$name)
#writeFasta(output_fwd.1[,c("name","seq")],"4_g4_lgr7.fasta")

output_fwd.1<-output_fwd.1%>% filter(length > 7)

##load g4_mat,test,train
# load("/bio/home/goonerrn/g_quadruplex/annotation_database/G_quad_clustering/g4_mat_test_train.RData")


set.seed(101) # Set Seed so that same sample can be reproduced in future also
# Now Selecting 75% of data as sample from total 'n' rows of the data  
sample <- sample.int(n = nrow(output_fwd.1), size = floor(partition_per*nrow( output_fwd.1)), replace = F)

train <-  as.data.frame(output_fwd.1[sample, ])
test  <-  output_fwd.1[-sample, ]


#name<-x2$name
g4_mat<-g4_df(as.data.frame(train))
 


if (!is.null(g4_mat$l3f_V5)) {
temp<-prop.table(table(is.na(g4_mat$l3f_V5)))%>% as.data.frame()
result<-temp%>% filter(Var1=="FALSE") %>% summarise(n=ifelse(Freq>0.2,1,0))
if (result==0){
g4_mat$l3f_V5=NULL}
}


if (!is.null(g4_mat$l3f_V6)) {
temp<-prop.table(table(is.na(g4_mat$l3f_V6)))%>% as.data.frame()
result<-temp%>% filter(Var1=="FALSE") %>% summarise(n=ifelse(Freq>0.2,1,0))
if (result==0){
g4_mat$l3f_V6=NULL}
}

if (!is.null(g4_mat$l3f_V7)) {
temp<-prop.table(table(is.na(g4_mat$l3f_V7)))%>% as.data.frame()
result<-temp%>% filter(Var1=="FALSE") %>% summarise(n=ifelse(Freq>0.2,1,0))
if (result==0){
g4_mat$l3f_V8=NULL}
}



g4_mat<-cbind(as.data.frame(train[,c("name","sequence","nggg","GC","GGr","AGr","GAr","GCr","GTr","TGr","CGr","AT.GCr","MFEadj.GC","DiffMFE.EFE","ED","MEA")]),g4_mat)
colnames(g4_mat)[2]<-"sequence"
g4_mat[is.na(g4_mat)] <- 0
g4_mat$length<-nchar(g4_mat$sequence)
g4_mat<-g4_mat%>% filter(length > 7)
#using k=6 for nggg==5 | nggg==6
kmers<-kmer_from_granges(g4_mat[,c("name","sequence")],k=3)
#kmers<-kmer_from_granges(g4_mat[,c("name","sequence")],k=6)  #20% data
#kmers<-kmer_from_granges(g4_mat[,c("name","sequence")],k=6)  #nggg=5 | 6

#kmers<-kmer_from_granges(g4_mat[,c("name","sequence")],k=7)


g4_mat.1<-cbind(g4_mat,kmers)
nc<-ncol(g4_mat.1)
rownames(g4_mat.1)<-g4_mat.1$name
g4_mat.1$nggg<-as.factor(g4_mat.1$nggg)
cols<-colnames(g4_mat.1)[-c(1,2,3)]

g4_mat.2<-g4_mat.1%>% mutate_at(cols,~(scale(.,center=TRUE)%>% as.vector))

g4_mat.3<-g4_mat.2[,-c(1,2,3)]

ind<- g4_mat.3[, colSums(g4_mat.3 != 0) > 0]
ind <- sapply(g4_mat.3, function(x) sum(x==0)) != nrow(g4_mat.3)

g4_mat.4<-as.data.frame(g4_mat.3)[,ind]
g4_mat.2$id<-rownames(g4_mat.4)

file1=paste0(added_name_for_output_file,"_initial_files",".RData")
save(g4_mat.4,g4_mat.2,train,test,file=file1)

 dat%>% as.table() %>% as.data.frame() %>%     
   subset(Var1 != Var2 & abs(Freq)>0.5) %>% 
   filter(!duplicated(paste0(pmax(as.character(Var1), as.character(Var2)), 
pmin(as.character(Var1),  as.character(Var2)))) %>% arrange(desc(Freq))


##for k=5, method=centroid, metric with euclidean 

#k=7 running with train

###currently, k=6 is running with 20% of all datasets. H
##currently,we are trying this,clara with k=6, for all train.
# Use map_dbl to run many models with varying value of k
if (nrow(g4_mat.4)> 3000) {
ncls<-3
fcls<-ncls + 15

}else if (nrow(g4_mat.4)> 10000) {
ncls<-4
fcls<-ncls + 10 
} else {
ncls<-3
fcls<-ncls + 30}


sil_width <- future_map_dbl(ncls:fcls,  function(k){

start.time <- Sys.time()

model<-clara(g4_mat.4,k,metric="euclidean",samples=500)

save(model,file=paste0(model_Dir,"clara_ml_",k,".RData"))
return(model$silinfo$avg.width)
end.time <- Sys.time()
time.taken <- end.time - start.time
print(time.taken)
})

##run CLARA at the beginning:

# Generate a data frame containing both k and sil_width
sil_df <- data.frame(
  k = ncls:fcls,
  sil_width = sil_width
)


optimal_cluster<-sil_df[which.max((sil_df$sil_width)),]$k
print_seq<-"the optimal cluster was found to be"
print_seq2<-"for cluster number"

model_name<-paste0(model_Dir,"clara_ml_",optimal_cluster,".RData")
load(model_name)
print(paste0("the model used is:",optimal_cluster))
#model<-clara(g4_mat.4,k = optimal_cluster,samples=2000)
model$clustering<-paste0("a_",model$clustering)

clus_1<-cbind(as.data.frame(g4_mat.2[,c("name","sequence")]),as.data.frame(model$clustering))
index<-which(colnames(clus_1)=="model$clustering")
colnames(clus_1)[index]<-"cluster"


assign("clus_1", clus_1, envir = .GlobalEnv)


new_df<-clus_1%>%group_by(cluster)%>%
summarise(Freq=n(),.groups = 'drop')

new_df<-as.data.frame(new_df)
a<-as.data.frame(table(new_df$cluster))
print("STEP 1: total should be equal to above as clus_1")

if (sum(new_df$Freq)==nrow(clus_1)){
  y<-paste("and total for clusters is ",sum(new_df$Freq))
z<-("After 1st cluster number of rows here  in clus_1 ")
print(paste(z,"is",nrow(clus_1),y)) 
print("TRUE step1")
}else {
    print ("FALSE step 1")}




sub_cluster_large<-as.character(new_df[new_df[,"Freq"] > 1000,]$cluster)

final_sub_cluster<-as.character(new_df[new_df[,"Freq"]<=40,]$cluster)

sub_cluster<-as.character(new_df[(new_df[,"Freq"] > 40) & (new_df[,"Freq"] <= 1000),]$cluster)



gc()

if (length(sub_cluster_large)>0) {

#clus_level_large <- future_map(list(sub_cluster_large),~level_3(.x))%>% bind_rows()


#clus_level_large <- mapply(level_3,sub_cluster_large,SIMPLIFY = FALSE) %>% unlist(recursive = FALSE) %>%    bind_rows()
clus_level_large <-mclapply(list(sub_cluster_large),FUN=function(X)level_3(X),mc.cores=8, mc.set.seed = TRUE) %>%   bind_rows()%>%
 as.data.frame()

colnames(clus_level_large)<-c("name","sequence","cluster")
clus_1.1<-merge(clus_1[,c("cluster","sequence","name")],clus_level_large,by="name",all.x=TRUE)
clus_1.1$cluster<-ifelse(is.na(clus_1.1$cluster.y)==TRUE,as.character(clus_1.1$cluster.x),as.character(clus_1.1$cluster.y))

drops <- c("cluster.x", "sequence.y","cluster.y")
clus_1.1<-clus_1.1[,!(names(clus_1.1) %in% drops)]
print(colnames(clus_1.1))
#clus_1.1$sequence.y<-NULL
#clus_1.1$cluster.x<-NULL
clus_1<-clus_1.1
colnames(clus_1)<-c("name","sequence","cluster")

#assign("clus_1", clus_1, envir = .GlobalEnv)

new_df.1<-clus_1.1 %>% group_by(cluster)%>%
summarise(Freq=n(),.groups = 'drop')

new_df.1<-as.data.frame(new_df.1)
new_df<-new_df.1

sub_cluster_large<-as.character(new_df[new_df[,"Freq"] > 1000,]$cluster)

final_sub_cluster<-as.character(new_df[new_df[,"Freq"]<=40,]$cluster)

sub_cluster<-as.character(new_df[(new_df[,"Freq"] > 40)  & (new_df[,"Freq"] <= 1000 ),]$cluster)


if (length(sub_cluster_large)>0) {
clus_level_large <-mclapply(list(sub_cluster_large),FUN=function(X)level_3(X),mc.cores=8, mc.set.seed = TRUE) %>%   bind_rows()%>%
 as.data.frame()
print("breaking large clusters")
colnames(clus_level_large)<-c("name","sequence","cluster")
clus_1.1<-merge(clus_1[,c("cluster","sequence","name")],clus_level_large,by="name",all.x=TRUE)
clus_1.1$cluster<-ifelse(is.na(clus_1.1$cluster.y)==TRUE,as.character(clus_1.1$cluster.x),as.character(clus_1.1$cluster.y))

drops <- c("cluster.x", "sequence.y","cluster.y")
clus_1.1<-clus_1.1[,!(names(clus_1.1) %in% drops)]
print(colnames(clus_1.1))
clus_1<-clus_1.1
colnames(clus_1)<-c("name","sequence","cluster")
#assign("clus_1", clus_1, envir = .GlobalEnv)

new_df.1<-clus_1.1 %>% group_by(cluster)%>%
summarise(Freq=n(),.groups = 'drop')

new_df.1<-as.data.frame(new_df.1)
new_df<-new_df.1

sub_cluster<-as.character(new_df[new_df[,"Freq"] > 40,]$cluster)

final_sub_cluster<-as.character(new_df[new_df[,"Freq"]<=40,]$cluster)
warning<-as.character(new_df[new_df[,"Freq"] > 2000,]$cluster)
if((length(sub_cluster_large)>0))
{
print(paste0("warning",":",length(warning)))
}

}
}

gc()

if (length(sub_cluster)>0) {

#clus_level_2 <- future_map(list(sub_cluster),~level_2(.x))%>% bind_rows()
#merge(

clus_level_2 <-mclapply(list(sub_cluster),FUN=function(X)level_2(X),mc.cores=4) %>%  bind_rows()
#clus_level_2 <-(mapply(level_2,sub_cluster,SIMPLIFY = FALSE)) %>%   bind_rows()
}

if (length(final_sub_cluster)>0){
#clus_level_1 <-mcmapply(final_sub_cluster,FUN=function(X)level_1(X),mc.cores=8) %>%unlist(recursive = FALSE) %>%    bind_rows()
clus_level_1<-(mapply(level_1,final_sub_cluster,SIMPLIFY = FALSE))  %>%  bind_rows()
}


if ((exists("clus_level_1")==TRUE & exists("clus_level_2")==FALSE)){
  final_df<-clus_level_1
	clus_2<-final_df
  z<-("number of rows here  in clus_1 ")
  print(paste(z,"is",nrow(clus_2),"output row is" ,(nrow(final_df))))
  
  #rint(colnames(clus_level_1))

}else if ((exists("clus_level_2")==TRUE & exists("clus_level_1")==FALSE)){
  final_df<-clus_level_2
	clus_2<-final_df
    print("Clus_level_2 number of rows at the end")
      z<-("number of rows here  in clus_2 ")
  print(paste(z,"is",nrow(clus_2),"output row is" ,(nrow(final_df))))  #print(colnames(clus_level_2))
}else if (exists("clus_level_1")==TRUE & exists("clus_level_2")==TRUE){
  clus_level_2<-clus_level_2[names(clus_level_1)]
final_df<-rbind(clus_level_1,clus_level_2)
clus_2<-final_df
print(paste("level 1 is",nrow(clus_level_1),"level 2 is",nrow(clus_level_2)))
#print(colnames(clus_level_1))
print("Now final step")
#print(paste(z,"is",nrow(clus_2),"output row is" ,(nrow(final_df))))
#print(colnames(clus_level_2))}
}

if (exists("final_df")){
rm(clus_2)}
if (exists("clus_level_2")){
rm(clus_level_2)}
if (exists("clus_level_1")){
rm(clus_level_2)}


total.end.time <- Sys.time()
time.taken_total<-total.end.time-total.start.time
time.print<-paste("TOTAL TIME TAKEN IS\n",time.taken_total)
print(time.print)
print(time.taken_total)
name<-paste0(SubDir,"/",added_name_for_output_file,".tsv")
fwrite(final_df,sep="\t",file=name)
