

#!/bin/bash

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
library(stringr)
library(purrr)

library(data.table)
#library(tidyverse)
library(stringdist)

source("/bio/home/goonerrn/g_quadruplex/annotation_database/G_quad_clustering/starcode/codes/starcode_functions_cdhit.R")
source("/bio/home/goonerrn/g_quadruplex/annotation_database/G_quad_clustering/clustering_output/cluster_functions.R")
source("/bio/home/goonerrn/g_quadruplex/annotation_database/G_quad_clustering/starcode/codes/cdhit_parser.R")
source("/bio/home/goonerrn/g_quadruplex/annotation_database/G_quad_clustering/starcode/codes/G4regex_oneintomultiple.R")

args <- commandArgs(trailingOnly = TRUE)

total.start.time <- Sys.time()

#metrictouse="G4score"
metrictouse="evalue"

added_name_for_output_file<-args[1]
partition_per<-1
partition_per<-as.numeric(partition_per)/100




output_fwd<-fread("/bio/home/goonerrn/g_quadruplex/annotation_database/query/output_fwd_withCHR.txt",header=FALSE)

chromosome=paste0("chr",c(1:22,"X","Y"))

output_fwd<-output_fwd[output_fwd$V1 %in% chromosome]
output_fwd$name<-paste0(output_fwd$V1,":",output_fwd$V2)
output_fwd<-output_fwd%>%separate(V3,c("nggg","nquad","fquad"),":")
output_fwd$nggg<-as.numeric(as.character(output_fwd$nggg))


output_fwd_g4<-rbindlist(apply(output_fwd, 1, function(x) DNAregexfinder(seq=x["V4"],name=x["name"])))

output_fwd.1<-merge(output_fwd_g4,output_fwd)

colnames(output_fwd.1)[6]<-"sequence"
output_fwd.1$length<-nchar(output_fwd.1$sequence)

nggg<-"all"


o<-output_fwd.1%>% filter(sequence!=0)%>% group_by(sequence)%>% 
summarise(name=paste0(unique(one_in_many),collapse="|"),name_sub=paste0(unique(one_in_many),collapse="|"), ggg=paste0(unique(nggg),collapse="|"),
length_trimmed=nchar(sequence)) %>% filter(length_trimmed>5)%>% distinct()%>% as.data.frame()
min_len<-min(o$length_trimmed)

o$id<-paste0("g4_",rownames(o))

#pairwise_dist_folder<-mat_Dir

 output_fwd_g4%>% filter(sequence!=0)%>%
group_by(sequence)%>% mutate(n=n())%>% arrange(desc(n))%>%  
 select(sequence, n)%>%
 distinct()%>% head(20)


mainDir<<-"/bio/home/goonerrn/g_quadruplex/annotation_database/G_quad_clustering/starcode/"
SubDir0 <<- paste0(mainDir,added_name_for_output_file,"/")

SubDir <<- paste0(mainDir,added_name_for_output_file,"/",nggg,"/")
#print(SubDir)
model_Dir <<- paste0(SubDir,"models_",added_name_for_output_file,"/")
cluster_fasta_Dir<<-paste0(SubDir,"fasta_",added_name_for_output_file,"/")


mat_Dir<<-paste0(SubDir,"G4dist_",added_name_for_output_file,"/")
path_out_mat<-mat_Dir
dir.create(file.path( SubDir0 ), showWarnings = FALSE,recursive=TRUE)
dir.create(file.path( SubDir ), showWarnings = FALSE,recursive=TRUE)
dir.create(file.path( model_Dir ), showWarnings = FALSE,recursive=TRUE)
dir.create(file.path( mat_Dir ), showWarnings = FALSE,recursive=TRUE)
dir.create(file.path( cluster_fasta_Dir ), showWarnings = FALSE,recursive=TRUE)

pairwise_Dir <<- paste0(SubDir,"SECOND_cluster_",added_name_for_output_file,"/")
pairwise_dist_folder <<- paste0(pairwise_Dir,"SECOND_cluster","/")



setwd(file.path( SubDir ))
#writeFasta(o[,c("id","sequence")],"g4_notduplicated.fasta")
print("working0")

print(paste0("partition percent used is",partition_per))
set.seed(121) # Set Seed so that same sample can be reproduced in future also
# Now Selecting 75% of data as sample from total 'n' rows of the data  
sample <- sample.int(n = nrow(o), size = floor(partition_per*nrow(o)), replace = F)

train <-  as.data.frame(o[sample, ])
test  <-  o[-sample, ]



fwrite(train,"G4_training.tsv",sep="\t",col.names=TRUE,quote=FALSE)
fwrite(test,"G4_test.tsv",sep="\t",col.names=TRUE,quote=FALSE)


train.original<-train

writeFasta(train[,c("id","sequence")],"g4_training.fasta")
#writeFasta(test[,c("id","sequence")],"g4_test.fasta")


#----------------------------------------------------------


latermatrix="cdhit"
if (latermatrix=="cdhit"){

commandoptions<-" -n 5 -aS 0.8 -b 11 -t 5 -g 1  "
#file_path<-"g4_training.fasta"

command1<-paste0("cd-hit  -i ./g4_training.fasta -o ./cdhit_fullfastafile.result ",commandoptions) 

system(command1)

data<-fread(file = 'cdhit_fullfastafile.result', sep = '\t', header = F)

clusters<-parsed_result('cdhit_fullfastafile.result.clstr')
clusters$sequence<-gsub(" ","",clusters$sequence)
clusters$sequence<-gsub("\\...*","",clusters$sequence)
clusters$stat<-gsub("%","",clusters$stat)
clusters$stat<-gsub(" ","",clusters$stat)

table(clusters$cluster)%>% as.data.frame()%>% arrange(desc(Freq))%>% head(20)%>% print()
colnames(clusters)<-c("cluster","nt","id","stat","count")
data.2<-left_join(clusters,train,by=c("id"="id"))%>% as.data.frame()
final.1<-as.data.frame(data.2)
}

final.1%>% filter(cluster=="Cluster 479")

print(paste0(rep("_",100),collapse=""))
print(paste0(rep("_",100),collapse=""))
print(paste0(rep("_",100),collapse=""))
print(paste0(rep("_",100),collapse=""))
print(paste0(rep("_",100),collapse=""))



pairwise_Dir <- paste0(SubDir,"SECOND_cluster_",added_name_for_output_file,"/")
pairwise_dist_folder <<- paste0(pairwise_Dir,"SECOND_cluster","/")

dir.create(file.path(pairwise_Dir), showWarnings = FALSE,recursive=TRUE)
dir.create(file.path(pairwise_dist_folder), showWarnings = FALSE,recursive=TRUE)

setwd(file.path(pairwise_Dir))



clusters_topair<-final.1%>%  select(cluster)%>%
 unlist() %>% unique()

#sequences_nottopair<-final.1%>% group_by(diana_clusters)%>% mutate(n=n())%>% filter(!(n>1 & !is.na(diana_clusters)))%>% select(sequence)%>%
# unlist()%>% unique()
#notpairedseq<-final.1[final.1$sequence %in% sequences_nottopair]
#notpairedseq$cluster<-c("bl_singleton",pairseq$diana_clusters,"a")

#notpairedseq<-notpairedseq[,c("sequence", "name" , "ggg", "length_trimmed","id","cluster")]

pairseq<-final.1[final.1$cluster %in% clusters_topair,]




X1<-pairseq%>% select("cluster","id","sequence")
print("X1")
print(head(X1))




colnames(X1)<-c("cluster","id","sequence")
X<-split(X1,X1$cluster)
#X<-purrr::imap(X1, ~mutate(.x, cluster = .y))

n<-length(X)
#print(head(X[[1]]))
print(paste0("length is",n))
print(head(pairseq))

n=3
final_result<-lapply(1:n,function(i){

 #tryCatch({
dat_frame<-as.data.frame(X[[i]])

result<-second_clustering(dat_frame)

return(result)  #},error=function(e) return(NULL))
})


temp<-rbindlist(final_result)
temp$final_c<-paste0(temp$cluster,"_",temp$hclusters2,"a")
print(head(final_result))




table(temp$cluster)%>% as.data.frame()%>% arrange(desc(Freq))%>% head(10)


#---------------------------------------------------------------------

save(final_result,file="FINAL.RData")



#---------------------------------------------------------------------------
cluster_fasta_Dir<-paste0(SubDir,"fasta_",added_name_for_output_file,"/")

#dir.create(file.path( cluster_fasta_Dir ), showWarnings = FALSE,recursive=TRUE)


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


final_result2<-final_result

library(purrr)



result3<-lapply(1:length(final_result2),function(y){

results3b<-map(final_result2,~.x %>% as.data.frame()%>% dplyr::group_by(hclusters2)%>% mutate(n=n())%>% filter(n<10)%>% ungroup()%>%
mutate(cluster2=ifelse(n >= 13,hclusters2,"small"))%>% mutate(cluster3=paste0(cluster,"_",cluster2,"b"))%>%
select(id,cluster3,sequence))

result4<-results3b %>% bind_rows(.id="source")
names(result4)[grepl("source", names(result4))] <- "osource"

print(head(result4))
print(paste0("NROW==",nrow(result4)))
return(results3b)})


result4<-bind_rows(result3,.id="source")
result4$cluster<-paste0("ggg_",result4$source,"_",result4$cluster3)
data1<-split(result4,result4$cluster)

data2<-map(data1,~writeFasta1(.x,cluster_fasta_Dir))



print('completed')


#})
#print(MSA_Dir)
print("hurray")


#})




