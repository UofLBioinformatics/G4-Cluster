

#!/bin/bash

library(data.table)
library(purrr)
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
library(dplyr)


source("/bio/home/goonerrn/g_quadruplex/annotation_database/G_quad_clustering/clustering_output/cluster_functions.R")


args <- commandArgs(trailingOnly = TRUE)

total.start.time <- Sys.time()

#metrictouse="G4score"
metrictouse="evalue"
added_name_for_output_file<-args[1]
latermatrix="blat"

partition_per<-args[2]
partition_per<-as.numeric(partition_per)/100


source("/bio/home/goonerrn/g_quadruplex/annotation_database/G_quad_clustering/starcode/codes/starcode_functions_blat.R")
#source("/bio/home/goonerrn/g_quadruplex/annotation_database/G_quad_clustering/starcode/codes/recovery_graphcluster.R")

output_fwd<-fread("/bio/home/goonerrn/g_quadruplex/annotation_database/query/output_fwd_withCHR_2.txt",header=FALSE)


chromosome=paste0("chr",c(1:22,"X","Y"))

output_fwd<-output_fwd[output_fwd$V1 %in% chromosome]
output_fwd$name<-paste0(output_fwd$V1,":",output_fwd$V2)
output_fwd<-output_fwd%>%separate(V3,c("nggg","nquad","fquad"),":")
output_fwd$nggg<-as.numeric(as.character(output_fwd$nggg))
output_fwd<-as.data.frame(output_fwd)
output_fwd.3<- output_fwd %>% as.data.frame() %>% dplyr::filter(nggg==4 | nggg==5 | nggg==6 | nggg==7 | nggg==8)
output_fwd.3$nggg_g<-ifelse(output_fwd.3$nggg<6,output_fwd.3$nggg,"6_edit")
colnames(output_fwd.3)[8]<-"sequence"
output_fwd.3$length<-nchar(output_fwd.3$sequence)

output_fwd.2<-split(output_fwd.3,output_fwd.3$nggg_g)
#------------------------------------------------------------------------

n=length(output_fwd.2)


final_result2<-invisible(lapply(1:n,function(i){
#output_fwd.4<-output_fwd.3%>% filter(nggg_g=="6_edit")
#output_fwd.4<-output_fwd.3%>% filter(nggg_g=="4")


output_fwd.4<-as.data.frame(output_fwd.2[[i]])

nggg_g<-as.character(output_fwd.4$nggg_g[1])
print(nggg_g)
print("Here we go")
########the function starts here
result<-nggg_function(output_fwd.4,nggg_g,latermatrix)
print(paste0(rep("-",50),collapse=""))
print(paste0(i,"            completed         ",nggg_g))
print(paste0(rep("-",50),collapse=""))
filename=paste0("final_result_",nggg_g,".RData")
print(filename)
save(result,file=filename)
return(result)}))#,mc.cores=8))


#---------------------------------------------------------------------

save(final_result2,file="FINAL.RData")

#---------------------------


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




library(purrr)



result3<-lapply(1:length(final_result2),function(y){

results3b<-map(final_result2[[y]],~.x %>% group_by(hclusters2)%>% mutate(n=n())%>% ungroup()%>%
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




