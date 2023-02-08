#!/bin/bash

library(data.table)
library(purrr)
library(dplyr)
library(tidyr)

output_fwd<-fread("/bio/home/goonerrn/g_quadruplex/annotation_database/query/output_fwd_withCHR.txt",header=FALSE)

chromosome=paste0("chr",c(1:22,"X","Y"))

output_fwd<-output_fwd[output_fwd$V1 %in% chromosome]
output_fwd$name<-paste0(output_fwd$V1,":",output_fwd$V2)
output_fwd<-output_fwd%>%separate(V3,c("nggg","nquad","fquad"),":")
output_fwd$nggg<-as.numeric(output_fwd$nggg)
#output_fwd$grp1<-ifelse(as.numeric(output_fwd$nggg>7),1,as.numeric(output_fwd$nggg))
#output_fwd$grp1<-as.factor(output_fwd$grp1)

output_fwd.1<-output_fwd%>% filter(nggg==4 | nggg==5 | nggg==6 | nggg==7)

colnames(output_fwd.1)[6]<-"sequence"

#sequence.1<-gsub('(^[G])\\1+','\\1',output_fwd.1$sequence)
#output_fwd.1$sequence.1<-gsub('([G])\\1+$','\\1',sequence.1)

o<-output_fwd.1%>% group_by(sequence)%>% 
summarise(name=paste0(unique(name),collapse="|"), ggg=paste0(unique(nggg),collapse="|"),
length_trimmed=nchar(sequence)) %>% filter(length_trimmed>5)%>% as.data.frame()
min_len<-min(o$length_trimmed)

add=list()
for (i in ((min_len):(min_len+10))){
add=append(add,(paste0(rep("G",i),collapse="")))
}

add<-as.data.frame(t(as.data.frame(add)))
row.names(add)<-NULL
add$V2<-"added"
add$V3<-0
add$V1<-as.character(add$V1)
add$V4<-nchar(as.character(add$V1))


colnames(add)<-colnames(o)
o<-rbind(add,o)


x<-o[,1]
y<-o[,1:ncol(o)]

#o[,3]$length<-nchar(o[,1])
clusters<-list()
for (i in c(1:4)){

filename.g4=paste0("G4_seq",i,".txt")
filename.name=paste0("name_seq",i,".txt")

fwrite(as.data.frame(x),file=filename.g4,sep="\t", col.names=FALSE,quote=FALSE)
fwrite(as.data.frame(y),file=filename.name,sep="\t", col.names=TRUE,quote=FALSE)

timestamp=system(" (date +%d-%m-%Y_%H-%M-%S)")
dvalue=i
file_out=paste0("Cluster.d",dvalue,"___",as.character(timestamp))
code=paste0("/bio/home/goonerrn/tools/starcode-1.3/starcode  -t8 -d",dvalue," -s --input ", filename.g4, " -o ", file_out, " --print-clusters")
system(code)

cls<-fread(file_out)

cls<-cls%>% filter(V2>1)

cls$id<-row.names(cls)
if (nrow(cls)>0){
cls.1<-cls%>% separate_rows(V3)
}

cls.1$V1len<-nchar(cls.1$V1)

cls.1$V3len<-nchar(cls.1$V3)
table(cls.1$V1len==cls.1$V3len)
#cls.2<-cls.1%>% mutate(div=ifelse(V1len==V3len,"same","notsame"))
#cls.2<-cls.1%>% mutate(div=ifelse(abs(V1len-V3len)<2,"2same",div))
#cls.1<-cls.1%>% filter(V1!=V3)
diff<-ifelse(i<=4,2,4)
cls.2<-cls.1%>% mutate(div=ifelse(abs(V1len-V3len)<diff,"same","notsame"))
table(cls.2$div)
cls.2%>% group_by(id)%>% summarise(n=n(),v1n=n_distinct(V1),v3n=n_distinct(V3))
cls.3<-cls.2%>% filter(div=="same")%>% group_by(id)%>% summarise(n=n())%>% filter(n>30)

if (nrow(cls.3)>0){
count<-unique(cls.3$id)

withcluster<-map(count,~cls.2%>% filter(id==.x & div=="same")%>% select(V3,V1)%>% unlist()%>% 
unique()%>% as.data.frame()) 
#withcluster<-data.table(withcluster)
withcluster<-rbindlist(withcluster,idcol=TRUE)
colnames(withcluster)<-c("clus","sequence")
toremove<-withcluster$sequence %>% unlist()
clusters[[length(clusters)+1]] <- list(withcluster)
x<-x[!x %in% toremove]
y<-y[!x %in% toremove,]

} else {(count<-NULL)
withcluster<-data.frame(clus=NA,sequence=NA)

clusters[[length(clusters)+1]] <- list(withcluster)
}

}




#clusters_1<-lapply(rbindlist(clusters,idcol=TRUE)

clusters_1<-clusters%>% map(as.data.frame)%>% bind_rows(.id="source")
clusters_1%>% filter(!is.na(clus))%>%group_by(source)%>% summarise(n=n())%>% arrange(desc(n))

clusters_1%>% filter(source==3)%>% head()
g4name<-fread("name_seq1.txt",sep="\t")


clusters_2<-left_join(g4name,clusters_1)%>% filter(name!="added")







writeFasta<-function(data, filename){
  fastaLines = c()
  for (rowNum in 1:nrow(data)){
    fastaLines = c(fastaLines, as.character(paste(">", data[rowNum,"name"], sep = "")))
    fastaLines = c(fastaLines,as.character(data[rowNum,"seq"]))
  }
  fileConn<-file(filename)
  writeLines(fastaLines, fileConn)
  close(fileConn)
}