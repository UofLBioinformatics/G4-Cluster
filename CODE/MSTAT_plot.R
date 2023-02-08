#!/bin/bash




library(ggplot2)
library(dplyr)
library(data.table)
library(tidyr)
library(purrr)

list.dirs <- function(path=".", pattern=NULL, all.dirs=FALSE,
  full.names=FALSE, ignore.case=FALSE) {
  # use full.names=TRUE to pass to file.info
  all <- list.files(path, pattern, all.dirs,
           full.names=TRUE, recursive=FALSE, ignore.case)
  dirs <- all[file.info(all)$isdir]
  # determine whether to return full names or just dir names
  if(isTRUE(full.names))
    return(dirs)
  else
    return(basename(dirs))
} 



whichfile <- grep(
  x = list.files(),
  pattern = "result.txt",
  value = TRUE
)


#args <- commandArgs(trailingOnly = TRUE)


#edit6_path<-args[1]
#path_to_use<-edit6_path
listscorepath<-paste0("MSTAT_out_trident_",times)
(pre <- list.files(path_to_use, listscorepath, recursive=TRUE, full.names=TRUE, include.dirs=TRUE))

if (length(pre)>1){
details <- file.info(list.files(path_to_use, listscorepath, recursive=TRUE, full.names=TRUE, include.dirs=TRUE))
details <- details[with(details, order(as.POSIXct(mtime))), ]
pre = rownames(details)
pre<-pre[length(pre)]
}


print(pre)
MSTAT_plot<-function(mstat_path,times){
foldername<-dirname(mstat_path)
base_file<-basename(foldername)
infile<-paste0((mstat_path),"/","output_screen.txt")

rl <- readLines(infile)
rl1<-gsub("Multiple alignment : ","",rl[grep('^Multiple*', rl)])
rl1<-as.data.frame(rl1)
rl1<-separate(rl1,col="rl1",into=c("seq","col"),sep=",")
rl1$seq<-gsub("nb seq = " ,"",rl1$seq)
rl1$col<-gsub("nb col = " ,"",rl1$col)

infile_filename<-fread(paste0((mstat_path),"/","filename.txt"),header=FALSE)
infile_result<-fread(paste0((mstat_path),"/","result.txt"),header=FALSE)


x<-cbind(infile_filename,infile_result,rl1)
x$id<-rownames(x)

colnames(x)<-c("file","score","seq","col","id")
x1<-x %>% as.data.frame()%>%
  filter_all(~ !is.na(.))
x1$basenames <- sub("\\.[[:alnum:]]+$", "", basename(as.character(x$file)))

filename<-paste0(foldername,"/",base_file,"_stat_bubble_gap_try_",times,".pdf")
pdf(filename)
print("filename")
p<-x1 %>%
  arrange(desc(as.numeric(col))) %>%
  mutate(id = factor(id)) %>%
  ggplot(aes(x=as.numeric(col), y=as.numeric(score), size = as.numeric(seq))) +
    geom_point(alpha=0.5) +
    scale_size(range = c(.1, 10), name="Number of sequences in family")+xlab("TOTAL COLUMNS IN MSA")+
ylab("gaps")+
  geom_point(data = x1[which.min(x1$score), ], color="blue", 
               size=5) +
    geom_point(data = x1[which.max(x1$score), ], color="red", 
               size=5) +
    geom_text(data = rbind(x1[which.min(x1$score), ], x1[which.max(x1$score),]), aes(as.numeric(seq)+10,score, label=basenames), size=2)

print(p)

dev.off()

}



MSTAT_plot_filter<-function(mstat_path,times){
foldername<-(dirname(mstat_path))
base_file<-basename(foldername)
infile<-paste0((mstat_path),"/","output_screen.txt")

rl <- readLines(infile)
rl1<-gsub("Multiple alignment : ","",rl[grep('^Multiple*', rl)])
rl1<-as.data.frame(rl1)
rl1<-separate(rl1,col="rl1",into=c("seq","col"),sep=",")
rl1$seq<-gsub("nb seq = " ,"",rl1$seq)
rl1$col<-gsub("nb col = " ,"",rl1$col)

infile_filename<-fread(paste0((mstat_path),"/","filename.txt"),header=FALSE)
infile_result<-fread(paste0((mstat_path),"/","result.txt"),header=FALSE)


x<-cbind(infile_filename,infile_result,rl1)
x$id<-rownames(x)

colnames(x)<-c("file","score","seq","col","id")
x1<-x %>% as.data.frame()%>%
  filter_all(~ !is.na(.))
x1$basenames <- sub("\\.[[:alnum:]]+$", "", basename(as.character(x$file)))
filename_filtered_df<-paste0(foldername,"/",base_file,"MSTAT_final.tsv")

if (times>1){
threshold_gapscore<-0.6
} else if (times==1) {
threshold_gapscore<-0.6
}

x2<-x1%>% as.data.frame()%>% filter((as.numeric(score)<threshold_gapscore) & (as.numeric(seq)>4))%>%
 arrange(desc(as.numeric(score)))
fwrite(x2, file=filename_filtered_df,sep="\t",col.names=TRUE,quote=FALSE)

filename<-paste0(foldername,"/",base_file,"_stat_bubble_gap_filtered_",times,".pdf")
pdf(filename)
print("filename")
p<-x2 %>%
  arrange(desc(as.numeric(col))) %>%
  mutate(id = factor(id)) %>%
  ggplot(aes(x=as.numeric(col), y=as.numeric(score), size = as.numeric(seq))) +
    geom_point(alpha=0.5) +
    scale_size(range = c(.1, 10), name="Number of sequences in Family")+xlab("TOTAL COLUMNS IN MSA")+
ylab("gaps")+
  geom_point(data = x2[which.min(x2$score), ], color="blue", 
               size=5) +
    geom_point(data = x2[which.max(x2$score), ], color="red", 
               size=5) +
    geom_text(data = rbind(x2[which.min(x2$score), ], x2[which.max(x2$score),]), aes(as.numeric(seq)+10,score, label=basenames), size=2)

print(p)

dev.off()

}




purrr::map(1:length(pre),~MSTAT_plot(pre[.x],times))

purrr::map(1:length(pre),~MSTAT_plot_filter(pre[.x],times))
