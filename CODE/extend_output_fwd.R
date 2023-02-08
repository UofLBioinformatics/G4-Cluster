require(stringr) 

require(Biostrings) 
library(BSgenome.Hsapiens.UCSC.hg38)
library(data.table)
genome<-BSgenome.Hsapiens.UCSC.hg38

output_fwd<-fread("/bio/home/goonerrn/g_quadruplex/annotation_database/query/output_fwd_withCHR.txt",header=FALSE)

output_fwd<-separate(output_fwd,col="V2",into=c("start","end"),sep="-",remove=FALSE)
output_fwd$start1<-as.numeric(output_fwd$start)-1

output_fwd$end1<-as.numeric(output_fwd$end)+3


seq_g4_full <- getSeq(genome, names =output_fwd$V1,
                               start = (as.numeric(output_fwd$start1)),
                               end = (as.numeric(output_fwd$end1)))
  

bich1<-cbind(output_fwd,as.data.frame(seq_g4_full))
bich1$start1<-NULL
bich1$end1<-NULL
bich1$start<-NULL
bich1$end<-NULL
fwrite(bich1,"/bio/home/goonerrn/g_quadruplex/annotation_database/query/output_fwd_withCHR_2.txt",sep="\t",quote=FALSE,col.names=FALSE)