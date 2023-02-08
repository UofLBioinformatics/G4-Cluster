#bed_to_granges
library("BiocParallel")
 library(igraph)
   library("GenomicRanges")
library(dplyr)
library(Matrix)
library(data.table)
library(tidyr)
library(GenomicRanges)
library(gtools)
library(cluster)
library(purrr)

library(Biostrings)

setwd("/bio/home/goonerrn/g_quadruplex/annotation_database/G_quad_clustering/")

bed_to_granges <- function(file,bed,sample){
   df1 <- fread(file,
                    header=F)
     x <- c("alt", "random","Un_KI")
     df1<-df1[!grepl(paste(x, collapse = "|"), df1$V1)]
     if (sample==TRUE){
     df2<-df1[sample(nrow(df1), 1000), ]}
     
     else (df2<-df1)
   if (bed==TRUE){
 ##remove it for all sequences
  
   
   df2$name<-paste0(df2$V1,":",df2$V2,"-",df2$V3)
      setkey(df2, name)
   setkey(cluster_file, V1) 
   
   df<-cluster_file[df2]
   if("cluster" %in% colnames(df))
{
  cat("Yep, there is a cluster column in there!\n");
   }
      
   if(length(df) > 7){
      df <- df[,-c(8:length(df))]
   }
 
   if(length(df)<3){
      stop("File has less than 3 columns")
   }
 
    header <- c('name','cluster','chr','start','end','strand','id')
   names(df) <- header[1:length(names(df))]
 
   if('strand' %in% colnames(df)){
      df$strand <- gsub(pattern="[^+-]+", replacement = '*', x = df$strand)
   }
 
 
   if(length(df)==3){
      gr <- with(df, GRanges(chr, IRanges(start, end)))
   } else if (length(df)==4 ){
      gr <- with(df, GRanges(chr, IRanges(start, end), id=id))
   } else if (length(df)==5 ){
      gr <- with(df, GRanges(chr, IRanges(start, end), id=id, score=score))
   } else if (length(df)==6 ){
      gr <- with(df, GRanges(chr, IRanges(start, end), id=id, score=score, strand=strand))
   } else if (length(df)==7 ){
      gr <- with(df, GRanges(chr, IRanges(start, end),name=name, id=id,strand=strand,cluster=cluster))
   }
   
   return(gr)
   }

  else if (bed==FALSE){
  #df1 <- fread("C:/Users/goone/OneDrive - University of Louisville/Research/Clusters/files/output_fwd_withCHR.txt",  header=F)
   
   
   df2$name<-paste0(df2$V1,":",df2$V2)
   df2<-separate(df2,V2,into=c("start","end"),sep="-")
  setkey(df2, name)
   setkey(cluster_file, V1) 
   
   df3<-cluster_file[df2]
   if("cluster" %in% colnames(df3))
{
  cat("Yep, there is a cluster column in there!\n
      Removing rows with NA cluster");
    df<- df3%>% drop_na("cluster")
   }
   
      
 
   if(length(df)<3){
      stop("File has less than 3 columns")
   }
 
    header <- c('name','cluster','chr','start','end','type_g4','sequence','strand')
   names(df) <- header[1:length(names(df))]
    if('strand' %in% colnames(df)){
      df$strand <- gsub(pattern="[^+-]+", replacement = '*', x = df$strand)
    }
cols.num <- c('start','end')
df<-as.data.frame(df)
df[cols.num] <- lapply(df[cols.num],as.numeric)
sapply(df, class)
 
    if(length(df)==3){
      gr <- with(df, GRanges(chr, IRanges(start, end)))
   } else if (length(df)==4 ){
      gr <- with(df, GRanges(chr, IRanges(start, end), id=id))
   } else if (length(df)==5 ){
      gr <- with(df, GRanges(chr, IRanges(start, end), id=id, score=score))
   } else if (length(df)==6 ){
      gr <- with(df, GRanges(chr, IRanges(start, end), id=id, score=score, strand=strand))
   } else if (length(df)==8 ){
      gr <- with(df, GRanges(chr, IRanges(start, end),name=name, sequence=sequence,strand=strand,cluster=cluster,type_g4=type_g4))
   }
   
   return(gr)

   }
   
 
}



#*#*#*#*#*#*#*#*#*#
#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*

kmer_from_granges<-function(seq_2_level2){
  print("reached kmer_from_granges")
new_seq<-DNAStringSet(seq_2_level2$sequence)
names(new_seq)<-(seq_2_level2$name)
#metadata(new_seq)$cluster <- seq_2_level2$sub_cluster
#metadata(new_seq)$type <- seq_2_level2$type_g4

#####WILL DO PAIRWISE IN PLACE OF THIS
G4.kmers<-kmer_by_gene(new_seq, 4)
G4.kmers[,-1]<-G4.kmers[,-1]+0.000001
df <- normalize_function(G4.kmers[,-1])
row.names(df)<-(G4.kmers$name)
####
return(df)

}


  #*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#
kmer_by_gene <- function(DNAStringSet, k){
   print("counting kmers...", quote = FALSE)
   kmer_counts <- Biostrings::oligonucleotideFrequency(DNAStringSet, width = 4) %>% 
     dplyr::as_tibble() %>% 
     dplyr::mutate(name = names(DNAStringSet)) %>% 
     dplyr::select(name, dplyr::everything())
   
   print("counting complete.", quote= FALSE)
   
   return(kmer_counts)
}
#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*
level_1<-function(final_sub_cluster){
z<-("working for level 1")
print(paste(z,(final_sub_cluster)))
g4_names<- clus_2 %>% filter(sub_cluster==final_sub_cluster) %>% select(name)
#g4_names<- clus_2[,"g4"][clus_2[,"sub_cluster"]==final_sub_cluster]
seq_2_level1<-as.data.frame(my_granges[my_granges$name %in% g4_names$name])
seq_2_level1$sub_cluster<-as.character(final_sub_cluster)
seq_2_level1$subcluster_level_2<-as.character(final_sub_cluster)

#print(colnames(seq_2_level1))

return(seq_2_level1)
}


#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*

level_2<-function(sub_cluster1){
    z<-("reached level_2")
print(paste(z,(sub_cluster1)))

 g4_names<- clus_2 %>% filter(as.character(clus_2$sub_cluster)==as.character(sub_cluster1)) %>% select(name)
   #g4_names<- clus_2[,"g4"][clus_2[,"sub_cluster"]==sub_cluster]
  seq_2_level2<-(my_granges[my_granges$name %in% g4_names$name])
  seq_2_level2$sub_cluster<-as.character(sub_cluster1)
  ###just for now, later must use pairwise distance
  
new_seq<-DNAStringSet(seq_2_level2$sequence)
names(new_seq)<-(seq_2_level2$name)


g4_kmer<-kmer_from_granges(seq_2_level2)
#G4.kmers<-kmer_by_gene(new_seq, 4)
print("round 1, NORMALIZE")
name<-row.names(g4_kmer)
#df <- normalize_function(G4.kmers[,-1])
#row.names(df)<-(G4.kmers$g4)

print("round 1, GET DISTANCE")

G4_dist_mat<-(dist(g4_kmer, method = "manhattan"))

# Use map_dbl to run many models with varying value of k
sil_width <- map_dbl(2:7,  function(k){
 model<-pam(G4_dist_mat,k = k,diss = TRUE)
  model$silinfo$avg.width
})

# Generate a data frame containing both k and sil_width
sil_df <- data.frame(
  k = 2:7,
  sil_width = sil_width
)


optimal_cluster<-sil_df[which.max((sil_df$sil_width)),]$k
print_seq<-"the optimal cluster was found to be"
print_seq2<-"for cluster number"
print(paste(print_seq,optimal_cluster,print_seq2,sub_cluster1))
model<-pam(G4_dist_mat,k = optimal_cluster,diss = TRUE)
model$clustering<-paste0(sub_cluster1,"_",model$clustering,"b")
clus_1<-cbind(as.data.frame(name),model$clustering)
clus_3<-merge(clus_2,clus_1,by="name")
colnames(clus_3)[11]<-"subcluster_level_2"





#x<-cbind(as.data.frame(seq_2_level2),subcluster_level_2)
x<-clus_3
#if (length(seq_2_level2)==length(x)){
#  print("level 2 same row number in input and output TRUE")
#}else {
#  print(paste(length(seq_2_level2)),"and",length(x),"the number of columns, and FALSE")
#}

return(x)
}





#######*#*#*#*#*#*#*#*#*#*#**#


pairwise_dist<-function(seq_2_level2){
  
sub_mat<-nucleotideSubstitutionMatrix(match = 1, mismatch = 0, baseOnly = TRUE, type = "DNA")
sub_mat['A','A']<-30
sub_mat['G','G']<-10
sub_mat['C','C']<-30
sub_mat['T','T']<-30
sub_mat['A','G']<--40
sub_mat['C','G']<--40
sub_mat['G','T']<--40
sub_mat['A','C']<--40
sub_mat['A','T']<--40
sub_mat['C','T']<--40

sub_mat<-Matrix::forceSymmetric(sub_mat)
sub_mat<-as.matrix(sub_mat)
#my_param <- McParam(workers = 8)
filename<-seq_2_level2$cluster[1] #selecting just one because all are same
seq_2<-c(seq_2_level2$sequence)
name<-c(seq_2_level2$name)
#my_comb <- combn(x = 1:length(seq_2), m = 2)
#my_comb_name <- combn(x = name, m = 2)
z<-"DISTANCE MATRIX FOR cluster"
print(paste(z,filename,sep=" "))
#my_comb_name<-combinations(length(name),2, name,repeats.allowed = TRUE)

print("GOT INSIDE DISTANCE MATRIX")


# 

my_comb <- combinations(length(seq_2),2, 1:length(seq_2),repeats.allowed = TRUE)
#my_comb <- combinations(5,2, 1:5,repeats.allowed = TRUE)

my_comb_name<-combinations(length(name),2, name,repeats.allowed = TRUE)
#my_comb_seq<-combinations(length(seq_2),2, seq_2,repeats.allowed = TRUE)

#as.data.frame(cbind(my_comb,my_comb_name))
#my_comb <- combn(x = 1:length(seq_2), m = 2)

my_comb_list <- lapply(apply(X = my_comb, MARGIN = 1, FUN = list), unlist)
length(my_comb_list)

# it should be list of pairs of seqs
#my_comb_list <- my_comb

my_seq<-c(seq_2_level2$sequence)
print("RUNNING BEFORE ALIGNMENT")
 align_wrapper <- function(x){
  pairwiseAlignment(pattern = my_seq[x[1]], subject = my_seq[x[2]], type = "global", substitutionMatrix=sub_mat,scoreOnly=TRUE)
}

system.time(
  my_alignments <- bplapply(my_comb_list, align_wrapper, BPPARAM = safeBPParam(8))
)

df<-as.data.frame((my_comb_name),row.names = NULL)
print("done alignment")
df1<-cbind(df,unlist((my_alignments)))
colnames(df1)<-c("From","To","Weight")
mygraph <- graph.data.frame(df1)


pair_mat<-get.adjacency(mygraph,sparse=FALSE,attr='Weight',type="upper")
pair_mat<-Matrix::forceSymmetric(pair_mat)
pair_mat<-as.matrix(pair_mat)
pair_mat<-pair_mat
if (isTRUE(all.equal( max(pair_mat) ,min(pair_mat)))==FALSE){
pair_mat_norm<-normalize_function(pair_mat)
pair_norm_dist<-(1-pair_mat_norm)
  path_out_mat<-"/bio/home/goonerrn/g_quadruplex/annotation_database/G_quad_clustering/sim_matrix/"
  filename2=paste(path_out_mat,"cluster_",filename,"sim_matrix.RData",sep="")
save(pair_mat,file=filename2)
  
  
#sim2dist <- function(mx) as.dist(sqrt(outer(diag(mx), diag(mx), "+") - 2*mx))

#d.mx<-as.matrix(sim2dist(pair_mat))
return(pair_norm_dist)
 

} else if ((isTRUE(all.equal( max(pair_mat) ,min(pair_mat)))==TRUE)){
  
 pair_mat_norm<-(pair_mat)
pair_norm_dist<-(max(pair_mat_norm))-pair_mat_norm
  path_out_mat<-"/bio/home/goonerrn/g_quadruplex/annotation_database/G_quad_clustering/sim_matrix/"
  filename2=paste(path_out_mat,"cluster_",filename,"sim_matrix.RData",sep="")
save(pair_mat,file=filename2)
  

return(pair_norm_dist)

  
}

 

}
#*#*#*#*#*#*#*#*#*#*#**#*#*#*#*#*#*#*#*#**#*#*#*#*#*#*#*#**#*#**#*
 

#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#**#*#*#*#*
 normalize_function <- function(x){(x-min(x))/(max(x)-min(x))}
##############################
 
 
 
 safeBPParam <- function(nworkers) {
    if (.Platform$OS.type=="windows") {
        BiocParallel::SerialParam()
    } else {
        BiocParallel::MulticoreParam(nworkers)
    }
 }
 
 
#################################################################################################################################################

##FILES TO READ


###read cluster files
cluster_file<-fread("/bio/home/goonerrn/g_quadruplex/annotation_database/G_quad_clustering/cluster_MCL_I16.I20.txt")

#bed file
#convert to g_ranges
#my_granges <- bed_to_granges("C:/Users/goone/OneDrive - University of Louisville/Research/annotation_G_quad_28feb2020/4_bed.bed",bed=TRUE)
my_granges <- bed_to_granges("/bio/home/goonerrn/g_quadruplex/annotation_database/query/output_fwd_withCHR.txt",bed=FALSE,sample=FALSE)
#knitr::opts_chunk$set(echo = TRUE)
#knitr::opts_knit$set(root.dir = "")
#getwd()
#getwd()



















#################################################
final_df<-function(cluster_name){
seq<-my_granges[my_granges$cluster %in% cluster_name]
new_seq<-DNAStringSet(seq$sequence)
names(new_seq)<-(seq$name)

if (length(new_seq)>7){

metadata(new_seq)$cluster <- seq$cluster
metadata(new_seq)$type <- seq$type_g4

print("round 1, Distance matrix ")

########################
dat_frame<-as.data.frame(seq)
dist_mat<-pairwise_dist(dat_frame)

if (all(dist_mat==0)==FALSE){
print("GOT DISTANCE MATRIX")
G4_dist_mat_norm<- as.dist(dist_mat)
 print("now for hclust")
hc<-hclust(G4_dist_mat_norm,method="ward.D2")
 hc$height <- round(hc$height,4)
mycl <- cutree(hc, h=max(hc$height/1.2))
#plot(hc)
idx <- order(names(mycl))
clust.cutree <- as.data.frame(mycl[idx])
clust.cutree$name <- rownames(clust.cutree)
df.merge <- merge(as.data.frame(seq),clust.cutree,by.x='name')
#df.merge.sorted <- df.merge[order(df.merge$y),]
df.merge$subcluster<-paste0(df.merge$cluster,"_",df.merge$mycl[idx],"a")
cols.dont.want<-"mycl[idx]"
df.merge <- df.merge[, !names(df.merge) %in% cols.dont.want, drop = F]
#df.merge
#######################################
print("DONE FIRST CLUSTERING")


#clus_2<-cbind(G4.kmers[,1],model$clustering)

clus_2<-df.merge
colnames(clus_2)[10]<-"sub_cluster"
assign("clus_2", clus_2, envir = .GlobalEnv)


new_df<-clus_2%>%group_by(sub_cluster)%>%
summarise(Freq=n(),.groups = 'drop')

new_df<-as.data.frame(new_df)
#a<-as.data.frame(table(new_df$sub_cluster))
print("STEP 1: total should be equal to above as clus_2")
#print(a)

if (sum(new_df$Freq)==nrow(clus_2)){
  y<-paste("and total for subclusters is ",sum(a$Freq))
z<-("After 1st cluster number of rows here  in clus_2 ")
print(paste(z,"is",nrow(clus_2),y)) 
print("TRUE step1")
}else {
    print ("FALSE step 1")}

final_sub_cluster<-new_df[new_df[,"Freq"]<=9,]$sub_cluster
sub_cluster<-as.character(new_df[new_df[,"Freq"] > 9,]$sub_cluster)
print("clear")



if (length(sub_cluster)>0) {
clus_level_2 <-(mapply(level_2,sub_cluster,SIMPLIFY = FALSE)) %>%   bind_rows()
}
if (length(final_sub_cluster)>0){
clus_level_1<-(mapply(level_1,final_sub_cluster,SIMPLIFY = FALSE))  %>% bind_rows()
}


if ((exists("clus_level_1")==TRUE & exists("clus_level_2")==FALSE)){
  final_df<-clus_level_1
  z<-("number of rows here  in clus_1 ")
  print(paste(z,"is",nrow(clus_2),"output row is" ,(nrow(final_df))))
  
  #rint(colnames(clus_level_1))

}else if ((exists("clus_level_2")==TRUE & exists("clus_level_1")==FALSE)){
  final_df<-clus_level_2
    print("Clus_level_2 number of rows at the end")
      z<-("number of rows here  in clus_2 ")
  print(paste(z,"is",nrow(clus_2),"output row is" ,(nrow(final_df))))  #print(colnames(clus_level_2))
}else if (exists("clus_level_1")==TRUE & exists("clus_level_2")==TRUE){
final_df<-rbind(clus_level_1,clus_level_2)
print(paste("level 1 is",nrow(clus_level_1),"level 2 is",nrow(clus_level_2)))
#print(colnames(clus_level_1))
print("Now final step")
print(paste(z,"is",nrow(clus_2),"output row is" ,(nrow(final_df))))
#print(colnames(clus_level_2))}
}

if (length(clus_2)==nrow(final_df)){
print("the length of clus 2 is equal to final df output.TRUE")
}else {z<-("the length of clus 2 is NOT equal to final df output.FALSE")
  print(paste(z,length(clus_2),nrow(final_df)))
  print("the two columns now are clus_2 and final_df")
  #print(head(clus_2))
  print("*************")
  #print(head(final_df))}
}
rm(clus_2)
rm(clus_level_2,clus_level_1)

return(final_df)

} else if (all(dist_mat==0)==TRUE) {
seq$sub_cluster<-as.character(cluster_name)
seq$subcluster_level_2<-as.character(cluster_name) 
seq<-as.data.frame(seq)
#print(paste(cluster_name," length is",nrow(seq)))
return(seq)
  
}
} else{
  seq$sub_cluster<-as.character(cluster_name)
seq$subcluster_level_2<-as.character(cluster_name) 
seq<-as.data.frame(seq)
#print(paste(cluster_name," length is",nrow(seq)))
return(seq)
}
}

######################################################################################################################



afastafile<-DNAStringSet(my_granges$sequence)
names(afastafile) <- my_granges$name
metadata(afastafile)$cluster <- my_granges$cluster
metadata(afastafile)$type <- my_granges$type_g4



cluster<-levels(as.factor(my_granges$cluster))
a<-my_granges[(elementMetadata(my_granges)[, "cluster"] %in% cluster)]

final_final_df <-(mapply(final_df,cluster,SIMPLIFY = FALSE))  %>%   bind_rows()

saveRDS(final_final_df,"final_clusters.rds")
