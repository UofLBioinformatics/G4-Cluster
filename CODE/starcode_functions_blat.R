#!/bin/bash

#--------------------------------------------------------------------------------------------
#--------------------------------------------------------------------------------------------
#--------------------------------------------------------------------------------------------

writeFasta<-function(data, filename){
  fastaLines = c()
  for (rowNum in 1:nrow(data)){
    fastaLines = c(fastaLines, as.character(paste(">", data[rowNum,"id"], sep = "")))
    fastaLines = c(fastaLines,as.character(data[rowNum,"sequence"]))
  }
  fileConn<-file(filename)
  writeLines(fastaLines, fileConn)
  close(fileConn)
}



#--------------------------------------------------------------------------------------------
#--------------------------------------------------------------------------------------------
#--------------------------------------------------------------------------------------------

best.cutree <-
  function(hc, min = 3, max = 20, loss = FALSE, graph = FALSE, ...) {
    if (!inherits(hc, "hclust")) hc <- as.hclust(hc)
    max <- min(max, length(hc$height))
    inert.gain <- rev(hc$height)
    intra <- rev(cumsum(rev(inert.gain)))
    relative.loss <- intra[min:(max)] / intra[(min - 1):(max - 1)]
    best <- which.min(relative.loss)
    names(relative.loss) <- min:max
    if (graph) {
      temp <- relative.loss
      temp[best] <- NA
      best2 <- which.min(temp)
      pch <- rep(1, max - min + 1)
      pch[best] <- 16
      pch[best2] <- 21
      plot(min:max, relative.loss, pch = pch, bg = "grey60", ...)
    }
    if (loss) {
      relative.loss
    } else {
      best + min - 1
    }
  }


#--------------------------------------------------------------------------------------------
#--------------------------------------------------------------------------------------------
#--------------------------------------------------------------------------------------------
ReadFasta<-function(file) {
   # Read the file line by line
   fasta<-readLines(file)
   # Identify header lines
   ind<-grep(">", fasta)
   # Identify the sequence lines
   s<-data.frame(ind=ind, from=ind+1, to=c((ind-1)[-1], length(fasta)))
   # Process sequence lines
   seqs<-rep(NA, length(ind))
   for(i in 1:length(ind)) {
      seqs[i]<-paste(fasta[s$from[i]:s$to[i]], collapse="")
   }
   # Create a data frame 
   DF<-data.frame(name=gsub(">", "", fasta[ind]), sequence=seqs)
   # Return the data frame as a result object from the function
   return(DF)
}




#--------------------------------------------------------------------------------------------
#--------------------------------------------------------------------------------------------

#--------------------------------------------------------------------------------------------
#--------------------------------------------------------------------------------------------




edited_pairwise_dist<-function(pairwise_dat_frame){
 print("reached pairwise distance")
sub_mat<-nucleotideSubstitutionMatrix(match = 15, mismatch = 2, baseOnly = TRUE, type = "DNA")
sub_mat['A','A']<-15
sub_mat['G','G']<-7
sub_mat['C','C']<-15
sub_mat['T','T']<-15
sub_mat['A','G']<- -12
sub_mat['C','G']<- -12
sub_mat['G','T']<- -8
sub_mat['A','C']<- 4
sub_mat['A','T']<- 6
sub_mat['C','T']<- 5


sub_mat<-Matrix::forceSymmetric(sub_mat)
sub_mat<-as.matrix(sub_mat)

#my_param <- McParam(workers = 8)
filename<-pairwise_dat_frame$cluster[1] #selecting just one because all are same
seq_2<-c(pairwise_dat_frame$sequence)
pairwise_dat_frame$name<-pairwise_dat_frame$id
name<-c(pairwise_dat_frame$name)
##my_comb <- combn(x = 1:length(seq_2), m = 2)
##my_comb_name <- combn(x = name, m = 2)
z<-"DISTANCE MATRIX FOR cluster"
print(paste(z,filename,sep=" "))
#my_comb_name<-combinations(length(name),2, name,repeats.allowed = TRUE)

print("GOT INSIDE DISTANCE MATRIX")
if (length(seq_2) >= 15){
#t(combn(length(seq_2), 2))
##alternate way to make combinations for n greater than 1200

dta<-expand.grid(1:length(seq_2),1:length(seq_2))
#######dta<-expand.grid(1:5,1:5)
dta_self<-dta[which(dta[,1] == dta[,2]), ]
dta_two<-t(combn(length(seq_2), 2))
colnames(dta_two)<-c("one","two")
colnames(dta_self)<-c("one","two")


my_comb<-rbind(dta_two,dta_self)

dta<-expand.grid(name,name)
dta_self<-dta[which(dta[,1] == dta[,2]), ]
dta_two<-t(combn((name), 2))
colnames(dta_two)<-c("one","two")
colnames(dta_self)<-c("one","two")

my_comb_name<-rbind(dta_two,dta_self)

}else {
my_comb <- combinations(length(seq_2),2, 1:length(seq_2),repeats.allowed = TRUE)
my_comb_name<-combinations(length(name),2, name,repeats.allowed = TRUE)

#my_comb <- combinations(5,2, 1:5,repeats.allowed = TRUE)
#my_comb_seq<-combinations(length(seq_2),2, seq_2,repeats.allowed = TRUE)
}
#as.data.frame(cbind(my_comb,my_comb_name))
#my_comb <- combn(x = 1:length(seq_2), m = 2)

my_comb_list <- lapply(apply(X = my_comb, MARGIN = 1, FUN = list), unlist)
length(my_comb_list)
print(paste0("total sequences in this cluster is:",length(seq_2)))
print(paste0("total number of pairwise alignment should be ",length(my_comb_list)))
# it should be list of pairs of seqs
#my_comb_list <- my_comb

my_seq<-c(pairwise_dat_frame$sequence)
my_name<-c(pairwise_dat_frame$name)



align_wrapper <- function(x){
answer<-paste0(my_name[x[1]],"-",my_name[x[2]])
#print(answer)

value<-  pairwiseAlignment(pattern = my_seq[x[1]], subject = my_seq[x[2]], type = "global", substitutionMatrix=sub_mat,gapOpening=-5, gapExtension=-10)
#print(value)

#g4score<-g4_similarity(as.character(my_seq[x[1]]),as.character(my_seq[x[2]]))
score<-as.numeric(value@score)

return(as.data.frame(cbind(my_name[x[1]],my_name[x[2]],score)))
}


print("RUNNING BEFORE ALIGNMENT")
system.time(
  my_alignments <- lapply(my_comb_list, align_wrapper)#, BPPARAM =  SerialParam())
)

df1<-rbindlist(my_alignments)%>% as.data.frame()
colnames(df1)<-c("From","To","Weight")
df1$Weight<-as.numeric(as.character(df1$Weight))
df1$Weight<- ifelse(df1$Weight< 1, 0, df1$Weight)
#print(head(df1))
print("is infinite")
table(is.infinite(df1$Weight))

mygraph <- graph.data.frame(df1)


pair_mat<-get.adjacency(mygraph,sparse=FALSE,attr='Weight',type="upper")
#print(row.names(pair_mat))

A<-pair_mat
k=which(diag(A)!=0,arr.ind=TRUE)
A.1<-A[k,k]


A[A < 0]<-0

pair_mat<-A
A.1 <- A/rowMaxs(A)

A<-Matrix::forceSymmetric(A.1)
A.1<-1-A
pair_norm_dist<-A.1

#print(A)
#print(A.1)
print("let's see")
filename<-as.character(filename)
  print("1111")
if ((isTRUE(all.equal(max(pair_mat) ,min(pair_mat)))==FALSE)){
  filename2=paste0(pairwise_dist_folder,"cluster_",filename,"sim_matrix.RData")
  filename1=paste0(pairwise_dist_folder,"cluster_",filename,"dist_normalized_matrix.RData")
  save(pair_mat,file=filename2)
  save(pair_norm_dist,file=filename1)
  print("let's see 3")
#rm(pair_mat,pair_norm_dist)

return(pair_norm_dist)
 
} else if ((isTRUE(all.equal(max(pair_mat) ,min(pair_mat)))==TRUE)){
  print("2222")
 pair_mat_norm<-(pair_mat)
pair_norm_dist<-(max(pair_mat_norm))-pair_mat_norm

  filename2=paste0(pairwise_dist_folder,"cluster_",filename,"sim_matrixequal.RData")
filename1=paste0(pairwise_dist_folder,"cluster_",filename,"dist_normalized_matrixequal.RData")
  
#write.table(x=pair_mat,file=filename2,append=FALSE)
 # write.table(x=pair_norm_dist,file=filename1,append=FALSE)
save(pair_mat_norm,file=filename2)
save(pair_norm_dist,file=filename1)

return(pair_norm_dist)

}
gc()

}



#------------------------------------------------------------------
#------------------------------------------------------------------
#------------------------------------------------------------------


combinations<-function (n, r, v = 1:n, set = TRUE, repeats.allowed = FALSE) 
{
  if (mode(n) != "numeric" || length(n) != 1 || n < 1 || (n%%1) != 
      0) 
    stop("bad value of n")
  if (mode(r) != "numeric" || length(r) != 1 || r < 1 || (r%%1) != 
      0) 
    stop("bad value of r")
  if (!is.atomic(v) || length(v) < n) 
    stop("v is either non-atomic or too short")
  if ((r > n) & repeats.allowed == FALSE) 
    stop("r > n and repeats.allowed=FALSE")
  if (set) {
    v <- unique(sort(v))
    if (length(v) < n) 
      stop("too few different elements")
  }
  v0 <- vector(mode(v), 0)
  if (repeats.allowed) 
    sub <- function(n, r, v) {
      if (r == 0) 
        v0
      else if (r == 1) 
        matrix(v, n, 1)
      else if (n == 1) 
        matrix(v, 1, r)
      else rbind(cbind(v[1], Recall(n, r - 1, v)), Recall(n - 
                                                            1, r, v[-1]))
    }
  else sub <- function(n, r, v) {
    if (r == 0) 
      v0
    else if (r == 1) 
      matrix(v, n, 1)
    else if (r == n) 
      matrix(v, 1, n)
    else rbind(cbind(v[1], Recall(n - 1, r - 1, v[-1])), 
               Recall(n - 1, r, v[-1]))
  }
  sub(n, r, v[1:n])
}


#----------------------------------------------------------------------------------------
#----------------------------------------------------------------------------------------
#----------------------------------------------------------------------------------------

euclidean <- function(a, b) sqrt(sum((a - b)^2))
 

#----------------------------------------------------------------------------------------
add.uneven <- function(x, y) {
    l <- max(length(x), length(y))
    length(x) <- l
    length(y) <- l
    x[is.na(x)] <- 0
    y[is.na(y)] <- 0
    return(abs(x - y))
}


#----------------------------------------------------------------------------------------
#can remove
#library(stringi)
#x1<-"GGGAGGGGGGTGGGGTGGGG"
#x2<-"GGGTGGGTGGGTGGGG"
#x2<-"GGGGGGGGGGGATTTTTTTTGGGGGGGGGGGGG"

g4_similarity<-function(a,b){
x1<-as.character(a)
x2<-as.character(b)

nquad.3<-lengths(str_match_all(x1, "GGG"))
nquad.4<-lengths(str_match_all(x1, "GGGG"))
nquad.5<-lengths(str_match_all(x1, "GGGGGG*"))
#nquad.6<-lengths(str_match_all(x1, "GGGGGG"))
g1<-lengths( str_match_all(x1, "-"))
g2<-lengths( str_match_all(x1, "-"))
array.3<-length(unlist(strsplit(x1, "GGG")))
loops.3<-nchar(unlist(strsplit(x1, "GGG")))
array.4<-length(unlist(strsplit(x1, "GGGG")))
loops.4<-nchar(unlist(strsplit(x1, "GGGG")))

#length(unlist(strsplit(x1, "GGGGG")))
array.5<-length(unlist(strsplit(x1, "GGGGGG*")))
loops.5<-nchar(unlist(strsplit(x1, "GGGGGG*")))



nquad.3.2<-lengths(str_match_all(x2, "GGG"))

nquad.4.2<-lengths(str_match_all(x2, "GGGG"))


nquad.5.2<-lengths(str_match_all(x2, "GGGGGG*"))
#nquad.6.2<-lengths(str_match_all(x2, "GGGGG"))

array.3.2<-length(unlist(strsplit(x2, "GGG")))
loops.3.2<-nchar(unlist(strsplit(x2, "GGG")))

#euclidean(

array.4.2<-length(unlist(strsplit(x2, "GGGG")))
loops.4.2<-nchar(unlist(strsplit(x2, "GGGG")))

#length(unlist(strsplit(x2, "GGGGG")))
loops.5.2<-nchar(unlist(strsplit(x2, "GGGGGG*")))


##this term will be minimal if the sequence are similar

a1<-median(sum(add.uneven(loops.3,loops.3.2)),sum(add.uneven(loops.4,loops.4.2)),sum(add.uneven(loops.5,loops.5.2)))
b1<-mean(abs(nquad.3-nquad.3.2),abs(nquad.4-nquad.4.2),abs(nquad.5-nquad.5.2))

lx1<-nchar(x1)
lx2<-nchar(x2)
ldiff<-abs(lx1-g1-lx2-g2)
lsum<-abs(lx1-g1+lx2-g2)

#the logic behind ldiff/lsum is that for a specific length of sequences, i.e x1, and x2, 
#suppose sum of length of sequences is 100.With that in mind, the difference in the length
#of either sequence will only contribute to make the whole term larger while being normalized by the 
#total length of the sequence being compared.

g4sc<-((a1*b1*10)^2+1)*((1+ldiff)/(lsum+0.0001))

return(g4sc)
}

#----------------------------------------------------------------------------------------

#----------------------------------------------------------------------------------------
#----------------------------------------------------------------------------------------


second_clustering<-function(dat_frame){
filename<-dat_frame$cluster[1] #selecting just one because all are same
print(filename)
dist_mat<-edited_pairwise_dist(dat_frame)


g4matrix<-as.matrix(dist_mat)


print("---------------PAIRWISE------------------------------------------------")

print("-------------------------------------------------------------------------")

maxnc<-as.numeric(nrow(g4matrix)-1)
calc_temp<-fastcluster::hclust(as.dist(g4matrix), method = "complete")
pdf(paste0(filename,"_bestcuttree_pairwise.pdf"))
x<-best.cutree(calc_temp,min=2,max=maxnc,graph=TRUE)
dev.off()
a<-x[which.max(x)]
nc<-as.numeric(a)
print(paste0("optimal clusters are: ", nc))
print("*****************************************************************")

res<-list()
res_df_nc<-list()
res_df_nc1<-list()
res_df_nc2<-list()


kmin<-ifelse((nc-10)<2,2,(nc-10))
kmax<-ifelse(((nc+5)>=(nrow(g4matrix)-2)),(nrow(g4matrix)-2),(nc+5))
kmax<- ifelse(((kmax-kmin)< 2),(nc+10),kmax)
kmin<- ifelse(((kmax-kmin)< 2),2,kmin)

if ((nc<10) & (nrow(g4matrix)<10)){
kmin<-2
kmax<-nc+5
}
if ((nc+5)>=(nrow(g4matrix)-10)){
kmin<-nrow(g4matrix)-15
kmax<-nrow(g4matrix)-7
print(kmax)
print(kmin)
}
method_type<-c( "ward.D", "ward.D2", "single", "complete", "average",
          "mcquitty", "median", "centroid")
metric2<-list()
metric<-list()
metric<-list()


optimal_clusters2<-function(method_type){
metrics<-c("mcclain","silhouette","dunn","cindex")
res<-map(metrics,~NbClust(diss=as.dist(g4matrix), distance = NULL, min.nc=kmin, max.nc=kmax,
method = method_type,index=.x))

metric<-map(res,~.x$Best.nc)%>% bind_rows()%>% as.data.frame()
row.names(metric)<-c("mcclain","silhouette","dunn","cindex")
metric1<-map(res,~.x$All.index)%>% bind_rows()%>% as.data.frame()
row.names(metric1)<-c("mcclain","silhouette","dunn","cindex")
metric2<-map(res,~.x$Best.partition )%>% bind_rows()%>% as.data.frame()%>%t()

colnames(metric2)<-c("mcclain","silhouette","dunn","cindex")
if (!is.null(metric1) && nrow(metric1)>0) res_df_nc1[[length(res_df_nc1)+1]] <- metric1
if (!is.null(metric2) && nrow(metric2)>0) res_df_nc2[[length(res_df_nc2)+1]] <- metric2
if (!is.null(metric) && nrow(metric)>0) res_df_nc[[length(res_df_nc)+1]] <- metric

return(list(res_df_nc,res_df_nc1,res_df_nc2))
}

#--------------------------------------------------------------------------------------------
#--------------------------------------------------------------------------------------------

hclustK<-function(x, k){
  x.hclust = fastcluster::hclust(as.dist(x),method = "complete")
  x.cluster = list(cluster=cutree(x.hclust,k=k))
  plot(x.hclust, cex = 0.6, hang = -1,labels=FALSE)
	
 #plot(as.dendrogram(res.diana), cex = 0.6,  horiz = TRUE)
  return(x.cluster)
}

agnesK<-function(x, k){
# x: Data matrix or frame, k: Number of clusters
  x.agnes = agnes(x=x,diss=TRUE,method="ward")
  nodes<-cutree(as.hclust(x.agnes),k=k)
 x.cluster = list(cluster=nodes)
    
plot(as.hclust(x.agnes), cex = 0.6, hang = -1,labels=FALSE)
  plot(x.agnes, sub = paste("Agglomerative Coefficient = ",round(x.agnes$ac, digits = 2)))
  return(x.cluster)
}


diana.reformat<-function(x, k){
# x: Data matrix or frame, k: Number of clusters
  x.diana = diana(x,diss=TRUE)
  x.cluster = list(cluster=cutree(x.diana,k=k))
  plot(as.hclust(x.diana), cex = 0.6, hang = -1,labels=FALSE)
  plot(x.diana,sub = paste("Divisive Coefficient = ", round(x.diana$dc, digits = 2)))
  return(x.cluster)
}

#--------------------------------------------------------------------------------------------
#--------------------------------------------------------------------------------------------
Mode <- function(x) {
   if ( length(x) <= 2 ) return(x[1])
   if ( anyNA(x) ) x = x[!is.na(x)]
   ux <- unique(x)
   ux[which.max(tabulate(match(x, ux)))]
}


#plan(multisession, workers = 4)
metricofmetrics<-map(method_type,~optimal_clusters2(.x))
result<-list() 
result2<-list()
for (i in 1:length(metricofmetrics)){
 #result<-metricofmetrics[[i]][1]
 #if (!is.null(metricofmetrics[[i]][3]) && nrow(metricofmetrics[[i]][3])>0) 
result[[i]] <- as.data.frame(metricofmetrics[[i]][[1]])
result2[[i]] <- as.data.frame(metricofmetrics[[i]][[2]])

}


resultdf<-result%>% bind_rows(.id="source")

resultdf$metric<-row.names(resultdf)

resultdf$method<-map(method_type,~rep(.x,4))%>% unlist()

resultdf<-as.data.frame(resultdf)
resultdf$metric<-gsub("\\...*","",resultdf$metric)

resultdf2<-result2%>% bind_rows(.id="source")
resultdf2$metric<-row.names(resultdf2)
resultdf2$method<-map(method_type,~rep(.x,4))%>% unlist()
resultdf2<-as.data.frame(resultdf2)
resultdf2$metric<-gsub("\\...*","",resultdf2$metric)

print(resultdf2)
nc<-Mode(resultdf$Number_clusters)
#epsc<-kNNdist(as.dist(g4_mat.4),k=nc)
nc<-ifelse(is.null(nc),1,nc)
#pdf("knnplot.pdf")
#kNNdistplot(as.dist(g4_mat.4),k=nc)
#dev.off()

#dbscanresult<-clusterDBSCAN1(g4matrix,eps=1,minPts=2)
#print("DBSCAN result pairwise")
#print(dbscanresult)
#dbscan_result<-cbind(rownames(g4matrix),dbscanresult$cluster)
#colnames(dbscan_result)<-c("id","dbscanclusters")


pdf(paste0(filename,"_hierachical_clustering_pairwise.pdf"))
clustering.diana<-diana.reformat(as.dist(g4matrix),k=nc)
clustering.agnes<-agnesK(as.dist(g4matrix),k=nc)
clustering.hclustK <-hclustK(as.dist(g4matrix),k=nc)
dev.off()

print(paste0("Optimal clusters are ",nc))
agnesclustering<-as.data.frame(cbind(id=names(clustering.hclustK$cluster),clustering.agnes$cluster))
dianaclustering<-as.data.frame(cbind(id=names(clustering.hclustK$cluster),clustering.diana$cluster))
hclustering<-as.data.frame(cbind(id=names(clustering.hclustK$cluster),clustering.hclustK$cluster))


colnames(agnesclustering)<-c("id","agnesclusters2")
colnames(dianaclustering)<-c("id","dianaclusters2")
colnames(hclustering)<-c("id","hclusters2")

final.agnes<-merge(dat_frame,agnesclustering,by="id")
final.diana<-merge(final.agnes,dianaclustering,by="id")
final_pairwise<-merge(final.diana,hclustering,by="id")
#final_pair<-merge(final_pairwise,dbscan_result,by="id")

final_pair<-final_pairwise
colnames(final_pair)[4:6]<-c("agnes_clusters2","diana_clusters2","hclust_clusters2")



#final.2<-left_join(final.1,final_pair,by="id")
file <- paste0(filename,"_pairwise_clusters.RData")
save(final_pairwise,resultdf2,resultdf,file=file)

print("about complete")
print("...")
print("plotting")


total.end.time<- Sys.time()
time_taken<-(total.end.time-total.start.time) 
print(paste0("total time taken", time_taken))
print("completed")
return(final_pairwise)
}


#---------------------------------------------------------------
#---------------------------------------------------------------
#---------------------------------------------------------------




optimal_clusters<-function(method_type){
metrics<-c("mcclain","silhouette","dunn","cindex")#"frey"
res<-map(metrics,~NbClust(diss=as.dist(g4_mat.4), distance = NULL, min.nc=kmin, max.nc=kmax,
method = method_type,index=.x))

metric<-map(res,~.x$Best.nc)%>% bind_rows()%>% as.data.frame()
row.names(metric)<-c("mcclain","silhouette","dunn","cindex")
metric1<-map(res,~.x$All.index)%>% bind_rows()%>% as.data.frame()
row.names(metric1)<-c("mcclain","silhouette","dunn","cindex")
metric2<-map(res,~.x$Best.partition )%>% bind_rows()%>% as.data.frame()%>%t()
print(metric1)
print(metric)
colnames(metric2)<-c("mcclain","silhouette","dunn","cindex")
if (!is.null(metric1) && nrow(metric1)>0) res_df_nc1[[length(res_df_nc1)+1]] <- metric1
if (!is.null(metric2) && nrow(metric2)>0) res_df_nc2[[length(res_df_nc2)+1]] <- metric2
if (!is.null(metric) && nrow(metric)>0) res_df_nc[[length(res_df_nc)+1]] <- metric

return(list(res_df_nc,res_df_nc1,res_df_nc2))
}




nggg_function<-function(output_fwd.4,nggg,train,test,o,SubDir){

print("working1")
min_len<-min(o$length_trimmed)
add=list()
for (i in ((min_len-2):(min_len+5))){
add=append(add,(paste0(rep("G",i),collapse="")))
}


lateruse<-train$sequence

add<-as.data.frame(t(as.data.frame(add)))
row.names(add)<-NULL
add$V2<-"added"
add$V3<-"added"

add$V4<-0
add$V1<-as.character(add$V1)
add$V5<-nchar(as.character(add$V1))
add$V6<-"added"


sequence.1<-gsub('(^[G])\\1+','\\1',train$sequence)
train$sequence<-gsub('([G])\\1+$','\\1',sequence.1)
colnames(add)<-colnames(train)

print("working")
train<-rbind(add,train)


x<-train[,1]
y<-train[,1:ncol(train)]

#o[,3]$length<-nchar(o[,1])
clusters<-list()

for (i in c(1:3)){

filename.g4=paste0("G4_seq","_",i,".txt")
filename.name=paste0("name_seq","_",i,".txt")

fwrite(as.data.frame(x),file=filename.g4,sep="\t", col.names=FALSE,quote=FALSE)
fwrite(as.data.frame(y),file=filename.name,sep="\t", col.names=TRUE,quote=FALSE)

timestamp=system(" (date +%d_%m_%Y_%H_%M_%S)")
dvalue=i
file_out=paste0("Cluster.d",dvalue,"___",(timestamp))
code=paste0("/bio/home/goonerrn/tools/starcode-1.3/starcode  -t8 -d",dvalue," -s --input ", filename.g4, " -o ", file_out, " --print-clusters")
system(code)
print("working2")

cls<-fread(file_out)

cls<-cls%>% filter(V2>1)

cls$id<-row.names(cls)
if (nrow(cls)>0){
cls.1<-cls%>% separate_rows(V3)
}

cls.1$V1len<-nchar(cls.1$V1)

cls.1$V3len<-nchar(cls.1$V3)
table(cls.1$V1len==cls.1$V3len)
#cls.1<-cls.1%>% filter(V1!=V3)


diff<-ifelse(((cls.1$V3len>15) & (cls.1$V1len>15)),4,3)
#diff<-ifelse(i<=5,2,3)
cls.2<-cls.1%>% mutate(div=ifelse(abs(V1len-V3len)<diff,"same","notsame"))
table(cls.2$div)

cls.2 %>% group_by(id) %>% summarise(n=n(),v1n=n_distinct(V1),v3n=n_distinct(V3))

cls.3<-cls.2%>% filter(div=="same")%>% group_by(id)%>% mutate(n=n())%>% filter(n>=13)%>% ungroup()

cls.3$n<-NULL
if (nrow(cls.3)>0){
count<-unique(cls.3$id)

withcluster<-map(count,~cls.3%>% filter(id==.x & div=="same")%>% select(V3,V1)%>% unlist()%>% 
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



clusters_1<-clusters%>% map(as.data.frame)%>% bind_rows(.id="source")
clusters_1%>% filter(!is.na(clus))%>%group_by(source,clus)%>% summarise(n=n())%>% arrange(desc(n))

#clusters_1%>% filter(source==3)%>% head()
g4name<-fread("name_seq_1.txt",sep="\t")


clusters_2<-left_join(g4name,clusters_1)%>% filter(name!="added")
print(head(clusters_2))
clusters_2$sequence1<-lateruse
print(head(clusters_2))
clusters_2$sequence<-clusters_2$sequence1
clusters_2$sequence1<-NULL

final.1<-clusters_2%>% filter(!is.na(source))

seq_blast<-clusters_2%>% filter(is.na(source))
seq_blast$source<-NULL
seq_blast$clus<-NULL 
train.original<-train


train<-seq_blast

fwrite(seq_blast,"G4_seq_blast.tsv",sep="\t",col.names=TRUE,quote=FALSE)
fwrite(final.1,"starcode_clusters.tsv",sep="\t",col.names=TRUE,quote=FALSE)

writeFasta(train[,c("id","sequence")],"g4_training_minusstarcode.fasta")
#writeFasta(test[,c("id","sequence")],"g4_test.fasta")


bf<-fread("starcode_clusters.tsv")
#----------------------------------------------------------



#writeFasta(test[,c("id","sequence")],"g4_test.fasta")

if (latermatrix=="blast"){

command1<-paste0("makeblastdb -in ./g4_training_minusstarcode.fasta -title G4DNA -dbtype nucl -parse_seqids -out ./blastdb_training")
system(command1)
command2<-paste0("/bio/home/goonerrn/g_quadruplex/annotation_database/G_quad_clustering/starcode/blast_everything/BLAST_g4.sh ./g4_training_minusstarcode.fasta ./blastdb_training g4_training_blast_results")
system(command2)

#R
###### Blast result remove duplicates and with specific score threshold
#train<-fread("G4_seq_blast.tsv")

data<-fread(file = 'g4_training_blast_results_a_blast', sep = '\t', header = F)

#remove the gaps from BLAST result to select hits with runs at least 4 runs of "GGG"
data$V3_temp<-str_replace_all(data$V3, "-", "")
data$V4_temp<-str_replace_all(data$V4, "-", "")

data.1 <- data%>% rowwise()%>% mutate(number.of.GGG=min(str_count(V3_temp, "GGG"),str_count(V4_temp, "GGG")),
number.of.GGGG=min(str_count(V3_temp, "GGGG"),str_count(V4_temp, "GGGG")))

data.1$gggtouse<-ifelse(data.1$number.of.GGG >= data.1$number.of.GGGG,data.1$number.of.GGG,data.1$number.of.GGGG)

data.2<-left_join(data.1,train,by=c("V1"="id"))
data.2<-left_join(data.2,train,by=c("V2"="id"))

data.3 <- data.2%>% rowwise()%>% mutate(GGGref=max(str_count(sequence.x, "GGG"),str_count(sequence.y, "GGG")),
GGGGref=max(str_count(sequence.x, "GGGG"),str_count(sequence.y, "GGGG")))

data.3$GGGreftouse<-ifelse(data.3$GGGref >= data.3$GGGGref,data.3$GGGref ,data.3$GGGGref)

#((as.numeric(data.3$ggg.x)+as.numeric(data.3$ggg.y))/2))

head(data.3)
head(data.3)%>% as.data.frame()
data.3$finalscore=data.3$V10*(data.3$V5/100)*(data.3$gggtouse/data.3$GGGreftouse)*(1/data.3$V12)*
(nchar(data.3$V3_temp)*nchar(data.3$V4_temp))/(data.3$length_trimmed.x*data.3$length_trimmed.y)
summary(data.3$finalscore)

#*((data.3$length_trimmed.x*data.3$length_trimmed.y)))
#(1/(log(data.2$length_trimmed.x)*log(data.2$length_trimmed.y)))



data.4<-data.2%>% mutate(finalscore=V12)

temp<-data.4%>% select(V1,V2,finalscore)%>% summarise(From=V1,To=V2,Weight=finalscore)

} else if (latermatrix=="blat"){
blatoptions<-"-stepSize=3 -minScore=14 minMatch=2 -oneOff=9  -minIdentity=50 -tileSize=7 -repMatch=2147483647 -fine " 

command1<-paste0("blat ",blatoptions," ./g4_training_minusstarcode.fasta ./g4_training_minusstarcode.fasta  out=blast8 ./g4_training_blat_results")
system(command1)

data<-fread(file = 'g4_training_blat_results', sep = '\t', header = F)

#remove the gaps from BLAST result to select hits with runs at least 4 runs of "GGG"

data.1<-data
data.2<-left_join(data.1,train,by=c("V1"="id"))
data.2<-left_join(data.2,train,by=c("V2"="id"))
data.3<-data.2
data.4<-setDT(data.3)[, .SD[which.max(V12)], by=list(V1,V2)]

temp<-data.4%>% select(V1,V2,V11)%>% summarise(From=V1,To=V2,Weight=V11)
}

#-------------------------------------------------------------------------------



metrictouse<-"evalue"
if (metrictouse=="G4score"){

data.4<-setDT(data.3)[, .SD[which.max(finalscore)], by=list(V1,V2)]

temp<-data.4%>% select(V1,V2,finalscore)%>% summarise(From=V1,To=V2,Weight=finalscore)

temp$ordered_vertices <- apply(temp[1:2], 1, FUN = function(x) paste0(sort(x), collapse = ""))
temp%>% group_by(ordered_vertices)%>% summarise(n=n())%>% arrange(desc(n))%>% head()

df <-temp %>% 
   group_by(ordered_vertices) %>% 
  arrange(desc(Weight)) %>% 
  dplyr::slice(1) %>% 
  ungroup()
temp<-df[,c(1,2,3)]

#data.4<-data.4%>% filter(finalscore)


toplot.1<-temp$Weight


save(temp,file="blast_edgelist.RData")


gex <- graph.data.frame(temp, directed=FALSE)
matex<-as_adjacency_matrix(gex, attr="Weight")


A<-matex
k=which(diag(A)!=0,arr.ind=TRUE)
A.1<-as.matrix(A[k,k])

#A<-t(apply(A.1, 1, function(x) replace(x, x== 0, min(x[x>0],na.rm=TRUE))))
A<-as.matrix(A.1)

##here the pairwise matrix gets transposed...
#A@x <-  (A@x / rep.int(diag(A), diff(A@p)))

A.1 <-  A/rowMaxs(A)
##converting to distance matrix from similarity
A<-(1/A.1+0.000000001)

g4_mat.4<-as.matrix(A)

D<-graph.adjacency(g4_mat.4,weighted=TRUE)

undirected_graph <- as.undirected(D,
                                  mode = "collapse", edge.attr.comb = list("min"))

E<-get.data.frame(graph.adjacency(as.matrix(A),weighted=TRUE))

} else if (metrictouse=="evalue"){


temp$ordered_vertices <- apply(temp[1:2], 1, FUN = function(x) paste0(sort(x), collapse = ""))
temp%>% group_by(ordered_vertices)%>% summarise(n=n())%>% arrange(desc(n))%>% head()

df <-temp %>% 
   group_by(ordered_vertices) %>% filter(Weight<1)%>%
  arrange((Weight)) %>%
  dplyr::slice(1) %>% 
  ungroup()
temp<-as.data.frame(df[,c(1,2,3)])

toplot.1<-temp$Weight

save(temp,file="evalue_edgelist.RData")
#---------------------------------------------------------------

 nnodes<-length(unique(c(temp$From,temp$To)))
 nodes_names<-unique(c(temp$From,temp$To))


gex <- graph.data.frame(temp, directed=FALSE)
matex<-as_adjacency_matrix(gex, attr="Weight")

A<-matex
k=which(diag(A)!=0,arr.ind=TRUE)
A.1<-as.matrix(A[k,k])
A<- -log10(A.1)
print(A[1:10,1:10])
A[!is.finite(A)] <- 0.0001

#min(x[x>0],na.rm=TRUE)

#a


#A <- as.matrix(a)/rowMaxs(as.matrix(a))
##here we make symmetric matric selecting > higher of the two values. this is for similarity score.

countnodes<-nrow(A)
b = matrix(ifelse(A>t(A),A,t(A)),nrow=countnodes,ncol=countnodes)

rownames(b)<-rownames(A)
colnames(b)<-colnames(A)

b<-b[rowSums(b==0.0001 ) < (nnodes-2),colSums(b==0.0001) < (nnodes-2) ]

print(b[1:10,1:10])
#a<-Matrix::forceSymmetric(a)
b<-as.matrix(b)
A.1<-1/b

print(A.1[1:10,1:10])

#D<-graph.adjacency(A.1,weighted=TRUE,mode='undirected',diag=TRUE)

#undirected_graph <- as.undirected(D,mode = "collapse",
# 				edge.attr.comb = list("max"))

#D<-get.adjacency(undirected_graph,sparse=FALSE,attr='weight')

#a<-D[rowSums(D == 0) < (nnodes),colSums(D == 0) < (nnodes) ]


#A<-t(apply(A.1, 2, function(x) replace(x, x== 0, max(x[x>0],na.rm=TRUE))))
D<-graph.adjacency(A.1,weighted=TRUE,mode='undirected',diag=TRUE)
undirected_graph <- as.undirected(D,mode = "collapse",	edge.attr.comb = list("min"))

E<-get.data.frame(undirected_graph)

#D<-get.adjacency(undirected_graph,sparse=TRUE,attr='weight')


A.1<-as.matrix(A.1)
 

A<- A.1
g4_mat.4<-A
}
#end of if



rm(data.2,data.3,data.1,D)


pdf("hist_weight_transformed.pdf")
temp1<-temp%>% filter(From!=To)
hist(temp1$Weight,main="This is similarity",breaks=10,
xlab="E value scores",
col="darkmagenta",freq=TRUE)
E1<-E%>% filter(from!=to)%>% filter(as.numeric(weight)!=10000)

hist(E1$weight,main="This is distance", breaks=50,
xlab="Distance scores transformed",col="magenta",freq=TRUE)

#plot(toplot.1, E$weight, main="Scatterplot of Scores transformed",
   #xlab="Similarity score", ylab="Distance score", pch=5)

dev.off()


best.cutree <-
  function(hc, min = 3, max = 150, loss = FALSE, graph = TRUE, ...) {
    if (!inherits(hc, "hclust")) hc <- as.hclust(hc)
    max <- min(max, length(hc$height))
    inert.gain <- rev(hc$height)
    intra <- rev(cumsum(rev(inert.gain)))
    relative.loss <- intra[min:(max)] / intra[(min - 1):(max - 1)]
    best <- which.min(relative.loss)
    names(relative.loss) <- min:max
    if (graph) {
      temp <- relative.loss
      temp[best] <- NA
      best2 <- which.min(temp)
      pch <- rep(1, max - min + 1)
      pch[best] <- 16
      pch[best2] <- 21
      plot(min:max, relative.loss, pch = pch, bg = "grey60", ...)
    }
    if (loss) {
      relative.loss
    } else {
      best + min - 1
    }
  }

 



idx <- ! rowSums(g4_mat.4) <=0

g4matrix<-as.matrix(g4_mat.4)

d <- as.dist(g4_mat.4)

mds.coor <- cmdscale(d)

#pdf("mds_plot1.pdf")
#plot(mds.coor[,1], mds.coor[,2], type="n", xlab="", ylab="",xlim=c(-1,1),ylim=c(-1,1))
#text(jitter(mds.coor[,1]), jitter(mds.coor[,2]),
#     rownames(mds.coor), cex=0.8)
#abline(h=0,v=0,col="gray75")
#dev.off()

hclust_K<-function(x, k){
  x.hclust = fastcluster::hclust(as.dist(x),method = "ward.D2")
  nodes<-cutree(x.hclust,k=k)
x.cluster = list(names=names(nodes),cluster=nodes)

#    x.cluster<-data.frame(cbind(cluster=nodes),id=names(nodes))

  plot(x.hclust, cex = 0.6, hang = -1,labels=FALSE)
	
 #plot(as.dendrogram(x.hclust), cex = 0.6,  horiz = TRUE)
  return(x.cluster)
}

agnes_K<-function(x, k){
  x.agnes = agnes(x,diss=TRUE,method="complete")
  nodes<-cutree(as.hclust(x.agnes),k=k)
  x.cluster = list(names=names(nodes),cluster=nodes)

#  x.cluster<-data.frame(cbind(cluster=nodes),id=names(nodes))
plot(x.agnes, sub = paste("Agglomerative Coefficient = ",round(x.agnes$ac, digits = 2)))
  plot(as.hclust(x.agnes), cex = 0.6, hang = -1,labels=FALSE)
  return(x.cluster)
}


diana_reformat<-function(x, k){
# x: Data matrix or frame, k: Number of clusters
  x.diana = diana(x=x,diss=TRUE)
    nodes<-cutree(as.hclust(x.diana),k=k)
  x.cluster = list(names=nodes)
  #x.cluster<-as.data.frame(cbind(cluster=nodes),id=names(nodes))
  plot(x.diana,sub = paste("Divisive Coefficient = ", round(x.diana$dc, digits = 2)))
    plot(as.hclust(x.diana), cex = 0.6, hang = -1,labels=FALSE)
	rownames(x.cluster)<-NULL
  return(x.cluster)
}



library(NbClust)

print(dput(g4_mat.4[1:10,1:10]))

calc_temp<-fastcluster::hclust(as.dist(g4_mat.4), method = "ward.D2")
pdf("bestcuttree.pdf")
x<-best.cutree(calc_temp,loss=FALSE,min=2,max=100,graph=TRUE)
dev.off()
print("best cutree result")
a<-x[which.max(x)]
print(a)
nc<-as.numeric(a)
print(nc)
print("THe optimal clusters are")

pdf("hierachical_clustering.pdf")

clustering.diana<-diana_reformat(as.dist(g4_mat.4),k=nc)
clustering.agnes<-agnes_K(as.dist(g4_mat.4),k=nc)
clustering.hclustK <-hclust_K((g4_mat.4),k=nc)

dev.off()

print(paste0("Optimal clusters are ",nc))

print(rownames(g4_mat.4)[1:10])
print(clustering.agnes)
agnesclustering<-clustering.agnes

hclustering<-clustering.hclustK
print("two done")
#print(head(agnesclustering))
#print(head(clustering.diana))
namesclus<-as.matrix(rownames(as.matrix(g4_mat.4)),ncol=1)
print(length(namesclus))

agnesclustering<-cbind(id=namesclus,clustering.agnes$cluster)
print("agnes done")
dianaclustering<-cbind(id=names(clustering.hclustK$cluster),clustering.diana$names)
print("dianadone")
hclustering<-cbind(id=names(clustering.hclustK$cluster),clustering.hclustK$cluster)
print("hclust done so columns issue maybe")

colnames(agnesclustering)<-c("id","agnesclusters")
colnames(dianaclustering)<-c("id","dianaclusters")
colnames(hclustering)<-c("id","hclusters")

final.agnes<-merge(train,agnesclustering,by="id")
final.diana<-merge(final.agnes,dianaclustering,by="id")
final<-merge(final.diana,hclustering,by="id")
print("4  done")

colnames(final)[7:9]<-c("agnes_clusters","diana_clusters","hclust_clusters")
print("5  done")


#dimnames(g4_mat.4)[[1]]

#clustering<-clustermethod(g4_mat.4,type="agnes",numberCluster=nc)
#cols=sample(rainbow(10),nc,replace=T)

mds.plot<-as.data.frame(mds.coor)
mds.plot$id<-row.names(mds.plot)
mds.plot<-left_join(mds.plot,final,by="id")

tSNE <- Rtsne(as.dist(g4_mat.4), is_distance = TRUE, pca=FALSE, 
verbose = TRUE, max_iter = 1000,perplexity=35)
tsne_df <- data.frame(x = tSNE$Y[,1], y = tSNE$Y[,2], id=dimnames(g4_mat.4)[[1]])
final_TSNE<-merge(final,tsne_df,by="id")


plotmds<-function(mdsmatrix,column_cluster){
 column_cluster<- rlang::sym(column_cluster)
p<-ggplot(mdsmatrix)+geom_point(aes(x = V1, y =  V2, color = {{column_cluster}},alpha = 0.1)) + theme_bw()+ theme(legend.position = "none")
print(p)
}

plottSNE<-function(tSNEdf,column_cluster,col){
column_cluster<-rlang::sym(col)
p<-ggplot(tSNEdf) + geom_point(aes(x=x, y=y, color={{column_cluster}},alpha = 0.01))+theme_bw()+ theme(legend.position = "none")
print(p)
}

train<-train%>% filter(id!="added")
final.1<-left_join(train,final)

save(final.1,file = "BLAT_clusters.RData")##resultdf2,resultdf

print("about complete")
print("...")
print("plotting")

#print(head(final.1))
pdf("mds_plot_clusters1.pdf")
plotmds(mds.plot,column_cluster="agnes_clusters")
plotmds(mds.plot,column_cluster="diana_clusters")
plotmds(mds.plot,column_cluster="hclust_clusters")

plottSNE(final_TSNE,col="agnes_clusters")
plottSNE(final_TSNE,col="diana_clusters")
plottSNE(final_TSNE,col="hclust_clusters")
print(fviz_dist(as.dist(g4_mat.4),show_labels=FALSE))
dev.off()

total.end.time<- Sys.time()
time_taken<-(total.end.time-total.start.time) 
print(paste0("total time taken", time_taken))
print("completed")
#save(clustering,file = "BLAST_clusters_agnes.RData")



#source("/bio/home/goonerrn/g_quadruplex/annotation_database/G_quad_clustering/starcode/pairwise.R")
print(paste0(rep("_",100),collapse=""))
print(paste0(rep("_",100),collapse=""))
print(paste0(rep("_",100),collapse=""))
print(paste0(rep("_",100),collapse=""))
print(paste0(rep("_",100),collapse=""))

#-----------------------------------------------------------------------------------------------------------
#-----------------------------------------------------------------------------------------------------------
#-----------------------------------------------------------------------------------------------------------
#-----------------------------------------------------------------------------------------------------------
#-----------------------------------------------------------------------------------------------------------

#print(SubDir)
#------------------------------------------------------------------------------------------------------------------------------
#---------------------------------------------------------------
#---------------------------------------------------------------
#---------------------------------------------------------------
print("-------------------------------------------------------------------------")

#bf<-fread("starcode_clusters.tsv")
SubDir<-SubDir 
#load("../BLAST_clusters.RData")
pairwise_Dir <- paste0(SubDir,"SECOND_cluster_",added_name_for_output_file,"/")
pairwise_dist_folder <<- paste0(pairwise_Dir,"SECOND_cluster","/")

dir.create(file.path(pairwise_Dir), showWarnings = FALSE,recursive=TRUE)
dir.create(file.path(pairwise_dist_folder), showWarnings = FALSE,recursive=TRUE)

setwd(file.path(pairwise_Dir))



if (nggg=="6_edit"){
final.1%>% group_by(hclust_clusters)%>% summarise(n=n())%>% arrange(desc(n))%>% head()%>% print()
print(head(final.1))
}

clusters_topair<-final.1%>% group_by(hclust_clusters)%>% mutate(n=n())%>% filter((n>=15) & (!is.na(hclust_clusters)))%>% select(hclust_clusters)%>%
 unlist() %>% unique()

#sequences_nottopair<-final.1%>% group_by(diana_clusters)%>% mutate(n=n())%>% filter(!(n>1 & !is.na(diana_clusters)))%>% select(sequence)%>%
# unlist()%>% unique()
#notpairedseq<-final.1[final.1$sequence %in% sequences_nottopair]
#notpairedseq$cluster<-c("bl_singleton",pairseq$diana_clusters,"a")

#notpairedseq<-notpairedseq[,c("sequence", "name" , "ggg", "length_trimmed","id","cluster")]

pairseq<-final.1[final.1$hclust_clusters %in% clusters_topair]



bf$cluster<-paste0("sc_",bf$source,"_",bf$clus,"a")
pairseq$cluster<-paste0("bl_",pairseq$hclusters,"a")
#print(head(bf))

a<-bf[,c("sequence", "name" , "ggg", "length_trimmed","id","cluster")]
b<-pairseq[,c("sequence", "name" , "ggg", "length_trimmed","id","cluster")]

final<-rbind(a,b)

X1<-final%>% select("cluster","id","sequence")
print("X1")
#print(head(X1))

colnames(X1)<-c("cluster","id","sequence")
X<-split(X1,X1$cluster)
#X<-purrr::imap(X1, ~mutate(.x, cluster = .y))

n<-length(X)
#print(head(X[[1]]))
print(paste0("length is",n))
#print(head(final))
final_result<-lapply(1:n,function(i){

 #tryCatch({
dat_frame<-as.data.frame(X[[i]])

result<-second_clustering(dat_frame)

return(result)  #},error=function(e) return(NULL))
})

#final_result<-map(X,~second_clustering(.x))

filename=paste0("final_result_",nggg,".RData")

save(final_result,file=filename)


return(final_result)
}


