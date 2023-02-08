


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
      plot(min:max, relative.loss, pch = pch, bg = "grey75", ...)
    }
    if (loss) {
      relative.loss
    } else {
      best + min - 1
    }
  }


edited_pairwise_dist<-function(pairwise_dat_frame){
  
sub_mat<-nucleotideSubstitutionMatrix(match = 1, mismatch = 0, baseOnly = FALSE, type = "DNA")
sub_mat['A','A']<-15
sub_mat['G','G']<-15
sub_mat['C','C']<-15
sub_mat['T','T']<-15
sub_mat['A','G']<- -8
sub_mat['C','G']<- -8
sub_mat['G','T']<- -4
sub_mat['A','C']<- -6
sub_mat['A','T']<- -6
sub_mat['C','T']<- -6

sub_mat<-Matrix::forceSymmetric(sub_mat)
sub_mat<-as.matrix(sub_mat)
#my_param <- McParam(workers = 8)
filename<-pairwise_dat_frame$cluster[1] #selecting just one because all are same
seq_2<-c(pairwise_dat_frame$sequence)
name<-c(pairwise_dat_frame$name)
##my_comb <- combn(x = 1:length(seq_2), m = 2)
##my_comb_name <- combn(x = name, m = 2)
z<-"DISTANCE MATRIX FOR cluster"
print(paste(z,filename,sep=" "))
#my_comb_name<-combinations(length(name),2, name,repeats.allowed = TRUE)

print("GOT INSIDE DISTANCE MATRIX")
if (length(seq_2) >= 20){
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
#answer<-paste0(my_seq[x[1]],"-",my_seq[x[2]]))
#print(answer)

value<-  pairwiseAlignment(pattern = my_seq[x[1]], subject = my_seq[x[2]], type = "global", substitutionMatrix=sub_mat,gapOpening=-7, gapExtension=-15)
print(value)
#g4score<-g4_similarity(as.character(my_seq[x[1]]),as.character(my_seq[x[2]]))
score<-as.numeric(value@score)
#score[!is.finite(score)] <- 1
#f=g4score/score
#f=1/score


#return(g4score/score)

#i=i+1
return(score)
}


print("RUNNING BEFORE ALIGNMENT")
system.time(
  my_alignments <- lapply(my_comb_list, align_wrapper)#, BPPARAM =  SerialParam())
)


df<-as.data.frame((my_comb_name),row.names = NULL)
print("done alignment")
df1<-cbind(df,unlist((my_alignments)))
colnames(df1)<-c("From","To","Weight")
df1$Weight<- ifelse(df1$Weight<= 0, 0.1, df1$Weight)
print(head(df1))
print("is infinite")
table(is.infinite(df1$Weight))

mygraph <- graph.data.frame(df1)


pair_mat<-get.adjacency(mygraph,sparse=FALSE,attr='Weight',type="upper")

pair_mat<-Matrix::forceSymmetric(pair_mat)

A<-pair_mat
k=which(diag(A)!=0,arr.ind=TRUE)
A.1<-A[k,k]

A<-as.matrix(t(apply(A.1, 1, function(x) replace(x, x <= 0, min(x[x>0],na.rm=TRUE)))))
pair_mat<-A
A.1 <-  A/rowMaxs(A)
A<- -log(A.1)
pair_norm_dist<-A
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










 
 
 
 safeBPParam <- function(nworkers) {
    if (.Platform$OS.type=="windows") {
        BiocParallel::SerialParam()
    } else {
        BiocParallel::MulticoreParam(nworkers)
    }
 }





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






second_clustering<-function(dat_frame){
filename<-dat_frame$cluster[1] #selecting just one because all are same

dist_mat<-edited_pairwise_dist(dat_frame)


g4matrix<-as.matrix(dist_mat)
################---------part2


print("---------------PAIRWISE------------------------------------------------")

print("-------------------------------------------------------------------------")



calc_temp<-hclust(as.dist(g4matrix), method = "average")
pdf(paste0(filename,"_bestcuttree_pairwise.pdf"))
x<-best.cutree(calc_temp,loss=TRUE,min=2,max=1000,graph=TRUE)
dev.off()
a<-x[which.max(x)]
nc<-as.numeric(names(a))
print(paste0("optimal clusters are: ", nc))
print("*****************************************************************")

res<-list()
res_df_nc<-list()
res_df_nc1<-list()
res_df_nc2<-list()
kmin<-ifelse((nc-10)<2,2,(nc-10))
kmax<-ifelse(((nc+5)>(nrow(g4matrix)-2)),(nrow(g4matrix)-2),(nc+5))
method_type<-c( "ward.D", "ward.D2", "single", "complete", "average",
          "mcquitty", "median", "centroid")
metric2<-list()
metric<-list()
metric<-list()


optimal_clusters2<-function(method_type){
metrics<-c("frey","mcclain","silhouette","dunn","cindex")
res<-map(metrics,~NbClust(diss=as.dist(g4matrix), distance = NULL, min.nc=kmin, max.nc=kmax,
method = method_type,index=.x))

metric<-map(res,~.x$Best.nc)%>% bind_rows()%>% as.data.frame()
row.names(metric)<-c("frey","mcclain","silhouette","dunn","cindex")
metric1<-map(res,~.x$All.index)%>% bind_rows()%>% as.data.frame()
row.names(metric1)<-c("frey","mcclain","silhouette","dunn","cindex")
metric2<-map(res,~.x$Best.partition )%>% bind_rows()%>% as.data.frame()%>%t()

colnames(metric2)<-c("frey","mcclain","silhouette","dunn","cindex")
if (!is.null(metric1) && nrow(metric1)>0) res_df_nc1[[length(res_df_nc1)+1]] <- metric1
if (!is.null(metric2) && nrow(metric2)>0) res_df_nc2[[length(res_df_nc2)+1]] <- metric2
if (!is.null(metric) && nrow(metric)>0) res_df_nc[[length(res_df_nc)+1]] <- metric

return(list(res_df_nc,res_df_nc1,res_df_nc2))
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

resultdf$method<-map(method_type,~rep(.x,5))%>% unlist()

resultdf<-as.data.frame(resultdf)
resultdf$metric<-gsub("\\...*","",resultdf$metric)

resultdf2<-result2%>% bind_rows(.id="source")
resultdf2$metric<-row.names(resultdf2)
resultdf2$method<-map(method_type,~rep(.x,5))%>% unlist()
resultdf2<-as.data.frame(resultdf2)
resultdf2$metric<-gsub("\\...*","",resultdf2$metric)

print(resultdf2)
nc<-Mode(resultdf$Number_clusters)
#epsc<-kNNdist(as.dist(g4_mat.4),k=nc)

#pdf("knnplot.pdf")
#kNNdistplot(as.dist(g4_mat.4),k=nc)
#dev.off()

dbscanresult<-clusterDBSCAN(g4matrix,eps=1,minPts=1)
print("DBSCAN result pairwise")
print(dbscanresult)
dbscan_result<-cbind(rownames(g4matrix),dbscanresult$cluster)
colnames(dbscan_result)<-c("id","dbscanclusters")


pdf(paste0(filename,"_hierachical_clustering_pairwise.pdf"))
clustering.diana<-diana.reformat(as.dist(g4matrix),k=nc)
clustering.agnes<-agnesK(as.dist(g4matrix),k=nc)
clustering.hclustK <-hclustK(as.dist(g4matrix),k=nc)
dev.off()

print(paste0("Optimal clusters are ",nc))
agnesclustering<-cbind(names(clustering.agnes$cluster),clustering.agnes$cluster)
dianaclustering<-cbind(names(clustering.hclustK$cluster),clustering.hclustK$cluster)
hclustering<-cbind(names(clustering.hclustK$cluster),clustering.hclustK$cluster)

colnames(agnesclustering)<-c("id","agnesclusters2")
colnames(dianaclustering)<-c("id","dianaclusters2")
colnames(hclustering)<-c("id","hclusters2")

final.agnes<-merge(final,agnesclustering,by="id")
final.diana<-merge(final.agnes,dianaclustering,by="id")
final_pairwise<-merge(final.diana,hclustering,by="id")
final_pair<-merge(final_pairwise,dbscan_result,by="id")


colnames(final_pair)[6:9]<-c("agnes_clusters2","diana_clusters2","hclust_clusters2","dbscan_result2")



final.2<-left_join(final.1,final_pair,by="id")
file <- paste0(filename,"_BLAST_clusters.RData")
save(final.2,resultdf2,resultdf,file=file)

print("about complete")
print("...")
print("plotting")


total.end.time<- Sys.time()
time_taken<-(total.end.time-total.start.time) 
print(paste0("total time taken", time_taken))
print("completed")
return(final.2)
}




#---------------------------------------------------------------
plan(multisession)

#cores=detectCores()
#cl <- makeCluster(cores[1]/2) #not to overload
#registerDoParallel(cl)
bf<-fread("starcode_clusters.tsv")

#load("../BLAST_clusters.RData")
pairwise_Dir <- paste0(SubDir,"SECOND_cluster_",added_name_for_output_file,"/")
dir.create(file.path(pairwise_Dir), showWarnings = FALSE,recursive=TRUE)
pairwise_dist_folder <- paste0(pairwise_Dir,"SECOND_cluster","/")

dir.create(file.path(pairwise_dist_folder), showWarnings = FALSE,recursive=TRUE)

setwd(file.path(pairwise_Dir))


clusters_topair<-final.1%>% group_by(diana_clusters)%>% mutate(n=n())%>% filter(n>19 & !is.na(diana_clusters))%>% select(diana_clusters)%>%
 unlist()%>% unique()

pairseq<-final.1[final.1$diana_clusters %in% clusters_topair]


bf$cluster<-paste0("sc_",bf$source,"_",bf$clus,"a")
pairseq$cluster<-paste0("bl_",pairseq$diana_clusters,"a")
a<-bf[,c("sequence", "name" , "ggg", "length_trimmed","id","cluster")]
b<-pairseq[,c("sequence", "name" , "ggg", "length_trimmed","id","cluster")]
final<-rbind(a,b)
X1<-final%>% select("cluster","id","sequence")
colnames(X1)<-c("cluster","name","sequence")
X<-split(X1,X1$cluster)
#X<-purrr::imap(X1, ~mutate(.x, cluster = .y))

final_result<-map(X,~second_clustering(dat_frame))

save(final_result,file="final_result.RData")









