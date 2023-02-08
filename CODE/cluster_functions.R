




level_1<-function(final_sub_cluster){
z<-("reached level_1")
print(paste(z,(final_sub_cluster)))

seq<-clus_1[clus_1$cluster %in% final_sub_cluster,]
new_seq<-DNAStringSet(seq$sequence)
names(new_seq)<-(seq$name)

metadata(new_seq)$cluster <- seq$cluster
#metadata(new_seq)$type <- seq$type_g4

seq_2_level1<-as.data.frame(seq)
seq_2_level1$cluster<-as.character(final_sub_cluster)
seq_2_level1$subcluster_level_2<-as.character(final_sub_cluster)

#print(colnames(seq_2_level1))
return(as.data.frame(seq_2_level1))
}


level_3<-function(sub_cluster_large){
#z<-("working for level 3")
print(paste(z,(sub_cluster_large)))
print("reached level_3")

seq<-clus_1[clus_1$cluster %in% sub_cluster_large,]
new_seq<-DNAStringSet(seq$sequence)
names(new_seq)<-(seq$name)

metadata(new_seq)$cluster <- seq$cluster
#metadata(new_seq)$type <- seq$type_g4

df<-g4_mat.2[g4_mat.2$name %in% seq$name]
name<-df$name
id<-df$id
df$id<-NULL
k<-as.numeric(ncol(df))


mat.3<-df%>% select(-c(name,sequence))

ind<- mat.3[, colSums(mat.3 != 0) > 0]
ind <- sapply(mat.3, function(x) sum(x==0)) != nrow(mat.3)

mat.4<-as.data.frame(mat.3)[,ind]
df$id<-rownames(mat.4)


hc<-hclust.vector(mat.4,method='ward',metric="euclidean")
hc$height<-round(hc$height,5)

mycl <- cutree(hc, k=2)

#plot(hc)
df1<-cbind(df[,c("name","sequence")],mycl)


df1$cluster<-paste0(seq$cluster,"-",mycl)
df1$mycl<-NULL

#print(colnames(seq_2_level1))

return(as.data.frame(df1))
}





#---------------------------------
level_2<-function(sub_cluster1){
z<-("reached level_2")
print(paste(z,(sub_cluster1)))
seq<-clus_1[clus_1$cluster %in% sub_cluster1,]
new_seq<-DNAStringSet(seq$sequence)
names(new_seq)<-(seq$name)


metadata(new_seq)$cluster <- seq$cluster
#metadata(new_seq)$type <- seq$type_g4

print("round 1, Distance matrix ")

########################
dat_frame<-as.data.frame(seq)

#here we edit the sequence to remove the front and end G ends to make it only one.
edit_dat_frame<-dat_frame
sequence<-dat_frame$sequence
sequence.1<-gsub('(^[G])\\1+','\\1',sequence)
sequence.2<-gsub('([G])\\1+$','\\1',sequence.1)
edit_dat_frame$sequence<-sequence.2

dist_mat<-edited_pairwise_dist(edit_dat_frame)


sil_width <- future_map_dbl(2:5,  function(k){
start.time <- Sys.time()
model<-pam(dist_mat,k = k,diss = TRUE)
model$silinfo$avg.width
save(model,file=paste0(model_Dir,"pairwise_Level2PAM_",k,".RData"))
end.time <- Sys.time()
time.taken <- end.time - start.time
print(time.taken)
})

# Generate a data frame containing both k and sil_width
sil_df <- data.frame(
  k = 2:5,
  sil_width = sil_width
)

optimal_cluster<-sil_df[which.max((sil_df$sil_width)),]$k
print_seq<-"the optimal cluster was found to be"
print_seq2<-"for cluster number"
#print(paste(print_seq,optimal_cluster,print_seq2,sub_cluster1))
model<-pam(dist_mat,k = optimal_cluster,diss = TRUE)
model$clustering<-paste0(as.character(sub_cluster1),"_",model$clustering,"b")

clus_2<-cbind(as.data.frame(seq$name),model$clustering)
#clus_2<-merge(clus_2,clus_1,by="name")
#colnames(clus_3)
index<-which(colnames(clus_2)=="model$clustering")
colnames(clus_2)[index]<-"subcluster_level_2"
colnames(clus_2)[1]<-"name"
clus_2<-as.data.frame(clus_2)
return(clus_2)
}




















euclidean <- function(a, b) sqrt(sum((a - b)^2))
 

add.uneven <- function(x, y) {
    l <- max(length(x), length(y))
    length(x) <- l
    length(y) <- l
    x[is.na(x)] <- 0
    y[is.na(y)] <- 0
    return(abs(x - y))
}
#x1<-"GGGGGTGGGTGGGTGAAGGGGG"
#x2<-"GGGTGGGTGCTTGGGATGGGG"
#x1<-"GGGAGGGGTGGGGCGGGGTGGGG"

g4_similarity<-function(a,b){
x1<-as.character(a)
x2<-as.character(b)
print(x1)
print(x2)
nquad.3<-lengths(regmatches(x1, gregexpr("GGG", x1)))
nquad.4<-lengths(regmatches(x1, gregexpr("GGGG",x1)))

nquad.5<-lengths(regmatches(x1, gregexpr("GGGGG", x1)))
nquad.6<-lengths(regmatches(x1, gregexpr("GGGGGG", x1)))
array.3<-length(unlist(strsplit(x1, "GGG")))
loops.3<-nchar(unlist(strsplit(x1, "GGG")))
g1<-lengths(regmatches(x1, gregexpr("-", x1)))
g2<-lengths(regmatches(x2, gregexpr("-", x1)))

array.4<-length(unlist(strsplit(x1, "GGGG")))
loops.4<-nchar(unlist(strsplit(x1, "GGGG")))

length(unlist(strsplit(x1, "GGGGG")))
loops.5<-nchar(unlist(strsplit(x1, "GGGGG")))

nquad.3.2<-lengths(regmatches(x2, gregexpr("GGG", x2)))
nquad.4.2<-lengths(regmatches(x2, gregexpr("GGGG", x2)))

nquad.5.2<-lengths(regmatches(x2, gregexpr("GGGGG", x2)))
nquad.6.2<-lengths(regmatches(x2, gregexpr("GGGGGG", x2)))
array.3.2<-length(unlist(strsplit(x2, "GGG")))
loops.3.2<-nchar(unlist(strsplit(x2, "GGG")))

#euclidean(

array.4.2<-length(unlist(strsplit(x2, "GGGG")))
loops.4.2<-nchar(unlist(strsplit(x2, "GGGG")))

length(unlist(strsplit(x2, "GGGGG")))
loops.5.2<-nchar(unlist(strsplit(x2, "GGGGG")))


##this term will be minimal if the sequence are similar

a1<-median(sum(add.uneven(loops.3,loops.3.2)),sum(add.uneven(loops.4,loops.4.2)),sum(add.uneven(loops.5,loops.5.2)))
b1<-median(abs(nquad.3-nquad.3.2),abs(nquad.4-nquad.4.2),abs(nquad.5-nquad.5.2),abs(nquad.6-nquad.6.2))

lx1<-nchar(x1)
lx2<-nchar(x2)
ldiff<-abs(lx1+g1-lx2+g2)
lsum<-abs(lx1-g1+lx2-g2)

#the logic behind ldiff/lsum is that for a specific length of sequences, i.e x1, and x2, 
#suppose sum of length of sequences is 100.With that in mind, the difference in the length
#of either sequence will only contribute to make the whole term larger while being normalized by the 
#total length of the sequence being compared.

g4sc<-((a1*b1)^2+1)*((1+ldiff)/(lsum+0.0001))
print("g4sc:")
print(g4sc)
return(g4sc)
}



###_------------------------------------
kmer_from_granges<-function(seq_2_level2,k){
print("reached kmer_from_granges")
x2<-seq_2_level2$sequence
name<-seq_2_level2$name
sequence.1<-gsub('(^[G])\\1+','\\1',x2)
x1<-gsub('([G])\\1+$','\\1',sequence.1)
#x1<-gsub("G{3}","U",x1)
#x1<-output_fwd.1$V4[1:20]

#new_seq<-as.character(x1)
new_seq<-DNAStringSet(x1)

names(new_seq)<-(name)

G4.kmers<-kmer_by_gene(new_seq,k)
#G4.kmers[,'GGGGG']<-G4.kmers[,'GGGGG']/str_count(new_seq)
G4.kmers[,-1]<-G4.kmers[,-1]+0.000000000000001
df <- normalize_function(G4.kmers[,-1])
row.names(df)<-(G4.kmers$name)
####
return(df)

}

#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#**#*#*#*#*
 normalize_function <- function(x){(x-min(x))/(max(x)-min(x))}

  #*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#
kmer_by_gene <- function(DNAStringSet, k){
   print("counting kmers...", quote = FALSE)
   kmer_counts <- Biostrings::oligonucleotideFrequency(DNAStringSet, width = k) %>% 
     dplyr::as_tibble() %>% 
     dplyr::mutate(name = names(DNAStringSet)) %>% 
     dplyr::select(name, dplyr::everything())
   
   print("counting complete.", quote= FALSE)
   
   return(kmer_counts)
}
#--------------------------------------------------














###make attributes for G4 clustering
g4_df<-function(x2){
nquad.3f<-list()
nquad.4f<-list()
nquad.5f<-list()
nquad.6f<-list()
loops.3f<-list()
loops.4f<-list()
loops.5f<-list()
array.3f<-list()
array.4f<-list()
array.5f<-list()
array.6f<-list()

for(i in 1:nrow(x2)) {       # for-loop over rows
 x1<-x2[i,"sequence"]
sequence.1<-gsub('(^[G])\\1+','\\1',x1)
x1<-gsub('([G])\\1+$','\\1',sequence.1)
 
nquad.3<-lengths(regmatches(x1, gregexpr("GGG", x1)))
nquad.4<-lengths(regmatches(x1, gregexpr("GGGG",x1)))

nquad.5<-lengths(regmatches(x1, gregexpr("GGGGG", x1)))
nquad.6<-lengths(regmatches(x1, gregexpr("GGGGGG", x1)))
array.3<-length(unlist(strsplit(x1, "GGG")))
loops.3<-list(nchar(unlist(strsplit(x1, "GGG"))))


array.4<-length(unlist(strsplit(x1, "GGGG")))
loops.4<-list(nchar(unlist(strsplit(x1, "GGGG"))))

array.5<-length(unlist(strsplit(x1, "GGGGG")))
loops.5<-list(nchar(unlist(strsplit(x1, "GGGGG"))))
array.6<-length(unlist(strsplit(x1, "GGGGGG")))
loops.6<-list(nchar(unlist(strsplit(x1, "GGGGGG"))))

#print(nquad.3)
#print(loops.3)
nquad.3f<-append(nquad.3f,nquad.3)
nquad.4f<-append(nquad.4f,nquad.4)
nquad.5f<-append(nquad.5f,nquad.5)
nquad.6f<-append(nquad.6f,nquad.6)
array.3f<-append(array.3f,array.3)
array.4f<-append(array.4f,array.4)
array.5f<-append(array.5f,array.5)
array.6f<-append(array.6f,array.6)

loops.3f<-append(loops.3f,(loops.3))

loops.4f<-append(loops.4f,(loops.4))
loops.5f<-append(loops.5f,(loops.5))
print(x1)
#print(loops.3)
#print(loops.3f)
}
x3 <- map_dfr(loops.3f, ~dplyr::as_data_frame(t(.)))
colnames(x3) <- paste("l3f",colnames(x3),sep="_")

x4 <- map_dfr(loops.4f, ~dplyr::as_data_frame(t(.)))
colnames(x4) <- paste("l4f",colnames(x4),sep="_")

x5 <- map_dfr(loops.5f, ~dplyr::as_data_frame(t(.)))
colnames(x5) <- paste("l5f",colnames(x5),sep="_")

print(x3)
print(x2)
x6<-as.data.frame(do.call(rbind, nquad.3f))
x7<-as.data.frame(do.call(rbind, nquad.4f))
x8<-as.data.frame(do.call(rbind, nquad.5f))
x9<-as.data.frame(do.call(rbind, nquad.6f))


x10<-as.data.frame(do.call(rbind, array.3f))
x11<-as.data.frame(do.call(rbind, array.4f))
x12<-as.data.frame(do.call(rbind, array.5f))
x13<-as.data.frame(do.call(rbind, array.6f))

#nquads<- (cbind(x6,x7,x8,x9,x10,x11,x12,x13,(x3),(x4),(x5)))
nquads<-as.data.frame(do.call(cbind,c(x6,x7,x8,x9,x10,x11,x12,x13,(x3),(x4),(x5))))
colnames(nquads)[1:8]<-c("nquad.3f","nquad.4f","nquad.5f","nquad.6f","nloops.3","nloops.4",
"nloops.5","nloops.6")

return(as.data.frame(nquads))
}



edited_pairwise_dist<-function(pairwise_dat_frame){
  
sub_mat<-nucleotideSubstitutionMatrix(match = 1, mismatch = 0, baseOnly = FALSE, type = "DNA")
sub_mat['A','A']<-5
sub_mat['G','G']<-3
sub_mat['C','C']<-5
sub_mat['T','T']<-5
sub_mat['A','G']<- -4
sub_mat['C','G']<--4
sub_mat['G','T']<--4
sub_mat['A','C']<--3
sub_mat['A','T']<--3
sub_mat['C','T']<--3

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
if (length(seq_2) > 500){
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
#my_comb <- combinations(5,2, 1:5,repeats.allowed = TRUE)
my_comb_name<-combinations(length(name),2, name,repeats.allowed = TRUE)
#my_comb_seq<-combinations(length(seq_2),2, seq_2,repeats.allowed = TRUE)
}
#as.data.frame(cbind(my_comb,my_comb_name))
#my_comb <- combn(x = 1:length(seq_2), m = 2)

my_comb_list <- lapply(apply(X = my_comb, MARGIN = 1, FUN = list), unlist)
length(my_comb_list)

# it should be list of pairs of seqs
#my_comb_list <- my_comb

my_seq<-c(pairwise_dat_frame$sequence)
my_name<-c(pairwise_dat_frame$name)

print("RUNNING BEFORE ALIGNMENT")
 align_wrapper <- function(x){
value<-  pairwiseAlignment(pattern = my_seq[x[1]], subject = my_seq[x[2]], type = "global", substitutionMatrix=sub_mat,gapOpening=-2, gapExtension=-3)

g4score<-g4_similarity(as.character(my_seq[x[1]]),as.character(my_seq[x[2]]))
score<-as.numeric(value@score)
score<-ifelse(score<=0,1,score)
#print("score")
#print(score)
f=g4score/score

#c0<-paste(as.character(my_name[x[1]]),"\t",as.character(my_name[x[2]]))
#c1<-paste(as.character(value@pattern),"\n",as.character(value@subject),"\n",f)

#catf(filename,'\n',c0,c1,'\n')
print(f)

print("final_score")
return(g4score/score)
}

system.time(
  my_alignments <- bplapply(my_comb_list, align_wrapper, BPPARAM =  SerialParam())
)

df<-as.data.frame((my_comb_name),row.names = NULL)
print("done alignment")
df1<-cbind(df,unlist((my_alignments)))
colnames(df1)<-c("From","To","Weight")
#df1$Weight<- ifelse(df1$Weight< 0, 0, df1$Weight)
mygraph <- graph.data.frame(df1)


pair_mat<-get.adjacency(mygraph,sparse=FALSE,attr='Weight',type="upper")
pair_mat<-Matrix::forceSymmetric(pair_mat)
pair_mat<-as.matrix(pair_mat)
#pair_mat<-t(apply(pair_mat, 1, function(x)(x-min(x))/(max(x)-min(x))))
pair_mat<-pair_mat


print("let's see")
if (isTRUE(all.equal(max(pair_mat) ,min(pair_mat)))==FALSE){
 pair_mat_norm<-t(apply(pair_mat, 1, function(x)(x-min(x))/(max(x)-min(x))))
 

pair_norm_dist<-(1-pair_mat_norm)
  
filename<-as.character(filename)
  print("1111")

  filename2=paste0(path_out_mat,"cluster_",filename,"sim_matrix.RData")
  filename1=paste0(path_out_mat,"cluster_",filename,"dist_normalized_matrix.RData")
#write.table(pair_mat,file=filename2,col.names=NA,append=FALSE)
 # write.table(pair_norm_dist,file=filename1,col.names=NA,append=FALSE)
save(pair_mat,file=filename2)
save(pair_norm_dist,file=filename1)
  print("let's see 3")

#sim2dist <- function(mx) as.dist(sqrt(outer(diag(mx), diag(mx), "+") - 2*mx))

#d.mx<-as.matrix(sim2dist(pair_mat))
return(pair_mat)
 

} else if ((isTRUE(all.equal(max(pair_mat) ,min(pair_mat)))==TRUE)){
  print("2222")
 pair_mat_norm<-(pair_mat)
pair_norm_dist<-(max(pair_mat_norm))-pair_mat_norm

  filename2=paste0(path_out_mat,"cluster_",filename,"sim_matrix.RData")
filename1=paste0(path_out_mat,"cluster_",filename,"dist_normalized_sim_matrix.RData")
  
#write.table(x=pair_mat,file=filename2,append=FALSE)
 # write.table(x=pair_norm_dist,file=filename1,append=FALSE)
save(pair_mat,file=filename2)
save(pair_norm_dist,file=filename1)


return(pair_mat)

  
}
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

