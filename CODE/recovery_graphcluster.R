
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

repeat_blastblat<-function(clust.g,ax_cl,times,train){
	new_train<-clust.g[(clust.g$clust.louisvain %in% ax_cl),]
	outfile_blast<-paste0("G4_seq_blast_",times,".tsv")
	fwrite(new_train,outfile_blast,sep="\t",col.names=TRUE,quote=FALSE)
	fasta_file_name<-paste0("g4_training_",times,".fasta")
	writeFasta(new_train[,c("id","sequence")],fasta_file_name)


	if (latermatrix=="blast"){
		fasta_file_path<-paste0("./",fasta_file_name)
		outfile_path<-paste0("./g4_training_blast_results_",times)
		db_training<-paste0("./blastdb_training_",times)
		temp_2<-blast_evalue(train,fasta_file_path,outfile_path,db_training)
	} else if (latermatrix=="blat"){
		fasta_file_path<-paste0("./",fasta_file_name)
		outfile_path<-paste0("./g4_training_blat_results_",times)
		temp_2<-blat_evalue(train,fasta_file_path,outfile_path)

	}
	if (length(unique(clust.g$clust.louisvain))>1){
		cl_1<-table(clust.g$clust.louisvain)%>% as.data.frame()%>% arrange(desc(Freq))%>% dplyr::filter(Freq>5)%>%tail(-1)
		rem_train<-clust.g[(clust.g$clust.louisvain%in% cl_1$Var1),]
		rem_train$clust.louisvain<-paste0("bl_rem_",rem_train$clust.louisvain)
	}


		clust_2<-graph_cluster(temp_2)
		node.degree2<-clust_2[["node.degree"]]
		clust_2<-clust_2[["clust"]]

		# Get cluster highest degree(top nodes) 
		clust.repr = merge(node.degree2, clust_2,by="id") %>%
  			group_by(clust.louisvain) %>%
  			top_n(1,degree) %>%
  			ungroup()

		######################### Save clusters ###################################
		message("# Saving sequence cluster mapping and cluster representatives...")

		clust.size = clust_2 %>% as.data.frame%>%
  			group_by(clust.louisvain) %>%
  			summarise(count = length(id), id = min(id)) %>%
  			ungroup()

		clust.size%>% arrange(desc(count))%>% head()


		clust.g2<-merge(clust_2,train)
		 table(clust.g2$clust.louisvain)%>% as.data.frame()%>% arrange(desc(Freq))%>% head()

		clust.g2$clust.louisvain<-paste0("bl_2_",clust.g2$clust.louisvain)
		if (exists("rem_train")){
			clust.gt<-rbind(clust.g2,rem_train)
		}else {
			clust.gt<-clust.g2
		}

		# Look at the distribution of cluster sizes
		p = ggplot(clust.size, aes(x = count)) +
  		geom_histogram(binwidth = 1, boundary = 0, alpha = 0.5, color = "black") +
  		#geom_line(stat = "bin", binwidth = 1, boundary = 0, alpha = 0.5) +
  		#scale_y_log10() +
  		xlim(0, max(clust.size$count)) +
  		xlab("Cluster size") +
	  	ylab("Number of clusters")


		pdf(sprintf("%s_cluster-size.pdf", "a"), 5, 5)
		print(p)
		log = dev.off()
		
return(clust.gt)
	}
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

graph_cluster<-function(edgelist){
temp<-edgelist
temp1<-temp %>% dplyr::filter(From!=To)


med1<-median(temp1$weight)
m1<-(mean(temp1$weight))
q1<-summary(temp1$weight)[["1st Qu."]]
q3<-summary(temp1$weight)[["3rd Qu."]]


blast.filter<- temp1  %>% dplyr::filter(weight>q3) #mutate(weight=ifelse(weight < q3, 0.001, weight))
print(head(blast.filter))

# Construct a graph and do the clustering with connected components
graph.filter = graph_from_data_frame(
  blast.filter,
  directed = F
)

graph.cc = components(graph.filter, mode ="strong")
graph.louvain<- cluster_walktrap(graph.filter)
#print(modularity(graph.filter,graph.cc))

print(head(graph.louvain))

# Extract the degree of nodes
node.degree = data.frame(
  id = names(degree(graph.filter)),
  degree = degree(graph.filter)
)

# Convert clustering results into a dataframe
clust = data.frame(
  id = names(graph.cc$membership),
  clust.id = graph.cc$membership,
clust.louisvain = graph.louvain$membership,
  stringsAsFactors = F)

return(list(clust=as.data.frame(clust),node.degree=as.data.frame(node.degree)))
}

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


blast_evalue<-function(train,fasta_file_path,outfile_path,db_training){

#-------------------------------------

command1<-paste0("makeblastdb -in ",fasta_file_path," -title G4DNA -dbtype nucl -parse_seqids -out ",db_training)
system(command1)
command2<-paste0("/bio/home/goonerrn/g_quadruplex/annotation_database/G_quad_clustering/starcode/blast_everything/BLAST_g4.sh ",fasta_file_path, "  ",db_training,"  ",outfile_path)
system(command2)

#------------------------------------

#R
###### Blast result remove duplicates and with specific score threshold
#train<-fread("G4_seq_blast.tsv")
in_out<-paste0(outfile_path,"_a_blast")
data<-fread(file = in_out, sep = '\t', header = F)

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


data.4<-setDT(data.3)[, .SD[which.max(V12)], by=list(V1,V2)]


data.4<-data.2%>% mutate(finalscore=V12/10)

temp<-data.4%>% select(V1,V2,finalscore)%>% summarise(From=V1,To=V2,weight=finalscore)

return(temp)
}



blat_evalue<-function(train,fasta_file_path,outfile_path){
blatoptions<-"-stepSize=3 -minScore=14 minMatch=2 -oneOff=9  -minIdentity=50 -tileSize=7 -repMatch=2147483647 -fine "
command1<-paste0("blat ",blatoptions,fasta_file_path ,"  ",fasta_file_path, "   out=blast8 ",outfile_path)
system(command1)

data<-fread(file = outfile_path, sep = '\t', header = F)

#remove the gaps from BLAST result to select hits with runs at least 4 runs of "GGG"

data.1<-data
data.2<-left_join(data.1,train,by=c("V1"="id"))
data.2<-left_join(data.2,train,by=c("V2"="id"))
data.3<-data.2
data.4<-setDT(data.3)[, .SD[which.max(V12)], by=list(V1,V2)]
data.4<-data.4%>% mutate(finalscore=V11/10)

temp<-data.4%>% select(V1,V2,finalscore)%>% summarise(From=V1,To=V2,weight=finalscore)
return(temp)
}







edited_pairwise_dist<-function(pairwise_dat_frame,pairwise_dist_folder){
 print("reached pairwise distance")

sub_mat<-nucleotideSubstitutionMatrix(match = 15, mismatch = 2, baseOnly = TRUE, type = "DNA")
sub_mat['A','A']<-3
sub_mat['G','G']<-2
sub_mat['C','C']<-3
sub_mat['T','T']<-3
sub_mat['A','G']<- -4
sub_mat['C','G']<- -4
sub_mat['G','T']<- -3
sub_mat['A','C']<- -2
sub_mat['A','T']<-  -2
sub_mat['C','T']<- -2

sub_mat<-Matrix::forceSymmetric(sub_mat)
sub_mat<-as.matrix(sub_mat)
pairwise_dat_frame$sequence<-as.character(pairwise_dat_frame$sequence)

#my_param <- McParam(workers = 8)
filename<-pairwise_dat_frame$cluster[1] #selecting just one because all are same
seq_2<-c(as.character(pairwise_dat_frame$sequence))
pairwise_dat_frame$name<-c(as.character(pairwise_dat_frame$id))
name<-c(pairwise_dat_frame$name)
##my_comb <- combn(x = 1:length(seq_2), m = 2)
##my_comb_name <- combn(x = name, m = 2)
z<-"DISTANCE MATRIX FOR cluster"
print(paste(z,filename,sep=" "))
#my_comb_name<-combinations(length(name),2, name,repeats.allowed = TRUE)

print("GOT INSIDE DISTANCE MATRIX")
if (length(seq_2) >= 6){
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
colnames(df1)<-c("From","To","weight")
df1$weight<-as.numeric(as.character(df1$weight))
df1$weight<- ifelse(df1$weight< 0, 0, df1$weight)
#print(head(df1))
print("is infinite")
table(is.infinite(df1$weight))

mygraph <- graph.data.frame(df1)


pair_mat<-get.adjacency(mygraph,sparse=FALSE,attr='weight',type="upper")
#print(row.names(pair_mat))

A<-pair_mat
k=which(diag(A)!=0,arr.ind=TRUE)
A.1<-A[k,k]


A.1[A.1 < 0]<-0

pair_mat<-A.1
A <- A.1/rowMaxs(A.1)

A.1<-Matrix::forceSymmetric(A)
A<-1-A.1
pair_norm_dist<-A

#print(A)
#print(A.1)
print("let's see")
filename<-as.character(filename)
  print("1111")
if ((isTRUE(all.equal(max(pair_mat) ,min(pair_mat)))==FALSE)){
  filename2=paste0(pairwise_dist_folder,"c_",filename,"sim_mat.RData")
  filename1=paste0(pairwise_dist_folder,"c_",filename,"dist_norm_mat.RData")
  save(pair_mat,file=filename2)
  save(pair_norm_dist,file=filename1)
  print("let's see 3")
#rm(pair_mat,pair_norm_dist)

return(pair_norm_dist)

} else if ((isTRUE(all.equal(max(pair_mat) ,min(pair_mat)))==TRUE)){
  print("2222")
 pair_mat_norm<-(pair_mat)
pair_norm_dist<-(max(pair_mat_norm))-pair_mat_norm

  filename2=paste0(pairwise_dist_folder,"c_",filename,"sim_matequal.RData")
filename1=paste0(pairwise_dist_folder,"c_",filename,"dist_norm_matequal.RData")

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


#----------------------------------------------------------------------------------------
#----------------------------------------------------------------------------------------



#--------------------------------------------------------------------------------------------

hclustK<-function(x, k){
  x.hclust = hclust(as.dist(x),method = "ward")
  x.cluster = list(cluster=cutree(x.hclust,k=k))
  #plot(x.hclust, cex = 0.6, hang = -1,labels=FALSE)

  return(x.cluster)
}

agnesK<-function(x, k){
# x: Data matrix or frame, k: Number of clusters
  x.agnes = agnes(x=x,diss=TRUE,method="complete")
  nodes<-cutree(x.agnes,k=k)
  x.cluster = list(id=names(nodes),cluster=nodes)
  #plot(as.hclust(x.agnes), cex = 0.6, hang = -1,labels=FALSE)
  #plot(x.agnes, sub = paste("Agglomerative Coefficient = ",round(x.agnes$ac, digits = 2)))
  return(x.cluster)
}


diana.reformat<-function(x, k){
# x: Data matrix or frame, k: Number of clusters
  x.diana = diana(x,diss=TRUE)
  x.cluster = list(cluster=cutree(x.diana,k=k))
  #plot(as.hclust(x.diana), cex = 0.6, hang = -1,labels=FALSE)
  #plot(x.diana,sub = paste("Divisive Coefficient = ", round(x.diana$dc, digits = 2)))
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


second_clustering<-function(dat_frame,pairwise_dist_folder){

total.start.time<- Sys.time()
filename<-(dat_frame$cluster[1]) #selecting just one because all are same
#print(filename)
dist_mat<-edited_pairwise_dist(dat_frame,pairwise_dist_folder)


g4matrix<-as.matrix(dist_mat)


print("---------------PAIRWISE------------------------------------------------")

print("-------------------------------------------------------------------------")



calc_temp<-hclust(as.dist(g4matrix), method = "complete")
#pdf(file=paste0(filename,"_bct_pair.pdf"))
x<-best.cutree(calc_temp,loss=TRUE,min=2,max=1000,graph=FALSE)
#dev.off()
a<-x[which.max(x)]
nc<-as.numeric(names(a))
print(paste0("optimal clusters are: ", nc))
print("*****************************************************************")
#res<-list()

if (nrow(g4matrix)>1500 & nrow(g4matrix)<1501){
	res<-list()
	res_df_nc<-list()
	res_df_nc1<-list()
	res_df_nc2<-list()

	kmin<-ifelse((nc-5)<2,2,(nc-3))
	kmax<-ifelse(((nc+5)>=(nrow(g4matrix)-2)),(nrow(g4matrix)-3),(nc+5))
	kmin<- ifelse(((kmax-kmin)< 2),2,kmin)
	kmax<- ifelse(((kmax-kmin)< 2),(kmin+5),kmax)


	if ((kmax)>=(nrow(g4matrix)-1)){
		kmax<-nrow(g4matrix)-2
		print(kmax)
		print(kmin)
	}






#--------------------------------------------------------------------------------------------
optimal_clusters2<-function(method_type,kmax,kmin){
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

print(kmax)
print(kmin)


	method_type<-c( "ward.D", "ward.D2", "single", "complete", "average",
          "mcquitty", "median", "centroid")
	metric2<-list()
	metric<-list()
	metric<-list()



#plan(multisession, workers = 4)
metricofmetrics<-map(method_type,~optimal_clusters2(.x,kmax,kmin))
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

}
else {
nc<-nc}


#pdf(paste0(filename,"_hc_pairwise.pdf"))
clustering.diana<-diana.reformat(as.dist(g4matrix),k=nc)
clustering.agnes<-agnesK(as.dist(g4matrix),k=nc)
clustering.hclustK <-hclustK(as.dist(g4matrix),k=nc)
#dev.off()

print(paste0("Optimal clusters are ",nc))

agnesclustering<-as.data.frame(cbind(id=dimnames(g4matrix)[[1]],clustering.agnes$cluster))
dianaclustering<-as.data.frame(cbind(id=names(clustering.hclustK$cluster),clustering.hclustK$cluster))
hclustering<-as.data.frame(cbind(id=names(clustering.hclustK$cluster),clustering.hclustK$cluster))


colnames(agnesclustering)<-c("id","agnesclusters2")
colnames(dianaclustering)<-c("id","dianaclusters2")
colnames(hclustering)<-c("id","hclusters2")

final.agnes<-merge(dat_frame,agnesclustering,by="id")
final.diana<-merge(final.agnes,dianaclustering,by="id")
final_pairwise<-merge(final.diana,hclustering,by="id")

final_pair<-final_pairwise
colnames(final_pair)[4:6]<-c("agnes_clusters2","diana_clusters2","hclust_clusters2")



file <- paste0(filename,"_pairwise_clusters.RData")
save(final_pairwise,file=file)

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




optimal_clusters<-function(method_type,kmax,kmin){
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









nggg_function<-function(output_fwd.4,nggg,latermatrix){
#nggg<-nggg_g
output_fwd.1<-output_fwd.4

o<-output_fwd.1%>% group_by(sequence)%>%
summarise(name=paste0(unique(name),collapse="|"), ggg=paste0(unique(nggg),collapse="|"),realseq=paste0(unique(V4),collapse="|"),
length_trimmed=nchar(sequence)) %>% dplyr::filter(length_trimmed>5)%>% distinct()%>% as.data.frame()
min_len<-min(o$length_trimmed)

o$id<-paste0("g4_",rownames(o))

#pairwise_dist_folder<-mat_Dir



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

if (nggg=="4"){
partition_per=partition_per
}else if (nggg=="5" ){
partition_per=ifelse(((partition_per)>.1),partition_per,(partition_per*3))
}else if ( nggg=="6_edit"){
partition_per=ifelse(((partition_per)>.1),partition_per,(partition_per*3))
}

print(paste0("partition percent used is",partition_per))
set.seed(17) # Set Seed so that same sample can be reproduced in future also
# Now Selecting 75% of data as sample from total 'n' rows of the data
sample <- sample.int(n = nrow(o), size = floor(partition_per*nrow(o)), replace = F)

train <-  as.data.frame(o[sample, ])
test  <-  o[-sample, ]

fwrite(train,"G4_training.tsv",sep="\t",col.names=TRUE,quote=FALSE)
fwrite(test,"G4_test.tsv",sep="\t",col.names=TRUE,quote=FALSE)

#test and train not found
#test<-ReadFasta(testf)
#testf="g4_test.fasta"
#test$name<-test$name.y
#test$name.x<-NULL
#test$name.y<-NULL

print("working1")

add=list()
for (i in ((min_len-2):(min_len+5))){
add=append(add,(paste0(rep("G",i),collapse="")))
}


lateruse<-train$realseq

add<-as.data.frame(t(as.data.frame(add)))
row.names(add)<-NULL
add$V2<-"added"
add$V3<-0
add$realseq<-add$V1

add$V1<-as.character(add$V1)
add$V4<-nchar(as.character(add$V1))
add$V5<-"added"


sequence.1<-gsub('(^[G])\\1+','\\1',train$realseq)
train$realseq<-gsub('([G])\\1+$','\\1',sequence.1)
colnames(add)<-colnames(train)


train<-rbind(add,train)


x<-as.character(train[,4])
y<-train[,c(1,2,3,4,5,6)]

#o[,3]$length<-nchar(o[,1])
clusters<-list()

for (i in c(2:4)){

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

cls<-cls%>% dplyr::filter(V2>1)

cls$id<-row.names(cls)
if (nrow(cls)>0){
cls.1<-cls%>% separate_rows(V3)
}

cls.1$V1len<-nchar(cls.1$V1)

cls.1$V3len<-nchar(cls.1$V3)
table(cls.1$V1len==cls.1$V3len)
#cls.1<-cls.1%>% dplyr::filter(V1!=V3)


diff<-ifelse(((cls.1$V3len>17) & (cls.1$V1len>17)),4,3)
#diff<-ifelse(i<=5,2,3)
cls.2<-cls.1%>% mutate(div=ifelse(abs(V1len-V3len)<diff,"same","notsame"))
table(cls.2$div)

cls.2 %>% group_by(id) %>% summarise(n=n(),v1n=n_distinct(V1),v3n=n_distinct(V3))

cls.3<-cls.2%>% dplyr::filter(div=="same")%>% group_by(id)%>% summarise(n=n())%>% dplyr::filter(n>=5)

if (nrow(cls.3)>0){
count<-unique(cls.3$id)

withcluster<-map(count,~cls.2%>% dplyr::filter(id==.x & div=="same")%>% select(V3,V1)%>% unlist()%>%
unique()%>% as.data.frame())
#withcluster<-data.table(withcluster)
withcluster<-rbindlist(withcluster,idcol=TRUE)
colnames(withcluster)<-c("clus","realseq")
toremove<-withcluster$realseq%>% unlist()
clusters[[length(clusters)+1]] <- list(withcluster)
x<-x[!x %in% toremove]
y<-y[!x %in% toremove,]

} else {(count<-NULL)
withcluster<-data.frame(clus=NA,realseq=NA)

clusters[[length(clusters)+1]] <- list(withcluster)
}

}



clusters_1<-clusters%>% map(as.data.frame)%>% bind_rows(.id="source")
clusters_1%>% dplyr::filter(!is.na(clus))%>%group_by(source,clus)%>% summarise(n=n())%>% arrange(desc(n))

#clusters_1%>% dplyr::filter(source==3)%>% head()
g4name<-fread("name_seq_2.txt",sep="\t")


clusters_2<-left_join(g4name,clusters_1,by=c("realseq"="realseq"))%>% dplyr::filter(name!="added")
clusters_2$realseq<-lateruse

final.1<-clusters_2%>% dplyr::filter(!is.na(source))

seq_blast<-clusters_2%>% dplyr::filter(is.na(source))
seq_blast$source<-NULL
seq_blast$clus<-NULL
train.original<-train



fwrite(seq_blast,"G4_seq_blast.tsv",sep="\t",col.names=TRUE,quote=FALSE)
fwrite(final.1,"starcode_clusters.tsv",sep="\t",col.names=TRUE,quote=FALSE)

writeFasta(train[,c("id","sequence")],"g4_training_minusstarcode.fasta")
#writeFasta(test[,c("id","sequence")],"g4_test.fasta")


bf<-fread("starcode_clusters.tsv")
#----------------------------------------------------------




train<-seq_blast


fwrite(seq_blast,"G4_seq_blast.tsv",sep="\t",col.names=TRUE,quote=FALSE)
fwrite(final.1,"starcode_clusters.tsv",sep="\t",col.names=TRUE,quote=FALSE)

writeFasta(train[,c("id","sequence")],"g4_training_minusstarcode.fasta")
#writeFasta(test[,c("id","sequence")],"g4_test.fasta")
#writeFasta(test[,c("id","realseq")],"g4_test.fasta")


if (latermatrix=="custom"){
temp<-train%>% select(id,sequence)

matrix<-similarity_selfmade(temp,filename="distance")
}


#-------------------------------------------------------------------------------


if (latermatrix=="blast"){
	times<-1

	fasta_file_path<-"./g4_training_minusstarcode.fasta"
	outfile_path<-paste0("./g4_training_blast_results_",times)
	db_training<-paste0("./blastdb_training_",times)

	temp<-blast_evalue(train,fasta_file_path,outfile_path,db_training)

} else if (latermatrix=="blat"){
	times<-1

	fasta_file_path<-"./g4_training_minusstarcode.fasta"
	outfile_path<-paste0("./g4_training_blast_results",times)

	temp<-blat_evalue(train,fasta_file_path,outfile_path)
}

metrictouse<-"evalue"


#if (metrictouse=="evalue"){

	temp$ordered_vertices <- apply(temp[1:2], 1, FUN = function(x) paste0(sort(x), collapse = ""))
	temp%>% group_by(ordered_vertices)%>% summarise(n=n())%>% arrange(desc(n))%>% head()
	
	df <-temp %>%
	   group_by(ordered_vertices) %>% dplyr::filter(weight<1)%>% mutate(weight=-log10(weight))%>%
	  arrange((weight)) %>%
	  dplyr::slice(1) %>%
 	 ungroup()
	temp<-as.data.frame(df[,c(1,2,3)])

	temp1<-temp %>% dplyr::filter(From!=To)


	med1<-median(temp1$weight)
	m1<-(mean(temp1$weight))
	q1<-summary(temp1$weight)[["1st Qu."]]
	q3<-summary(temp1$weight)[["3rd Qu."]]


  	temp<- temp1  %>% dplyr::filter(weight>m1) #

	toplot.1<-temp$weight

	save(temp,file="evalue_edgelist.RData")
#---------------------------------------------------------------

	clust<-graph_cluster(temp)
	node.degree<-clust[["node.degree"]]
	clust<-clust[["clust"]]

	# Extract cluster representatives
	clust.repr = merge(node.degree, clust,by="id") %>%
  	group_by(clust.louisvain) %>%
  	top_n(1,degree) %>%
  	ungroup()

######################### Save clusters ###################################
	message("# Saving sequence cluster mapping and cluster representatives...")

	clust.size = clust %>% as.data.frame%>%
  	group_by(clust.louisvain) %>%
  	summarise(count = length(id), id = min(id)) %>%
  	ungroup()

	clust.size%>% arrange(desc(count))%>% head(20)


	clust.g<-merge(clust,train)

	clust.g%>% dplyr::filter(clust.louisvain==39)%>% head()

	 max_cl<-table(clust.g$clust.louisvain)%>% as.data.frame()%>% arrange(desc(Freq))%>% dplyr::filter(Freq>1000)%>% select(Var1)%>% mutate(Var1=as.character(Var1))%>% unlist()
	
	if (length(max_cl)>0){
	times<-times+1
	clust.gt<-purrr::map(max_cl,~repeat_blastblat(clust.g,ax_cl=.x,times,train) %>% as.data.table)	%>% rbindlist(idcol="source")	
}  else { clust.gt<-clust.g}

setDT(clust.gt)
clust.gt1<-clust.gt[ , if (.N >= 6) .SD, by = .(clust.louisvain)]

temp<-as.data.frame(df[,c(1,2,3)]) %>% dplyr::filter(From!=To)


med1<-median(temp$weight)
m1<-(mean(temp$weight))
q1<-summary(temp$weight)[["1st Qu."]]
q3<-summary(temp$weight)[["3rd Qu."]]


temp<- temp  %>% dplyr::filter(weight>q3)
blast.filter<-temp
nnodes<-length(unique(c(temp$From,temp$To)))
nodes_names<-unique(c(temp$From,temp$To))

gex <- graph.data.frame(temp, directed=FALSE)
matex<-as_adjacency_matrix(gex, attr="weight")

A<-matex
k=which(diag(A)!=0,arr.ind=TRUE)
A.1<-as.matrix(A[k,k])



##here we make symmetric matric selecting > higher of the two values. this is for similarity score.

countnodes<-nrow(A)
b = matrix(ifelse(A>t(A),A,t(A)),nrow=countnodes,ncol=countnodes)

rownames(b)<-rownames(A)
colnames(b)<-colnames(A)

 nnodes<-nrow(b)

#b<-b[rowSums(b==value) < (nnodes-2),colSums(b==value) < (nnodes-2) ]

print(b[1:10,1:10])
#a<-Matrix::forceSymmetric(a)
b<-as.matrix(b)
A.1<-1/b
value<-max(A.1[is.finite(A.1)])
A.1[!is.finite(A.1)] <- value*2

print(A.1[1:10,1:10])
print(paste0("Number of rows:",nrow(A.1)))

D<-graph.adjacency(A.1,weighted=TRUE,mode='undirected',diag=TRUE)
undirected_graph <- as.undirected(D,mode = "collapse",  edge.attr.comb = list("min"))

E<-get.data.frame(undirected_graph)

#D<-get.adjacency(undirected_graph,sparse=TRUE,attr='weight')


A.1<-as.matrix(A.1)


A<- A.1
g4_mat.4<-A



#end of if



rm(data.2,data.3,data.1,D)


pdf("hist_weight_transformed.pdf")
temp1<-temp%>% dplyr::filter(From!=To)
hist(temp1$weight,main="This is similarity",breaks=10,
xlab="E value scores",
col="darkmagenta",freq=TRUE)
E1<-E%>% dplyr::filter(from!=to)%>% dplyr::filter(as.numeric(weight)!=1/value)

hist(E1$weight,main="This is distance", breaks=50,
xlab="Distance scores transformed",col="magenta",freq=TRUE)

#plot(toplot.1, E$weight, main="Scatterplot of Scores transformed",
   #xlab="Similarity score", ylab="Distance score", pch=5)

dev.off()






best.cutree <-
  function(hc, min = 3, max = 100, loss = FALSE, graph = FALSE, ...) {
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
  x.hclust = hclust(as.dist(x),method = "ward.D2")
  nodes<-cutree(x.hclust,k=k)
  x.cluster<-data.frame(cbind(id=names(nodes),cluster=nodes))
  #plot(x.hclust, cex = 0.6, hang = -1,labels=FALSE)

  return(x.cluster)
}

agnes_K<-function(x, k){
  x.agnes = agnes(x,diss=TRUE,method="ward")
  nodes<-cutree(as.hclust(x.agnes),k=k)
  x.cluster<-data.frame(cbind(id=names(nodes),cluster=nodes))

 # plot(x.agnes, sub = paste("Agglomerative Coefficient = ",round(x.agnes$ac, digits = 2)))
 # plot(as.hclust(x.agnes), cex = 0.6, hang = -1,labels=FALSE)
  return(x.cluster)
}


diana_reformat<-function(x, k){
# x: Data matrix or frame, k: Number of clusters
  x.diana = diana(x=x,diss=TRUE)
  nodes<-cutree((x.diana),k=k)
  x.cluster<-data.frame(cbind(id=x.diana$order.lab,cluster=nodes))
  #plot(x.diana,sub = paste("Divisive Coefficient = ", round(x.diana$dc, digits = 2)))
   # plot(as.hclust(x.diana), cex = 0.6, hang = -1,labels=FALSE)

  return(x.cluster)
}



library(NbClust)
print(dput(g4_mat.4[1:10,1:10]))


clust.gt$clust.louisvain<-as.character(clust.gt$clust.louisvain)
final<-clust.gt

mds.plot<-as.data.frame(mds.coor)
mds.plot$id<-row.names(mds.plot)
mds.plot<-left_join(mds.plot,final,by="id")

perplexitytsne<-ifelse((nrow(g4_mat.4)<1000),12,31)

tSNE <- Rtsne(as.dist(g4_mat.4), is_distance = TRUE,theta=0.2,pca_center=FALSE,
verbose = TRUE, max_iter = 500,perplexity=perplexitytsne)
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

train<-train%>% dplyr::filter(id!="added")
#tricky part, because we have extended sequences, we replace those with normal G4 for pairwise
#so we remove sequence and replace with realseq in train, already done for final above
head(train)
head(final)
final.1<-left_join(train,final,by=c("name","sequence","id","length_trimmed","ggg","realseq"))%>% dplyr::filter(!is.na(clust.louisvain))
final.1$sequence<-final.1$realseq
final.1$realseq<-NULL

save(final.1,file = "BLAST_clusters.RData")##resultdf2,resultdf

print("about complete")
print("...")
print("plotting")

#print(head(final.1))
pdf("mds_plot_clusters1.pdf")
plotmds(mds.plot,column_cluster="clust.louisvain")
dev.off()
pdf("tsne_plot_clusters1.pdf")
plottSNE(final_TSNE,col="clust.louisvain")
dev.off()

pdf("dist_matrix_clusters1.pdf")
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


#------------------------------------------------------------------------------------------------------------------------------
#---------------------------------------------------------------
#---------------------------------------------------------------
#---------------------------------------------------------------
print("-------------------------------------------------------------------------")

#bf<-fread("starcode_clusters.tsv")

#load("../BLAST_clusters.RData")
pairwise_Dir <- paste0(SubDir,"SECOND_cluster_",added_name_for_output_file,"/")
pairwise_dist_folder <<- paste0(pairwise_Dir,"SECOND_cluster","/")

dir.create(file.path(pairwise_Dir), showWarnings = FALSE,recursive=TRUE)
dir.create(file.path(pairwise_dist_folder), showWarnings = FALSE,recursive=TRUE)

setwd(file.path(pairwise_Dir))



if (nggg=="6_edit"){
final.1%>% group_by(clust.louisvain)%>% summarise(n=n())%>% arrange(desc(n))%>% head()%>% print()
print(head(final.1))
}

clusters_topair<-final.1%>% group_by(clust.louisvain)%>% mutate(n=n())%>% dplyr::filter((n>=6) & (!is.na(clust.louisvain)))%>% select(clust.louisvain)%>%
 unlist() %>% unique()

pairseq<-final.1[final.1$clust.louisvain %in% clusters_topair]


bf$cluster<-paste0("sc_",bf$source,"_",bf$clus,"a")
pairseq$cluster<-paste0(pairseq$clust.louisvain ,"a")
print(head(bf))
bf$sequence<-bf$realseq
bf$realseq<-NULL
a<-bf[,c("sequence", "name" , "ggg", "length_trimmed","id","cluster")]
b<-pairseq[,c("sequence", "name" , "ggg", "length_trimmed","id","cluster")]

final<-rbind(a,b)

if (nrow(final) == 0) {
    return (NULL)
} else {
   
X1<-final%>% select("cluster","id","sequence")
print("X1")
print(head(X1))
colnames(X1)<-c("cluster","id","sequence")

X<-split(X1,X1$cluster)
#X<-purrr::imap(X1, ~mutate(.x, cluster = .y))

n<-length(X)
#print(head(X[[1]]))
print(paste0("length is",n))
print(head(final))
 


final_result<-lapply(1:n,function(i){

	 #tryCatch({
	dat_frame<-as.data.frame(X[[i]])

	result<-second_clustering(dat_frame,pairwise_dist_folder)
	result[] <- lapply(result, as.character)
	return(result)  #},error=function(e) return(NULL))
})



#final_result<-map(X,~second_clustering(.x,pairwise_dist_folder))
filename=paste0("final_result_",nggg,".RData")
save(final_result,file=filename)

return(final_result)
}
}

