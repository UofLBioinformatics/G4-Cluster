require(tidyr)
library(stringr)
library(tidyr)
library(stringr)
library(dplyr)
library(data.table)

library(dplyr)

parsed_result<-function(clstr.result){
clstr <- read.csv(clstr.result, sep = "\t", row.names = NULL, header = FALSE, stringsAsFactors = FALSE)

clstr2 <- clstr
n = nrow(clstr)
x = 0
numbers_only <- function(x) !grepl("\\D", x)
for (row in c(1:n)) {
  if (numbers_only(clstr2[row,1]) == TRUE) {
    clstr2[row,1] <- x}
  else {NULL}
  x <- clstr2[row,1]
}
#kable(head(clstr2))
#library(knitr)
clstr.sums <- data.frame(dplyr::count(clstr2,V1))
#kable(head(clstr.sums))

clstr_m3<- clstr.sums[(clstr.sums[,2]>4),]
#Yes, the switch properly took place. Now we need to get rid of rows that have empty V2
switch <- clstr.sums[1,2]
clstr3 <- cbind(clstr2[1], clstr)
#kable(clstr3[c((switch-2):(switch+2)),])

clstr4 <- clstr2[-which(clstr2$V2 == ""), ]
#kable(clstr4[c(1:5,(switch-2):(switch+2)),])

clstr5 <- clstr4
clstr5[] <- lapply(clstr5, gsub, pattern='>', replacement='')
clstr5.2 <- data.frame(str_split_fixed(clstr5$V2, "nt,", 2))
clstr5.3 <- data.frame(str_split_fixed(clstr5.2$X2, "... ", 2))
clstr6 <- cbind(clstr5[1],clstr5.2[1],clstr5.3[1:2])
colnames(clstr6) <- c("cluster","nt","sequence","stat")
clstr<-data.table(clstr6)
clstr_1<-setDT(clstr)[, if (.N > 4) .SD, by = cluster]
return(clstr_1)
}}