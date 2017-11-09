library(stringr)
library(gplots) 
library(ggplot2)
library(magrittr)
args <- commandArgs(TRUE)
file<-args[1]
st_file = args[2]
for_vplot<-read.csv(file,header=T)
st<-read.table(st_file,header=T)

#create new assembly column of for_vplot
tmp<-str_split_fixed(for_vplot$isolate, "_", 3)
for_vplot["assembly"]<-as.numeric(tmp[,3])

#Merge for_vplot with st info
merged<-merge(for_vplot,st,by.x="assembly",by.y="Assembly_ID",all.x=TRUE)


#trim merged
df<-merged[,colnames(merged) %in% c("k1_k2","insertion_sequence_length","MLST")]

#make violin plots stratified by mlst for each k1k2 pair
for (var in unique(df$k1_k2)) {
    pdf(paste(var,"_violin.pdf",sep=""))
    sub_mat<-df[df$k1_k2 == var,]
    p<-ggplot(sub_mat[sub_mat$MLST %in% c("1","3","110","3","28","42","54","8"),],aes(x=as.factor(MLST),y=insertion_sequence_length)) + geom_violin() + ggtitle(paste(var," insertion size distribtion\n per ST",sep=""))
   print(p)
   dev.off()
}

