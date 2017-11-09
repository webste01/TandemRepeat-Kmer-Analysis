library(grid)
library(ggplot2)
library(reshape2)
library(ggplot2)
library(magrittr)
library(gplots)
library(RColorBrewer)
library(plyr)

args <- commandArgs(TRUE)
fn = args[1]
kmers = args[2]
kmers_of_interest <- read.table(kmers,header=F)
STs<-read.table("../mlst_file.txt",header=T)

ST_freq<-data.frame(table(STs$MLST))
colnames(ST_freq)<-c("MLST","Freq")

tab<-read.csv(fn,header=T)
tmp<-read.table(text = as.character(tab$isolate), sep = "_", colClasses = "character")
tmp$V3<-as.numeric(tmp$V3)
tab["isolate"]<-tmp$V3

tab_ST<-merge(tab,STs,by.x="isolate",by.y="Assembly_ID",all.x=T)
tab_ST$MLST<-as.character(tab_ST$MLST)
y <- which(is.na(tab_ST$MLST)==TRUE)
tab_ST$MLST[y] <- 0

kmer_tab_sub<-subset(tab_ST,tab_ST$k1_k2 %in% kmers_of_interest$V1)
deduped.data <- unique(kmer_tab_sub)
kmer_tab_sub<-deduped.data[ , -which(names(deduped.data) %in% c("k1_start","k2_end","insertion_sequenceuence_length","orientation","isolate"))]



for (kmer in unique(kmer_tab_sub$k1_k2)){

#Get counts of each MLST type for each kmer
kmer_tab_sub2<-kmer_tab_sub[kmer_tab_sub$k1_k2 == kmer,]

CT<-dcast(kmer_tab_sub2,k1_k2*insertion_sequence~MLST,fill=0,length)

#Write out mapping of insertion sequence to allele ID
write.table(CT[,c("insertion_sequence")],paste(kmer,"_allele_mapping.csv",sep=""),row.names=T,quote=F)

#Remove insertion sequence and columns for mlst "0"
CT<-CT[,-c(1,2,3)]

#Make rowsums a column
CT["row_total"]<-rowSums(CT)

#Make allele column
CT["Allele"]<-row.names(CT)
#Write out contingency table
write.csv(CT,paste(kmer,"_contingency_table.csv",sep=""),row.names=F,quote=F)

#Make colsums  a row
cols<-ncol(CT)
colsum_row<-c(colSums(CT[,1:(cols-1)]))
#CT_tmp<-rbind(CT,colsum_row)

#Make heatmap
CT.m<-melt(CT,id.vars = "Allele")
colnames(CT.m)<-c("X1","X2","value")
xint = length(unique(CT.m$X2)) - 0.5
yint = 0.9 
p2<-ggplot(CT.m, aes(x=X2,y=X1))+
geom_tile(aes(fill=value))+
scale_fill_gradient(low = "white",high = "steelblue",name="Isolate Count",trans = "log") +
#geom_text(aes(label=value), data=cbind(aggregate(value~X1, CT.m, sum), X2="total"),size=6) +
scale_x_discrete(position = "top") +
scale_y_discrete(breaks=c(seq(1:nrow(CT))),limits=c(CT$X1),expand = c(0,0)) +
geom_text(aes(X2, label= value),size=6) +
ggtitle("MLST counts per allele") +
theme_bw()  +
theme(axis.text.x=element_text(angle=45,hjust=-.25,size=16,face="bold")) +
theme(axis.text.y=element_text(size=16,face="bold")) +
theme(panel.border = element_blank(), panel.grid.major = element_blank(),panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"),axis.title=element_text(size=18,face="bold"),plot.title = element_text(size=20,face="bold")) +
geom_vline(xintercept = xint) +
geom_hline(yintercept = yint) + 
labs(y="Allele",x="MLST") +
annotate("text", x = 1:length(colsum_row), y =.7, label = colsum_row,size = 6) 
ggsave(filename=paste(kmer,"_cont_tab_logscale.pdf",sep=""), width = 34, height = 28, plot=p2)
}



