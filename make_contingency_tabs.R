library(reshape2)
library(ggplot2)
library(magrittr)
library(gplots)
library(RColorBrewer)
args <- commandArgs(TRUE)
fn <- args[1]

STs<-read.table("mlst_file.txt",header=F)
colnames(STs)<-c("Isolate","ST")
ST_freq<-data.frame(table(STs$ST))
colnames(ST_freq)<-c("ST","Freq")

fn1<-read.table(fn,header=T)
colnames(fn1)<-c("K1_K2","Insertion_Seq")
name = fn1[1,1]


K1K2_ins_ST<-read.table("K1K2_ins_ST.txt",header=T)

colnames(K1K2_ins_ST) <- c("Isolate","K1_K2","Insertion_Seq","ST")

#Merge again with ST information to make sure its consistent
tmp<-merge(K1K2_ins_ST,STs,by="Isolate",all.x=TRUE)
tmp<-as.data.frame(tmp)
tmp<-tmp[,-c(4)] 
tmp["ST"]<-as.factor(tmp$ST.y)
tmp<-tmp[,-c(4)] 
#Remove NA
tmp[is.na(tmp)] <- "none"

K1K2_ins_ST<-tmp

rm(tmp)

test<-merge(fn1,K1K2_ins_ST,by=c("K1_K2","Insertion_Seq"),all.x=TRUE)
deduped.test <- unique( test[ , 1:ncol(test) ] )
CT<-dcast(deduped.test,K1_K2*Insertion_Seq~ST,fill=0,length)

write.table(CT[,-1],paste(name,"_contingency.txt",sep=""),quote=FALSE,row.names=FALSE)

write.table(deduped.test[,-1],paste(name,"_ST_count_iso.txt",sep=""),quote=FALSE,row.names=FALSE)

CT["Allele"]<-unclass(CT$Insertion_Seq) %>% as.numeric

#Write out Allele and Number mapping
write.table(CT[,c("Allele","Insertion_Seq")],paste(name,"_allele_map.txt",sep=""),quote=FALSE,row.names=FALSE)

#Remove insertion sequence and K1K2 (for simplicity)
CT<-CT[,-c(1,2)]

#Move Allele Collumn to the front
row.names(CT)<-CT$Allele
lastcol<-ncol(CT)
CT<-CT[,-c(lastcol)]
CT<-as.matrix(CT)

#Make heatmap
CT.m<-melt(CT)
colnames(CT.m)<-c("X1","X2","value")
xint = length(unique(CT.m$X2)) + 0.5
p<-ggplot(CT.m, aes(x=X2,y=X1))+
geom_tile(aes(fill=value))+
scale_fill_gradient(low = "white",high = "steelblue",name="Isolate Count") +  
geom_text(aes(label=value), data=cbind(aggregate(value~X1, CT.m, sum), X2="total"),size=6) + 
scale_x_discrete(position = "top") +
scale_y_discrete(breaks=c(seq(1:nrow(CT.m))),limits=c(CT.m$X1),expand = c(0,0)) +
geom_text(aes(X2, label= value),size=6) +
ggtitle("MLST counts per allele") +
theme_bw()  +
theme(axis.text.x=element_text(angle=45,hjust=-.25,size=18,face="bold")) +
theme(axis.text.y=element_text(size=18,face="bold")) +
theme(panel.border = element_blank(), panel.grid.major = element_blank(),panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"),axis.title=element_text(size=18,face="bold"),plot.title = element_text(size=20,face="bold")) + 
geom_vline(xintercept = xint) +
labs(y="Allele",x="MLST") 


#ggsave(filename=paste(name,"_cont_tab.pdf",sep=""), width = 24, height = 20, plot=p) 
ggsave(filename=paste(name,"_cont_tab.pdf",sep=""), width = 34, height = 28, plot=p) 




