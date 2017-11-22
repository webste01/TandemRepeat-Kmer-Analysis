library(grid)
library(ggplot2)
library(reshape2)
library(ggplot2)
library(magrittr)
library(gplots)
library(RColorBrewer)
library(plyr)

tab<-read.csv("table4plotting.csv",header=T)

#Remove first column (artifact from writing out or python pandas tab)
tab<-tab[,-1]

for (kmer in unique(tab$k1_k2)){

#Subset df for a kmer to generate individual plots
tab_sub<-tab[tab$k1_k2 == kmer,]

#Get counts of each MLST type for each kmer
CT<-dcast(tab_sub,k1_k2*insertion_sequence~MLST,fill=0,length)

#Make allele column
CT["Allele"]<-paste(nchar(as.character(CT$insertion_sequence)),as.character(row.names(CT)),sep="_")

#Write out contingency table
write.csv(CT,paste(kmer,"_contingency_table.csv",sep=""),row.names=F,quote=F)

#Make heatmap
CT.m<-melt(CT,id.vars = "Allele")
nskip<- (length(unique(CT.m$Allele))*2)+1

p<-ggplot(CT.m[nskip:nrow(CT.m),], aes(x=variable,y=Allele))+
geom_tile(aes(fill=as.numeric(value)))+
scale_fill_gradient(low = "white",high = "steelblue",name="Isolate Count",trans = "log") +
scale_x_discrete(position = "top") +
geom_text(aes(variable, label= value),size=6) +
ggtitle("MLST counts per allele") +
theme_bw()  +
theme(axis.text.x=element_text(angle=45,hjust=-.25,size=16,face="bold")) +
theme(axis.text.y=element_text(size=16,face="bold")) +
theme(panel.border = element_blank(), panel.grid.major = element_blank(),panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"),axis.title=element_text(size=18,face="bold"),plot.title = element_text(size=20,face="bold")) +
labs(y="Allele",x="MLST") 

tmp<-CT.m[nskip:nrow(CT.m),]
rowsums<-aggregate(as.numeric(tmp$value), by=list(Category=tmp$Allele), FUN=sum)
colsums<-aggregate(as.numeric(tmp$value), by=list(Category=tmp$variable), FUN=sum)
total=sum(colsums$x)
xlabels<-c(as.character(colsums$x),sum(colsums$x))

n_mlst<-length(CT[,3:ncol(CT)])-1
p+annotate("rect", xmin=n_mlst+1, xmax=n_mlst+2, ymin=0, ymax=length(unique(CT.m$Allele)),fill="white") + 
annotate("text",x=n_mlst+1,y=c(seq(1:length(rowsums))),label=as.character(rowsums$x)) +
annotate("rect", xmin=0, xmax=n_mlst+1, ymin=0, ymax=0.5,fill="white") +
annotate("text",x=c(seq(1:(n_mlst+1))),y=0.25,label=as.character(xlabels))

ggsave(filename=paste(kmer,"_cont_tab_logscale.pdf",sep=""), width = 34, height = 28, plot=p2)
}

