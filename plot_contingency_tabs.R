library(grid)
library(stringr)
library(gridExtra)
library(ggdendro)
library(ggplot2)
library(reshape2)
library(ape)
library(ggplot2)
library(magrittr)
library(gplots)
library(RColorBrewer)
library(plyr)
require(gridBase)
require(gtable)
args <- commandArgs(TRUE)

#Read in parameters from the command line
allele_cutoff<-args[1]
percentage_cutoff <- args[2]
kmers_of_interest<-args[3]

#allele_cutoff<-20
#percentage_cutoff<-0.9


mers_of_interest<-read.table(kmers_of_interest,header=F)
table4plotting<-args[4]
tab<-read.csv(table4plotting,header=T)

#Remove first column (artifact from writing out or python pandas tab)
tab<-tab[,-1]

#Remove duplicates
tab<-tab[!duplicated(tab), ]

#Subset based on number of alleles
tab<-tab[tab$alleles>allele_cutoff,]

for (kmer in mers_of_interest$V1){
#kmer<-"gatggagtagcaggtccaac_caacaggtgctagtgcaata"

print(kmer)
#Subset df for a kmer to generate individual plots
tab_sub<-tab[tab$k1_k2 == kmer,]


if (nrow(tab_sub) >10){ 

print("more than 10 rows")

#Get counts of each MLST type for each kmer
CT<-dcast(tab_sub,k1_k2*insertion_sequence~MLST,fill=0,length)

#Make allele column and name it by sorting by length
CT["length"] <- nchar(as.character(CT$insertion_sequence)) + nchar(as.character(CT$k1_k2))-1

#Order by length
CT<-CT[order(CT$length),]

#Re-do rownames
rownames(CT)<-1:nrow(CT) 

CT["Allele"]<-paste((nchar(as.character(CT$insertion_sequence)) + nchar(as.character(CT$k1_k2))-1),as.character(row.names(CT)),sep="_")
CT$Allele<-as.factor(CT$Allele)

#Remove length column
CT<-CT[ , -which(names(CT) %in% c("length"))] 

#Write out contingency table
write.csv(CT,paste(kmer,"_contingency_table.csv",sep=""),row.names=F,quote=F)


tmp<-merge(tab,CT,by.x="insertion_sequence",by.y="insertion_sequence")
tmp2<-tmp[,c("isolate","length_allele_mapping","Allele","allele","k1_k2.x")]
colnames(tmp2)<-c("isolate","original_allele_mapping","updated_allele_mapping","allele","kmer")
write.csv(tmp2,paste(kmer,"_ins_seq_mapping.csv",sep=""),row.names=F,quote=F)
rm(tmp)
rm(tmp2)

#Remove insertion sequence column
CT<-CT[,-which(names(CT) %in% c("insertion_sequence","k1_k2"))]


#Make heatmap
CT.m<-melt(CT,id.vars = "Allele")

print("successful melt")

#Create not in function
'%!in%' <- function(x,y)!('%in%'(x,y))


order<-read.table("order_of_mlst.txt",header=F)
x_order<-c(order$V2)
missing_mlsts<-x_order[x_order %!in% CT.m$variable]

if (length(missing_mlsts)!=0) {
	n_alleles<-length(unique(as.factor(CT.m$Allele)))
	n_rows2add<-n_alleles*length(missing_mlsts)
	unique_alleles<-as.factor(as.character(unique(as.factor(CT.m$Allele))))
	print("length missing_MLSTS")
	print(length(missing_mlsts))


	#Format data frame so that mlsts not represented by the kmers are present in the contingency tables
	datalist = list()
	for (i in (1:length(missing_mlsts))){
		tmp<-data.frame(unique_alleles,rep(missing_mlsts[i],n_rows2add),rep(0,n_rows2add))
		datalist[[i]] <- tmp
	}

	mlsts2add = do.call(rbind, datalist)
	colnames(mlsts2add)<-c("Allele","variable","value")
	mlsts2add$variable<-as.factor(mlsts2add$variable)

	CT.m<-rbind(CT.m,mlsts2add)

}
else if (length(missing_mlsts)==0){
	n_alleles<-length(unique(as.factor(CT.m$Allele)))
}

CT.m<-CT.m[order(match(CT.m$variable,x_order)),]
CT.m$variable<-factor(CT.m$variable,levels=CT.m$variable)



#get isolates that are not in the table
mlst_file<-read.table("mlst_file.txt",header=T)
dfs<-read.table("data_freeze_samples.txt")
missing<-dfs[dfs$V1 %!in% tab_sub$isolate,]
missing_samples<-tab[tab$isolate %in% missing,]
missing_samples<-missing_samples[,c("isolate","MLST")]
missing_samples<-missing_samples[!duplicated(missing_samples), ]

#Create vector from missing samples
missing_samples_tab<-table(missing_samples)
sums_missing_samples<-data.frame(colSums(missing_samples_tab))
sums_missing_samples["MLST"]<-row.names(sums_missing_samples)
colnames(sums_missing_samples)<-c("count","MLST")
x_order<-data.frame(x_order)
missing_samples_merged<-merge(x_order, sums_missing_samples, by.y="MLST",by.x="x_order",all.x=TRUE)
missing_samples_merged[is.na(missing_samples_merged)] <- 0
missing_samples_merged<-missing_samples_merged[order(match(missing_samples_merged$x_order,x_order$x_order)),]




#Get row and column sums for contingency matrix
colsums<-aggregate(as.numeric(CT.m$value), by=list(Category=CT.m$Allele), FUN=sum)
rowsums<-aggregate(as.numeric(CT.m$value), by=list(Category=CT.m$variable), FUN=sum)
total=sum(colsums$x)
print(total)

#Merge in eburst data
eburst<-read.table("eburst_groups.txt",header=T)
tmp<-merge(CT.m,eburst,by.x="variable",by.y="ST",all.x=T)
CT.m<-tmp
CT.m[is.na(CT.m)] <- 0

#Make labels to annotate table with eburst groups
eburst_labs<-CT.m[,c("variable","eburst_group")]
eburst_labs<-eburst_labs[!duplicated(eburst_labs), ]

#Make the contingencty tables
n<-as.numeric(percentage_cutoff)*(length(unique(tab$isolate)))
print(n)

if (total > n){
	xlabels<-c(as.character(colsums$x),sum(colsums$x))
	n_mlst<-length(unique(CT.m$variable))
	p<-ggplot(CT.m, aes(x=Allele,y=variable))+
	geom_tile(aes(fill=as.numeric(value)))+
	scale_fill_gradient(low = "white",high = "steelblue",name="Isolate Count",trans = "log") +
	scale_x_discrete(position = "top",expand=c(0,0)) +
	coord_cartesian(xlim = c(-2,length(unique(CT.m$Allele))+4),ylim=c(-1,n_mlst+0.5)) +
	geom_text(aes(Allele, label= value),size=8) +
	ggtitle(paste("MLST counts per allele\n",total," total isolates\n",kmer,sep=" ")) +
	theme_bw()  +
	theme(axis.text.x=element_text(angle=45,hjust=-.25,size=18)) +
	theme(axis.text.y=element_text(size=22)) +
	theme(panel.border = element_blank(), panel.grid.major = element_blank(),panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"),axis.title=element_text(size=18,face="bold"),plot.title = element_text(size=20,face="bold")) +
	labs(x="Allele",y=" ") +
	annotate("rect", xmin=n_alleles+1.2, xmax=n_alleles+2, ymin=-1, ymax=length(unique(CT.m$Allele))+1,fill="white") + 
	annotate("text",x=n_alleles+1.2,y=c(seq(1:(n_mlst))),label=as.character(rowsums$x),size=8) +
	annotate("rect", xmin=0, xmax=nrow(colsums)+1, ymin=-1.5, ymax=0.5,fill="white") +
	annotate("text",x=c(seq(1:nrow(colsums))),y=-0.5,label=as.character(colsums$x),size=8) +
        annotate("rect", xmin=-1, xmax=0, ymin=0, ymax=length(eburst_labs$eburst_group),fill="white") +
	annotate("text",x=-0.48,y=c(seq(1:length(eburst_labs$eburst_group))),label=eburst_labs$eburst_group,size=8) +
	annotate("rect", xmin=n_alleles+2, xmax=n_alleles+2.8, ymin=-1, ymax=length(unique(CT.m$Allele))+2,fill="white") +
        annotate("text",x=n_alleles+2,y=c(seq(1:(n_mlst))),label=as.character(missing_samples_merged$count),size=8)
      	ggsave(filename=paste(kmer,"_cont_tab_logscale.pdf",sep=""), width = 38, height = 28, plot=p)
      	ggsave(filename=paste(kmer,"_cont_tab_logscale.svg",sep=""), width = 38, height = 28, plot=p)
	
}
}
}
