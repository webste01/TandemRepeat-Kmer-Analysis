library(gplots) 
library(ggplot2)
library(magrittr)
for_vplot<-read.table("all_pairs_vplot.txt",header=T)
vp<-ggplot(for_vplot,aes(x=K1_K2,y=Length)) + geom_violin() + facet_wrap(K1_K2 ~ K1_K2,nrow=1) 

for_vplot["Pair"]<-unclass(for_vplot$K1_K2) %>% as.numeric 


pdf("ViolinPlot_Length2.pdf",height=8, width=24)
ggplot(for_vplot,aes(factor(K1_K2),Length ))+
geom_violin(aes(fill=factor(K1_K2)))+
geom_boxplot(alpha=0.3, color="black", width=.1)+
labs(x = "", y = "")+
theme_bw()+
theme(legend.title = element_blank())+
facet_wrap(~K1_K2, scales="free")
dev.off()


