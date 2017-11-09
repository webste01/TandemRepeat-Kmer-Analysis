tab<-read.csv("all_insertion_seqs.csv",header=T)
test<- tab %>% group_by(k1_k2) %>% summarise(count = n_distinct(insertion_sequence_length))
#find alleles with more than 10 variants
test2<-test[test$count > 10,]
write.table(test2$k1_k2,"k1k2_of_interest.txt",quote=F,row.names=F,col.names=F)
