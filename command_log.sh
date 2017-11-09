#Run code to get kmers
sh run_kmer_indentify.sh

#Move files to results directory
mv /sc/orga/scratch/webste01/final_TRF_pipeline/data/*txt* /sc/orga/scratch/webste01/final_TRF_pipeline/results/kmer_results/

#Change directories
cd /sc/orga/scratch/webste01/final_TRF_pipeline/results/kmer_results/

#Concatenate and find unique  (run as job on cluster:get_kmer_frequencies.sh)
cat u*freq1 | awk '{print $1}' | sort | uniq -c  > total-freq.txt

#Get kmers that occur in >95% number of genomes: singleton kmers in X genomes (rename outfile)
awk '$1>236' total-freq.txt > singleton_kmers.txt

#Change directory back to scripts 
cd /sc/orga/scratch/webste01/final_TRF_pipeline/scripts

#Extract k1k2 pairs
sh get_k1k2_pairs.sh

#Move files to results directory
mv *csv /sc/orga/scratch/webste01/final_TRF_pipeline/results/kmer_results/
cd /sc/orga/scratch/webste01/final_TRF_pipeline/results/kmer_results/

#Get all k2k1 pairs
cat *csv | awk -F, '$1>0' | awk -F, '$2>0' | awk -F, '{print $1,$2}' | sort | uniq > all_k1k2_pairs.txt

#Get the insertion sequences between these pairs
cd /sc/orga/scratch/webste01/final_TRF_pipeline/scripts
sh get_insertion_seqs.sh

#Cat files and only select those rows which have insertion sequences
mv *ins_seqs.csv /sc/orga/scratch/webste01/final_TRF_pipeline/results/insertion_sequence_results/
cd  /sc/orga/scratch/webste01/final_TRF_pipeline/results/insertion_sequence_results/
cat *ins_seqs.csv | awk -F, '$5>1' > all_insertion_seqs.csv

#Remove duplicate lines (aka header line that remained from cat command)
awk '!x[$0]++' all_insertion_seqs.csv > tmp
mv tmp all_insertion_seqs.csv

#perform analysis to get the kmers of interest
cat  all_insertion_seqs.csv | awk -F, '{out[$1 "\t" $4]++} END{for(var in out){print var, out[var]}}' > pairs-sequence-counts.txt

#Find insertion sequences greater than X
awk -F, '$5>300' all_insertion_seqs.csv  > ins_longer_300.csv

#Make violin plots stratified by mlst
cd ../plots
Rscript ../scripts/InsertionSize_violinPlot.R results/insertion_sequence_results/ins_longer_300.csv ../mlst_file.txt 
/sc/orga/scratch/webste01/final_TRF_pipeline/results/insertion_sequence_results/kmers_of_interest.txt

