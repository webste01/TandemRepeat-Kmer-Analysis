#!/bin/bash
#Get location of fasta files, and desired kmer size, and working directory (i.e. directory holding scripts)
DATA_DIR=$1
KMER_SIZE=$2
WD=$3
TRF_dir=$4 #Directory with the TRF files
FLANKING=$5 #Size of flanking region to TR allowable
INS_SIZE=$6 #Maximum insertion size limit

echo ${DATA_DIR}

#Make directory to hold batch scripts
#mkdir scripts

#Make directory to hold scratch files
#mkdir scratch

#Make directory for fasta files
#mkdir fasta_dir

#Make results directory
#mkdir results
#mkdir results/kmer_results

#Copy all fasta files into separate fasta directory
#cp ${DATA_DIR}/u0*fa fasta_dir/

#CD into fasta directory, run script to identify kmers
#cd fasta_dir
#sh ${WD}/identify_kmers.sh ${KMER_SIZE} ${WD}


#Move the kmer frequency files into the results directory
#mv ${WD}/fasta_dir/*freq1 ${WD}/results/kmer_results/

#Change to kmer_results dir
#cd ${WD}/results/kmer_results

#Concatenate and find unique kmers
#cat *freq1 | awk '{print $1}' | sort | uniq -c  > total-freq.txt

#Get kmers that occur in > 95% of genomes *will have to change*
#awk '$1>236' total-freq.txt > singleton_kmers.txt

#Change directory back to working directory
#cd ${WD}

#Get a "fofn" of fasta files
#cd ${DATA_DIR}
#ls -1 *fa | sed -e 's/\.fa$//' > ${WD}/samples.txt
#cd ${WD}

#Run script to get k1_k2
#sh get_k1k2_pairs.sh ${DATA_DIR} ${TRF_dir} ${WD} ${KMER_SIZE} ${FLANKING} 

#Move the files to results directory
#mv ${WD}/*csv ${WD}/results/kmer_results/

#Move to that directory
#cd ${WD}/results/kmer_results/

#Get all k1_k2 pairs
#cat *csv | awk -F, '$1>0' | awk -F, '$2>0' | awk -F, '{print $1,$2}' | sort | uniq > all_k1k2_pairs.txt

#Replace space with underbar in all_k1k2_pairs.txt
#sed 's/ /_/g' all_k1k2_pairs.txt > tmp
#mv tmp all_k1k2_pairs.txt

#Get the insertion sequence between the k1_k2 pairs
#cd ${WD}
#sh get_insertion_sequence_btw_kmers.sh ${DATA_DIR} ${WD} ${INS_SIZE}

#Move insertion files to results directory, get table for insertion sequence across all isolates
mkdir ${WD}/results/insertion_sequences

mv *ins_seqs.csv ${WD}/results/insertion_sequences
cd ${WD}/results/insertion_sequences
cat *ins_seqs.csv | awk -F, '$5>1' > all_insertion_seqs.csv

#Remove duplicate lines (aka header line that remained from cat command)
awk '!x[$0]++' all_insertion_seqs.csv > tmp
mv tmp all_insertion_seqs.csv

#Format the insertion table for plotting
cd ${WD}
module purge
module load python py_packages
python format_table_for_plotting.py ${WD}/results/insertion_sequences/all_insertion_seqs.csv table4plotting.csv






