DATA_DIR=$1
WD=$2
FLANKING=$3

while read name;
do
echo "
#!/bin/bash
#BSUB -J ${name}
#BSUB -P acc_PacbioGenomes
#BSUB -q private
#BSUB -n 1
#BSUB -W 6:00
#BSUB -m bashia02
#BSUB -o ${WD}/scratch/get_ins_seq_${name}.stdout
#BSUB -eo ${WD}/scratch/get_ins_seq_${name}.stderr
#BSUB -L /bin/bash

cd ${WD}
module load python py_packages

python ${WD}/get_ins_seq_per_k1k2.py ${DATA_DIR}/${name}.fa ${WD}/results/kmer_results/all_k1k2_pairs.txt ${FLANKING} > ${WD}/scratch/${name}_other_alt_kmer_counts.txt " > ${WD}/scripts/${name}.get_ins_seqs_other_alt_kmers.sh


bsub < ${WD}/scripts/${name}.get_ins_seqs_other_alt_kmers.sh

done <  samples.txt


