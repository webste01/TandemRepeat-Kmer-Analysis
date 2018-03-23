DATA_DIR=$1
TRF_DIR=$2
WD=$3
KMER_size=$4
FLANKING=$5

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
#BSUB -o ${WD}/scratch/findk1k2_${name}.stdout
#BSUB -eo ${WD}/scratch/findk1k2_${name}.stderr
#BSUB -L /bin/bash

cd ${WD}
module load python py_packages

python get_k1k2_pairs.py ${DATA_DIR}/${name}.fa ${TRF_DIR}/${name}.fa_out ${WD}/results/kmer_results/singleton_kmers.txt ${FLANKING} ${KMER_size}" > ${WD}/scripts/${name}.getk1k2.sh

sh ${WD}/scripts/${name}.getk1k2.sh

done <  samples.txt


