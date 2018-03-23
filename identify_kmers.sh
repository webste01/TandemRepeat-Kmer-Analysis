KMER_SIZE=$1
WD=$2
for f in *fa*;
do
echo "#!/bin/bash
#BSUB -J ${f}_identify_kmers
#BSUB -P acc_PacbioGenomes
#BSUB -q private
#BSUB -n 1
#BSUB -W 24:00
#BSUB -m bashia02
#BSUB -o ${WD}/scratch/${f}_kmer_identify.stdout
#BSUB -eo ${WD}/scratch/${f}_kmer_identify.stderr
#BSUB -L /bin/bash

module load python py_packages

cd ${WD}/fasta_dir

python ${WD}/identify_kmers.py $f ${KMER_SIZE} " > ${WD}/scripts/${f}_kmer_identify_${KMER_SIZE}.sh

#sh ${WD}/scripts/${f}_kmer_identify_${KMER_SIZE}.sh
bsub <  ${WD}/scripts/${f}_kmer_identify_${KMER_SIZE}.sh

done 

