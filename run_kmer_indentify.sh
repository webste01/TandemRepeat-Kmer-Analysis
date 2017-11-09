cd /sc/orga/scratch/webste01/final_TRF_pipeline/data

for f in *fa;
do
echo "#!/bin/bash
#BSUB -J ${f}_identify_kmers
#BSUB -P acc_PacbioGenomes
#BSUB -q private
#BSUB -n 1
#BSUB -W 24:00
#BSUB -m bashia02
#BSUB -o /sc/orga/scratch/webste01/final_TRF_pipeline/scratch/${f}_kmer_identify.stdout
#BSUB -eo /sc/orga/scratch/webste01/final_TRF_pipeline/scratch/${f}_kmer_identify.stderr
#BSUB -L /bin/bash

module load python py_packages

cd /sc/orga/scratch/webste01/final_TRF_pipeline/data

python /sc/orga/scratch/webste01/final_TRF_pipeline/scripts/identify_kmers.py $f 15

#Get unique kmers aka with frequency 1
#Now we have the singleton kmers per genome
awk '{print \$1}' ${f}_kmers.txt | sort | uniq -u > tmp_${f}

#Remove lines that are not 15 characters long (artifact of alis code)
sed -r '/^.{,14}$/d' tmp_${f} > ${f}_kmers.txt_freq1


rm tmp_${f}" > /sc/orga/scratch/webste01/final_TRF_pipeline/scripts/batch_scripts/${f}_kmer_identify.sh

bsub < /sc/orga/scratch/webste01/final_TRF_pipeline/scripts/batch_scripts/${f}_kmer_identify.sh

done 

