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
#BSUB -o /sc/orga/scratch/webste01/final_TRF_pipeline/scratch/get_ins_seq_${name}.stdout
#BSUB -eo /sc/orga/scratch/webste01/final_TRF_pipeline/scratch/get_ins_seq_${name}.stderr
#BSUB -L /bin/bash

cd /sc/orga/scratch/webste01/final_TRF_pipeline/scripts
module load python py_packages

python get_ins_seq_per_K1K2.py /sc/orga/scratch/webste01/final_TRF_pipeline/data/${name}.fa /sc/orga/scratch/webste01/final_TRF_pipeline/results/kmer_results/all_k1k2_pairs.txt" > batch_scripts/${name}.get_ins_seqs.sh

bsub < batch_scripts/${name}.get_ins_seqs.sh

done <  data_freeze_samples.txt


