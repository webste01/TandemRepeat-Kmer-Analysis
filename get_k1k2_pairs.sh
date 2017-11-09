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
#BSUB -o /sc/orga/scratch/webste01/final_TRF_pipeline/scratch/findk1k2_${name}.stdout
#BSUB -eo /sc/orga/scratch/webste01/final_TRF_pipeline/scratch/findk1k2_${name}.stderr
#BSUB -L /bin/bash

cd /sc/orga/scratch/webste01/final_TRF_pipeline/scripts
module load python py_packages


python get_k1k2_pairs.py /sc/orga/scratch/webste01/final_TRF_pipeline/data/${name}.fa /sc/orga/scratch/webste01/final_TRF_pipeline/results/trf_results/${name}.fa_out /sc/orga/scratch/webste01/final_TRF_pipeline/results/kmer_results/singleton_kmers.txt 20000" > batch_scripts/${name}.getk1k2.sh

bsub < batch_scripts/${name}.getk1k2.sh

done <  data_freeze_samples.txt


