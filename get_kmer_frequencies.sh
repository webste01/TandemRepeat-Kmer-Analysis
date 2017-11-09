#!/bin/bash
#BSUB -J kmer_freqs
#BSUB -P acc_PacbioGenomes
#BSUB -q private
#BSUB -n 1
#BSUB -W 24:00
#BSUB -m bashia02
#BSUB -o /sc/orga/scratch/webste01/TRF_analysis_2/scratch/kmer_freqs.stdout
#BSUB -eo /sc/orga/scratch/webste01/TRF_analysis_2/scratch/kmer_freqs.stderr
#BSUB -L /bin/bash

cd /sc/orga/scratch/webste01/TRF_analysis_2/results/kmer_results
cat u*freq1 | awk '{print $1}' | sort | uniq -c  > total-freq.txt
