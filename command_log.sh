#Identify kmers and get insertion sequences
#Data directory with fasta files, kmer size, working directory, directory with tandem repeat files (per fasta), max flanking size to look for kmers (bp), max insertion size (bp)
sh command.sh /sc/orga/projects/InfectiousDisease/studies/CDiff_first_paper/final_TRF_pipeline/data 20 /sc/orga/projects/InfectiousDisease/studies/CDiff_first_paper/pipeline_test /sc/orga/projects/InfectiousDisease/studies/CDiff_first_paper/final_TRF_pipeline/results/trf_results 20000 2000

#Perform machine learning to get the most predictive kmer sets
#python get_top_scoring_kmers.py [table4plotting] [score_cutoff] 
python get_top_scoring_kmers.py table4plotting.csv 0.9 0.7 

#Make the contingency tables
#Rscript plot_contingency_tabs.R [allele_cutoff] [percentage_cutoff] [kmers_of_interest] [table4plotting]
sh make_plots.sh 20 0.9 kmers_of_interest.txt table4plotting.csv
