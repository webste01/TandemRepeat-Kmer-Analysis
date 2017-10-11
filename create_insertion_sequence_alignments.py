import operator
import os
import sys
import csv
from Bio import SeqIO
from Bio.Seq import Seq
#Script to parse the output of *csv files with kmer,TRF and insertion sequence data and create fasta files for alignment
#Fasta files for each unique pair of unique kmers  are produced
#Outputs multi entry fasta files for each k1k2 pair, with the insertion sequence between k1 and k2 (inclusive of k1k2) per isolate

#infile     = sys.argv[1] #contig_outputs.csv: upstream kmer, downsteam kmer, insertion sequence, length of insertion sequence, isolate name
infile = 'contig_outputs.csv'


with open(infile) as  in_file:
	for l in in_file:
		l=l.strip().split(',')
		if l[3]!= '0' and l[3]!='none':
			out_fn = str(l[0]) + "_" + str(l[1]) + ".fasta"	
			iso_name = str(l[5]).split("_")[2]
			header = ">" + str(iso_name) 
			with open(out_fn,'a') as of:
				of.write("\n")
				of.write(header)
				of.write("\n")
				of.write(str(l[2]))
			of.close()



