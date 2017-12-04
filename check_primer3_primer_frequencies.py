import operator
import os
import sys
import csv
from Bio import SeqIO
from Bio.Seq import Seq

#Check is the primer3 predicted primers occur in the set of singleton kmers
all_kmers=sys.argv[1]
primer_outfile=sys.argv[2]

mers = []

with open(all_kmers,'r') as kmer_file:
	for l in kmer_file:
		l=l.strip().split()
		mers.append(str(l[1]).lower())

with open(primer_outfile) as p_of:
	for l in p_of:
		if "LEFT" in l:
			ll = l.strip().split("PRIMER")[1]	
			kmer=ll.strip().split(" ")[14]
			if kmer in mers:
				print l
		if "RIGHT" in l:
			ll = l.strip().split("PRIMER")[1]
                        kmer=ll.strip().split(" ")[14]
			revcomp_seq = str(Seq(str(kmer)).reverse_complement())
                        if revcomp_seq in mers:
                                print l

		
