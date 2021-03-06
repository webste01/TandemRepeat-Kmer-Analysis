import operator
import os
import sys
import csv
from Bio import SeqIO
from Bio.Seq import Seq

#Get the insertion sequence plus ~500 bp flanking region to put into a file that can be run through primer3

#Read in the fasta file, insertion sequence file, desired primer length, and the set of candidate kmers to find alternative kmers for
in_fa = sys.argv[1]
ins_seq_file = sys.argv[2]
primer_length = sys.argv[3]
candidate_kmers = sys.argv[4]

#Parse fasta file
full_fasta_seq = SeqIO.parse(open(in_fa),'fasta')
for fasta in full_fasta_seq:
                fa_string = str(fasta.seq)

tmp= os.path.basename(os.path.normpath(in_fa))
isolate = str(tmp).split('.')[0]
i=0

#Read in candidate kmers
kmers = []
with open(candidate_kmers,'r') as kmer_file:
	for l in kmer_file:
		l = l.strip().split()
		kmers.append(l[0])
kmer_file.close()

#Only get insertion sequences flanked by the candidate kmers, and get coordinates and sequence for 500bp flanking regions
with open(ins_seq_file,'r') as ins_seq:
	next(ins_seq)
        for l in ins_seq:
		l=l.strip().split(",")
                start_ins_seq = int(l[1])
		start_primer = int(int(l[1]) - 500)
                end_ins_seq = int(l[2])
		end_primer = int(int(l[2]) + 500)
		kmer_pair = str(l[0])
		orientation=str(l[6])
		if kmer_pair in kmers:
			length_insertion_seq = end_ins_seq - start_ins_seq
			length_primer_seq = end_primer - start_primer
			if length_insertion_seq > 200:
				i=i+1
				out_name = str(isolate) +"_"+ str(i) + "_primer3_infile" + primer_length +".txt"
				with open(out_name, 'w') as out_file:
					out_file.write("PRIMER_TASK=generic")
		        	        out_file.write("\n")
                			out_file.write("PRIMER_OPT_SIZE="+primer_length)
                			out_file.write("\n")
              				out_file.write("PRIMER_MIN_SIZE="+primer_length)
            				out_file.write("\n")
          				out_file.write("PRIMER_MAX_SIZE="+primer_length)
        				out_file.write("\n")
					out_file.write("SEQUENCE_ID="+kmer_pair)
					out_file.write("\n")
					out_file.write("PRIMER_PRODUCT_SIZE_RANGE="+str(length_insertion_seq)+"-"+str(length_primer_seq))
                                	out_file.write("\n")
					end_window_start = int(length_primer_seq - 500)
                                	out_file.write("SEQUENCE_PRIMER_PAIR_OK_REGION_LIST=1,500,"+str(end_window_start)+",500;")
					out_file.write("\n")
					if orientation == "forward":	
	                                	seq = fa_string[start_primer:end_primer]
						out_file.write("SEQUENCE_TEMPLATE="+seq)
						out_file.write("\n")
						out_file.write("=")
						out_file.close()
					else:
						seq = fa_string[start_primer:end_primer]
						rc_seq = str(Seq(seq).reverse_complement())
                                	        out_file.write("SEQUENCE_TEMPLATE="+rc_seq)
                                	        out_file.write("\n")
                                	        out_file.write("=")
                                	        out_file.close()					
ins_seq.close()
