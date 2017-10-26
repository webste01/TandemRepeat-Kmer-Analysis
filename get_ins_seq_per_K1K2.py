#Get the insertion sequence from each k1k2
import operator
import os
import sys
import bisect
from Bio import SeqIO
from Bio.Seq import Seq

in_fa = sys.argv[1] #input fasta sequence
#k1k2   = "/sc/orga/scratch/webste01/TRF_analysis/results/insertion_sequences/all_k1k2_pairs.txt"
k1k2 = sys.argv[2]
k=15

#Limit to the insertion size
x=2000

#Extract isolate ID from fastafile name
tmp= os.path.basename(os.path.normpath(in_fa))
isolate = str(tmp).split('.')[0]

#Fasta reader, read in fasta and then get the reverse complement
def fasta_reader(in_fa,rc):
	full_fasta_seq = SeqIO.parse(open(in_fa),'fasta')
	for fasta in full_fasta_seq:
		fa_string = str(fasta.seq)
	if str(rc) == "yes":
		rc_fa_string = str(fasta.seq.reverse_complement())
		return rc_fa_string
	else:	
		return fa_string

#Checks if a Kmer is in a string sequence, returns a list of start positions if the kmer is in the sequence, returns false is the kmer is not in the sequence
def kmer_check(seq,mer):
	index = 0
	positions = []
	while index < len(seq):
        	index = seq.find(mer, index)
        	if index == -1:
            		break
     		index += len(mer)
		positions.append(index)
	if len(positions) == 1:
		return positions[0]	
	elif len(positions) > 1:
		return None
	else:
		return None

#Get sequence of fasta (as string) between two coordinates
def get_seq_btw_2coords(fasta_seq,start_coord,end_coord):
	subset_fa =fasta_seq[int(start_coord):int(end_coord)]
        return subset_fa

out_file = open(str(isolate) + "ins_seqs.csv",'w')
out_file.write("k1_k2,k1_start,k2_end,insertion_sequence,insertion_sequence_length,isolate,orientation \n")

with open(k1k2) as of:
        for l in of:
                l = l.strip().split("_")
                k1= l[0].lower() 
                k2= l[1].lower()
		fasta    = fasta_reader(in_fa,"no")
                rc_fasta = fasta_reader(in_fa,"yes") #Get the reverse complement of the sequence
		if kmer_check(fasta,k1): #Check to see if K1 is in the fasta
			orientation = "forward"
			k1_start = kmer_check(fasta,k1) 
			k1_start_inclusive = int(k1_start - k)
			if kmer_check(fasta,k2): #Check to see if K2 is in the fasta
                        	k2_start = kmer_check(fasta,k2)
				k2_end = int(k2_start)
				ins_seq = get_seq_btw_2coords(fasta,k1_start_inclusive,k2_end)
				ins_seq_length = len(ins_seq)
				if ins_seq_length < x:
					out_file.write("%s,%s,%s,%s,%s,%s,%s \n" % (str(k1) + "_" + str(k2),k1_start_inclusive,k2_end,ins_seq,ins_seq_length,isolate,orientation))
				
			else:
				out_file.write("%s,%s,%s,%s,%s,%s,%s \n" % (str(k1) + "_" + str(k2),k1_start,"0","0","0",isolate,orientation))	
		elif kmer_check(rc_fasta,k1): #Check to see if K1 is in the reverse complement of the sequence if it is not in the forward orientation
			orientation = "reverse"
			k1_rc_start = kmer_check(rc_fasta,k1)
			k1_rc_start_inclusive = int(k1_rc_start - k)
			if kmer_check(rc_fasta,k2): #Check to see if K2 is in the reverse complement of the fasta
				k2_rc_start = kmer_check(rc_fasta,k2)
				k2_rc_end = k2_rc_start 
				ins_seq = get_seq_btw_2coords(rc_fasta,k1_rc_start_inclusive,k2_rc_end)
				ins_seq_length = len(ins_seq)
				if ins_seq_length < x:
					out_file.write("%s,%s,%s,%s,%s,%s,%s \n" % (str(k1) + "_" + str(k2),k1_rc_start,k2_rc_end,ins_seq,ins_seq_length,isolate,orientation))
			else:
                                out_file.write("%s,%s,%s,%s,%s,%s,%s \n" % (str(k1) + "_" + str(k2),k1_rc_start,"0","0","0",isolate,orientation))          


                        		

of.close()
