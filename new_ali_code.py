import operator
import os
import sys
import bisect
from Bio import SeqIO
from Bio.Seq import Seq

in_fa          = sys.argv[1]
in_trf         = sys.argv[2]
in_jf          = sys.argv[3]
flanking_bp    = int(20000)
k = int(15)

#in_fa          = "/sc/orga/scratch/webste01/TRF_analysis/data/u00012crpx_c_025914.fa"
#in_trf         = "/sc/orga/scratch/webste01/TRF_analysis/results/trf_results/u00012crpx_c_025914.fa_out"
#in_jf          = "/sc/orga/scratch/webste01/TRF_analysis/results/jellyfish_results/high_freq_kmers.txt"



#Get the isolate name and set the name of the outfile
tmp= os.path.basename(os.path.normpath(in_fa))
isolate = str(tmp).split('.')[0]
out_name = str(isolate) + "paired_K1K2.csv"
out_name_2 = str(isolate) + "unpaired_K1K2.csv"

#Initialize a set of high frequency kmers
kmers = set()
with open(in_jf) as jellyfish:
	for l in jellyfish:
		l = l.strip().split()
		kmers.add(l[0].lower())
	jellyfish.close()

def get_seq_btw_2coords(start_coord,end_coord):
#'''Given a fasta sequence, returns the sequence between indices start_coord,end_coord'''
	full_fasta_seq = SeqIO.parse(open(in_fa),'fasta')
	for fasta in full_fasta_seq:
		fa_string = str(fasta.seq)
		subset_fa =fa_string[int(start_coord):int(end_coord)]
		return subset_fa


#'''Given a fasta sequence, a kmer-size (k), and a set of kmers (NOT A LIST) returns a dict mapping the position to the kmer'''
def build_kmer_dict(in_fa, k, kmers):
	full_fasta_seq = SeqIO.parse(open(in_fa),'fasta')
	f_count = 0
	pos_kmer_dict = {}
	for fasta in full_fasta_seq:
		f_count += 1
		fa_string = str(fasta.seq).lower()
		for i in range(len(fa_string)-k):
			kmer = fa_string[i:i+k]
			if kmer in kmers:
				pos_kmer_dict[i] = kmer
		if f_count > 1:
			print "Too many fasta sequences!!"
		return pos_kmer_dict

#'''Returns a file with the format k1,k2,start_pos_k1,end_pos_k2,len_of_insertion_seq,insertion_seq,trf_repeat'''
def read_trf_and_find_kmer (in_trf, pos_kmer_dict, flank):
	kmer_pos_sorted = pos_kmer_dict.keys()
	kmer_pos_sorted.sort()
	with open(in_trf,'r') as trf:
		next(trf)
		for l in trf:
			if not l.lstrip().startswith('@'):
				l=l.strip().split()
				start = int(l[0])
				#print "start="
				#print start
				end = int(l[1])
				#print "end="
				#print end
				repeat = str(l[13]).lower()
				# find closest upstream kmer
				start_index = bisect.bisect_left(kmer_pos_sorted,start-1)
				k1_pos = kmer_pos_sorted[start_index-1]
				#print "k1 pos (closest upstream kmer)="
				#print k1_pos
				k1 = pos_kmer_dict[k1_pos]
				#print  k1
				# find closest downstream kmer
				end_index = bisect.bisect_right(kmer_pos_sorted,end+1)
				k2_pos = kmer_pos_sorted[end_index]
				k2 = pos_kmer_dict[k2_pos]
				#print "k2 pos (closest downstream kmer)="
                                #print k2_pos
				#print k2
				ins_length = abs(k2_pos - k1_pos)
				if ins_length > 0:
				# Check to see if fullfills flank criteria i.e. the repeat must fall within flank bp of start/end of repeats
					if abs(start-k1_pos) > flank:
						#print "larger region"
						#print abs(start-k1_pos)
						b_file.write("%s,%s,%d,%d,%d,%s,%s \n" % (k1,k2,k1_pos,k2_pos,ins_length,repeat,isolate))
					elif abs(end-k2_pos) > flank:
						#print "larger region"
						#print abs(end-k2_pos)
						b_file.write("%s,%s,%d,%d,%d,%s,%s \n" % (k1,k2,k1_pos,k2_pos,ins_length,repeat,isolate))
					else:
						ins_seq = get_seq_btw_2coords(k1_pos,k2_pos)
						a_file.write("%s,%s,%d,%d,%d,%s,%s,%s \n" % (k1,k2,k1_pos,k2_pos,ins_length,ins_seq,repeat,isolate))


pos_kmer_dict = build_kmer_dict(in_fa,k,kmers)
			
#output file:
with open(out_name, 'w') as a_file:
	with open(out_name_2, 'w') as b_file:
		a_file.write("K1,K2,start_pos_k1,end_pos_k2,len_of_insertion_seq,insertion_seq,trf_repeat")
	        b_file.write("K1,K2,start_pos_k1,end_pos_k2,len_of_insertion_seq,trf_repeat")
		a_file.write('\n')
		b_file.write('\n')
		read_trf_and_find_kmer(in_trf, pos_kmer_dict, flanking_bp)

