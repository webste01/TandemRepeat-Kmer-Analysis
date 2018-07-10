#Get the insertion sequence from each k1k2
import operator
import os
import sys
import bisect
from Bio import SeqIO
from Bio.Seq import Seq
import re

in_fa = sys.argv[1] #input fasta sequence
k1k2 = sys.argv[2]


#Limit to the insertion size
x=int(sys.argv[3])

#Extract isolate ID from fastafile name
tmp= os.path.basename(os.path.normpath(in_fa))
isolate = str(tmp).split('.')[0]

def fasta_reader(in_fa):
        full_fasta_seq = SeqIO.parse(open(in_fa),'fasta')
        for fasta in full_fasta_seq:
                fa_string = str(fasta.seq)
        return fa_string

#def kmer_check(seq,mer):
#        index = 0
#	print mer
#        positions = []
#        while index < len(seq):
#                index = seq.find(mer, index)
#                if index == -1:
#                        break
#                index += len(mer)
#                positions.append(index)
#        if len(positions) == 1:
#                return positions[0]
#        elif len(positions) > 1:
#                return None
#        else:
#                return None

def kmer_check(seq,mer):
	regex = re.compile(mer)
	match = regex.search(seq)
	if match:
		return int(match.start())
	else:
		return None


def revcomp(kmer):
	'''get reverse complement sequence of sequence'''
	revcomp_seq = str(Seq(kmer).reverse_complement())
	return revcomp_seq


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
    print k1
    k2= l[1].lower()
    print k2
    k1_start = None
    k2_end = None
    k1_start_inclusive = None
    fasta    = fasta_reader(in_fa)
    orientation = 'None'
    
    # check for the forward or reverse k 
    forward_k1_check = kmer_check(fasta,k1)
    reverse_k1_check = kmer_check(fasta,revcomp(k1))
    forward_k2_check = kmer_check(fasta,k2)
    reverse_k2_check = kmer_check(fasta,revcomp(k2))
    
    # update output line for sanity checking (analagous to harms)
    output_line = [None, None, None, None]
    if forward_k1_check: output_line[0] = 1
    if reverse_k1_check: output_line[1] = 1
    if forward_k2_check: output_line[2] = 1
    if reverse_k2_check: output_line[3] = 1
      
    found = None
    # identify forward pair
    if forward_k1_check and forward_k2_check:
      k1_start = forward_k1_check
      k2_end = forward_k2_check
      if k1_start < k2_end:
        orientation = 'forward'
        k1_start_inclusive = int(k1_start - len(str(k1)))
        ins_seq = get_seq_btw_2coords(fasta,k1_start_inclusive,k2_end)
      else:
        orientation = 'forward_k2k1'
	print >>sys.stderr, 'forward_k2k1'
	print >>sys.stderr, "k1= " + str(k1) + " start= " + str(k1_start)
        print >>sys.stderr, "k2= " + str(k2) + " end= " + str(k2_end)

    # identify reverse pair
    if reverse_k1_check and reverse_k2_check:
      k1_start = reverse_k1_check
      k2_end = reverse_k2_check
      if k2_end < k1_start:
        orientation = 'reverse'
        k1_start_inclusive = int(k2_end - len(str(k2)))
	k2_end = k1_start
        ins_seq = revcomp(get_seq_btw_2coords(fasta,k1_start_inclusive,k2_end))
      else:
        orientation = 'reverse_k1k2'
	print >>sys.stderr, "reverse_k1k2"
        print >>sys.stderr, "k1= " + str(k1) + " start= " + str(k1_start)
	print >>sys.stderr, "k2= " + str(k2) + " end= " + str(k2_end)

        
    # if present on the same orientation, print out the k1, k2 pair 
    if orientation in ('forward','reverse'):
      ins_seq_length = len(ins_seq)
      if ins_seq_length < x:
        out_file.write("%s,%s,%s,%s,%s,%s,%s \n" % (str(k1) + "_" + str(k2), k1_start_inclusive, k2_end, ins_seq, ins_seq_length, isolate, orientation))
      else:
        print >>sys.stderr, "too long!"
        print >>sys.stderr, ins_seq_length
        
    print "%s_%s\t%s" %(k1, k2, "\t".join(map(str, output_line)))

