import operator
import os
import sys
import csv
from Bio import SeqIO
from Bio.Seq import Seq
#Script to get 200bp flanking region of each tandem repeat, and step through each repeat for the closest kmer
in_fa          = sys.argv[1]
in_trf         = sys.argv[2]
in_kmers          = sys.argv[3]
flanking_bp = int(sys.argv[4])
k=int(sys.argv[5])

tmp= os.path.basename(os.path.normpath(in_fa))
isolate = str(tmp).split('.')[0]
out_name = str(isolate) + ".csv"

full_fasta_seq = SeqIO.parse(open(in_fa),'fasta')
for fasta in full_fasta_seq:
                fa_string = str(fasta.seq)

#List of high frequency kmers
kmers=[]

#idictionary of Tuple that is keyed upon [leftmer,rightmer,tandem repeat] with the corresponding isolate name as the value
d={}

with open(in_kmers) as mers:
        for l in mers:
                l = l.strip().split()
                kmers.append(l[1].lower())

kmers_set = set(kmers)

def find_kmer_in_seq(repeat,start_coord,end_coord,isolate):
#'''given two sequences that flank the region with a TR and the kmers from the list of all high frequency kmers, get the kmer closest to the end of the upstream flanking sequence and the kmer closest to the start of the downstream flanking sequence'''
        downstream_kmers = {}
        upstream_kmers ={}
	upstream_start = int(end_coord)
	upstream_end = int(end_coord + flanking_bp)
	downstream_end = int(start_coord)
	downstream_start = int(start_coord - flanking_bp)
	if downstream_start < 1:
		downstream_end = 1
        for i in range (upstream_start,upstream_end):
                kmer = fa_string[i:i+k]
                if (kmer in kmers_set):
                        max_upmer = kmer
			max_upmer_end = i + k
                        finish="TRUE"
                        break
                else:
                        max_upmer = "0"
			max_upmer_end = "0"
                        finish = "FALSE"
                if finish == "TRUE":
                        break
        for i in range (1,flanking_bp):
		kmer = fa_string[downstream_end - (k+i):downstream_end - i]
                if (kmer in kmers_set):
                        min_downmer = kmer
			min_downmer_start = downstream_end - (k+i)
                        finish="TRUE"
                        break
                else:
                        min_downmer = "0"
			min_downmer_start = "0"
                        finish = "FALSE"
                if finish == "TRUE":
                        break
        d[min_downmer,max_upmer,min_downmer_start,max_upmer_end,repeat]=isolate


with open(in_trf,'r') as trf:
        next(trf)
        for l in trf:
                if not l.lstrip().startswith('@'):
                        l=l.strip().split()
                        start = int(l[0])
                        end = int(l[1])
                        repeat = str(l[13])
                        find_kmer_in_seq(repeat,start,end,isolate)
trf.close()
t=d.items()

with open(out_name, 'w') as a_file:
    for result in t:
        dtr = str(result[0][0]) #Kmer closest upstream to TR
        utr = str(result[0][1]) #Kmer closest downstream to TR
	min_downmer_start = str(result[0][2])
	max_upmer_end = str(result[0][3])
        repeat = str(result[0][4]) #TR
        iso = result[1] #Isolate
        final = ','.join([dtr,utr,min_downmer_start,max_upmer_end,repeat,iso])
        a_file.write(final + '\n')
a_file.close()

