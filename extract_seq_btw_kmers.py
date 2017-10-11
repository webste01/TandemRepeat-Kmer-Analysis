import operator
import os
import sys
import csv
from Bio import SeqIO
from Bio.Seq import Seq
#Script to get 200bp flanking region of each tandem repeat, and step through each repeat for the closest kmer
#Also extract the sequence between the kmers (to later be compared for similarity) 
in_fa          = sys.argv[1]
in_trf         = sys.argv[2]
in_jf          = 'high_freq_kmers.txt'
flanking_bp    = int(20000)

isolate = str(in_fa).split('/')[1].split('.')[0]

out_name = str(isolate) + ".csv"


full_fasta_seq = SeqIO.parse(open(in_fa),'fasta')

#List of high frequency kmers
kmers=[]

#idictionary of Tuple that is keyed upon [leftmer,rightmer,tandem repeat] with the corresponding isolate name as the value
d={}


with open(in_jf) as jellyfish:
	for l in jellyfish:
		l = l.strip().split()
		kmers.append(l[0].lower())

def get_seq_btw_2coords(start_coord,end_coord):
#'''Given a fasta sequence, returns the sequence between indices start_coord,end_coord'''
	full_fasta_seq = SeqIO.parse(open(in_fa),'fasta')
	for fasta in full_fasta_seq:
		fa_string = str(fasta.seq)
		subset_fa =fa_string[int(start_coord):int(end_coord)]
		return subset_fa

def get_seq_btw_2kmers(K1,K2):
#'''Given a fasta sequence and two kmers, returns the full insertion sequence between the end of K1 and the start of K2'''
	full_fasta_seq = SeqIO.parse(open(in_fa),'fasta')
	k1=K1.lower()
	k2=K2.lower()
	for fasta in full_fasta_seq:
		fa_string = str(fasta.seq)
		start = fa_string.find(k1)
		end = fa_string.find(k2)
		end_coord = end + len(k2)
		length = end_coord - start
		if length < flanking_bp:
			subset_fa =fa_string[int(start):int(end_coord)]
		else:
			subset_fa = "too long"
	return subset_fa, length

def get_flanking(start_coord,end_coord,flanking_bp):
#'''Given the start and end coordinates of region, find the flanking region with length  X base pairs upstream and downstream'''
	start_flank              = start_coord - flanking_bp
	if start_flank < 1:
		start_flank = 1
	end_flank                = end_coord + flanking_bp
	upstream_flank_seq	 = get_seq_btw_2coords(start_flank,start_coord)
	downstream_flank_seq     = get_seq_btw_2coords(end_coord,end_flank)
	return downstream_flank_seq,upstream_flank_seq

def find_kmer_in_seq(repeat,upstream_flank_seq,downstream_flank_seq,isolate):
#'''given two sequences that flank the region with a TR and the kmers from the list of all high frequency kmers, get the kmer closest to the end of the upstream flanking sequence and the kmer closest to the start of the downstream flanking sequence'''
	downstream_kmers = {}
	upstream_kmers ={}
	for kmer in kmers:
		if upstream_flank_seq.find(kmer) > 0:
			upstream_kmers[kmer]=upstream_flank_seq.find(kmer)
		if downstream_flank_seq.find(kmer) > 0:
			downstream_kmers[kmer]=downstream_flank_seq.find(kmer)
	if bool(upstream_kmers):
		max_upmer = max(upstream_kmers.iteritems(), key=operator.itemgetter(1))[0]
	else:
		max_upmer = "0"
	if bool(downstream_kmers):
		min_downmer = min(downstream_kmers.iteritems(), key=operator.itemgetter(1))[0]
	else: min_downmer = "0"
	d[max_upmer,min_downmer,repeat]=isolate


with open(in_trf,'r') as trf:
        next(trf)
        for l in trf:
		if not l.lstrip().startswith('@'):
                	l=l.strip().split()
                	start = int(l[0])
                	end = int(l[1])
                	repeat = str(l[13])
                	flanks = get_flanking(start,end,flanking_bp)
                	downstream_flank = flanks[0]
                	upstream_flank   = flanks[1]
                	find_kmer_in_seq(repeat,upstream_flank,downstream_flank,isolate)
			my_seq = Seq(repeat)
trf.close()
t=d.items()


with open(out_name, 'w') as a_file:
    for result in t:
	utr = str(result[0][0]) #Kmer closest upstream to TR
	dtr = str(result[0][1]) #Kmer closest downstream to TR
	if utr!="0" and dtr!="0":
		between_2_kmers = get_seq_btw_2kmers(utr,dtr)
		ins_seq = between_2_kmers[0]
		length  = str(between_2_kmers[1])
	else:
		ins_seq = "none"
		length = "0"
	repeat = str(result[0][2]) #TR
	iso = result[1] #Isolate
        final = ','.join([utr,dtr,ins_seq,length,repeat,iso])
        a_file.write(final + '\n')

a_file.close()
