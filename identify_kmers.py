import sys
from collections import Counter

fn = sys.argv[1]
k  = int(sys.argv[2])
outfile_name = str(fn) + "_kmers.txt"
outfile = open(outfile_name,"w")

# hack to read in file                                                                                                                                                    
seq = ""
with open(fn) as f:
	seq = f.read()

seq = "".join(seq.split("\n")[1:])
seq = seq.upper()

kmer_counter = Counter()

# get forward kmers                                                                                                                                                       
for i in range(len(seq)):
	kmer = seq[i:i+k]
	outfile.write(kmer.lower())
	outfile.write('\n')
	kmer_counter[kmer] += 1

# reverse and add reverse kmers                                                                                                                                           
rc_map = dict(zip('ACGT', 'TGCA'))

def rc (s):
	return ''.join(map(lambda x: rc_map[x], s[::-1]))

rc_seq = rc(seq)
for i in range(len(rc_seq)):
	kmer = rc_seq[i:i+k]
	outfile.write(kmer.lower())
	outfile.write('\n')
	kmer_counter[kmer] += 1

print "number of kmers identified:"
print len(kmer_counter)
