from collections import Counter
import sys
import pandas as pd

table = 'table4plotting.csv'
df = pd.read_csv(table)
 df.filter(items=['one', 'three'])

sub_df = df.filter(items=['k1_k2','insertion_sequence','insertion_sequence_length','isolate','MLST','allele_mapping'])

def create_allele_mapping_frac (df, total_genomes=None, total_mlst_types=None):
  if total_genomes is None:
    total_genomes = len(df)
  if total_mlst_types is None:
    total_mlst_types = len(df.MLST.unique())
  
  # build counters
  allele_mapping_dict = {}
  for index, row in df.iterrows():
    allele = row['insertion_sequence']
    mlst = row['MLST']
    allele_mapping_dict.setdefault(allele, Counter())[mlst] += 1
  
  # Get max MLST
  allele_to_mlst = {}
  for allele, mlst_counter in allele_mapping_dict.items():
    mlst = mlst_counter.most_common(1)[0][0]
    allele_to_mlst[allele] = mlst
    
  # get fraction
  allele_mapping_dict = {}
  correct = 0
  for index, row in df.iterrows():
    allele = row['insertion_sequence']
    mlst = row['MLST']
    if allele_to_mlst[allele] == mlst:
      correct += 1
  return float(correct) / float(total_genomes)


pair_scores = []
for mer_pair in df['k1_k2'].unique():
  tmp_df = df[df['k1_k2'] == mer_pair]
  pair_scores.append((mer_pair, create_allele_mapping_frac(tmp_df, total_genomes=248)))
  
pair_scores.sort(lambda x, y: cmp(y[1], x[1]))

for pair_score in pair_scores[0:50]:
  print pair_score
