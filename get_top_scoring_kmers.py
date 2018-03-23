import pandas as pd
import numpy as np
import sklearn
from sklearn import linear_model, datasets
from sklearn import ensemble
from sklearn.model_selection import cross_val_score
from sklearn.model_selection import ShuffleSplit
from sklearn.model_selection import train_test_split
import matplotlib.pyplot as plt
import operator
from sklearn.tree import DecisionTreeClassifier
from sklearn.metrics import confusion_matrix
from sklearn.metrics import roc_curve, auc
from sklearn.model_selection import ShuffleSplit
from sklearn import preprocessing
from sklearn import ensemble
import os
import sys
import bisect
from Bio import SeqIO
from Bio.Seq import Seq
import math

table = sys.argv[1]
percent_threshold = float(sys.argv[2])
score_threshold = float(sys.argv[3])
df = pd.read_csv(table)

header = 'row_number,k1_k2,k1_start,k2_end,allele,length,genome,orientation,num_of_alleles,assembly_id,mlst,k1_k2_mlst_alleles,allele,allele_mapping,length_with_allele_mapping'
header_list = header.split(',')
df.columns = header_list

#Get number of samples that are needed to hit cutoff threshold
n_total = len(df['assembly_id'].unique())
n = int(math.ceil(n_total*percent_threshold))


#Run a decision tree classifier to rank the top kmer candidates
dt_scores = {} 
for pair in df.k1_k2.unique():
  df_sub = df.loc[df['k1_k2'] == pair]
  if len(df_sub['assembly_id'].unique()) > n:
    if len(df_sub['mlst'].unique()) >1:
      dtree = DecisionTreeClassifier()
      X = df_sub[['allele_mapping']]
      y = df_sub.mlst
      dtree.fit(X,y)
      dt_scores[pair]=dtree.score(X,y)


sorted_x = sorted(dt_scores.items(), key=operator.itemgetter(1))
sorted_scores = pd.DataFrame.from_records(sorted_x,columns=['Kmer','score'])
sorted_scores.to_csv('all_scores_for_kmers.csv',sep=',',index=False)

sorted_scores_top = sorted_scores.loc[sorted_scores['score'] > score_threshold]
sorted_scores_top.to_csv('top_scoring_kmers.csv',sep=',',index=False)

sorted_scores_top.Kmer.to_csv('kmers_of_interest.txt',index=False)



