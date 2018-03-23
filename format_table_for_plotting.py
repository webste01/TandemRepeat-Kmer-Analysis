import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
import sys
import os

file = sys.argv[1]
trf=pd.read_csv(file)
out_name=sys.argv[2]

k1k2_alleles=trf.groupby(['k1_k2'])['insertion_sequence'].nunique()
k1k2_alleles.name = 'alleles'

k1k2_alleles = pd.DataFrame(k1k2_alleles)
k1k2_alleles['k1_k2'] = k1k2_alleles.index
trf_w_alleles = pd.merge(trf, k1k2_alleles, on='k1_k2')

#Read in mlst file
path = '/sc/orga/projects/InfectiousDisease/studies/CDiff_first_paper/final_TRF_pipeline/scripts/mlst_file.txt'
st_iso_df = pd.read_table(path, dtype={'Assembly_ID': str}, sep='\t',header='infer')


def isolate_to_assembly_id (isolate):
  return isolate.split("_")[-1]


trf_w_alleles['Assembly_ID'] = trf_w_alleles.isolate.apply(isolate_to_assembly_id)


trf_w_alleles_st = pd.merge(trf_w_alleles, st_iso_df, on='Assembly_ID')


#Group by allele instead of insertion_sequence_length
trf_mlst_alleles = trf.groupby([trf_w_alleles_st.k1_k2, 
                                trf_w_alleles_st.insertion_sequence,
                                trf_w_alleles_st.MLST], as_index=False)
trf_mlst_alleles_size = pd.DataFrame({'k1_k2_mlst_alleles': trf_mlst_alleles.size()}).reset_index()
trf_w_alleles_st_count = pd.merge(trf_w_alleles_st, trf_mlst_alleles_size, on=['k1_k2', 'insertion_sequence', 'MLST'])

trf_w_alleles_st_count['allele'] = trf_w_alleles_st_count['insertion_sequence']
trf_w_alleles_st_count.allele = pd.Categorical(trf_w_alleles_st_count.allele)
trf_w_alleles_st_count['allele_mapping'] = trf_w_alleles_st_count.allele.cat.codes

#Add length to allele_mapping
trf_w_alleles_st_count["length_allele_mapping"] = trf_w_alleles_st_count["insertion_sequence_length"].map(str) + "_" + trf_w_alleles_st_count["allele_mapping"].map(str)

trf_w_alleles_st_count.to_csv(out_name, sep=',')

