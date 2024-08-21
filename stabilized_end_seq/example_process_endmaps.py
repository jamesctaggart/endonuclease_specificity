'''
This script takes EndMap objects uses the, to generate
a summary DataFrame that can easily be filtered downstream to identify cleavage sites. 
These DataFrames are pickled for use in call_sites_endseq_RR_analysis.py

This script does this for multiple reference backgrounds and associated knockout datasets
'''

import pickle
import pandas as pd
from Bio import SeqIO
from endmap.calc_endseq_columns import build_endseq_table

# region Import references (genome seq, mochiview)
# Import genome/annotation
BS168_genome = SeqIO.to_dict(SeqIO.parse("reference_files/B_subtilis_168_NC_000964.3.fa", "fasta"))['NC_000964.3'].seq
mochi_file = 'reference_files/CDS_168.mochiview_v13.txt'
# endregion

# Define common parameters to use throughout analysis
data_column_to_use = 'cds_normalized'
skip_thresh = 1
group_thresh = 7.5
pseudocount_val = 0.001

# region build the DataFrames
print('Building table...')
# Each em corresponds to a different EndMap object generated from a pair of .wig files.
# Ref refers to a strain with no endonucleases knocked out, ko refers to dataset where a single endoribonuclease (e.g. RNase Y) is ablated.
# A single reference might be associated with multiple knockouts, in which case those EndMaps can be given as a list here, each with a different "ko_id" associated.
# The final table will include an entry for each position/strand for each knockout
ses_df = build_endseq_table(ref_5end_em=ref_5end_em, ko_5end_ems=ko_5end_ems, ref_rend5_em=ref_rend5_em, \
     ko_rend5_ems=ko_rend_ems, ref_rend3_em=ref_rend3_em, ref_3end_em=em_3endseq, ko_ids=ko_ids, genome_seq=BS168_genome,group_metric='ref_bg_ratio',\
        group_thresh=group_thresh, data_column=data_column_to_use, skip_thresh=skip_thresh, pseudocount_val=pseudocount_val)
# endregion

with open(pickle/file/path,'wb') as f:
    pickle.dump(ses_df,f)
