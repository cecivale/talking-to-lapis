# ------------------------------------------------------------------------------
#          ---        
#        / o o \    Snakemake workflow talking to LAPIS
#        V\ Y /V    Script to query lapis sequences inside snakemake workflow
#    (\   / - \     
#     )) /    |     
#     ((/__) ||     Code by Ceci VA 
# ------------------------------------------------------------------------------

import urllib, json, requests, os, sys, ast
import numpy as np
import pandas as pd

from Bio import SeqIO
from io import StringIO


from lapis_functions import query_lapis

if __name__ == '__main__':
    database = snakemake.params["database"] 
    endpoint = "fasta-aligned"
    # access_key = ""
    df_ids = pd.read_csv(snakemake.input["ids"], sep='\t')
    batch_size = snakemake.params["batch_size"]
    output_file = snakemake.output["alignment"]
    seq_id = snakemake.params["seq_id"] 

    i = 0

    while i < len(df_ids):
        n = i + batch_size if i + batch_size < len(df_ids) else len(df_ids) 
        ids = df_ids["seq_id"][i:n] 
        seq_names = df_ids[i:n].set_index("strain", append = False) 
        attributes = dict() 
        attributes[seq_id] = ids.to_string(index = False)
        data = query_lapis(database, endpoint, attributes) #, accessKey = access_key)
        temp = output_file.replace(".fasta", "_temp.fasta")
        SeqIO.write(data, temp, "fasta")

        with open(temp) as original, open(output_file, 'a') as corrected:
            records = SeqIO.parse(original, 'fasta')
            for record in records:
                # TODO QC in another step
                if snakemake.params["drop_incomplete"] and len(record.seq.replace("-","")) < 29000: continue # Write only full SARS-CoV-2 genomes
                record.id = seq_names.loc[record.id]['sample_id']
                record.description = ""
                SeqIO.write(record, corrected, 'fasta')
        print("Saved " + str(i) + " to " + str(n-1) + " sequences")
        i = n
    os.remove(temp)

