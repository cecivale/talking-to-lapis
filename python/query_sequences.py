# ------------------------------------------------------------------------------
#          ---        
#        / o o \    Project:  cov-euphydyn
#        V\ Y /V    Script to query lapis sequences inside snakemake workflow
#    (\   / - \     29 August 2022
#     )) /    |     
#     ((/__) ||     Code by Ceci VA 
# ------------------------------------------------------------------------------

import urllib, json, requests, os, sys, ast
import argparse
import numpy as np
import pandas as pd

from Bio import SeqIO
from io import StringIO


from lapis_functions import query_lapis

if __name__ == '__main__':
    parser = argparse.ArgumentParser(
        description="query lapis sequences",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter
    )
    parser.add_argument("--config", type=ast.literal_eval, required=True, help="Config yaml")
    parser.add_argument("--ids_file", type=str, required=True, help="File with sequences id to download")
    parser.add_argument("--output_file", type=str, required=True,  help="Output file for sequences")
    parser.add_argument("--drop_incomplete", default=True, help="Drop <29000 length sequences")
    parser.add_argument("--seq_identifier", default="gisaidEpiIsl", help="Sequence identifier used in ids file")
    parser.add_argument("--batch_size", default=50, help="Lapis donwload batch size")
    args = parser.parse_args()
    
    database = args.config['lapis']["database"]
    endpoint = "fasta-aligned"
    access_key = args.config['lapis']["access_key"]
    df_ids = pd.read_csv(args.ids_file, sep='\t')
    batch_size = args.batch_size
    output_file = args.output_file

    i = 0
    
    while i < len(df_ids):
        n = i+batch_size if i+batch_size < len(df_ids) else len(df_ids) 
        ids = df_ids["gisaidEpiIsl"][i:n] # TODO use other seq identifier
        seq_names = df_ids[i:n].set_index("gisaidEpiIsl", append = False) 
        attributes = dict(gisaidEpiIsl = ids.to_string(index = False)) 
        data = query_lapis(database, endpoint, attributes, accessKey = access_key)
        temp = output_file.replace(".fasta", "_temp.fasta")
        SeqIO.write(data, temp, "fasta")
        with open(temp) as original, open(output_file, 'a') as corrected:
            records = SeqIO.parse(original, 'fasta')
            for record in records:
                if args.drop_incomplete and len(record.seq.replace("-","")) < 29000: continue # Write only full SARS-CoV-2 genomes
                record.id = seq_names.loc[record.id]['seq_name']
                record.description = ""
                SeqIO.write(record, corrected, 'fasta')
        print("Saved " + str(i) + " to " + str(n-1) + " sequences")
        i = n
    os.remove(temp)


