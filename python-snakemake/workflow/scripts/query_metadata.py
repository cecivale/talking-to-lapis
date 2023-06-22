# ------------------------------------------------------------------------------
#          ---        
#        / o o \    Snakemake workflow talking to LAPIS
#        V\ Y /V    Script to query lapis metadata inside snakemake workflow
#    (\   / - \     
#     )) /    |     
#     ((/__) ||     Code by Ceci VA 
# ------------------------------------------------------------------------------

import urllib, json, requests, os, sys, ast
import argparse
import numpy as np
import pandas as pd


from lapis_functions import query_lapis

if __name__ == '__main__':
    parser = argparse.ArgumentParser(
        description="query lapis metadata",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter
    )
    parser.add_argument("--config", type=ast.literal_eval, required=True, help="Config yaml")
    parser.add_argument("--dataset", type=str, required=True, help="Dataset name")
    parser.add_argument("--deme", type=str, required=True, help="Deme name")
    parser.add_argument("--output", type=str, required=True,  help="Output file for metadata")
    args = parser.parse_args()
    

    database = args.config['lapis']["database"]
    endpoint = "details"
    access_key = args.config['lapis']["access_key"]
    
    if args.config['data'][args.dataset][args.deme].get('ids_file') is None:
        attributes = args.config['data'][args.dataset][args.deme]
        data = query_lapis(database, endpoint, attributes, accessKey = access_key)
        data['deme'] = args.deme

        if 'exclude_country' in attributes:
            data = data.query("country not in @attributes['exclude_country']").reset_index()
            
    else:
        df_ids = pd.read_csv(args.config['data'][args.dataset][args.deme]['ids_file'], sep='\t')
        batch_size = 100
        i = 0
        seq_id =  args.config['lapis']['seq_id']
        data = pd.DataFrame()
    
        while i < len(df_ids):
            n = i+batch_size if i+batch_size < len(df_ids) else len(df_ids) 
            ids = df_ids[seq_id][i:n] # TODO use other seq identifier
            seq_names = df_ids[i:n].set_index(seq_id, append = False) 
            ids_dict = dict() 
            ids_dict[seq_id] = ids.to_string(index = False)
            attributes = {**args.config['data'][args.dataset][args.deme],**ids_dict}
            data_query = query_lapis(database, endpoint, attributes, accessKey = access_key)
            data = pd.concat([data, data_query], ignore_index=True)
            i = n
            
        data['deme'] = args.deme
    data.to_csv(args.output, index = False, sep="\t")

