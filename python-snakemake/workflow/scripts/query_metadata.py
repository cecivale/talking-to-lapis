# ------------------------------------------------------------------------------
#          ---        
#        / o o \    Snakemake workflow talking to LAPIS
#        V\ Y /V    Script to query lapis metadata inside snakemake workflow
#    (\   / - \     
#     )) /    |     
#     ((/__) ||     Code by Ceci VA 
# ------------------------------------------------------------------------------

import urllib, json, requests, os, sys, ast
import numpy as np
import pandas as pd


from lapis_functions import query_lapis


if __name__ == '__main__':

    with open(snakemake.log[0], "w") as f:
        sys.stderr = sys.stdout = f
        database = snakemake.params["database"]
        endpoint = "details"
        access_key = snakemake.params["access_key"]


        dataset = snakemake.wildcards["dataset"]
        deme = snakemake.params["deme"]
        structured = True if deme is not None else False
        attributes = snakemake.config['datasets'][dataset]["structure"][deme]["select"] if structured else snakemake.config['datasets'][dataset]["select"]
        seq_id =  snakemake.params['seq_id']
        

        if attributes.get('ids_file') is None:
            attributes_ext = attributes.copy()
            if attributes.get('add_query') is not None:
                attributes_ext.pop("add_query")
            metadata = query_lapis(database, endpoint, attributes_ext, accessKey = access_key)
            if structured: metadata['deme'] = deme

        else:
            df_ids = pd.read_csv(attributes.get('ids_file'), sep='\t')

            batch_size = 50
            i = 0
            metadata = pd.DataFrame()

 
            while i < len(df_ids):

                n = i + batch_size if i + batch_size < len(df_ids) else len(df_ids) 

                ids = df_ids[seq_id][i:n] 

                seq_names = df_ids[i:n].set_index(seq_id, append = False) 
                ids_dict = dict() 
                ids_dict[seq_id] = ids.to_list()

                attributes_ext = {**attributes, **ids_dict}.copy()
                attributes_ext.pop("ids_file")
                if attributes.get('add_query') is not None:
                    attributes_ext.pop("add_query")

                data_query = query_lapis(database, endpoint, attributes_ext, accessKey = access_key)

                metadata = pd.concat([metadata, data_query], ignore_index=True)

                i = n

                
            if structured: metadata['deme'] = deme


        if attributes.get('add_query') is not None:
            metadata = metadata.query(attributes.get('add_query')).reset_index()
        metadata["seq_id"] = metadata[seq_id]

        # assign sample id
        ids = pd.DataFrame()
        if "sample_id" in metadata:
            ids["sample_id"] = metadata["sample_id"]
        elif "seq_id" in metadata and "date" in metadata and snakemake.params["deme"] is not None:
            ids["sample_id"] = metadata["seq_id"] + "|" + snakemake.params["deme"] + "|" + metadata["date"]
            metadata["sample_id"] = metadata["seq_id"] + "|" + snakemake.params["deme"] + "|" + metadata["date"]
        elif "seq_id" in metadata and "date" in metadata:
            ids["sample_id"] = metadata["seq_id"] + "|" + metadata["date"]
            metadata["sample_id"] = metadata["seq_id"] + "|" + metadata["date"]
        elif "seq_id" in metadata:
            ids["sample_id"] = metadata["seq_id"]
            metadata["sample_id"] = metadata["seq_id"]
        else:
            print("Error: No sequence or sample id in metadata.")
       
        if "date" in metadata:
            ids["date"] = metadata["date"]
        if "seq_id" in metadata:
            ids["seq_id"] = metadata["seq_id"]
        else:
            ids["seq_id"] = metadata[[snakemake.params["seq_id"]]]
        if "strain" in metadata:
            ids["strain"] = metadata["strain"]
        if snakemake.params["deme"] is not None: 
            ids[["deme"]] = snakemake.params["deme"] 
            metadata[["deme"]] = snakemake.params["deme"] 



        metadata.to_csv(snakemake.output["metadata"], index = False, sep="\t")
        ids.to_csv(snakemake.output["ids"], index = False, sep="\t")


