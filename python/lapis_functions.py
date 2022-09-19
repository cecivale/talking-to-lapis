# ------------------------------------------------------------------------------
#          ---        
#        / o o \    Project:  cov-euphydyn
#        V\ Y /V    Functions to query lapis metadata and sequences 
#    (\   / - \     29 July 2022
#     )) /    |     
#     ((/__) ||     Code by Ceci VA 
# ------------------------------------------------------------------------------

import urllib, json, requests, os, sys
from Bio import SeqIO
from io import StringIO
import pandas as pd


def attr_dict_to_lapis(attr_dict):
    s = str(attr_dict)
    lapis_s = s.replace(
        "{", "").replace(
        "}", "").replace(
        "'", "").replace(
        ",", "&").replace(
        ":", "=").replace(
        " ", "").replace(
        "\\n", ",").replace(
        "UnitedKingdom", "United Kingdom")
    
    return lapis_s


def lapis_build_query(database, endpoint, attributes, fields = None, accessKey = None, version = "v1"):
    url = "https://lapis.cov-spectrum.org/" + \
        database + "/" + \
        version + "/sample/" + \
        endpoint + "?" + \
        (("fields=" + fields) if fields is not None else "") + \
        attr_dict_to_lapis(attributes) + \
        (("&accessKey=" + accessKey) if accessKey is not None else "")
    print("LAPIS Query: ", url, "\n")    
    return url


def query_lapis(database, endpoint, attributes, **kwargs):

        print("Building query...\n")
        url = lapis_build_query(database, endpoint, attributes, **kwargs)
        print("Sending request to LAPIS...\n")
        req = requests.get(url, headers={'User-Agent': 'Mozilla/5.0'})  
        print("Response received from LAPIS!\n")
        # TODO manage LAPIS error
        
        if endpoint == "fasta-aligned":
            fasta_iterator = SeqIO.parse(StringIO(req.text), "fasta")
            return fasta_iterator
        else:
            json_data = req.json()
            data = pd.DataFrame.from_dict(json_data['data'])
            return data


def query_lapis_metadata(config, deme, output, log_file):
    database = config['lapis']["database"]
    endpoint = "details"
    attributes = config['data'][deme]
    access_key = config['lapis']["access_key"]
    
    sys.stdout = open(str(log_file), "w")
    data = query_lapis(database, endpoint, attributes, accessKey = access_key)
    data['deme'] = deme

    if 'exclude_country' in attributes:
        data = data.query("country not in @attributes['exclude_country']").reset_index()

    data.to_csv(output, index = False, sep="\t")


def query_lapis_aligned_seqs(config, ids_file, output_file, log_file, drop_incomplete = True, seq_identifier = "gisaidEpiIsl", batch_size = 50):
    database = config['lapis']["database"]
    endpoint = "fasta-aligned"
    access_key = config['lapis']["access_key"]
    df_ids = pd.read_csv(ids_file, sep='\t')
    i = 0
    
    sys.stdout = open(str(log_file), "w")
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
                if drop_incomplete and len(record.seq.replace("-","")) < 29000: continue # Write only full SARS-CoV-2 genomes
                record.id = seq_names.loc[record.id]['seq_name']
                record.description = ""
                SeqIO.write(record, corrected, 'fasta')
        print("Saved " + str(i) + " to " + str(n-1) + " sequences")
        i = n
    os.remove(temp)




