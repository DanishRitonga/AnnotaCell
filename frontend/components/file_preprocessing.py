import pandas as pd
import json
from components.ontology_translator import ontl

def check_columns(df: pd.DataFrame) -> None:
    df.columns = df.columns.str.lower()
    if "cluster" not in df.columns or "gene" not in df.columns:
        raise ValueError("The input File must contain 'Cluster' and 'Gene' columns.")

def extract_clusters(df: pd.DataFrame):
    # This will raise an error and stop execution if columns are missing.
    check_columns(df) 
    grouped = df.groupby("cluster")["gene"].apply(list)
    cluster_dict = grouped.to_dict()
    
    return cluster_dict

def translate_to_ensembl(cluster_dict):
    ontology_dict = {}
    for cluster_id, name_list in cluster_dict.items():
        ontology_id_list = []
        for name in name_list:
            ontology_id = ontl.gene(name, rev=True)
            ontology_id_list.append(ontology_id)
        
        ontology_dict[cluster_id] = ontology_id_list

    return ontology_dict    