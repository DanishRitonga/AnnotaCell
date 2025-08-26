import pandas as pd
import json

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

def translate_to_gene_ontology(cluster_dict, dictfile):
    with open(dictfile, 'r') as f:
        gene_ontology_dict = json.load(f)
    