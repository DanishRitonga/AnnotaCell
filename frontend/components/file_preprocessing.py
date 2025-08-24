import pandas as pd

def check_columns(df: pd.DataFrame) -> None:
    df.columns = df.columns.str.lower()
    if "cluster" not in df.columns or "gene" not in df.columns:
        raise ValueError("The input File must contain 'Cluster' and 'Gene' columns.")

def extract_clusters(df: pd.DataFrame):
    # This will raise an error and stop execution if columns are missing.
    check_columns(df) 
    grouped = df.groupby("Cluster")["Gene"].apply(list)
    cluster_dict = grouped.to_dict()
    
    return cluster_dict