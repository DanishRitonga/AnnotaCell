import pandas as pd
import csv

def csv_read(file):
    file.seek(0)
    sample_bytes = file.read(1024)
    sample_str = sample_bytes.decode('utf-8', errors='ignore')
    
    try:
        sep = csv.Sniffer().sniff(sample_str).delimiter
    except csv.Error:
        sep = ',' # Default to comma if sniffing fails

    file.seek(0) 
    df = pd.read_csv(file, sep=sep)
    
    return df