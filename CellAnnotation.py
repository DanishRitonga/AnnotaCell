import pandas as pd 
import requests, json, csv

#Input and Preprocess Data
def csv_read(file_path):
    with open(file_path, 'r') as file:
        sample = file.read(1024) 
        file.seek(0)
        detected_delimiter = csv.Sniffer().sniff(sample).delimiter    
    df = pd.read_csv(file_path, sep=detected_delimiter)
    return df

def extract_clusters(df):
    # Ensure the column names are correct
    if "Cluster" not in df.columns or "Gene" not in df.columns:
        raise ValueError("The file must contain 'Cluster' and 'Gene' columns")

    # Group the data by clusters and pad genes to ensure columns have equal lengths
    grouped = df.groupby("Cluster")["Gene"].apply(list)

    # Create a DataFrame where each column is a cluster and rows contain the genes
    max_length = max(grouped.apply(len))  # Determine the maximum number of genes in a cluster
    return pd.DataFrame({cluster: genes + [None] * (max_length - len(genes)) for cluster, genes in grouped.items()})

def translate_genes(clustered_df):
    # Load the Ensembl ID mapping from JSON
    with open('dictionary.json', 'r') as file:
        gene_to_ensembl = json.load(file)["ensembl_to_gene"]

    ensembl_to_gene = {v: k for k, v in gene_to_ensembl.items()}

    # Translate gene names to Ensembl IDs in the clustered DataFrame
    return clustered_df.apply(lambda col: col.map(lambda gene: ensembl_to_gene.get(gene, None) if pd.notna(gene) else None))

def data_preprocessing(df):
    extracted_df = extract_clusters(df)
    return translate_genes(extracted_df)

#Process the Data
def fetch_expression_data(cluster, ensembl_ids):
    payload = {
        "filter": {
            "dataset_ids": [],
            "development_stage_ontology_term_ids": [],
            "disease_ontology_term_ids": [],
            "gene_ontology_term_ids": ensembl_ids,
            "organism_ontology_term_id": "NCBITaxon:9606",
            "self_reported_ethnicity_ontology_term_ids": [],
            "sex_ontology_term_ids": [],
            "publication_citations": [],
        },
        "is_rollup": True
    }

    # URL for the POST request
    API_URL = "https://api.cellxgene.cziscience.com/wmg/v2/query"

    try:
        response = requests.post(API_URL, json=payload)
        response.raise_for_status()  # Raise an exception for HTTP errors
        return response.json()['expression_summary']
    except requests.RequestException as e:
        print(f"Error fetching data for Cluster {cluster}: {e}")
        return None
    
def expression_data_to_df(data):
    flattened_data = []
    for gene_id, anatomical_structures in data.items():
        for anatomical_id, cell_types in anatomical_structures.items():
            for cell_type_id, aggregated_data in cell_types.items():
                metrics = aggregated_data['aggregated']
                flattened_data.append({
                    'gene': gene_id,
                    'tissue': anatomical_id,
                    'cell': cell_type_id,
                    'expression': metrics.get('me', None),
                    'cell count': metrics.get('n', None),
                    'cell percentage': metrics.get('pc', None),
                    'tissue composition': metrics.get('tpc', None)
                })
    return pd.DataFrame(flattened_data)

def filter_expression_data(response_df, TARGET_UBERON_ID):
    tissue_df = response_df[response_df['tissue'] == TARGET_UBERON_ID].drop(columns=['tissue'])
    return tissue_df[~tissue_df['cell'].isin(['tissue_stats', 'CL:0000000'])]

def translate_ontology(df):
    with open('cell_dict.json', 'r') as cell_file, open('ensembl_to_gene.json', 'r') as gene_file:
        cl_to_cell = json.load(cell_file)
        ensembl_to_gene = json.load(gene_file)

    df['cell'] = df['cell'].map(cl_to_cell)
    df['gene'] = df['gene'].map(ensembl_to_gene).fillna(df['gene'])

    return df

def transform_results(results):
    # Count the number of entries for each cell
    cell_counts = results['cell'].value_counts().reset_index()
    cell_counts.columns = ['cell', 'count']

    # Merge the counts back to the original DataFrame
    merged_df = results.merge(cell_counts, on='cell')

    # Sort the DataFrame based on the count of entries in descending order
    return merged_df.groupby(['cell', 'gene']).first().sort_values(by=['count', 'cell', 'expression', 'cell percentage'], ascending=[False, True, False, False]).reset_index().drop(columns=['count'])

def calculate_cell_score(results):
    cell_df = results.copy()

    # Group by 'cell' and calculate the total cell count for each group
    sum_cell_count = cell_df.groupby('cell')['cell count'].sum().reset_index()
    sum_cell_count.columns = ['cell', 'total cell count']
    cell_df = cell_df.merge(sum_cell_count, on='cell')

    cell_df = cell_df[cell_df['total cell count'] >= 100]

    # Calculate scores
    cell_df['score'] = (1) * ((cell_df['expression'] - cell_df['expression'].min()) / (cell_df['expression'].max() - cell_df['expression'].min())) + (1.5) * (cell_df['cell percentage']) + (2.5) * ((cell_df['cell count'] - cell_df['cell count'].min()) / (cell_df['cell count'].max() - cell_df['cell count'].min()))

    # Group by 'cell' and calculate the standard deviation of 'score' for each group
    std_scores = cell_df.groupby('cell')['score'].std().reset_index()
    std_scores.columns = ['cell', 'std']
    cell_df = cell_df.merge(std_scores, on='cell')

    # Group by 'cell' and calculate the mean of 'score' for each group
    mean_scores = cell_df.groupby('cell')['score'].mean().reset_index().fillna(0)
    mean_scores.columns = ['cell', 'mean']
    cell_df = cell_df.merge(mean_scores, on='cell')

    # Calculate the coefficient of variation (CV) for each cell
    cell_df['CV'] = cell_df['std'] / cell_df['mean']

    # Group by 'cell' and calculate the kurtosis of 'score' for each group
    kurtosis_scores = cell_df.groupby('cell')['score'].apply(pd.Series.kurt).reset_index()
    kurtosis_scores.columns = ['cell', 'kurtosis']
    cell_df = cell_df.merge(kurtosis_scores, on='cell')

    # Normalize the scores
    cell_df['std norm'] = (cell_df['std'] - cell_df['std'].min()) / (cell_df['std'].max() - cell_df['std'].min())
    cell_df['mean norm'] = (cell_df['mean'] - cell_df['mean'].min()) / (cell_df['mean'].max() - cell_df['mean'].min())
    cell_df['CV norm'] = (cell_df['CV'] - cell_df['CV'].min()) / (cell_df['CV'].max() - cell_df['CV'].min())
    cell_df['kurtosis norm'] = ((cell_df['kurtosis'] - cell_df['kurtosis'].min()) / (cell_df['kurtosis'].max() - cell_df['kurtosis'].min()))

    # Count the number of entries for each cell
    cell_counts = cell_df['cell'].value_counts().reset_index()
    cell_counts.columns = ['cell', 'count']
    cell_df = cell_df.merge(cell_counts, on='cell')

    cell_df['cell score'] = (0.1 * (1 - cell_df['std norm']).fillna(0) + 0.7 * cell_df['mean norm'] + 0.1 * (1 - cell_df['CV norm']) + 0.1 * (1 - cell_df['kurtosis norm']).fillna(0)) * cell_df['count']
    return cell_df.sort_values(by=['count', 'cell score', 'cell', 'gene'], ascending=[False, False, False, True]).reset_index(drop=True)

def run_data_processing(cluster, ensembl_ids, target_uberon_id):
    try:
        # Step 1: Fetch expression data
        data = fetch_expression_data(cluster, ensembl_ids)

        # Step 2: Convert data to DataFrame
        response_df = expression_data_to_df(data)

        # Step 3: Filter data based on target UBERON ID
        filtered_df = filter_expression_data(response_df, target_uberon_id)

        # Step 4: Translate ontology terms
        translated_df = translate_ontology(filtered_df)

        # Step 5: Transform the results
        transformed_df = transform_results(translated_df)

        # Step 6: Rank the cells based on the calculated score
        ranked_df = calculate_cell_score(transformed_df)

        return ranked_df

    except Exception as e:
        print(f"An error occurred during the pipeline execution: {e}")
        return pd.DataFrame()
    
def main_data_analysis(ensembl_df, target_uberon_id):
    results = {}

    # Iterate over each cluster (column) in the DataFrame
    for cluster in ensembl_df.columns:
        # Extract Ensembl IDs for the current cluster, dropping any NaN values
        ensembl_ids = ensembl_df[cluster].dropna().tolist()
        
        if not ensembl_ids:
            continue
        
        # Run the expression pipeline for the current cluster
        result_df = run_data_processing(cluster, ensembl_ids, target_uberon_id)
        
        # Store the result in the dictionary
        results[cluster] = result_df.reset_index(drop=True)

    return results

def export_dataframe(df):
    unique_cell = df['cell'].unique()[:20]
    final_df = pd.DataFrame(unique_cell)
    final_df.columns = ['cell']
    genes, cell_score = [], []
    for cell in unique_cell:
        gene_list = df['gene'][df['cell'] == cell].unique()
        score_list = df['cell score'][df['cell'] == cell].unique()
        
        genes.append(", ".join(map(str, gene_list)))
        cell_score.append(", ".join(map(str, score_list)))
    final_df['gene'], final_df['cell score'] = genes, cell_score
    return final_df

def get_top3_cells(df):
    return df['cell'].unique()[:3]

def export_top3_cells(results):

    top3cells = pd.DataFrame()
    for cluster in results:
        top3cells[f'cluster {cluster}'] = get_top3_cells(results[cluster])

    top3cellsT = top3cells.T.reset_index()
    top3cellsT.columns = ['Cluster'] + list(top3cells.index)
    top3cellsT.columns = ['Cluster', 'Cell 1', 'Cell 2', 'Cell 3']

    return top3cellsT

# file_path = "10-top-cluster.csv"
# df = csv_read(file_path)
# clustered_df = extract_clusters(df)
# ensembl_df = translate_genes(clustered_df)

# # Specify the UBERON ID to filter
# TARGET_UBERON_ID = "UBERON:0002097"
# results = main_data_analysis(ensembl_df, TARGET_UBERON_ID)
# final_df = export_dataframe(calculate_cell_score(23, results))

# print(final_df)