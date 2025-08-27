import pandas as pd
import numpy as np
import os, requests, json, csv
import psycopg2
import os
from openai import OpenAI
from tavily import TavilyClient
from dotenv import load_dotenv

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
    df['cell'] = df['cell'].apply(lambda x: ontl.cell(x))
    df['gene'] = df['gene'].apply(lambda x: ontl.gene(x)).fillna(df['gene'])

    return df

def calculate_cell_score(
    results,
    w_e=1.0, w_p=1.5, w_ct=2.5,             # base score weights
    w_m=7.0, w_f=3.0,                       # mean vs flatness weights
    alpha=1.0, beta=1.0, gamma=1.0,         # flatness breakdown
    w_c=1.0                                 # count weight
):

    cell_df = results.copy()

    # Calculate total cell count per cell
    sum_cell_count = cell_df.groupby('cell')['cell count'].sum().reset_index()
    sum_cell_count.columns = ['cell', 'total cell count']
    cell_df = cell_df.merge(sum_cell_count, on='cell')

    # Adaptive cutoff: max(μ - 0.5σ, min(total count))
    mu = cell_df['total cell count'].mean()
    sigma = cell_df['total cell count'].std()
    minimum = cell_df['total cell count'].min()
    min_count = max(mu - 0.5 * sigma, minimum)
    cell_df = cell_df[cell_df['total cell count'] >= min_count]

    # Normalize expression and cell count
    norm_expr = (cell_df['expression'] - cell_df['expression'].min()) / (cell_df['expression'].max() - cell_df['expression'].min())
    norm_count = (cell_df['cell count'] - cell_df['cell count'].min()) / (cell_df['cell count'].max() - cell_df['cell count'].min())

    # Compute base score
    cell_df['score'] = (
        w_e * norm_expr +
        w_p * cell_df['cell percentage'] +
        w_ct * norm_count
    )

    # Compute flatness components
    std_scores = cell_df.groupby('cell')['score'].std().reset_index().rename(columns={'score': 'std'})
    mean_scores = cell_df.groupby('cell')['score'].mean().reset_index().rename(columns={'score': 'mean'})
    kurtosis_scores = cell_df.groupby('cell')['score'].apply(pd.Series.kurt).reset_index(name='kurtosis')

    cell_df = cell_df.merge(std_scores, on='cell')
    cell_df = cell_df.merge(mean_scores, on='cell')
    cell_df = cell_df.merge(kurtosis_scores, on='cell')

    cell_df['CV'] = cell_df['std'] / (cell_df['mean'].replace(0, np.nan))

    # Normalize flatness metrics
    for col in ['std', 'CV', 'kurtosis', 'mean']:
        col_min, col_max = cell_df[col].min(), cell_df[col].max()
        cell_df[f'{col} norm'] = (cell_df[col] - col_min) / (col_max - col_min + 1e-9)

    # Stability (flatness) score
    stability = (
        alpha * (1 - cell_df['std norm'].fillna(0)) +
        beta * (1 - cell_df['CV norm'].fillna(0)) +
        gamma * (1 - cell_df['kurtosis norm'].fillna(0))
    )
    flatness_weight = stability / (alpha + beta + gamma)

    # Count-based log-scaled soft weight
    cell_counts = cell_df['cell'].value_counts().reset_index()
    cell_counts.columns = ['cell', 'count']
    cell_df = cell_df.merge(cell_counts, on='cell')

    log_scaled_count = np.log1p(cell_df['count']) / np.log1p(cell_df['count'].max())

    # Final cell score
    cell_df['cell score'] = (
        (w_m * cell_df['mean norm'] + w_f * flatness_weight) *
        (1 + w_c * log_scaled_count)
    )

    return cell_df.sort_values(by=['count', 'cell score', 'cell', 'gene'], ascending=[False, False, False, True]).reset_index(drop=True)

def run_data_processing(cluster, ensembl_ids, target_uberon_id):
    try:
        data = fetch_expression_data(cluster, ensembl_ids)

        response_df = expression_data_to_df(data)

        filtered_df = filter_expression_data(response_df, target_uberon_id)

        translated_df = translate_ontology(filtered_df)

        ranked_df = calculate_cell_score(translated_df)

        return ranked_df

    except Exception as e:
        print(f"An error occurred during the pipeline execution: {e}")
        return pd.DataFrame()
    
def main_data_analysis(cluster_dict, target_uberon_id):
    results = {}

    # Iterate over each cluster (column) in the DataFrame
    for cluster, name_list in cluster_dict.items():
        # Extract Ensembl IDs for the current cluster, dropping any NaN values
        ensembl_ids = name_list.dropna()
        
        if not ensembl_ids:
            continue
        
        # Run the expression pipeline for the current cluster
        result_df = run_data_processing(cluster, ensembl_ids, target_uberon_id)
        
        # Store the result in the dictionary
        results[cluster] = result_df.reset_index(drop=True)

    return results

def export_analysis_dataframe(results):
    ex_df = results.copy()
    for cluster in ex_df:
        unique_cells = ex_df[cluster]['cell'].unique()[:10]
        ex_df[cluster] = ex_df[cluster][['gene', 'cell', 'cell score']][ex_df[cluster]['cell'].isin(unique_cells)]
    return ex_df


#Validation
def initialize_tavily():
    # Initialize Tavily and Gemini clients
    return TavilyClient(api_key=os.getenv("TAVILY_API_KEY"))

def initialize_openai():
    return OpenAI(api_key=os.getenv("OPENAI_API_KEY"))

# Connect to PostgreSQL
def get_db_connection():
    return psycopg2.connect(os.getenv("DATABASE_URL"))

def check_database(tissue_name, cell_name, gene_name):
    query = """
        SELECT relationcertainty FROM gene_expression
        WHERE uberon = (
            SELECT uberon 
            FROM tissue 
            WHERE tissuename = %s
        )
        AND cl = (
            SELECT cl 
            FROM cell 
            WHERE cellname = %s
        )
        AND ensembl = (
            SELECT ensembl 
            FROM gene 
            WHERE genename = %s
        )
    """
    try:
        with get_db_connection() as conn:
            with conn.cursor() as cur:
                cur.execute(query, (tissue_name, cell_name, gene_name))
                result = cur.fetchone()
                return result[0] if result else None
    except Exception as e:
        print(f"Database error: {e}")
        return None
    
# Function to query Tavily using SDK
def query_tavily(tissue_name, cell_name, gene_name, tavily_client = initialize_tavily()):
    query = f"Does {cell_name} express {gene_name} in {tissue_name}?"
    try:
        response = tavily_client.search(query= query, search_depth= "basic")
        sources = [result["url"] for result in response.get("results", [])]
        return sources if sources else None
    except Exception as e:
        print(f"Tavily Error: {e}")
        return None
    

# Function to summarize using OpenAI
def summarize_with_openai(tissue_name, cell_name, gene_name, sources, client = initialize_openai()):
    prompt = f"""
    Given the following sources, determine the confidence level that {cell_name} expresses {gene_name} in {tissue_name}.
    Return one of these values in valid JSON format:
    {{"status": 1.0}} → Yes
    {{"status": 0.67}} → Mostly Yes
    {{"status": 0.33}} → Mostly No
    {{"status": 0.0}} → No
    never return with markdown formatting.
    Sources: {sources}
    """

    try:
        response = client.chat.completions.create(
            model="gpt-4o-mini",
            messages=[{"role": "system", "content": "You are an expert in cellular biology."},
                      {"role": "user", "content": prompt}]
        )
        output = response.choices[0].message.content
        return json.loads(output)["status"]
    except Exception as e:
        print(f"OpenAI Error: {e}")
        return None
    
def insert_into_database(tissue_name, cell_name, gene_name, source_urls, certainty):
    conn = psycopg2.connect(os.getenv("DATABASE_URL"))
    cur = conn.cursor()

    try:
        # Lookup IDs
        cur.execute("SELECT uberon FROM tissue WHERE tissuename = %s", (tissue_name,))
        uberon = cur.fetchone()
        if not uberon:
            raise ValueError(f"Tissue '{tissue_name}' not found.")
        uberon = uberon[0]

        cur.execute("SELECT cl FROM cell WHERE cellname = %s", (cell_name,))
        cl = cur.fetchone()
        if not cl:
            raise ValueError(f"Cell '{cell_name}' not found.")
        cl = cl[0]

        cur.execute("SELECT ensembl FROM gene WHERE genename = %s", (gene_name,))
        ensembl = cur.fetchone()
        if not ensembl:
            raise ValueError(f"Gene '{gene_name}' not found.")
        ensembl = ensembl[0]

        # Insert or update gene_expression
        cur.execute("""
            INSERT INTO gene_expression (uberon, cl, ensembl, relationcertainty)
            VALUES (%s, %s, %s, %s)
            ON CONFLICT (uberon, cl, ensembl)
            DO UPDATE SET relationcertainty = EXCLUDED.relationcertainty
            RETURNING expressionid
        """, (uberon, cl, ensembl, certainty))
        expressionid = cur.fetchone()[0]

        # Insert sources and validate links
        for url in source_urls:
            # Check or insert source
            cur.execute("SELECT sourceid FROM source WHERE url = %s", (url,))
            row = cur.fetchone()
            if row:
                sourceid = row[0]
            else:
                cur.execute("INSERT INTO source (url) VALUES (%s) RETURNING sourceid", (url,))
                sourceid = cur.fetchone()[0]

            # Insert validate link if not exists
            cur.execute("""
                INSERT INTO validate (expressionid, sourceid)
                VALUES (%s, %s)
                ON CONFLICT DO NOTHING
            """, (expressionid, sourceid))

        conn.commit()
        print(f"✅ Inserted/Updated expressionid {expressionid} with sources.")

    except Exception as e:
        conn.rollback()
        print(f"❌ Error: {e}")

    finally:
        cur.close()
        conn.close()


def validate_cell_gene_relation(tissue_name, cell_name, gene_name):
    print(f"\n[Validation Start] Processing: {tissue_name} - {cell_name} - {gene_name}")

    # Step 2: Check database first
    stored_status = check_database(tissue_name, cell_name, gene_name)
    if stored_status is not None:
        print(f"[Validation] Relation already exists in DB with status {stored_status}. Skipping process.")
        return float(stored_status)  # Ensure status is returned as a float

    # Step 3: Query Tavily
    try:
        sources = query_tavily(tissue_name, cell_name, gene_name)
    except Exception as e:
        print(f"Tavily Error: {e}")
        return None

    # Step 3a: If no sources found, set status to 0.0 and insert into database
    if not sources:
        status = 0.0
        print(f"[Validation] No references found in Tavily. Setting status to {status}.")
        insert_into_database(tissue_name, cell_name, gene_name, [], status)
        return float(status)  # Ensure status is a float

    # Step 4: Summarize with OpenAI
    status = summarize_with_openai(tissue_name, cell_name, gene_name, sources)
    status = float(status)  # Ensure status is a float
    print(f"[Validation] OpenAI summary status: {status}")

    # Step 5: Insert into Database
    insert_into_database(tissue_name, cell_name, gene_name, sources, status)
    print(f"[Validation] Inserted/Updated in DB with status: {status}")

    return status  # Return the validation score as a float

def main_data_validation(results, target_uberon_id):
    ex_df = results.copy()
    for cluster in ex_df:
        unique_cells = ex_df[cluster]['cell'].unique()[:10]
        ex_df[cluster] = ex_df[cluster][['gene', 'cell', 'cell score']][ex_df[cluster]['cell'].isin(unique_cells)]
    tissue_name = ontl.tissue(target_uberon_id)

    for cluster in ex_df:
        ex_df[cluster]["certainty score"] = ex_df[cluster].apply(lambda row: validate_cell_gene_relation(tissue_name, row["cell"], row["gene"]), axis=1)

    return ex_df

def export_validation_dataframe(results):
    ex_df = results.copy()
    for cluster in ex_df:
        unique_cells = ex_df[cluster]['cell'].unique()[:10]
        ex_df[cluster] = ex_df[cluster][['gene', 'cell', 'cell score', 'certainty score']][ex_df[cluster]['cell'].isin(unique_cells)]
    return ex_df