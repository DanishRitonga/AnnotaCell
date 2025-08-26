import streamlit as st
import pandas as pd
import logging  # 1. Import the logging module
from components.file_input import csv_read
from components.file_preprocessing import extract_clusters

# 2. Configure logging at the start of your script
logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s - %(name)s - %(levelname)s - %(message)s',
    handlers=[
        logging.FileHandler("app.log"),  # Log to a file named app.log
        logging.StreamHandler()          # Log to the console
    ]
)

# Create a logger instance for this module
logger = logging.getLogger(__name__)

# --- Page Configuration ---
st.set_page_config(
    page_title="AnnotaCell",
    page_icon="🧬",
    layout="wide"
)

st.title("🧬 AnnotaCell")
logger.info("Application started and page loaded.") # Log app start

st.sidebar.header("Upload your CSV file")
uploaded_file = st.sidebar.file_uploader("Choose a file", type="csv")
df = None

# --- Main page logic and display ---
if uploaded_file is not None:
    logger.info(f"File uploaded: {uploaded_file.name}") # Log when a file is uploaded
    try:
        df = csv_read(uploaded_file)
        logger.info("CSV file read successfully into DataFrame.")
        
        cluster_dict = extract_clusters(df)
        logger.info(f"Successfully extracted {len(cluster_dict)} clusters.")
        
        display_df = pd.DataFrame(list(cluster_dict.items()), columns=['Cluster', 'Genes'])
        
        st.toast("File uploaded and processed successfully!", icon="✅")

        st.dataframe(display_df, 
                     column_config={
                         "Genes": st.column_config.Column(width=1500)
                     },
                     hide_index=True,
        )

    except Exception as e:
        # 3. Log the full error with traceback for debugging
        logger.error(f"An error occurred while processing {uploaded_file.name}.", exc_info=True)
        st.error(f"Error reading the file. Please check the separator and file format. Details: {e}")
else:
    st.info("Please upload a CSV file using the controls on the left.")