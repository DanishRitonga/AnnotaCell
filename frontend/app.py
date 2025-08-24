import streamlit as st
import pandas as pd
from components.file_input import csv_read
from components.file_preprocessing import extract_clusters

# --- Page Configuration ---
st.set_page_config(
    page_title="AnnotaCell",
    page_icon="🧬",
    layout="wide"
)

st.title("🧬 AnnotaCell")

st.sidebar.header("Upload your CSV file")
uploaded_file = st.sidebar.file_uploader("Choose a file", type="csv")
df = None

# --- Main page logic and display ---
if uploaded_file is not None:
    try:
        df = csv_read(uploaded_file)
        cluster_dict = extract_clusters(df)
        display_df = pd.DataFrame(list(cluster_dict.items()), columns=['Cluster', 'Genes'])
        
        st.toast("File uploaded and processed successfully!", icon="✅")

        st.dataframe(display_df, 
                     column_config={
                        "Genes": st.column_config.Column(width=1500)
                    },
                    hide_index=True,
        )

        

    except Exception as e:
        st.error(f"Error reading the file. Please check the separator and file format. Details: {e}")
else:
    st.info("Please upload a CSV file using the controls on the left.")