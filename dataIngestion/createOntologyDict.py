import requests
import json
import logging

# --- Configuration ---
# Configure logging for clear status updates
logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(levelname)s - %(message)s')

# The specific API endpoint to get all available tissues and cell types
CELLXGENE_FILTERS_URL = "https://api.cellxgene.cziscience.com/wmg/v2/filters"
OUTPUT_FILENAME = "ontology_index.json"

# --- Main Script Logic ---

def fetch_and_build_index():
    """
    Fetches the ontology index from the cellxgene API and saves it
    to a local JSON file.
    """
    try:
        logging.info("Attempting to fetch ontology data from cellxgene API...")
        response = requests.get(CELLXGENE_FILTERS_URL)
        
        # This will raise an HTTPError if the HTTP request returned an unsuccessful status code
        response.raise_for_status()
        
        api_data = response.json()
        logging.info("Successfully fetched data from API.")

        # Extract the relevant lists from the API response
        # We also filter to only include 'id' and 'name' to keep the file clean
        tissues = [
            {"id": t["id"], "name": t["name"]} 
            for t in api_data.get("snapshot", {}).get("tissues", [])
        ]
        
        cell_types = [
            {"id": c["id"], "name": c["name"]} 
            for c in api_data.get("snapshot", {}).get("cell_types", [])
        ]

        if not tissues or not cell_types:
            logging.warning("API response did not contain tissues or cell types. The output file will be empty.")

        # Structure the data into the final dictionary format
        ontology_index = {
            "tissues": tissues,
            "cell_types": cell_types
        }

        # Write the dictionary to the JSON file
        logging.info(f"Writing {len(tissues)} tissues and {len(cell_types)} cell types to '{OUTPUT_FILENAME}'...")
        with open(OUTPUT_FILENAME, 'w') as f:
            # indent=2 makes the JSON file human-readable
            json.dump(ontology_index, f, indent=2)
            
        logging.info(f"✅ Successfully created '{OUTPUT_FILENAME}'.")

    except requests.exceptions.RequestException as e:
        logging.error(f"❌ Failed to connect to the cellxgene API: {e}")
    except KeyError as e:
        logging.error(f"❌ Unexpected API response format. Missing key: {e}")
    except Exception as e:
        logging.error(f"❌ An unexpected error occurred: {e}")


# --- Run the script ---
if __name__ == "__main__":
    fetch_and_build_index()