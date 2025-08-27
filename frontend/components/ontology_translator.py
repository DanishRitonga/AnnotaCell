import psycopg2
import os
from dotenv import load_dotenv

def load_ontology() -> dict:
    script_dir = os.path.dirname(__file__)
    dotenv_path = os.path.join(script_dir, '..', '.env')
    load_dotenv(dotenv_path=dotenv_path)
    conn = psycopg2.connect(os.getenv("DATABASE_URL"))
    cur = conn.cursor()

    ontology_data = {
        "uberon-tissue": {},
        "cl-cell": {},
        "ensembl-gene": {}
    }

    table_mappings = [
        ("TISSUE", "UBERON", "tissueName", "uberon-tissue"),
        ("CELL", "CL", "cellName", "cl-cell"),
        ("GENE", "ENSEMBL", "geneName", "ensembl-gene")
    ]

    for table, key_col, value_col, json_key in table_mappings:
        cur.execute(f"SELECT {key_col}, {value_col} FROM {table}")
        rows = cur.fetchall()
        ontology_data[json_key] = {key: value for key, value in rows}

    cur.close()
    conn.close()

    return ontology_data

class ontl:
    _ontology_dict = None  # Class-level cache

    @classmethod
    def _load_ontology(cls):
        if cls._ontology_dict is None:
            cls._ontology_dict = load_ontology()

    @classmethod
    def tissue(cls, value: str, rev: bool = False) -> str | None:
        cls._load_ontology()
        table = cls._ontology_dict.get("uberon-tissue", {})
        return (
            next((k for k, v in table.items() if v.lower() == value.lower()), None)
            if rev else table.get(value)
        )

    @classmethod
    def cell(cls, value: str, rev: bool = False) -> str | None:
        cls._load_ontology()
        table = cls._ontology_dict.get("cl-cell", {})
        return (
            next((k for k, v in table.items() if v.lower() == value.lower()), None)
            if rev else table.get(value)
        ) 

    @classmethod
    def gene(cls, value: str, rev: bool = False) -> str | None:
        cls._load_ontology()
        table = cls._ontology_dict.get("ensembl-gene", {})
        return (
            next((k for k, v in table.items() if v.lower() == value.lower()), None)
            if rev else table.get(value)
        )