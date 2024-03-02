import sqlite3
import json
from ftplib import FTP
import os
import gzip
import tempfile
import io

from ..data_structure import *
from ..database_ops import DATABASE_VERSION, DATABASE_FILE
from .ensembl import download_ensembl_data, load_ensembl_jsonfile
from .uniprot import download_uniprot_data, load_uniprot_table








def create_sqlite_database(sqlite_file):
    conn = sqlite3.connect(sqlite_file)
    c = conn.cursor()

    c.execute(f"CREATE TABLE database_version (version_number TEXT)")
    c.execute(f"INSERT INTO database_version (version_number) VALUES (?)", (DATABASE_VERSION,))

    c.execute(f"CREATE TABLE entity_table ({', '.join(ENTITY_TABLE_COLS)})")
    c.execute(f"CREATE TABLE metadata_table ({', '.join(METADATA_COLS)})")
    c.execute(f"CREATE TABLE identifier_directory ({', '.join(DIRECTORY_COLS)})")
    c.execute(f"CREATE TABLE species_table ({', '.join(SPECIES_COLS)})")

    for gene_col, _ in GENE_COLS:
        c.execute(f"CREATE TABLE {gene_col} ({', '.join(IDENTIFIER_TABLE_COLS)})")
        c.execute(f"INSERT INTO metadata_table (identifier_type, entity_type) VALUES (?, ?)", (gene_col, 'gene'))
    for isoform_col, _ in ISOFORM_COLS:
        c.execute(f"CREATE TABLE {isoform_col} ({', '.join(IDENTIFIER_TABLE_COLS)})")
        c.execute(f"INSERT INTO metadata_table (identifier_type, entity_type) VALUES (?, ?)", (isoform_col, 'mrna'))
    for proteoform_col, _ in PROTEOFORM_COLS:
        c.execute(f"CREATE TABLE {proteoform_col} ({', '.join(IDENTIFIER_TABLE_COLS)})")
        c.execute(f"INSERT INTO metadata_table (identifier_type, entity_type) VALUES (?, ?)", (proteoform_col, 'prot'))

    conn.commit()
    conn.close()






def compile_full_database(sqlite_file = None, include_list=None, cache_dir=None):
    if sqlite_file is None:
        sqlite_file = DATABASE_FILE
        if not os.path.exists(os.path.dirname(sqlite_file)):
            os.makedirs(os.path.dirname(sqlite_file))
    if cache_dir is None:
        cache_dir = tempfile.mkdtemp()

    assert(os.path.exists(cache_dir))
    assert(not os.path.exists(sqlite_file))
   
    create_sqlite_database(sqlite_file)   

    for data_buffer in download_ensembl_data(include_list, [], cache_dir):
        load_ensembl_jsonfile(data_buffer, sqlite_file)
        del data_buffer





