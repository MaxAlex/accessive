import os
import sqlite3
import requests
import gzip
import io
import pandas as pd
import tempfile
from ..data_structure import *

def download_nextprot_map_files(cache_dir):
    enst_map_file = os.path.join(cache_dir, 'nextprot_enst.txt')
    if not os.path.exists(enst_map_file): 
        res = requests.get('https://download.nextprot.org/pub/current_release/mapping/nextprot_enst.txt')
        res.raise_for_status()
        with open(enst_map_file, 'w') as f:
            f.write(res.text)

    ensg_map_file = os.path.join(cache_dir, 'nextprot_ensg.txt')
    if not os.path.exists(ensg_map_file):
        res = requests.get('https://download.nextprot.org/pub/current_release/mapping/nextprot_ensg.txt')
        res.raise_for_status()
        with open(ensg_map_file, 'w') as f:
            f.write(res.text)
    return enst_map_file, ensg_map_file


def joined_table(conn, c, join_table, target_data, target_main_col, target_join_col, target_table_name):
    assert(all([x in target_data.columns for x in [target_main_col, target_join_col, 'taxon', 'is_canonical']]))
    target_data.to_sql('temp_table', conn, if_exists='replace', index=False)
    c.execute(f"DROP TABLE IF EXISTS {target_table_name}")
    c.execute(f"CREATE TABLE {target_table_name} ({', '.join(IDENTIFIER_TABLE_COLS)})")
    cmd = f"""INSERT INTO {target_table_name} (entity_index, identifier, taxon, is_canonical) 
              SELECT {join_table}.entity_index, temp_table.{target_main_col}, temp_table.taxon, temp_table.is_canonical 
              FROM temp_table JOIN {join_table} ON temp_table.{target_join_col} = {join_table}.identifier
           """
    c.execute(cmd)
    conn.commit()
    c.execute("DROP TABLE temp_table")
    conn.commit()



def load_nextprot_accessions(enst_map_file, ensg_map_file, sqlite_file):
    ensgs = pd.read_csv(ensg_map_file, sep='\t', header=None, names=['nextprot', 'ensg'])
    ensts = pd.read_csv(enst_map_file, sep='\t', header=None, names=['nextprot', 'enst'])

    ensgs['taxon'] = 9606 # Nextprot is only for human stuff!
    ensts['taxon'] = 9606
    ensgs['is_canonical'] = 1 # Nextprot is non-redundant
    ensts['is_canonical'] = 1   

    conn = sqlite3.connect(sqlite_file)
    c = conn.cursor()
    joined_table(conn, c, 'ensembl_gene', ensgs, 'nextprot', 'ensg', 'nextprot')
    joined_table(conn, c, 'ensembl_mrna', ensts, 'nextprot', 'enst', 'nextprot_isoform')
    c.execute("INSERT INTO metadata_table (identifier_type, entity_type) VALUES (?, ?)", ('nextprot', 'gene'))
    c.execute("INSERT INTO metadata_table (identifier_type, entity_type) VALUES (?, ?)", ('nextprot_isoform', 'mrna'))
    conn.commit()
    conn.close()


if __name__ == '__main__':
    load_nextprot_accessions('/data/biostuff/ensembl_data/accessive_standard.sqlite', '/data/biostuff/ensembl_data/cache')
