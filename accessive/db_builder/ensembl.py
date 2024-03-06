import sqlite3
import json
from ftplib import FTP
import os
import gzip
import io

from ..data_structure import *

def download_ensembl_data(include_list=None, already_loaded = [], data_save_dir = None):
    if include_list is not None:
        include_list = [x[1] for x in include_list] # Only the scientific name strings, which match Ensembl's directory names
    with FTP("ftp.ensembl.org") as ftp:
        ftp.login()
        ftp.cwd("pub/current_json")
        for dir_name in ftp.nlst():
            print(dir_name)
            if include_list is not None and dir_name not in include_list:
                print(f"Skipping {dir_name} because it is not in the include list.")
                continue
            elif dir_name in already_loaded:
                print(f"Skipping {dir_name} because it is already loaded.")
                continue
            ftp.cwd(dir_name)
            for item in ftp.nlst():
                if item.endswith('.json'):
                    if data_save_dir:
                        cache_file = os.path.join(data_save_dir, item)+'.gz'
                        if os.path.exists(cache_file):
                            print('Cached: %s' % cache_file)
                            yield gzip.open(cache_file, 'rb')  
                            break 
                        else:
                            print('Caching: %s' % cache_file)
                            data = gzip.open(cache_file, 'wb')
                    else:
                        data = io.BytesIO()

                    try:
                        print("Downloading " + item)
                        ftp.retrbinary('RETR ' + item, data.write)
                    except ConnectionResetError:
                        print("Connection reset error. Trying again.")
                        ftp.retrbinary('RETR ' + item, data.write)

                    if data_save_dir:
                        data.close()
                        data = gzip.open(os.path.join(data_save_dir, item)+'.gz', 'rb')  
                    else:
                        data.seek(0)

                    print("Downloaded " + item)
                    yield data 
                    break
            ftp.cwd('..')


def load_ensembl_jsonfile(json_file, sqlite_file):
    if isinstance(json_file, str):
        try:
            data = json.load(open(json_file, 'r')) # NB this is typically very large!
        except json.decoder.JSONDecodeError:
            import pickle
            data = pickle.load(open(json_file, 'rb'))
    else:
        data = json.load(json_file)

    conn = sqlite3.connect(sqlite_file)
    c = conn.cursor()


    taxon = data['organism']['taxonomy_id']
    c.execute("INSERT INTO species_table (taxon, name, common_name) VALUES (?, ?, ?)", (taxon, data['organism']['name'], data['organism']['display_name']))
    
    c.execute("SELECT MAX(MAX(gene_index), MAX(mrna_index), MAX(prot_index)) FROM entity_table;")
    next_index = c.fetchone()[0]
    if next_index is None:
        next_index = 0

    skipped_lrg = 0 
    for gene in data['genes']:
        gene_index = next_index
        next_index += 1

        if gene['id'][:3] == 'LRG':
            continue

        for db_name, json_name in GENE_COLS:
            for item in list_item(gene, json_name):
                c.execute(f"INSERT INTO {db_name} (entity_index, identifier, taxon, is_canonical) VALUES (?, ?, ?, ?)", (gene_index, item, taxon, 1))
                c.execute(f"INSERT INTO identifier_directory (identifier, identifier_type) VALUES (?, ?)", (item, db_name))

        if not gene.get('transcripts'):
            c.execute(f"INSERT INTO entity_table (taxon, gene_index) VALUES (?, ?)", (taxon, gene_index))
        for isoform in gene.get('transcripts', []):
            isoform_index = next_index
            next_index += 1

            for db_name, json_name in ISOFORM_COLS:
                for item in list_item(isoform, json_name):
                    c.execute(f"INSERT INTO {db_name} (entity_index, identifier, taxon, is_canonical) VALUES (?, ?, ?, ?)", (isoform_index, item, taxon, 1))
                    c.execute(f"INSERT INTO identifier_directory (identifier, identifier_type) VALUES (?, ?)", (item, db_name))
           
            if not isoform.get('translations'):
                c.execute(f"INSERT INTO entity_table (taxon, gene_index, mrna_index) VALUES (?, ?, ?)", (taxon, gene_index, isoform_index))
            for proteoform in isoform.get('translations', []):
                proteoform_index = next_index
                next_index += 1

                for db_name, json_name in PROTEOFORM_COLS:
                    for item in list_item(proteoform, json_name):
                        c.execute(f"INSERT INTO {db_name} (entity_index, identifier, taxon, is_canonical) VALUES (?, ?, ?, ?)", (proteoform_index, item, taxon, 1))
                        c.execute(f"INSERT INTO identifier_directory (identifier, identifier_type) VALUES (?, ?)", (item, db_name))
                
                c.execute(f"INSERT INTO entity_table (taxon, gene_index, mrna_index, prot_index) VALUES (?, ?, ?, ?)", (taxon, gene_index, isoform_index, proteoform_index))
        
        if gene_index % 1000 == 0:
            conn.commit()
            print(f"Processed {gene_index} genes.")


    conn.commit()
    conn.close()
    print(f"Skipped {skipped_lrg} LRG genes.")
    print(f"Finished loading {json_file}.")

