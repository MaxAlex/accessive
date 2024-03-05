import sqlite3
import json
from ftplib import FTP
import os
import gzip
import tempfile
import io
import datetime

from ..data_structure import *
from ..database_ops import DATABASE_VERSION, DATABASE_FILE
from .ensembl import download_ensembl_data, load_ensembl_jsonfile
from .nextprot import download_nextprot_map_files, load_nextprot_accessions
# from .uniprot import download_uniprot_data, load_uniprot_table








def create_sqlite_database(sqlite_file):
    conn = sqlite3.connect(sqlite_file)
    c = conn.cursor()

    c.execute(f"CREATE TABLE accessive_meta (key TEXT, val TEXT)")
    c.execute(f"INSERT INTO accessive_meta (key, val) VALUES (?, ?)", ("database_version", DATABASE_VERSION))
    c.execute(f"INSERT INTO accessive_meta (key, val) VALUES (?, ?)", ("creation_time", datetime.datetime.now().isoformat()))

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


def deprecate_trembl_accessions(sqlite_file):
    # This is special handling for Uniprot accessions, which are split between "SwissProt" and "TrEMBL" based on
    # whether they're reviewed. When a SwissProt accession exists for a given gene, it's generally better to
    # use that. So, here we're marking all TrEMBL accessions for genes with a SwissProt accession as non-canonical.

    conn = sqlite3.connect(sqlite_file)
    c = conn.cursor()
    c.execute("""
              UPDATE uniprot_gene
              SET is_canonical = 0
              WHERE identifier IN (
                  SELECT uniprot_trembl.identifier
                  FROM uniprot_trembl 
                  INNER JOIN uniprot_swissprot
                  WHERE uniprot_trembl.entity_index = uniprot_swissprot.entity_index
              )
              """)
    conn.commit()
    conn.close()

def remove_redundancies(sqlite_file):
    conn = sqlite3.connect(sqlite_file)
    c = conn.cursor()

    c.execute("DROP TABLE IF EXISTS dedup_directory")
    c.execute("CREATE TABLE dedup_directory AS SELECT * FROM identifier_directory GROUP BY identifier, identifier_type")
    c.execute("DROP TABLE identifier_directory")
    c.execute("ALTER TABLE dedup_directory RENAME TO identifier_directory")

    conn.commit()
    conn.close()


def build_indexes(sqlite_file):
    conn = sqlite3.connect(sqlite_file)
    c = conn.cursor()
    
    c.execute("CREATE INDEX IF NOT EXISTS gene_entity_index ON entity_table (gene_index)")
    c.execute("CREATE INDEX IF NOT EXISTS mrna_entity_index ON entity_table (mrna_index)")
    c.execute("CREATE INDEX IF NOT EXISTS prot_entity_index ON entity_table (prot_index)")

    for acc_table, _ in GENE_COLS + ISOFORM_COLS + PROTEOFORM_COLS:
        c.execute(f"CREATE INDEX IF NOT EXISTS {acc_table}_entity_index ON {acc_table} (entity_index)")
        print(acc_table)

    conn.commit()
    conn.close()
    print("Built indexes")


def compile_full_database(sqlite_file = None, include_list=None, cache_dir=None):
    if sqlite_file is None:
        sqlite_file = DATABASE_FILE
        if not os.path.exists(os.path.dirname(sqlite_file)):
            os.makedirs(os.path.dirname(sqlite_file))
    if cache_dir is None:
        cache_dir = tempfile.mkdtemp()

    assert(os.path.exists(cache_dir))
    assert(not os.path.exists(sqlite_file)), "Database file already exists: %s" % sqlite_file
  
    print("Initializing database...")
    create_sqlite_database(sqlite_file)   

    print("Downloading and loading Ensembl data...")
    for data_buffer in download_ensembl_data(include_list, [], cache_dir):
        load_ensembl_jsonfile(data_buffer, sqlite_file)
        del data_buffer

    print("Downloading and loading Nextprot data...")
    nextprot_ensts, nextprot_ensgs = download_nextprot_map_files(cache_dir)
    load_nextprot_accessions(nextprot_ensts, nextprot_ensgs, sqlite_file)
   
    print("Adjusting canonical accession designations...")
    deprecate_trembl_accessions(DATABASE_FILE)
    print("Removing redundant rows...")
    remove_redundancies(DATABASE_FILE)
    print("Building indexes...")
    build_indexes(DATABASE_FILE)

    print("Vacuuming database...")
    conn = sqlite3.connect(DATABASE_FILE)
    conn.execute("VACUUM")
    conn.commit()
    conn.close()

    print("Done")


if __name__ == '__main__':
    # This is chosen somewhat arbitrarily. Will be updated according to demand vs reasonable disk space usage.
    # Alternatelty, an enterprising user can change this list and build their own DB (although many
    # of the data sources don't cover all species!)
    species_manifest = [(9913, 'bos_taurus', 'Cow'),
                        (6239, 'caenorhabditis_elegans', 'Caenorhabditis elegans (Nematode, N2)'),
                        (9615, 'canis_lupus_familiaris', 'Dog'),
                        (7955, 'danio_rerio', 'Zebrafish'),
                        (7227, 'drosophila_melanogaster', 'Drosophila melanogaster (Fruit fly)'),
                        (9685, 'felis_catus', 'Cat'),
                        (9606, 'homo_sapiens', 'Human'),
                        (10090, 'mus_musculus', 'Mouse'),
                        (10116, 'rattus_norvegicus', 'Rat'),
                        (559292, 'saccharomyces_cerevisiae', 'Saccharomyces cerevisiae')]
    
    compile_full_database(sqlite_file = None, 
                          include_list = species_manifest, 
                          cache_dir = '/data/biostuff/ensembl_data/cache')


    
