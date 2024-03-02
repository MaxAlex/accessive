from .compile_database import download_ensembl_data 
from .data_structure import *
import sqlite3
import json


entity_table_columns = ['taxon INTEGER', 'gene_index INTEGER', 'mrna_index INTEGER', 'prot_index INTEGER']
identifier_table_colummns = ['entity_index INTEGER', 'identifier TEXT', 'taxon INTEGER', 'is_canonical INTEGER'] # Whether the index is _gene or etc is determined in metadata table
metadata_table_columns = ['identifier_type TEXT', 'entity_type TEXT']
identifier_directory_columns = ['identifier TEXT', 'identifier_type TEXT']


def create_sqlite_database(sqlite_file):
    conn = sqlite3.connect(sqlite_file)
    c = conn.cursor()

    c.execute(f"CREATE TABLE IF NOT EXISTS entity_table ({', '.join(entity_table_columns)})")
    c.execute(f"CREATE TABLE IF NOT EXISTS metadata_table ({', '.join(metadata_table_columns)})")
    c.execute(f"CREATE TABLE IF NOT EXISTS identifier_directory ({', '.join(identifier_directory_columns)})")

    conn.commit()
    conn.close()



# ENSEMBL_GENE_COLS = ['gene_description', 'gene_name', 'arrayexpress', 'biogrid', 'ens_lrg_gene', 'entrez_gene', 'genecards', 'hgnc', 'mim_gene', 'pfam', 'uniprot_gene', 'wikigene']
ENSEMBL_GENE_COLS = [('ensembl_gene', 'id'), ('gene_description', 'description'), ('gene_name', 'name'), ('arrayexpress', 'ArrayExpress'), ('biogrid', 'BioGRID'), 
                     ('ens_lrg_gene', 'ENS_LRG_gene'), ('entrez_gene', 'EntrezGene'), ('genecards', 'GeneCards'), ('hgnc', 'HGNC'), 
                     ('mim_gene', 'MIM_GENE'), ('pfam', 'Pfam'), ('uniprot_gene', 'Uniprot_gn'), ('wikigene', 'WikiGene')]
# ENSEMBL_ISOFORM_COLS = ['ccds', 'ens_lrg_transcript', 'refseq_mrna', 'refseq_ncrna', 'ucsc', 'isoform_biotype'] 
ENSEMBL_ISOFORM_COLS = [('ensembl_mrna', 'id'), ('ccds', 'CCDS'), ('ens_lrg_transcript', 'ENS_LRG_transcript'), ('refseq_mrna', 'RefSeq_mRNA'), ('refseq_ncrna', 'RefSeq_ncRNA'),
                        ('ucsc', 'UCSC'), ('isoform_biotype', 'biotype')]
# ENSEMBL_PROTEOFORM_COLS = ['uniparc', 'alphafold', 'uniprot_swissprot', 'uniprot_trembl', 'uniprot_isoform', 'refseq_peptide', 'embl', 'pdb']
ENSEMBL_PROTEOFORM_COLS = [('ensembl_prot', 'id'), ('uniparc', 'UniParc'), ('alphafold', 'alphafold'), ('uniprot_swissprot', 'Uniprot/SWISSPROT'), ('uniprot_trembl', 'Uniprot/SPTREMBL'),
                           ('uniprot_isoform', 'Uniprot_isoform'), ('refseq_peptide', 'RefSeq_peptide'), ('embl', 'EMBL'), ('pdb', 'PDB')]

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

    for gene_col, _ in ENSEMBL_GENE_COLS:
        c.execute(f"CREATE TABLE IF NOT EXISTS {gene_col} ({', '.join(identifier_table_colummns)})")
        c.execute(f"INSERT INTO metadata_table (identifier_type, entity_type) VALUES (?, ?)", (gene_col, 'gene'))
    for isoform_col, _ in ENSEMBL_ISOFORM_COLS:
        c.execute(f"CREATE TABLE IF NOT EXISTS {isoform_col} ({', '.join(identifier_table_colummns)})")
        c.execute(f"INSERT INTO metadata_table (identifier_type, entity_type) VALUES (?, ?)", (isoform_col, 'mrna'))
    for proteoform_col, _ in ENSEMBL_PROTEOFORM_COLS:
        c.execute(f"CREATE TABLE IF NOT EXISTS {proteoform_col} ({', '.join(identifier_table_colummns)})")
        c.execute(f"INSERT INTO metadata_table (identifier_type, entity_type) VALUES (?, ?)", (proteoform_col, 'prot'))

    taxon = data['organism']['taxonomy_id']

    c.execute("SELECT MAX(MAX(gene_index), MAX(mrna_index), MAX(prot_index)) FROM entity_table;")
    next_index = c.fetchone()[0]
    if next_index is None:
        next_index = 0

    skipped_lrg = 0 
    for gene in data['genes']:
        gene_index = next_index
        next_index += 1

        for db_name, json_name in ENSEMBL_GENE_COLS:
            for item in list_item(gene, json_name):
                c.execute(f"INSERT INTO {db_name} (entity_index, identifier, taxon, is_canonical) VALUES (?, ?, ?, ?)", (gene_index, item, taxon, 1))
                c.execute(f"INSERT INTO identifier_directory (identifier, identifier_type) VALUES (?, ?)", (item, db_name))

        if not gene.get('transcripts'):
            c.execute(f"INSERT INTO entity_table (taxon, gene_index) VALUES (?, ?)", (taxon, gene_index))
        for isoform in gene.get('transcripts', []):
            isoform_index = next_index
            next_index += 1

            for db_name, json_name in ENSEMBL_ISOFORM_COLS:
                for item in list_item(isoform, json_name):
                    c.execute(f"INSERT INTO {db_name} (entity_index, identifier, taxon, is_canonical) VALUES (?, ?, ?, ?)", (isoform_index, item, taxon, 1))
                    c.execute(f"INSERT INTO identifier_directory (identifier, identifier_type) VALUES (?, ?)", (item, db_name))
           
            if not isoform.get('translations'):
                c.execute(f"INSERT INTO entity_table (taxon, gene_index, mrna_index) VALUES (?, ?, ?)", (taxon, gene_index, isoform_index))
            for proteoform in isoform.get('translations', []):
                proteoform_index = next_index
                next_index += 1

                for db_name, json_name in ENSEMBL_PROTEOFORM_COLS:
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


import pandas as pd
class Interface():
    def __init__(self, sqlite_file):
        self.conn = sqlite3.connect(sqlite_file)
        self.c = self.conn.cursor()

    def _get_identifier_type(self, acc, taxon = None):
        if taxon:
            self.c.execute("SELECT identifier_type FROM identifier_directory WHERE identifier = ? AND taxon = ?", (acc, taxon))
        else:
            self.c.execute("SELECT identifier_type FROM identifier_directory WHERE identifier = ?", (acc,))
        return self.c.fetchone()[0]

    def _get_type_metadata(self, idtypes):
        self.c.execute("SELECT * FROM metadata_table WHERE identifier_type IN (%s)" % ','.join(['?']*len(idtypes)), idtypes)
        return dict(self.c.fetchall())


    def _query(self, accs, from_type, dest_types, taxon = None):
        if from_type not in dest_types:
            dest_types = [from_type] + dest_types

        type_meta = self._get_type_metadata(dest_types)
        assert(len(type_meta) == len(dest_types))

        # gene_types = [x for x in type_meta if type_meta[x] == 'gene']
        # isoform_types = [x for x in type_meta if type_meta[x] == 'isoform']
        # proteoform_types = [x for x in type_meta if type_meta[x] == 'proteoform']
        
        if taxon is None:
            self.c.execute(f"SELECT taxon, entity_index FROM {from_type} WHERE identifier IN (%s)" % ','.join(['?']*len(accs)), accs) 
            taxon, entity_indices = zip(*self.c.fetchall())
            assert(len(set(taxon)) == 1), "Multi-species lookup not currently supported (found taxons %s)" % ', '.join(set(taxon))
            taxon = taxon[0]
        else:
            self.c.execute(f"SELECT entity_index FROM {from_type} WHERE taxon = ? AND identifier IN (%s)" % ','.join(['?']*len(accs)), [taxon]+accs)
            entity_indices = [x[0] for x in self.c.fetchall()]

        self.c.execute(f"SELECT gene_index, mrna_index, prot_index FROM entity_table WHERE taxon = ? AND {type_meta[from_type]}_index IN ({','.join(['?']*len(entity_indices))})", [taxon]+entity_indices)
        entity_table_indices = self.c.fetchall()

        # Construct the base of the final query
        base_query = f"SELECT et.taxon, et.gene_index, et.mrna_index, et.prot_index"

        # Initialize a list to hold JOIN clauses
        join_clauses = []
        # Initialize a list to hold the columns to select from destination identifier tables
        select_columns = []
        column_names = ['taxon', 'gene_index', 'mrna_index', 'prot_index']
        # For each destination type, add the necessary JOIN clause and select column
        for dest_type in dest_types:
            # Determine the entity column name based on the entity type
            entity_col = f"{type_meta[dest_type]}_index"
            # Construct the JOIN clause for the current destination type table
            join_clause = f"LEFT JOIN {dest_type} ON et.{entity_col} = {dest_type}.entity_index AND et.taxon = {dest_type}.taxon"
            join_clauses.append(join_clause)
            # Add the identifier column from the destination type table to the select list
            select_columns.append(f"{dest_type}.identifier AS {dest_type}_identifier")
            column_names.append(dest_type)

        # Add the JOIN clauses and select columns to the base query
        final_query = base_query + ", " + ", ".join(select_columns) + " FROM entity_table et " + " ".join(join_clauses)

        # Add a WHERE clause to filter by taxon and entity indices, if necessary
        final_query += f" WHERE et.taxon = ? AND et.{type_meta[from_type]}_index IN ({','.join(['?']*len(entity_indices))})"

        # Execute the final query
        self.c.execute(final_query, [taxon]+entity_indices)
        results = self.c.fetchall()
        result_table = pd.DataFrame(results, columns=column_names)
        # Process results as needed...
        return result_table[column_names[4:]]

# Note: Ensure that the implementation of _get_type_metadata and other support functions correctly handle
# the mapping of identifier types to their corresponding columns and tables as expected.

        # Now this has to do a set of joins so that we get each of the destination types from the corresponding type table




















