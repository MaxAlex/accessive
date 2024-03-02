import os
import sqlite3
import pandas as pd

from .data_structure import *
from .database_ops import DATABASE_FILE

class Accessive():
    def __init__(self, sqlite_file = None):
        if sqlite_file is None:
            if not os.path.exists(DATABASE_FILE):
                raise RuntimeError(f"Database file not found in default location: {DATABASE_FILE} . Download the database or specify a different file.")
            sqlite_file = DATABASE_FILE
        self.conn = sqlite3.connect(sqlite_file)
        self.c = self.conn.cursor()

        database_ver = self._get_db_version()
        if database_ver != DATABASE_VERSION:
            print(f"WARNING: Database version {database_ver} does not match expected version {DATABASE_VERSION}. It may be incompatible with this version of Accessive.")
            print("You can download the correct database version by running the command 'python -m accessive.database_ops --download'")


    def _get_db_version(self):
        self.c.execute("SELECT version_number FROM database_version")
        return self.c.fetchone()[0]


    def _get_identifier_type(self, acc, taxon = None):
        if taxon:
            self.c.execute("SELECT identifier_type FROM identifier_directory WHERE identifier = ? AND taxon = ?", (acc, taxon))
        else:
            self.c.execute("SELECT identifier_type FROM identifier_directory WHERE identifier = ?", (acc,))
        return self.c.fetchone()[0]


    def _get_type_metadata(self, idtypes):
        self.c.execute("SELECT * FROM metadata_table WHERE identifier_type IN (%s)" % ','.join(['?']*len(idtypes)), idtypes)
        return dict(self.c.fetchall())

    
    def identify_taxon(self, taxon):
        self.c.execute("SELECT name, common_name FROM species_table WHERE taxon = ?", (taxon,))
        return self.c.fetchone()


    def available_taxons(self):
        self.c.execute("SELECT taxon, common_name FROM species_table")
        return [x for x in self.c.fetchall()]
    def available_taxa(self):
        return self.available_taxons()


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


    def map(self, ids, from_type = None, to_types = None, taxon = None, return_query_info = False, return_format=None):
        ids = [ids] if isinstance(ids, str) else ids 
        assert(to_types is None or isinstance(to_types, list))
        if not from_type:
            from_type = self._get_identifier_type(ids[0], taxon)

        if to_types is None:
            _, fromtype_entity = self._get_type_metadata([from_type])
            if fromtype_entity == 'gene':
                to_types = [x for x, _ in ENSEMBL_GENE_COLS]
            elif fromtype_entity == 'mrna':
                to_types = [x for x, _ in ENSEMBL_ISOFORM_COLS]
            elif fromtype_entity == 'prot':
                to_types = [x for x, _ in ENSEMBL_PROTEOFORM_COLS]
            else:
                raise Exception("Someone messed up the database and added an entity type: " + fromtype_entity)

        if len(to_types) < len(set(to_types)):
            to_types = list(set(to_types))

        try:
            from_type = TO_DATABASE_NAME[from_type]
        except KeyError:
            raise Exception(f"Source identifier type {from_type} is not recognized.")

        try:
            to_types = [TO_DATABASE_NAME[x] for x in to_types]
        except KeyError:
            raise Exception(f"Destination identifier type {[x for x in to_types if x not in TO_DATABASE_NAME]} is not recognized.")

        result = self._query(ids, from_type, to_types, taxon)

        dedup_ind = result.applymap(lambda x: x if not isinstance(x, list) else ','.join(x)).drop_duplicates().index
        result = result.loc[dedup_ind]
        result.rename(columns=TO_PRETTIER_NAME, inplace=True)
    
        if return_format == 'txt':
            result = result.applymap(lambda x: x if isinstance(x, str) else ','.join(x))
            result = result.to_csv(index=False, sep='\t') 
        elif return_format == 'json':
            result = result.to_json(orient='records')
        elif return_format == 'pandas':
            pass
        elif return_format is not None:
            raise Exception(f"Return format {return_format} is not recognized.")

        if return_query_info:
            return {'result': result, 'from_type': from_type, 'to_types': to_types, 'taxon': taxon}
        else:
            return result


















