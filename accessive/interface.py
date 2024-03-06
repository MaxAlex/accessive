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

        try:
            database_ver = self._get_db_version()
            if database_ver != DATABASE_VERSION:
                print(f"WARNING: Database version {database_ver} does not match expected version {DATABASE_VERSION}. It may be incompatible with this version of Accessive.")
                print("You can download the correct database version by running the command 'python -m accessive.database_ops --download'")
        except sqlite3.OperationalError:
            print("WARNING: Database version not found. This may be an old version of the database that does not include version information.")
            print("You can download the correct database version by running the command 'python -m accessive.database_ops --download'")


    def _get_db_version(self):
        self.c.execute("SELECT val FROM accessive_meta WHERE key = 'database_version'")
        return self.c.fetchone()[0]


    def _get_identifier_type(self, acc):
        self.c.execute("SELECT identifier_type FROM identifier_directory WHERE identifier = ?", (acc,))
        types = [x[0] for x in self.c.fetchall()]
        if len(set(types)) > 1:
            raise Exception(f"Identifier {acc} is associated with multiple types: {', '.join(types)}")
        else:
            return types[0]


    def _get_type_metadata(self, idtypes):
        self.c.execute("SELECT * FROM metadata_table WHERE identifier_type IN (%s)" % ','.join(['?']*len(idtypes)), idtypes)
        return dict(self.c.fetchall())

    
    def identify_taxon(self, taxon):
        """
        Utility function to identify what species a taxon number in the database corresponds to.
        """
        self.c.execute("SELECT name, common_name FROM species_table WHERE taxon = ?", (taxon,))
        return self.c.fetchone()


    def available_taxons(self):
        """
        Returns a list of all available taxa in the database.
        """
        self.c.execute("SELECT taxon, common_name FROM species_table")
        return [x for x in self.c.fetchall()]
    def available_taxa(self):
        """
        Alias for .available_taxons()
        """
        return self.available_taxons()


    def _query(self, accs, from_type, dest_types, taxon = None, require_canonical = False):
        if from_type not in dest_types:
            dest_types = [from_type] + dest_types

        type_meta = self._get_type_metadata(dest_types)
        assert(len(type_meta) == len(dest_types))
        
        if taxon is None:
            self.c.execute(f"SELECT taxon, entity_index FROM {from_type} WHERE identifier IN (%s)" % ','.join(['?']*len(accs)), accs) 
            taxon, entity_indices = zip(*self.c.fetchall())
            assert(len(set(taxon)) == 1), f"Multi-species lookup not currently supported (found taxons {', '.join(set(map(str, taxon)))}.) It is recommended to specify a taxon."
            taxon = taxon[0]
            entity_indices = list(entity_indices)
        else:
            self.c.execute(f"SELECT entity_index FROM {from_type} WHERE taxon = ? AND identifier IN (%s)" % ','.join(['?']*len(accs)), [taxon]+accs)
            entity_indices = [x[0] for x in self.c.fetchall()]

        self.c.execute(f"SELECT gene_index, mrna_index, prot_index FROM entity_table WHERE taxon = ? AND {type_meta[from_type]}_index IN ({','.join(['?']*len(entity_indices))})", 
                       [taxon]+entity_indices)

        base_query = f"SELECT et.taxon, et.gene_index, et.mrna_index, et.prot_index"

        join_clauses = []
        select_columns = []
        column_names = ['taxon', 'gene_index', 'mrna_index', 'prot_index']
        for dest_type in dest_types:
            entity_col = f"{type_meta[dest_type]}_index"
            join_clause = f"LEFT JOIN {dest_type} ON et.{entity_col} = {dest_type}.entity_index AND et.taxon = {dest_type}.taxon"
            join_clauses.append(join_clause)
            select_columns.append(f"{dest_type}.identifier AS {dest_type}_identifier")
            column_names.append(dest_type)

        final_query = base_query + ", " + ", ".join(select_columns) + " FROM entity_table et " + " ".join(join_clauses)

        final_query += f" WHERE et.taxon = ? AND et.{type_meta[from_type]}_index IN ({','.join(['?']*len(entity_indices))})"

        if require_canonical:
            # TODO test this more extensively!
            for to_type in dest_types:
                final_query += f" AND {to_type}.is_canonical = 1"

        self.c.execute(final_query, [taxon]+entity_indices)
        results = self.c.fetchall()
        result_table = pd.DataFrame(results, columns=column_names)
        return result_table[column_names[4:]]


    def map(self, ids, from_type = None, to_types = None, taxon = None, require_canonical = False,
            return_query_info = False, 
            return_format=None,
            extensive = False, 
            pretty_names = True):
        """
        Converts a set of biological identifiers from one type to another.

        Parameters:
        - ids (str or list of str): The accession identifiers to be converted. Can be a single ID as a string or a list of IDs.
        - from_type (str, optional): The type of the input identifiers. If not specified, Accessive will attempt to infer the type.
        - to_types (str or list of str, optional): The target identifier types to convert to. If not provided, defaults to all gene-level accession types.
        - taxon (str, optional): The taxonomic species identifier; this is recommended to avoid ambiguity
        - require_canonical (bool, optional): Only return canonical or 'recommended' identifiers (avoids less-common gene names, old versions of identifiers, etc.)
        - return_query_info (bool, optional): Return additional inforamtion about the query.
        - return_format (str, optional): The format of the returned data ('txt', 'json', 'pandas'). If not specified, returns a Pandas DataFrame.
        - extensive (bool, optional): Returns all relevant identifiers for the named genes/transcripts/proteins, including additional mappings back to the source accession type.
        - pretty_names (bool, optional): Renames the columns to more user-friendly names. Defaults to True.

        Returns:
        A table (in pandas Dataframe, JSON, or text TSV format) containing the requested identifiers.

        Raises:
        - Exception: If source or destination identifier types are not recognized.
        - Exception: If the return format is not recognized.

        Examples:
        Convert a single Ensembl Gene ID to UniProt and RefSeq peptide identifiers:
        >>> accessive.map(ids='ENSG00000139618', from_type='ensembl_gene', to_types=['uniprot_swissprot', 'refseq_peptide'])
        
        Convert a list of Gene Names to their corresponding HGNC identifiers without specifying source type:
        >>> accessive.map(ids=['BRCA1', 'TP53'], to_types=['hgnc'])
        """
        ids = [ids] if isinstance(ids, str) else ids 
        if isinstance(to_types, str):
            to_types = [to_types]
        assert(to_types is None or isinstance(to_types, list))
        if not from_type:
            from_type = self._get_identifier_type(ids[0])

        # if to_types is None:
        #     _, fromtype_entity = self._get_type_metadata([from_type])
        #     if fromtype_entity == 'gene':
        #         to_types = [x for x, _ in GENE_COLS]
        #     elif fromtype_entity == 'mrna':
        #         to_types = [x for x, _ in ISOFORM_COLS]
        #     elif fromtype_entity == 'prot':
        #         to_types = [x for x, _ in PROTEOFORM_COLS]
        #     else:
        #         raise Exception("Someone messed up the database and added an entity type: " + fromtype_entity)
        if to_types is None:
            to_types = [x[0] for x in GENE_COLS]

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

        result = self._query(ids, from_type, to_types, taxon, require_canonical)

        dedup_ind = result.applymap(lambda x: x if not isinstance(x, list) else ','.join(x)).drop_duplicates().index
        result = result.loc[dedup_ind]

        if not extensive:
            result = result[result[from_type].isin(ids)]

        if pretty_names:
            result.rename(columns=TO_PRETTIER_NAME, inplace=True)
    
        if return_format == 'txt':
            result = result.applymap(lambda x: x if isinstance(x, str) else ','.join(x))
            result = result.to_csv(index=False, sep='\t') 
        elif return_format == 'json':
            result = result.to_json(orient='records')
        elif return_format == 'pandas' or return_format == None:
            result = result.set_index(from_type, drop=(from_type not in to_types)) 
        else:
            raise Exception(f"Return format {return_format} is not recognized.")

        if return_query_info:
            return {'result': result, 'from_type': from_type, 'to_types': to_types, 'taxon': taxon}
        else:
            return result


    def get(self, accession, from_type, to_type, taxon = None):
        """
        Converts a single biological identifier from one type to another. Note that a list is returned to accomodate multiple mappings.

        Parameters:
        - accession (str): The accession identifier to be converted.
        - from_type (str): The type of the input identifier.
        - to_type (str): The target identifier type to convert to.
        - taxon (str, optional): The taxonomic species identifier; this is recommended to avoid ambiguity

        Returns:
        The requested identifier.

        Raises:
        - Exception: If source or destination identifier types are not recognized.

        Examples:
        Convert a single Ensembl Gene ID to UniProt and RefSeq peptide identifiers:
        >>> accessive.get(accession='ENSG00000139618', from_type='ensembl_gene', to_type='uniprot_swissprot')
        """
        return self.map(ids=accession, from_type=from_type, to_types=[to_type], taxon=taxon, return_format='pandas')[to_type].tolist()


    def make_lookup(self, from_type, to_type, taxon = None):
        """
        Creates a lookup table for converting identifiers from one type to another. Useful for scripts that require
        many repeated lookups of one accession at at time.

        Parameters:
        - from_type (str): The type of the input identifiers.
        - to_type (str): The target identifier type to convert to.
        - taxon (str, optional): The taxonomic species identifier; this is recommended to limit the size of the lookup table

        Returns:
            A dictionary mapping from_type accessions to to_type accessions

        """
        from_type = TO_DATABASE_NAME[from_type]
        to_type = TO_DATABASE_NAME[to_type]

        levels = self._get_type_metadata([from_type, to_type])

        if taxon is None:
            taxon_line = ''
        else:
            taxon_line = f"WHERE {from_type}.taxon = {taxon}"

        self.c.execute(f"""SELECT {from_type}.identifier, {to_type}.identifier 
                           FROM {from_type}
                           JOIN entity_table ON {from_type}.entity_index = entity_table.{levels[from_type]}_index
                           JOIN {to_type} ON entity_table.{levels[to_type]}_index = {to_type}.entity_index
                           {taxon_line}
                       """)
        lookup = dict(self.c.fetchall())
        return lookup

















