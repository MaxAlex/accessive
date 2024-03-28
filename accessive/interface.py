import os
import sqlite3
import pandas as pd

from .data_structure import *
from .database_ops import DATABASE_FILE


GENE_COLS = ['ensembl_gene', 'gene_description', 'gene_name', 'arrayexpress', 'biogrid', 'ens_lrg_gene', 'entrez_gene', 
             'genecards', 'hgnc', 'mim_gene', 'pfam', 'uniprot_gene', 'wikigene', 'nextprot']
ISOFORM_COLS = ['ensembl_mrna', 'ccds', 'ens_lrg_transcript', 'refseq_mrna', 'refseq_ncrna', 'ucsc', 'isoform_biotype', 'nextprot_isoform']
PROTEOFORM_COLS = ['ensembl_prot', 'uniparc', 'alphafold', 'uniprot_swissprot', 'uniprot_trembl', 'uniprot_isoform', 'refseq_peptide', 'embl', 'pdb']


KNOWN_IDENTIFIERS = set(GENE_COLS + ISOFORM_COLS + PROTEOFORM_COLS)

class Accessive():
    def __init__(self, sqlite_file = None, 
                 default_from_type = None,
                 default_to_types = None,
                 default_format = 'pandas', 
                 default_taxon = None, 
                 default_require_canonical = False):
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

        self.default_from_type = default_from_type
        self.default_to_types = default_to_types
        self.default_format = default_format
        self.default_taxon = default_taxon
        self.default_require_canonical = default_require_canonical

        assert(self.default_from_type is None or self.default_from_type in TO_DATABASE_NAME), f"default_from_type {self.default_from_type} is not recognized."
        assert(self.default_to_types is None or all(x in TO_DATABASE_NAME for x in self.default_to_types)), f"Some default to_types are not recognized: {[x for x in self.default_to_types if x not in TO_DATABASE_NAME]}"
        assert(self.default_format in ['pandas', 'json', 'txt']), "default_format must be one of 'pandas', 'json', or 'txt'."
        assert(self.default_taxon is None or isinstance(self.default_taxon, int)), "default_taxon must be an integer."
        assert(isinstance(self.default_require_canonical, bool)), "default_require_canonical must be a boolean."



    def _get_db_version(self):
        self.c.execute("SELECT val FROM accessive_meta WHERE key = 'database_version'")
        return self.c.fetchone()[0]


    def _get_identifier_type(self, acc, allow_multiple = False):
        self.c.execute("SELECT identifier_type FROM identifier_directory WHERE identifier = ?", (acc,))
        types = list(set([x[0] for x in self.c.fetchall()]))
        if not allow_multiple:
            if len(types) > 1:
                raise Exception(f"Identifier {acc} is associated with multiple types: {', '.join(types)}")
            else:
                return types[0]
        else:
            return types


    def _get_type_metadata(self, idtypes):
        self.c.execute("SELECT * FROM metadata_table WHERE identifier_type IN (%s)" % ','.join(['?']*len(idtypes)), idtypes)
        return dict(self.c.fetchall())

    
    def identify_taxon(self, taxon):
        """
        Utility function to identify what species a taxon number in the database corresponds to.
        """
        self.c.execute("SELECT name, common_name FROM species_table WHERE taxon = ?", (taxon,))
        return self.c.fetchone()

    
    def identify(self, acc):
        """
        Identifies the type(s) of an accession identifier.

        Parameters:
        - acc (str): The accession identifier to be identified.

        Returns:
        A list of strings representing the type(s) of the identifier, if any.
        """
        return self._get_identifier_type(acc)


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


    def map(self, ids, from_type = None, to_types = None, taxon = None, require_canonical = None,
            format=None,
            return_query_info = False, 
            extensive = False, 
            ):
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

        if to_types is None:
            if self.default_to_types is not None:
                to_types = self.default_to_types
            else:
                to_types = [x for x in GENE_COLS]

        if not from_type:
            if self.default_from_type is not None:
                from_type = self.default_from_type
            else:
                ## Automatic from_type inference is disabled for now; I have mixed feelings about making it easy
                ## to introduce effectively-nondeterministic behavior into a scientific workflows.
                # from_type = self._get_identifier_type(ids[0])
                raise ValueError("Missing from_type")

        if len(to_types) < len(set(to_types)):
            to_types = list(set(to_types))

        if taxon is None:
            taxon = self.default_taxon

        if require_canonical is None:
            require_canonical = self.default_require_canonical

        if from_type not in KNOWN_IDENTIFIERS:
            raise ValueError(f"Source identifier type {from_type} is not recognized.")
        if not all(x in KNOWN_IDENTIFIERS for x in to_types):
            raise ValueError(f"Some destination identifier types are not recognized: {[x for x in to_types if x not in KNOWN_IDENTIFIERS]}")


        result = self._query(ids, from_type, to_types, taxon, require_canonical)

        try:
            dedup_ind = result.map(lambda x: x if not isinstance(x, list) else ','.join(x)).drop_duplicates().index
        except AttributeError:
            # Pandas went from not having .map() at all to having an annoying deprecation message on .applymap(), in one version! 
            dedup_ind = result.applymap(lambda x: x if not isinstance(x, list) else ','.join(x)).drop_duplicates().index
        result = result.loc[dedup_ind]

        if not extensive:
            result = result[result[from_type].isin(ids)]

        result = result.set_index(from_type, drop=(from_type not in to_types)) 
        
        # Lots of queries will return all-None rows, for various complicated reasons, usually of the form 
        # "rows correspond to proteoforms since a proteoform accession was requested, but some genes/transcripts
        # in the result are non-coding or missing" 
        result = result[~result.isnull().all(axis=1)]

        if format == 'txt':
            result = result.applymap(lambda x: ', '.join(x) if isinstance(x, list) else x)
            result = result.to_csv(sep='\t') 
        elif format == 'json' or format == 'dict':
            # result = result.to_dict(orient='dict')
            ## None of the pandas to_dict options do quite what we want here.
            d_lookup = {ftype:{ttype:[] for ttype in result.columns} for ftype in set(result.index)}
            for acc, row in result.iterrows():
                for acctype in result.columns:
                    if row[acctype] is not None:
                        d_lookup[acc][acctype].append(row[acctype])
            result = d_lookup
        elif format == 'pandas' or format == None:
            pass
        else:
            raise Exception(f"Return format {format} is not recognized.")
        
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
        return self.map(ids=accession, from_type=from_type, to_types=[to_type], taxon=taxon, format='pandas')[to_type].tolist() # type: ignore













