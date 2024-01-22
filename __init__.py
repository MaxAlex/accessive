import psycopg2
from collections import defaultdict
import pandas as pd

from data_structure import *

class Accessive():
    def __init__(self, dbname='accessive', user='postgres', password=None, host='localhost', port=5432):
        self.conn = psycopg2.connect(dbname=dbname, user=user, password=password, host=host, port=port)
        self.c = self.conn.cursor()


    def _identifier_type(self, accs):
        accs = [accs] if isinstance(accs, str) else accs
        query = f"SELECT identifier, referent_type FROM identifiers WHERE identifier IN ({','.join(['%s']*len(accs))})"
        self.c.execute(query, accs)
        return self.c.fetchall()


    def _query(self, accs, acc_type, dest_types):
        source_table, source_map_key = (('proteoforms', 'ensembl_prot') if acc_type in PROTEOFORM_COLUMNS 
                                           else ('isoforms', 'ensembl_mrna') if acc_type in ISOFORM_COLUMNS 
                                           else ('genes', 'ensembl_gene'))
        dest_keys = defaultdict(list)
        query_items = []
        for dest_type in dest_types:
            if dest_type in GENE_COLUMNS:
                dest_keys['genes'].append(dest_type)
                query_items.append(f"genes.{dest_type}")
            elif dest_type in ISOFORM_COLUMNS:
                dest_keys['isoforms'].append(dest_type)
                query_items.append(f"isoforms.{dest_type}")
            elif dest_type in PROTEOFORM_COLUMNS:
                dest_keys['proteoforms'].append(dest_type)
                query_items.append(f"proteoforms.{dest_type}")
            else:
                raise Exception(f"Identifier type {dest_type} is not recognized.") 

        selection_subquery = f"SELECT {', '.join(query_items)} FROM entity_map"

        join_subqueries = []
        for join_table in set(list(dest_keys.keys()) + [source_table]):
            # join_table_key, join_dest_key = self.__table_join_keys(source_table, dest_table)
            table_key = 'ensembl_prot' if 'proteoform' in join_table else 'ensembl_mrna' if 'isoform' in join_table else 'ensembl_gene'
            join_subqueries.append(f"JOIN {join_table} ON entity_map.{table_key} = {join_table}.{table_key}") 

        query = f"{selection_subquery}\n{' '.join(join_subqueries)}\nWHERE\n{source_table}.{source_map_key} IN ({','.join(['%s']*len(accs))})"

        self.c.execute(query, accs)
        out = pd.DataFrame(self.c.fetchall(), columns=dest_types)
        return out


    def map_identifiers(self, accs, acc_type = None, destination_types = None):
        accs = [accs] if isinstance(accs, str) else accs
        assert(destination_types is None or isinstance(destination_types, list))
        if not acc_type:
            acc_types = self._identifier_type(accs)
            assert(len(set([x[1] for x in acc_types])) == 1), "All source identifiers must be of the same type. (Found " + ', '.join(set([x[1] for x in acc_types])) + ")"
            acc_type = acc_types[0][1]

        if destination_types is None:
            if acc_type in GENE_COLUMNS:
                destination_types = GENE_COLUMNS
            elif acc_type in ISOFORM_COLUMNS:
                destination_types = ['ensembl_gene']+ISOFORM_COLUMNS
            elif acc_type in PROTEOFORM_COLUMNS:
                destination_types = ['ensembl_gene', 'ensembl_mrna']+PROTEOFORM_COLUMNS
            else:
                raise Exception(f"Identifier type {acc_type} is not recognized.")
    
        return self._query(accs, acc_type, destination_types)



