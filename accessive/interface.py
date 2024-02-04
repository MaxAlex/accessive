import psycopg2
from collections import defaultdict
import pandas as pd

from .data_structure import *
from .database_ops import load_config

class Accessive():
    def __init__(self, dbname=None, user=None, password=None, host='localhost', port=5432):
        if not (dbname and user and password):
            config = load_config()
            dbname = config['database_name']
            user = config['username']
            password = config['password']
            host = config['host']
            port = config['port']

        try:
            self.conn = psycopg2.connect(dbname=dbname, user=user, password=password, host=host, port=port)
            self.c = self.conn.cursor()
        except psycopg2.OperationalError:
            raise Exception(f"Could not connect to database {dbname} with user {user}.")
    
        self.assert_database_exists()

        self.config = dbname, user, host, port


    def assert_database_exists(self):
        self.c.execute("SELECT 1 FROM pg_tables WHERE schemaname = 'public' AND tablename = 'genes'")
        if not self.c.fetchone():
            if self.config[2]=='localhost':
                raise Exception("Accessive database seems to be missing or empty. Try running 'python -m accessive.deploy_database' to download the database.")
            else:
                raise Exception(f"Accessive database seems to be missing or empty. Ensure that the database on {self.config[2]} has been properly configured.")
        self.conn.commit()


    def _identifier_type(self, accs, taxon, require_unambiguous = False):
        accs = [accs] if isinstance(accs, str) else accs
        query = f"SELECT identifier, identifier_type FROM identifiers WHERE taxon = %s AND identifier IN ({','.join(['%s']*len(accs))})" 
        self.c.execute(query, [taxon]+accs)
        agg = defaultdict(set)
        for id, idtype in self.c.fetchall():
            agg[id].add(idtype)
        
        common_idtype = set.intersection(*agg.values())
        if not agg:
            raise Exception(f"No identifiers found for any of [{', '.join(accs)}].")
        elif not common_idtype:
            raise Exception(f"No common identifier type found for any of [{', '.join(accs)}].")
        elif len(common_idtype) > 1:
            if require_unambiguous:
                raise Exception(f"Multiple common identifier types found for [{', '.join(accs)}]: {', '.join(common_idtype)}.")
            else:
                return sorted(common_idtype, key = lambda x: (x in PROTEOFORM_COLUMNS, x in ISOFORM_COLUMNS, x in GENE_COLUMNS))[-1]
        else:
            return common_idtype.pop()



    def _query(self, accs, from_type, dest_types, taxon = None):
        source_table, source_map_key = (('proteoforms', 'ensembl_prot') if from_type in PROTEOFORM_COLUMNS 
                                           else ('isoforms', 'ensembl_mrna') if from_type in ISOFORM_COLUMNS 
                                           else ('genes', 'ensembl_gene'))
        if from_type not in dest_types:
            dest_types = [from_type] + dest_types

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

        where_subqueries = []
        if taxon:
            for join_table in set(list(dest_keys.keys()) + [source_table]):
                where_subqueries.append(f"{join_table}.taxon = {taxon} AND ")
        
        if from_type in {'ensembl_gene', 'ensembl_mrna', 'ensembl_prot', 'taxon'}:
            query = f"""{selection_subquery}
                        {' '.join(join_subqueries)}
                        WHERE
                        {' '.join(where_subqueries)}
                        (
                        {source_table}.{from_type} IN ({','.join(['%s']*len(accs))})
                        )
                        """
        else:
            # Same but we're matching by a list-type column
            listmatch_subqueries = []
            for _ in accs:
                listmatch_subqueries.append(f"%s = ANY({source_table}.{from_type})")

            query = f"""{selection_subquery}
                        {' '.join(join_subqueries)}
                        WHERE
                        {' '.join(where_subqueries)}
                        (
                        {' OR '.join(listmatch_subqueries)}
                        )
                        """

        print(query % tuple(accs))
        self.c.execute(query, accs)
        self.conn.commit()
        out = pd.DataFrame(self.c.fetchall(), columns=dest_types)
        return out


    def map(self, ids, from_type = None, to_types = None, taxon = None, return_query_info = False, return_format=None):
        ids = [ids] if isinstance(ids, str) else ids 
        assert(to_types is None or isinstance(to_types, list))
        if not from_type:
            from_type = self._identifier_type(ids, taxon)

        if to_types is None:
            if from_type in GENE_COLUMNS:
                to_types = GENE_COLUMNS
            elif from_type in ISOFORM_COLUMNS:
                to_types = ['ensembl_gene']+ISOFORM_COLUMNS
            elif from_type in PROTEOFORM_COLUMNS:
                to_types = ['ensembl_gene', 'ensembl_mrna']+PROTEOFORM_COLUMNS
            else:
                raise Exception(f"Identifier type {from_type} is not recognized.")
        elif isinstance(to_types, str):
            to_types = [to_types]
        
        result = self._query(ids, from_type, to_types, taxon)

        dedup_ind = result.applymap(lambda x: x if not isinstance(x, list) else ','.join(x)).drop_duplicates().index
        result = result.loc[dedup_ind]

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


if __name__ == '__main__':
    foo = Accessive()
    print(foo.map(['ENSG00000139618', 'ENSG00000139618.20', ], return_query_info=True, return_format='txt', taxon=9606)) 
    print(foo.map(['VGF', 'JUN'], taxon=9606))
