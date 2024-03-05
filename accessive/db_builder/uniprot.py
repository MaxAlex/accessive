import sqlite3
import requests
import gzip
import io
import pandas as pd
import time

from ..data_structure import * 

# This takes awhile!
def download_uniprot_table(taxons):
    print("Downloading Uniprot table")
    pages = []
    res = requests.get('https://rest.uniprot.org/uniprotkb/search',
                       params = {'compressed': "true",
                                 'fields': ['accession', 'reviewed', 'organism_name', 'id', 'protein_name', 'gene_names', 'gene_oln', 
                                            'gene_primary', 'gene_synonym', 'xref_refseq', 'xref_ccds', 'xref_embl', 'xref_ensembl', 'xref_geneid',
                                            'organism_id'],
                                 'format': 'tsv',
                                 'query': '( ' + ' OR '.join(['(model_organism:%d)' % x for x in taxons]) + ' )',
                                 'size': 500})
    pages.append(gzip.decompress(res.content).decode('utf-8'))
    next_link = res.headers['Link'].split(';')[0].strip('<>')
    while next_link:
        try:
            res = requests.get(next_link)
        except:
            try:
                print("Failed to get next link, retrying in 30 seconds")
                time.sleep(30)
                res = requests.get(next_link)
            except:
                print("Failed to get next link")
                break
        pages.append(gzip.decompress(res.content).decode('utf-8'))
        assert(len(pages[-1]))
        try:
            next_link = res.headers['Link'].split(';')[0].strip('<>')
        except KeyError:
            next_link = None
        print('.', sep='', end='', flush=True)

    print("Done")

    table_name = 'uniprot_table_full.tsv'
    pd.concat(
            [pd.read_csv(io.StringIO(page), sep='\t') for page in pages]
    ).to_csv(table_name, sep='\t', index=False)

    return table_name


def add_to_table(conn, c, updating_identifier, joining_identifier, target_table, target_update_col, target_join_col):
    c.execute(f"SELECT index, identifier FROM {updating_identifier}")
    before_update = c.fetchall()
    

    target_table.to_sql('temp_table', conn, if_exists='replace', index=False)
    c.execute(f"""
        INSERT INTO {updating_identifier} ({', '.join(IDENTIFIER_TABLE_COLS)})
            SELECT {joining_identifier}.index, temp_table.{target_update_col}, temp_table.taxon, temp_table.is_canonical
            FROM temp_table
            JOIN {joining_identifier} ON temp_table.{target_join_col} = {joining_identifier}.identifier
            WHERE NOT EXISTS (
                SELECT 1 FROM {updating_identifier} WHERE {updating_identifier}.identifier = temp_table.{target_update_col}
                )
              """)
    conn.commit()
    c.execute(f"DROP TABLE temp_table")
    conn.commit()

    c.execute(f"SELECT index, identifier FROM {updating_identifier}")
    after_update = c.fetchall()

    print(f"Added {len(after_update) - len(before_update)} new entries to {updating_identifier}")
    print(f"From {len(set([x[1] for x in before_update]))} to {len(set([x[1] for x in after_update]))} unique entries (change of {len(set([x[1] for x in after_update])) - len(set([x[1] for x in before_update]))}"))




def load_uniprot_table(sqlite_file, uniprot_file_name):
    conn = sqlite3.connect(sqlite_file)
    c = conn.cursor()

    table = pd.read_csv(uniprot_file_name, sep='\t')
    gene_synonym_tab = table.apply(axis=1, func=lambda x: [(x['Entry'], gn, x['Organism (ID)']) for gn in
                                                            x['Gene Names (synonym)'].split()]
                                                            if isinstance(x['Gene Names (synonym)'], str) else [])
    gene_synonym_tab = pd.DataFrame(sum(gene_synonym_tab, []), columns=['uniport_gn', 'gene_name', 'taxon'])
    gene_synonym_tab['is_canonical'] = 0
    add_to_table(conn, c, 'gene_name', 'uniprot_gn', gene_synonym_tab, 'gene_name', 'uniport_gn')



# if __name__ == '__main__':
#     download_uniprot_table([9606, 10090, 10116, 559292])
