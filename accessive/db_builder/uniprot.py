import sqlite3
import requests
import gzip
import io
import pandas as pd


def download_uniprot_table(taxons):
    # taxons = [
    #           9606, # human 
    #           # 10090, # mouse
    #           # 10116, # rat
    #           # 559292] # yeast
    #             ]

    print("Downloading Uniprot table")
    pages = []
    res = requests.get('https://rest.uniprot.org/uniprotkb/search',
                       params = {'compressed': "true",
                                 'fields': ['accession', 'reviewed', 'organism_name', 'id', 'protein_name', 'gene_names', 'gene_oln', 
                                            'gene_primary', 'gene_synonym', 'xref_refseq', 'xref_ccds', 'xref_embl', 'xref_ensembl', 'xref_geneid'],
                                 'format': 'tsv',
                                 'query': '( ' + ' OR '.join(['(model_organism:%d)' % x for x in taxons]) + ' )',
                                 'size': 500})
    pages.append(gzip.decompress(res.content).decode('utf-8'))
    next_link = res.headers['Link'].split(';')[0].strip('<>')
    while next_link:
        res = requests.get(next_link)
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






def load_uniprot_table(sqlite_file, table_name):
    # TODO 
    pass
