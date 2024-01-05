import subprocess 
import requests 
import pandas as pd 
import gzip
import io
import redis
import redis.exceptions
import argparse
import pickle


## CONNECT TO REDIS
if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('--host', type=str, default='localhost')
    parser.add_argument('--port', type=int, default=6379)
    parser.add_argument('--db', type=int, default=0)
    args = parser.parse_args()

    HOST, PORT, DB = args.host, args.port, args.db
else:
    HOST, PORT, DB = 'localhost', 6379, 0



### BIOMART
def download_biomart_tables():
    biomart_datasets = [('human', 'hsapiens_gene_ensembl'), 
                        ('mouse', 'mmusculus_gene_ensembl'), 
                        ('rat', 'rnorvegicus_gene_ensembl'), 
                        ('yeast', 'scerevisiae_eg_gene')]
    
    biomart_request_template = """<?xml version="1.0" encoding="UTF-8"?>
    <!DOCTYPE Query>
    <Query  virtualSchemaName = "default" formatter = "TSV" header = "1" uniqueRows = "0" count = "" datasetConfigVersion = "0.6" >
    			
    	<Dataset name = "{DATASET}" interface = "default" >
    <Attribute name = "ensembl_gene_id_version" />
    <Attribute name = "ensembl_transcript_id_version" />
    <Attribute name = "ensembl_peptide_id_version" />
        {HUMAN_ATTRIBUTES}
    	</Dataset>
    
    </Query>""".replace('\n', '')
    
    
    main_table_atts = """
    <Attribute name = "ensembl_gene_id" />
    <Attribute name = "ensembl_transcript_id" />
    <Attribute name = "ensembl_peptide_id" />
    <Attribute name = "transcript_is_canonical" />
    <Attribute name = "uniprot_gn_symbol" />
    <Attribute name = "uniprotswissprot" />
    <Attribute name = "uniprotsptrembl" />""".replace('\n', '')
    
    
    # Biomart complains if we ask for too many attributes at once (which can be as few as two or three!)
    # so we do additional queries for the rest of the attributes in pairs.
    secondary_tables_atts = ['<Attribute name = "uniprot_gn_id" /> <Attribute name = "uniprot_isoform" /> <Attribute name = "entrezgene_id" />',
                             '<Attribute name = "hgnc_id" /> <Attribute name = "hgnc_symbol" />', 
                             '<Attribute name = "refseq_mrna" /> <Attribute name = "refseq_peptide" />']
    biomart_tables = []
    for species, dataset_name in biomart_datasets[1:]:
        biomart_table_main = f'BIOMART_{species}.tsv'
        print(f"Downloading main {species} table from Biomart...")
        biomart_request = biomart_request_template.format(HUMAN_ATTRIBUTES=main_table_atts, DATASET=dataset_name)
        biomart_response = requests.post('http://www.ensembl.org/biomart/martservice', data={'query': biomart_request}, stream=True)
        with open(biomart_table_main, 'wb') as f:
            for chunk in biomart_response.iter_content(chunk_size=1024):
                if chunk:
                    f.write(chunk)
        
        
        secondary_tables = []
        for i, atts in enumerate(secondary_tables_atts, start=1):
            biomart_table = f'BIOMART_{i}_{species}.tsv'
            secondary_tables.append(biomart_table)
            print("Downloading secondary table %d from Biomart..." % i)
            biomart_request = biomart_request_template.format(HUMAN_ATTRIBUTES=atts, DATASET=dataset_name)
            biomart_response = requests.post('http://www.ensembl.org/biomart/martservice', data={'query': biomart_request}, stream=True)
            with open(biomart_table, 'wb') as f:
                for chunk in biomart_response.iter_content(chunk_size=1024):
                    if chunk:
                        f.write(chunk)
    
    
        biomart_table = pd.read_csv(biomart_table_main, sep='\t')
        for subtab in secondary_tables:
            sub = pd.read_csv(subtab, sep='\t')
            sub = sub[sub['Gene stable ID version'].notnull() & 
                      sub['Transcript stable ID version'].notnull() &
                      sub['Protein stable ID version'].notnull()]
            biomart_table = biomart_table.merge(sub, 
                                on=['Gene stable ID version', 
                                    'Transcript stable ID version', 
                                    'Protein stable ID version'],
                                how='outer')
        # table.to_csv('BIOMART.tsv', sep='\t', index=False)
        
        assert((biomart_table['Gene stable ID version'].apply(lambda x: x.split('.')[0]) == biomart_table['Gene stable ID']).all())
        assert((biomart_table['Transcript stable ID version'].apply(lambda x: x.split('.')[0]) == biomart_table['Transcript stable ID']).all())
        assert((biomart_table['Protein stable ID version'].apply(lambda x: x.split('.')[0]) == biomart_table['Protein stable ID']).all())
        
        biomart_table = biomart_table.drop(['Gene stable ID version', 'Transcript stable ID version', 'Protein stable ID version'], axis=1)
        biomart_table = biomart_table.drop_duplicates()
        
        full_table = f'BIOMART_{species}_FULL.tsv'
        biomart_table.to_csv(full_table, sep='\t', index=False) # type: ignore
        biomart_tables.append(full_table)
    
    return biomart_tables



### UNIPROT

def download_uniprot_table():
    taxons = [9606, # human 
              # 10090, # mouse
              # 10116, # rat
              559292] # yeast

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


# def load_uniprot_into_redis(table_file):
#     table = pd.read_csv(table_file, sep='\t')

# def set_obj(red, key, obj):
#     for subkey, val in obj:
#         pre_val = red.hget(key, subkey)
#         assert(pre_val is None or pre_val == val), (pre_val, val)
#         red.hset(key, subkey, val)


class EntityObj:
    def pkl(self):
        return pickle.dumps(self)

    @classmethod
    def unpkl(cls, pickled):
        return pickle.loads(pickled)


class Gene(EntityObj):
    def __init__(self, identifier):
        self.identifier = identifier
        self.uniprot_kb_gene = None
        self.ensembl_gene = None
        self.geneid = None
        self.hgnc_id = None
        self.hgnc_symbol = None

        self.transcripts = []

class Transcript(EntityObj):
    def __init__(self, identifier):
        self.identifier = identifier
        self.ensembl_trans = None
        self.refseq_mrna = None
   
        self.parent_gene = None
        self.isoforms = []

class Isoform(EntityObj):
    def __init__(self, identifier):
        self.identifier = identifier
        self.refseq_mrna = None

        self.parent_gene = None
        self.parent_transcript = None
        self.protein = None

class Protein(EntityObj):
    def __init__(self, identifier):
        self.identifier = identifier
        self.ensembl_prot = None
        self.trembl_prot = None
        self.refseq_prot = None
        self.uniprot_swissprot = None
        self.uniprot_trembl = None
        self.uniprot_isoform = None

        self.gene = None
        self.parent_transcript = None
        self.parent_isoform = None
        self.proteoforms = []

class Proteoform(EntityObj):
    def __init__(self, identifier):
        self.identifier = identifier
        self.uniprot_isoform = None
        self.ensembl_prot = None
        self.refseq_prot = None

        self.parent_gene = None
        self.parent_transcript = None
        self.parent_isoform = None
        self.parent_protein = None

# TODO determine which of these values actually fits where- what level do given accessions refer to?

# def load_biomart_tables_into_redis(biomart_tables, red):
    for biomart_table in biomart_tables:
        print(f"Loading {biomart_table} into Redis")
        table = pd.read_csv(biomart_table, sep='\t')
        for _, row in table.iterrows():
            assert(row['Gene stable ID'].startswith('ENS'))  # type: ignore
            gene_entry = 'OBJ:GENE:' + row['Gene stable ID']
            trans_entry = 'OBJ:TRANS:' + row['Transcript stable ID']
            prot_entry = 'OBJ:PROT:' + row['Protein stable ID']

            set_obj(red, gene_entry, {'uniprot_kb_gene': row['UniProtKB Gene Name symbol'],
                                      'emsembl_gene': row['Gene stable ID'],})
            set_obj(red, trans_entry, {'ensembl_trans': row['Transcript stable ID']})
            set_obj(red, prot_entry, {'ensembl_prot': row['Protein stable ID'],
                                      'trembl_prot': row['UniProtKB/TrEMBL ID']})
    


if __name__ == '__main__':
    # try:
    #     red = redis.Redis(host=host, port=port, db=db)
    #     red.ping()
    # except redis.exceptions.ConnectionError as e:
    #     print("Error on attempting to connect to Redis. Please make sure that Redis is running.")
    #     raise e

    # uniprot_table = download_uniprot_table()
    biomart_tables = download_biomart_tables()

    load_biomart_tables_into_redis(biomart_tables)
    load_uniprot_into_redis(uniprot_table)
    
    raise Exception()



