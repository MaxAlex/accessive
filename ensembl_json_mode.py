import json
import sqlite3
import pandas as pd




SPECIES_COLUMNS = ['name', 'display_name', 'taxon', 'valid']
IDENTIFIER_COLUMNS = ['identifier', 'identifier_type', 'referent', 'referent_type']

GENE_COLUMNS = ['ensembl_gene', 'arrayexpress', 'biogrid', 'ens_lrg_gene', 'entrez_gene', 'genecards', 'hgnc', 'mim_gene', 'pfam', 'uniprot_gene', 'wikigene',
                'isoform_list', 'proteoform_list']
ISOFORM_COLUMNS = ['ensembl_mrna', 'ccds', 'ens_lrg_transcript', 'refseq_mrna', 'refseq_ncrna', 'ucsc', 'isoform_biotype',
                   'parent_gene', 'proteoform_list']
# PROTEIN_COLUMNS = ['uniprot_swissprot', 'pdb', 'parent_gene', 'parent_isoform', 'proteoform_list']
PROTEOFORM_COLUMNS = ['ensembl_prot', 'uniparc', 'alphafold', 'uniprot_swissprot', 'uniprot_trembl', 'uniprot_isoform', 'refseq_peptide', 'embl', 'pdb',
                      'parent_gene', 'parent_isoform']

AMBIGUOUS_IDENTIFIERS = ((set(GENE_COLUMNS) & set(ISOFORM_COLUMNS) | set(GENE_COLUMNS) & set(PROTEOFORM_COLUMNS) | set(ISOFORM_COLUMNS) & set(PROTEOFORM_COLUMNS)) 
                         - set(['parent_gene', 'proteoform_list']))
assert(len(AMBIGUOUS_IDENTIFIERS) == 0), AMBIGUOUS_IDENTIFIERS
# Protein and proteoform levels aren't represented directly in the JSON strucutre; each gene has a list of transcripts, each with
# precisely one translation listed. However, multiple transcripts can be the same protein, as going by uniprot_swissprot accession. (How? Why??)
# Group transcripts/translations by swissprot identifier to get the set of "protein", and then take each translation as a separate "proteoform."
# (Noting that not all transcripts have translations.)
# NB each translation has its own ENSP, so I'ms ure that ensembl_prot goes to proteoform. Ditto for refseq_peptide.
# NB NB Not all translations have uniprot_isoforms? What does that mean?
# ........and alphafold can show up without an isoform identifier, although they all have refseq_peptide identifiers. (They're also translation-unique.)

# PDB is *probably* at gene level, because of super high nonspecificity of how PDB works apparently? But this may not always be true! Have a check.
# Pfam is inconsistent; seems mostly same across gene, but also differs by isoform!


def _build_list(things):
    assert(all([isinstance(x, str) for x in things]))  # Else use pickle?
    return 'L#|'+'|'.join(things)

def _load_list(thing):
    assert(thing[:3]=='L#|')
    return thing[3:].split('|')

def _single_item(data, key):
    try:
        thing = data[key]
    except KeyError:
        return None
    if isinstance(thing, list):
        assert(len(thing) == 1), thing
        return thing[0]
    else:
        assert(isinstance(thing, str))
        return thing

def _list_item(data, key):
    try:
        thing = data[key]
    except KeyError:
        return None
    if isinstance(thing, list):
        return _build_list(thing)
    else:
        # raise Exception((key, thing))
        assert(isinstance(thing, str))
        return thing

def load_ensembl_jsonfile(json_file, db_file = 'accessive_sqlite.db'):
    try:
        data = json.load(open(json_file, 'r')) # NB this is typically very large!
    except json.decoder.JSONDecodeError:
        import pickle
        data = pickle.load(open(json_file, 'rb'))

    conn = sqlite3.connect(db_file) 
    c = conn.cursor()

    c.execute('CREATE TABLE IF NOT EXISTS species (' + ','.join(SPECIES_COLUMNS) + ')')
    c.execute('CREATE TABLE IF NOT EXISTS identifiers (' + ','.join(IDENTIFIER_COLUMNS) + ')')
    c.execute('CREATE TABLE IF NOT EXISTS genes (' + ','.join(GENE_COLUMNS) + ')')
    c.execute('CREATE TABLE IF NOT EXISTS isoforms (' + ','.join(ISOFORM_COLUMNS) + ')')
    c.execute('CREATE TABLE IF NOT EXISTS proteoforms (' + ','.join(PROTEOFORM_COLUMNS) + ')')
    c.execute('CREATE TABLE IF NOT EXISTS entity_map (ensembl_gene, ensembl_mrna, emsembl_prot)')

    print('Current size:')
    print(c.execute('SELECT COUNT(*) FROM species').fetchall())
    print(c.execute('SELECT COUNT(*) FROM identifiers').fetchall())
    print(c.execute('SELECT COUNT(*) FROM genes').fetchall())
    print(c.execute('SELECT COUNT(*) FROM isoforms').fetchall())
    print(c.execute('SELECT COUNT(*) FROM proteoforms').fetchall())

    # Check if species already exists
    c.execute(f"SELECT valid FROM species WHERE taxon = '{data['organism']['taxonomy_id']}'")
    previous_valid = c.fetchone()
    if previous_valid is not None and previous_valid[0] == 1:
        print(f"Species {data['organism']['display_name']} already exists in database. Skipping.")
        return
    elif previous_valid is not None and previous_valid[0] == 0:
        print(f"Species {data['organism']['display_name']} previously failed to load. Reloading.")

    # Add species
    c.execute(f"INSERT OR IGNORE INTO species VALUES ('{data['organism']['name']}', '{data['organism']['display_name']}', '{data['organism']['taxonomy_id']}', 0)")
    conn.commit()
    
    for row, gene in enumerate(data['genes']):
        gene_ensembl = _single_item(gene, 'id')
        gene_arrayexpress = _list_item(gene, 'ArrayExpress')
        gene_biogrid = _list_item(gene, 'BioGRID')
        gene_ens_lrg_gene = _list_item(gene, 'ENS_LRG_gene')
        gene_entrez = _list_item(gene, 'EntrezGene')
        gene_genecards = _list_item(gene, 'GeneCards')
        gene_hgnc = _list_item(gene, 'HGNC')
        gene_mim_gene = _list_item(gene, 'MIM_GENE')
        gene_pfam = _list_item(gene, 'Pfam')
        gene_uniprot_gene = _list_item(gene, 'Uniprot_gn')
        gene_wikigene = _list_item(gene, 'WikiGene')

        child_isoforms = _build_list([x['id'] for x in gene.get('transcripts', [])])
        child_proteins = _build_list([x['translations'][0]['id'] for x in gene.get('transcripts', []) if x.get('translations', [])])

        c.execute(f"INSERT OR IGNORE INTO genes VALUES (" + ', '.join(['?' for _ in GENE_COLUMNS]) + ')', 
                  (gene_ensembl, gene_arrayexpress, gene_biogrid, gene_ens_lrg_gene, gene_entrez, gene_genecards, gene_hgnc, gene_mim_gene,
                   gene_pfam, gene_uniprot_gene, gene_wikigene, child_isoforms, child_proteins))

        gene_acc_types = [(gene_ensembl, 'ensembl_gene'), (gene_entrez, 'entrez_gene'), (gene_hgnc, 'hgnc'), (gene_uniprot_gene, 'uniprot_gene'),
                          (gene_arrayexpress, 'arrayexpress'), (gene_biogrid, 'biogrid'), (gene_ens_lrg_gene, 'ens_lrg_gene'),
                          (gene_genecards, 'genecards'), (gene_mim_gene, 'mim_gene'), (gene_pfam, 'pfam'), (gene_wikigene, 'wikigene')]
        for acc, acc_type in gene_acc_types:
            c.execute(f"INSERT OR IGNORE INTO identifiers VALUES (" + ', '.join(['?' for _ in IDENTIFIER_COLUMNS]) + ')', 
                      (acc, acc_type, gene_ensembl, 'ensembl_gene'))
        
        ## Not necessary?
        # if 'transcripts' not in gene or len(gene['transcripts']) == 0:
        #     c.execute(f"INSERT OR IGNORE INTO entity_map VALUES (?, ?, ?)", (gene_ensembl, None, None))

        for isoform in gene.get('transcripts', []):
            iso_ensembl = _single_item(isoform, 'id')
            iso_ccds = _list_item(isoform, 'CCDS')
            iso_ens_lrg_transcript = _list_item(isoform, 'ENS_LRG_transcript')
            iso_refseq_mrna = _list_item(isoform, 'RefSeq_mRNA')
            iso_refseq_ncrna = _list_item(isoform, 'RefSeq_ncRNA')
            iso_ucsc = _list_item(isoform, 'UCSC')
            iso_biotype = _list_item(isoform, 'biotype')

            child_proteins = _build_list([x['id'] for x in isoform.get('translations', [])])
            parent_gene = gene_ensembl

            c.execute(f"INSERT OR IGNORE INTO isoforms VALUES (" + ', '.join(['?' for _ in ISOFORM_COLUMNS]) + ')',
                      (iso_ensembl, iso_ccds, iso_ens_lrg_transcript, iso_refseq_mrna, iso_refseq_ncrna, iso_ucsc, iso_biotype, parent_gene, child_proteins))

            iso_acc_types = [(iso_ensembl, 'ensembl_mrna'), (iso_ccds, 'ccds'), (iso_ens_lrg_transcript, 'ens_lrg_transcript'), (iso_refseq_mrna, 'refseq_mrna'),
                             (iso_refseq_ncrna, 'refseq_ncrna'), (iso_ucsc, 'ucsc'), (iso_biotype, 'isoform_biotype')]
            for acc, acc_type in iso_acc_types:
                c.execute(f"INSERT OR IGNORE INTO identifiers VALUES (" + ', '.join(['?' for _ in IDENTIFIER_COLUMNS]) + ')',
                          (acc, acc_type, iso_ensembl, 'ensembl_transcript'))

            if 'translations' not in isoform or len(isoform['translations']) == 0:
                c.execute(f"INSERT OR IGNORE INTO entity_map VALUES (?, ?, ?)", (gene_ensembl, iso_ensembl, None))

            for proteoform in isoform.get('translations', []):
                # TODO TODO check if uniprot-isoform-lacking translations are real proteforms... once ensembl is WORKING AGAIN >:(

                pform_ensembl = _single_item(proteoform, 'id')
                pform_uniparc = _list_item(proteoform, 'UniParc')
                pform_uniprot_isoform = _list_item(proteoform, 'Uniprot_isoform')
                pform_refseq_peptide = _list_item(proteoform, 'RefSeq_peptide')
                pform_alphafold = _list_item(proteoform, 'alphafold')
                pform_embl = _list_item(proteoform, 'EMBL')
                pform_uniprot_swissprot = _list_item(proteoform, 'Uniprot/SWISSPROT') 
                pform_uniprot_trembl = _list_item(proteoform, 'Uniprot/SPTREMBL')
                pform_pdb = _list_item(proteoform, 'PDB')

                parent_gene = gene_ensembl
                parent_isoform = iso_ensembl

                c.execute(f"INSERT OR IGNORE INTO proteoforms VALUES (" + ', '.join(['?' for _ in PROTEOFORM_COLUMNS]) + ')',
                          (pform_ensembl, pform_uniparc, pform_alphafold, pform_uniprot_swissprot, pform_uniprot_trembl, pform_uniprot_isoform, pform_refseq_peptide, pform_embl, pform_pdb, parent_gene, parent_isoform))

                pform_acc_types = [(pform_ensembl, 'ensembl_prot'), (pform_uniparc, 'uniparc'), (pform_uniprot_isoform, 'uniprot_isoform'), (pform_uniprot_swissprot, 'uniprot_swissprot'),
                                   (pform_uniprot_trembl, 'uniprot_trembl'), (pform_refseq_peptide, 'refseq_peptide'), (pform_alphafold, 'alphafold')]
                for acc, acc_type in pform_acc_types:
                    c.execute(f"INSERT OR IGNORE INTO identifiers VALUES (" + ', '.join(['?' for _ in IDENTIFIER_COLUMNS]) + ')',
                              (acc, acc_type, pform_ensembl, 'ensembl_prot'))

                c.execute(f"INSERT OR IGNORE INTO entity_map VALUES (?, ?, ?)", (gene_ensembl, iso_ensembl, pform_ensembl))


        if row % 1000 == 0:
            conn.commit()
            print(f"Processed {row} genes.")

    conn.commit()
    c.execute(f"UPDATE species SET valid = 1 WHERE taxon = '{data['organism']['taxonomy_id']}'")
    print(f"Done compiling {data['organism']['display_name']} database ({row} genes.)") # type: ignore



from collections import defaultdict

class Accessive():
    def __init__(self, db_file = 'accessive_sqlite.db'):
        self.conn = sqlite3.connect(db_file)
        self.c = self.conn.cursor()

    def _identifier_type(self, accs):
        accs = [accs] if isinstance(accs, str) else accs
        query = f"SELECT identifier, referent_type FROM identifiers WHERE identifier IN ({','.join(['?']*len(accs))})"
        self.c.execute(query, accs)
        return self.c.fetchall()


    def __table_join_keys(self, source_table, dest_table):
        assert(source_table != dest_table)
        source_table_key = 'proteoform_list' if 'proteoform' in dest_table else 'isoform_list' if 'isoform' in dest_table else 'ensembl_gene'
        dest_table_key = 'ensembl_prot' if 'proteoform' in dest_table else 'ensembl_mrna' if 'isoform' in dest_table else 'parent_gene' 
        return source_table_key, dest_table_key

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

        query = f"{selection_subquery}\n{' '.join(join_subqueries)}\nWHERE\n{source_table}.{source_map_key} IN ({','.join(['?']*len(accs))})"

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



if __name__ == '__main__':
    import sys
    # load_ensembl_jsonfile(sys.argv[1]) 
    foo = Accessive('/home/max/biostuff/accessive/accessive_sqlite.db')
    bar = foo.map_identifiers(['uc002iys'], 'ucsc', ['ensembl_mrna'])
    print(bar)





