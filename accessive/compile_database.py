import json
import pandas as pd
from ftplib import FTP
import os
import gzip
import psycopg2
import psycopg2.errors
import tempfile
import io

from .data_structure import *
from .database_ops import load_config





def download_ensembl_data(include_list=None, already_loaded = [], data_save_dir = None):
    with FTP("ftp.ensembl.org") as ftp:
        ftp.login()
        ftp.cwd("pub/current_json")
        for dir_name in ftp.nlst():
            print(dir_name)
            if include_list is not None and dir_name not in include_list:
                print(f"Skipping {dir_name} because it is not in the include list.")
                continue
            elif dir_name in already_loaded:
                print(f"Skipping {dir_name} because it is already loaded.")
                continue
            ftp.cwd(dir_name)
            for item in ftp.nlst():
                if item.endswith('.json'):
                    if data_save_dir:
                        cache_file = os.path.join(data_save_dir, item)+'.gz'
                        if os.path.exists(cache_file):
                            print('Cached: %s' % cache_file)
                            yield gzip.open(cache_file, 'rb')  
                            break 
                        else:
                            print('Caching: %s' % cache_file)
                            data = gzip.open(cache_file, 'wb')
                    else:
                        data = io.BytesIO()

                    try:
                        print("Downloading " + item)
                        ftp.retrbinary('RETR ' + item, data.write)
                    except ConnectionResetError:
                        print("Connection reset error. Trying again.")
                        ftp.retrbinary('RETR ' + item, data.write)

                    if data_save_dir:
                        data.close()
                        data = gzip.open(os.path.join(data_save_dir, item)+'.gz', 'rb')  
                    else:
                        data.seek(0)

                    print("Downloaded " + item)
                    yield data 
                    break
            ftp.cwd('..')


def load_ensembl_jsonfile(json_file):
    config = load_config()
    if isinstance(json_file, str):
        try:
            data = json.load(open(json_file, 'r')) # NB this is typically very large!
        except json.decoder.JSONDecodeError:
            import pickle
            data = pickle.load(open(json_file, 'rb'))
    else:
        data = json.load(json_file)

    # conn = sqlite3.connect(db_file) 
    # c = conn.cursor()
    print(f"Connecting to database... {config['database_name']}")
    conn = psycopg2.connect(dbname=config['database_name'], user=config['username'], password=config['password'], 
                            host=config['host'], port=config['port'])
    c = conn.cursor()

    taxon_id = int(data['organism']['taxonomy_id'])
    print("Compiling " + data['organism']['display_name'] + " database (" + str(len(data['genes'])) + " genes.)")
    print(data['organism'])

    # Check if species already exists
    c.execute(f"SELECT taxon, valid FROM species WHERE taxon = '{data['organism']['taxonomy_id']}'")
    previous_valid = c.fetchone()
    if previous_valid is not None and previous_valid[1] == True:
        print(f"Species {data['organism']['display_name']} already exists in database. Skipping.")
        return
    elif previous_valid is not None and previous_valid[1] == False:
        print(f"Species {data['organism']['display_name']} previously failed to load. Reloading.")

    c.execute(f'CREATE TABLE IF NOT EXISTS genes_{taxon_id} PARTITION OF genes FOR VALUES IN ({taxon_id})')
    c.execute(f'CREATE TABLE IF NOT EXISTS isoforms_{taxon_id} PARTITION OF isoforms FOR VALUES IN ({taxon_id})')
    c.execute(f'CREATE TABLE IF NOT EXISTS proteoforms_{taxon_id} PARTITION OF proteoforms FOR VALUES IN ({taxon_id})')
    c.execute(f'CREATE TABLE IF NOT EXISTS entity_map_{taxon_id} PARTITION OF entity_map FOR VALUES IN ({taxon_id})')

    # Add species
    singlequote = "'"
    c.execute(f"INSERT INTO species VALUES ('{data['organism']['name']}', '{data['organism']['display_name'].replace(singlequote, '')}', '{data['organism']['taxonomy_id']}', false) ON CONFLICT DO NOTHING")
    conn.commit()
   
    skipped_lrg = 0  # LRG sequences are causing problemes and not even standard? Why did ensembl include them like they did? Ignore for now.

    for row, gene in enumerate(data['genes']):
        gene_ensembl = single_item(gene, 'id')
        gene_name = list_item(gene, 'name')
        gene_description = list_item(gene, 'description')
        gene_arrayexpress = list_item(gene, 'ArrayExpress')
        gene_biogrid = list_item(gene, 'BioGRID')
        gene_ens_lrg_gene = list_item(gene, 'ENS_LRG_gene')
        gene_entrez = list_item(gene, 'EntrezGene')
        gene_genecards = list_item(gene, 'GeneCards')
        gene_hgnc = list_item(gene, 'HGNC')
        gene_mim_gene = list_item(gene, 'MIM_GENE')
        gene_pfam = list_item(gene, 'Pfam')
        gene_uniprot_gene = list_item(gene, 'Uniprot_gn')
        gene_wikigene = list_item(gene, 'WikiGene')

        if gene_ensembl is None:
            raise Exception(f"Gene {row} has no ensembl identifier.")
        elif gene_ensembl[:3] == 'LRG':
            skipped_lrg += 1
            continue

        c.execute(f"INSERT INTO genes ({', '.join(GENE_COLUMNS)}) VALUES (" + ', '.join(['%s' for _ in GENE_COLUMNS]) + ')', 
                  (gene_ensembl, taxon_id, gene_description, gene_name, gene_arrayexpress, gene_biogrid, gene_ens_lrg_gene, gene_entrez, gene_genecards, gene_hgnc, gene_mim_gene,
                   gene_pfam, gene_uniprot_gene, gene_wikigene))

        gene_acc_types = [(gene_ensembl, 'ensembl_gene'), (gene_name, 'gene_name'), (gene_entrez, 'entrez_gene'), (gene_hgnc, 'hgnc'), (gene_uniprot_gene, 'uniprot_gene'),
                          (gene_arrayexpress, 'arrayexpress'), (gene_biogrid, 'biogrid'), (gene_ens_lrg_gene, 'ens_lrg_gene'),
                          (gene_genecards, 'genecards'), (gene_mim_gene, 'mim_gene'), (gene_pfam, 'pfam'), (gene_wikigene, 'wikigene')]
        for acc, acc_type in gene_acc_types:
            if acc is None:
                continue
            elif isinstance(acc, str):
                c.execute(f"INSERT INTO identifiers VALUES (" + ', '.join(['%s' for _ in IDENTIFIER_COLUMNS]) + ') ON CONFLICT DO NOTHING', 
                          (acc, acc_type, gene_ensembl, 'ensembl_gene', taxon_id))
            elif isinstance(acc, list):
                for acc_ in acc:
                    c.execute(f"INSERT INTO identifiers VALUES (" + ', '.join(['%s' for _ in IDENTIFIER_COLUMNS]) + ') ON CONFLICT DO NOTHING', 
                              (acc_, acc_type, gene_ensembl, 'ensembl_gene', taxon_id))
        
        for isoform in gene.get('transcripts', []):
            iso_ensembl = single_item(isoform, 'id')
            iso_ccds = list_item(isoform, 'CCDS')
            iso_ens_lrg_transcript = list_item(isoform, 'ENS_LRG_transcript')
            iso_refseq_mrna = list_item(isoform, 'RefSeq_mRNA')
            iso_refseq_ncrna = list_item(isoform, 'RefSeq_ncRNA')
            iso_ucsc = list_item(isoform, 'UCSC')
            iso_biotype = list_item(isoform, 'biotype')

            c.execute(f"INSERT INTO isoforms ({', '.join(ISOFORM_COLUMNS)}) VALUES (" + ', '.join(['%s' for _ in ISOFORM_COLUMNS]) + ')',
                      (iso_ensembl, taxon_id, iso_ccds, iso_ens_lrg_transcript, iso_refseq_mrna, iso_refseq_ncrna, iso_ucsc, iso_biotype))

            iso_acc_types = [(iso_ensembl, 'ensembl_mrna'), (iso_ccds, 'ccds'), (iso_ens_lrg_transcript, 'ens_lrg_transcript'), (iso_refseq_mrna, 'refseq_mrna'),
                             (iso_refseq_ncrna, 'refseq_ncrna'), (iso_ucsc, 'ucsc'), (iso_biotype, 'isoform_biotype')]
            for acc, acc_type in iso_acc_types:
                if acc is None:
                    continue
                elif isinstance(acc, str):
                    c.execute(f"INSERT INTO identifiers VALUES (" + ', '.join(['%s' for _ in IDENTIFIER_COLUMNS]) + ') ON CONFLICT DO NOTHING', 
                              (acc, acc_type, iso_ensembl, 'ensembl_mrna', taxon_id))
                elif isinstance(acc, list):
                    for acc_ in acc:
                        c.execute(f"INSERT INTO identifiers VALUES (" + ', '.join(['%s' for _ in IDENTIFIER_COLUMNS]) + ') ON CONFLICT DO NOTHING', 
                                  (acc_, acc_type, iso_ensembl, 'ensembl_mrna', taxon_id))

            if 'translations' not in isoform or len(isoform['translations']) == 0:
                c.execute(f"INSERT INTO entity_map (ensembl_gene, ensembl_mrna, ensembl_prot, taxon) VALUES (%s, %s, %s, %s)", (gene_ensembl, iso_ensembl, None, taxon_id))

            for proteoform in isoform.get('translations', []):
                # TODO TODO check if uniprot-isoform-lacking translations are real proteforms... once ensembl is WORKING AGAIN >:(

                pform_ensembl = single_item(proteoform, 'id')
                pform_uniparc = list_item(proteoform, 'UniParc')
                pform_uniprot_isoform = list_item(proteoform, 'Uniprot_isoform')
                pform_refseq_peptide = list_item(proteoform, 'RefSeq_peptide')
                pform_alphafold = list_item(proteoform, 'alphafold')
                pform_embl = list_item(proteoform, 'EMBL')
                pform_uniprot_swissprot = list_item(proteoform, 'Uniprot/SWISSPROT') 
                pform_uniprot_trembl = list_item(proteoform, 'Uniprot/SPTREMBL')
                pform_pdb = list_item(proteoform, 'PDB')

                c.execute(f"INSERT INTO proteoforms ({', '.join(PROTEOFORM_COLUMNS)}) VALUES (" + ', '.join(['%s' for _ in PROTEOFORM_COLUMNS]) + ')',
                          (pform_ensembl, taxon_id, pform_uniparc, pform_alphafold, pform_uniprot_swissprot, pform_uniprot_trembl, 
                           pform_uniprot_isoform, pform_refseq_peptide, pform_embl, pform_pdb))

                pform_acc_types = [(pform_ensembl, 'ensembl_prot'), (pform_uniparc, 'uniparc'), (pform_uniprot_isoform, 'uniprot_isoform'), (pform_uniprot_swissprot, 'uniprot_swissprot'),
                                   (pform_uniprot_trembl, 'uniprot_trembl'), (pform_refseq_peptide, 'refseq_peptide'), (pform_alphafold, 'alphafold')]
                for acc, acc_type in pform_acc_types:
                    if acc is None:
                        continue
                    elif isinstance(acc, str):
                        c.execute(f"INSERT INTO identifiers VALUES (" + ', '.join(['%s' for _ in IDENTIFIER_COLUMNS]) + ') ON CONFLICT DO NOTHING', 
                                  (acc, acc_type, pform_ensembl, 'ensembl_prot', taxon_id))
                    elif isinstance(acc, list):
                        for acc_ in acc:
                            c.execute(f"INSERT INTO identifiers VALUES (" + ', '.join(['%s' for _ in IDENTIFIER_COLUMNS]) + ') ON CONFLICT DO NOTHING', 
                                      (acc_, acc_type, pform_ensembl, 'ensembl_prot', taxon_id))

                c.execute(f"INSERT INTO entity_map (ensembl_gene, ensembl_mrna, ensembl_prot, taxon) VALUES (%s, %s, %s, %s)", (gene_ensembl, iso_ensembl, pform_ensembl, taxon_id))


        if row % 1000 == 0:
            conn.commit()
            print(f"Processed {row} genes.")

    c.execute(f"UPDATE species SET valid = true WHERE taxon = '{data['organism']['taxonomy_id']}'")
    conn.commit()
    c.close()
    conn.close()
    print("Skipped " + str(skipped_lrg) + " LRG sequences.")
    print(f"Done compiling {data['organism']['display_name']} database ({row} genes.)") # type: ignore



def current_species_manifest():
    config = load_config()
    conn = psycopg2.connect(dbname=config['database_name'], user=config['username'], password=config['password'], 
                            host=config['host'], port=config['port'])
    c = conn.cursor()
    try:
        c.execute('SELECT * FROM species WHERE valid = true')
        species = [x[0] for x in c.fetchall()]
        c.close()
        conn.close()
        return species 
    except psycopg2.errors.UndefinedTable:
        c.close()
        conn.close()
        return []

# TODO TODO test this out!
def compile_full_database(include_list=None, cache_dir=None):
    assert(os.path.exists(cache_dir))
    already_loaded = current_species_manifest()

    for data_buffer in download_ensembl_data(include_list, already_loaded, cache_dir):
        load_ensembl_jsonfile(data_buffer)
        del data_buffer


from .database_ops import initialize_database, initialize_tables, clear_tables
if __name__ == '__main__':
    import sys
    # clear_tables()
    # load_ensembl_jsonfile(sys.argv[1]) 
    # foo = Accessive('/home/max/biostuff/accessive/accessive_sqlite.db')
    # bar = foo.map_identifiers(['uc002iys'], 'ucsc', ['ensembl_mrna'])
    # print(bar)
    species_list = open('/home/max/biostuff/accessive/shorter_species_list.txt').read().strip().split()
    print(species_list[:10])
    initialize_tables()
    compile_full_database(species_list, cache_dir = '/data/biostuff/ensembl_data/cache')
    # load_ensembl_jsonfile('/data/biostuff/ensembl_data/homo_sapiens_PARTIAL.json')   





