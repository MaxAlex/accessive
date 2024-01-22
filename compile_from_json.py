import json
import sqlite3
import pandas as pd

from data_structure import *

import psycopg2


def initialize_database(postgres_user=None, postgres_password=None):
    conn = psycopg2.connect(dbname='postgres', user=postgres_user, password=postgres_password, host='localhost', port='5432')
    c = conn.cursor()
    c.execute('CREATE DATABASE accessive')
    c.close()
    conn.close()


def initialize_tables(postgres_user=None, postgres_password=None):
    conn = psycopg2.connect(dbname='accessive', user=postgres_user, password=postgres_password, host='localhost', port='5432')
    c = conn.cursor()
    c.execute('CREATE TABLE IF NOT EXISTS species (' + ','.join(SPECIES_COLUMNS) + ')')
    c.execute('CREATE TABLE IF NOT EXISTS identifiers (' + ','.join(IDENTIFIER_COLUMNS) + ')')
    c.execute('CREATE TABLE IF NOT EXISTS genes (' + column_definitions(GENE_COLUMNS) + ')')
    c.execute('CREATE TABLE IF NOT EXISTS isoforms (' + column_definitions(ISOFORM_COLUMNS) + ')')
    c.execute('CREATE TABLE IF NOT EXISTS proteoforms (' + column_definitions(PROTEOFORM_COLUMNS) + ')')
    c.execute('CREATE TABLE IF NOT EXISTS entity_map (id SERIAL PRIMARY KEY, parent_gene TEXT, ensembl_mrna TEXT, ensembl_prot TEXT)')
    conn.commit()
    c.close()
    conn.close()


def clear_tables(postgres_user=None, postgres_password=None):
    conn = psycopg2.connect(dbname='accessive', user=postgres_user, password=postgres_password, host='localhost', port='5432')
    c = conn.cursor()
    c.execute('DROP TABLE IF EXISTS species')
    c.execute('DROP TABLE IF EXISTS identifiers')
    c.execute('DROP TABLE IF EXISTS genes')
    c.execute('DROP TABLE IF EXISTS isoforms')
    c.execute('DROP TABLE IF EXISTS proteoforms')
    c.execute('DROP TABLE IF EXISTS entity_map')
    conn.commit()
    c.close()
    conn.close()


def load_ensembl_jsonfile(json_file, db_file = 'accessive_sqlite.db'):
    try:
        data = json.load(open(json_file, 'r')) # NB this is typically very large!
    except json.decoder.JSONDecodeError:
        import pickle
        data = pickle.load(open(json_file, 'rb'))

    # conn = sqlite3.connect(db_file) 
    # c = conn.cursor()
    conn = psycopg2.connect(dbname='accessive', user='postgres', password='foobar', host='localhost', port='5432')
    c = conn.cursor()

    c.execute('CREATE TABLE IF NOT EXISTS species (' + ','.join(SPECIES_COLUMNS) + ')')
    c.execute('CREATE TABLE IF NOT EXISTS identifiers (' + ','.join(IDENTIFIER_COLUMNS) + ')')
    c.execute('CREATE TABLE IF NOT EXISTS genes (' + column_definitions(GENE_COLUMNS) + ')')
    c.execute('CREATE TABLE IF NOT EXISTS isoforms (' + column_definitions(ISOFORM_COLUMNS) + ')')
    c.execute('CREATE TABLE IF NOT EXISTS proteoforms (' + column_definitions(PROTEOFORM_COLUMNS) + ')')
    c.execute('CREATE TABLE IF NOT EXISTS entity_map (id SERIAL PRIMARY KEY, parent_gene TEXT, ensembl_mrna TEXT, ensembl_prot TEXT)')

    print('Current size:')
    c.execute('SELECT COUNT(*) FROM species')
    print(c.fetchone())
    c.execute('SELECT COUNT(*) FROM identifiers')
    print(c.fetchone())
    c.execute('SELECT COUNT(*) FROM genes')
    print(c.fetchone())
    c.execute('SELECT COUNT(*) FROM isoforms')
    print(c.fetchone())
    c.execute('SELECT COUNT(*) FROM proteoforms')
    print(c.fetchone())

    taxon_id = int(data['organism']['taxonomy_id'])

    # Check if species already exists
    c.execute(f"SELECT taxon, valid FROM species WHERE taxon = '{data['organism']['taxonomy_id']}'")
    previous_valid = c.fetchone()
    if previous_valid is not None and previous_valid[0] == 1:
        print(f"Species {data['organism']['display_name']} already exists in database. Skipping.")
        return
    elif previous_valid is not None and previous_valid[0] == 0:
        print(f"Species {data['organism']['display_name']} previously failed to load. Reloading.")

    # Add species
    c.execute(f"INSERT INTO species VALUES ('{data['organism']['name']}', '{data['organism']['display_name']}', '{data['organism']['taxonomy_id']}', false) ON CONFLICT DO NOTHING")
    conn.commit()
   
    skipped_lrg = 0  # LRG sequences are causing problemes and not even standard? Why did ensembl include them like they did? Ignore for now.

    for row, gene in enumerate(data['genes']):
        gene_ensembl = single_item(gene, 'id')
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

        c.execute(f"INSERT INTO genes VALUES (" + ', '.join(['%s' for _ in GENE_COLUMNS]) + ')', 
                  (gene_ensembl, taoxn_id, gene_arrayexpress, gene_biogrid, gene_ens_lrg_gene, gene_entrez, gene_genecards, gene_hgnc, gene_mim_gene,
                   gene_pfam, gene_uniprot_gene, gene_wikigene))

        gene_acc_types = [(gene_ensembl, 'ensembl_gene'), (gene_entrez, 'entrez_gene'), (gene_hgnc, 'hgnc'), (gene_uniprot_gene, 'uniprot_gene'),
                          (gene_arrayexpress, 'arrayexpress'), (gene_biogrid, 'biogrid'), (gene_ens_lrg_gene, 'ens_lrg_gene'),
                          (gene_genecards, 'genecards'), (gene_mim_gene, 'mim_gene'), (gene_pfam, 'pfam'), (gene_wikigene, 'wikigene')]
        for acc, acc_type in gene_acc_types:
            if acc is not None:
                c.execute(f"INSERT INTO identifiers VALUES (" + ', '.join(['%s' for _ in IDENTIFIER_COLUMNS]) + ') ON CONFLICT DO NOTHING', 
                          (acc, acc_type, gene_ensembl, 'ensembl_gene'))
        
        ## Not necessary?
        # if 'transcripts' not in gene or len(gene['transcripts']) == 0:
        #     c.execute(f"INSERT INTO entity_map VALUES (?, ?, ?)", (gene_ensembl, None, None))

        for isoform in gene.get('transcripts', []):
            iso_ensembl = single_item(isoform, 'id')
            iso_ccds = list_item(isoform, 'CCDS')
            iso_ens_lrg_transcript = list_item(isoform, 'ENS_LRG_transcript')
            iso_refseq_mrna = list_item(isoform, 'RefSeq_mRNA')
            iso_refseq_ncrna = list_item(isoform, 'RefSeq_ncRNA')
            iso_ucsc = list_item(isoform, 'UCSC')
            iso_biotype = list_item(isoform, 'biotype')

            c.execute(f"INSERT INTO isoforms VALUES (" + ', '.join(['%s' for _ in ISOFORM_COLUMNS]) + ')',
                      (iso_ensembl, taxon_id, iso_ccds, iso_ens_lrg_transcript, iso_refseq_mrna, iso_refseq_ncrna, iso_ucsc, iso_biotype))

            iso_acc_types = [(iso_ensembl, 'ensembl_mrna'), (iso_ccds, 'ccds'), (iso_ens_lrg_transcript, 'ens_lrg_transcript'), (iso_refseq_mrna, 'refseq_mrna'),
                             (iso_refseq_ncrna, 'refseq_ncrna'), (iso_ucsc, 'ucsc'), (iso_biotype, 'isoform_biotype')]
            for acc, acc_type in iso_acc_types:
                if acc is not None:
                    c.execute(f"INSERT INTO identifiers VALUES (" + ', '.join(['%s' for _ in IDENTIFIER_COLUMNS]) + ') ON CONFLICT DO NOTHING',
                              (acc, acc_type, iso_ensembl, 'ensembl_transcript'))

            if 'translations' not in isoform or len(isoform['translations']) == 0:
                c.execute(f"INSERT INTO entity_map (parent_gene, ensembl_mrna, ensembl_prot) VALUES (%s, %s, %s)", (gene_ensembl, iso_ensembl, None))

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

                c.execute(f"INSERT INTO proteoforms VALUES (" + ', '.join(['%s' for _ in PROTEOFORM_COLUMNS]) + ')',
                          (pform_ensembl, taxon_id, pform_uniparc, pform_alphafold, pform_uniprot_swissprot, pform_uniprot_trembl, 
                           pform_uniprot_isoform, pform_refseq_peptide, pform_embl, pform_pdb))

                pform_acc_types = [(pform_ensembl, 'ensembl_prot'), (pform_uniparc, 'uniparc'), (pform_uniprot_isoform, 'uniprot_isoform'), (pform_uniprot_swissprot, 'uniprot_swissprot'),
                                   (pform_uniprot_trembl, 'uniprot_trembl'), (pform_refseq_peptide, 'refseq_peptide'), (pform_alphafold, 'alphafold')]
                for acc, acc_type in pform_acc_types:
                    if acc is not None:
                        c.execute(f"INSERT INTO identifiers VALUES (" + ', '.join(['%s' for _ in IDENTIFIER_COLUMNS]) + ') ON CONFLICT DO NOTHING',
                                  (acc, acc_type, pform_ensembl, 'ensembl_prot'))

                c.execute(f"INSERT INTO entity_map (parent_gene, ensembl_mrna, ensembl_prot) VALUES (%s, %s, %s)", (gene_ensembl, iso_ensembl, pform_ensembl))


        if row % 1000 == 0:
            conn.commit()
            print(f"Processed {row} genes.")

    conn.commit()
    c.execute(f"UPDATE species SET valid = true WHERE taxon = '{data['organism']['taxonomy_id']}'")
    c.close()
    conn.close()
    print("Skipped " + str(skipped_lrg) + " LRG sequences.")
    print(f"Done compiling {data['organism']['display_name']} database ({row} genes.)") # type: ignore





if __name__ == '__main__':
    import sys
    clear_tables('postgres', 'foobar')
    load_ensembl_jsonfile(sys.argv[1]) 
    # foo = Accessive('/home/max/biostuff/accessive/accessive_sqlite.db')
    # bar = foo.map_identifiers(['uc002iys'], 'ucsc', ['ensembl_mrna'])
    # print(bar)





