DATABASE_VERSION = '0.1'

ENTITY_TABLE_COLS = ['taxon INTEGER', 'gene_index INTEGER', 'mrna_index INTEGER', 'prot_index INTEGER']
IDENTIFIER_TABLE_COLS = ['entity_index INTEGER', 'identifier TEXT', 'taxon INTEGER', 'is_canonical INTEGER'] # Whether the index is _gene or etc is determined in metadata table
METADATA_COLS = ['identifier_type TEXT', 'entity_type TEXT']
DIRECTORY_COLS = ['identifier TEXT', 'identifier_type TEXT']
SPECIES_COLS = ['taxon INTEGER', 'name TEXT', 'common_name TEXT']



GENE_COLS = ['ensembl_gene', 'gene_description', 'gene_name', 'arrayexpress', 'biogrid', 'ens_lrg_gene', 'entrez_gene', 
             'genecards', 'hgnc', 'mim_gene', 'pfam', 'uniprot_gene', 'wikigene', 'nextprot']
ISOFORM_COLS = ['ensembl_mrna', 'ccds', 'ens_lrg_transcript', 'refseq_mrna', 'refseq_ncrna', 'ucsc', 'isoform_biotype', 'nextprot_isoform']
PROTEOFORM_COLS = ['ensembl_prot', 'uniparc', 'alphafold', 'uniprot_swissprot', 'uniprot_trembl', 'uniprot_isoform', 'refseq_peptide', 'embl', 'pdb']


KNOWN_IDENTIFIERS = set(GENE_COLS + ISOFORM_COLS + PROTEOFORM_COLS)
# Some remaining notes re levels:
# PDB is *probably* always constant at gene level, because of super high nonspecificity of how PDB works- though it seems like it should be proteoform-specific.
# Pfam is inconsistent. I'm calling it gene-level for now based on how Ensembl handles it, but again, logically it ought to be proteoform!



