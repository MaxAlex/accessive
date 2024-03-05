DATABASE_VERSION = '0.1'

ENTITY_TABLE_COLS = ['taxon INTEGER', 'gene_index INTEGER', 'mrna_index INTEGER', 'prot_index INTEGER']
IDENTIFIER_TABLE_COLS = ['entity_index INTEGER', 'identifier TEXT', 'taxon INTEGER', 'is_canonical INTEGER'] # Whether the index is _gene or etc is determined in metadata table
METADATA_COLS = ['identifier_type TEXT', 'entity_type TEXT']
DIRECTORY_COLS = ['identifier TEXT', 'identifier_type TEXT']
SPECIES_COLS = ['taxon INTEGER', 'name TEXT', 'common_name TEXT']

GENE_COLS = [('ensembl_gene', 'id'), ('gene_description', 'description'), ('gene_name', 'name'), ('arrayexpress', 'ArrayExpress'), ('biogrid', 'BioGRID'), 
             ('ens_lrg_gene', 'ENS_LRG_gene'), ('entrez_gene', 'EntrezGene'), ('genecards', 'GeneCards'), ('hgnc', 'HGNC'), 
             ('mim_gene', 'MIM_GENE'), ('pfam', 'Pfam'), ('uniprot_gene', 'Uniprot_gn'), ('wikigene', 'WikiGene'),
             ('nextprot', 'Nextprot')]
ISOFORM_COLS = [('ensembl_mrna', 'id'), ('ccds', 'CCDS'), ('ens_lrg_transcript', 'ENS_LRG_transcript'), ('refseq_mrna', 'RefSeq_mRNA'), ('refseq_ncrna', 'RefSeq_ncRNA'),
                ('ucsc', 'UCSC'), ('isoform_biotype', 'biotype'), ('nextprot_isoform', 'NextProt_isoform')]
PROTEOFORM_COLS = [('ensembl_prot', 'id'), ('uniparc', 'UniParc'), ('alphafold', 'alphafold'), ('uniprot_swissprot', 'Uniprot/SWISSPROT'), ('uniprot_trembl', 'Uniprot/SPTREMBL'),
                   ('uniprot_isoform', 'Uniprot_isoform'), ('refseq_peptide', 'RefSeq_peptide'), ('embl', 'EMBL'), ('pdb', 'PDB')]

# Some remaining notes re levels:
# PDB is *probably* always constant at gene level, because of super high nonspecificity of how PDB works- though it seems like it should be proteoform-specific.
# Pfam is inconsistent. I'm calling it gene-level for now based on how Ensembl handles it, but again, logically it ought to be proteoform!

TO_DATABASE_NAME = dict(GENE_COLS + ISOFORM_COLS + PROTEOFORM_COLS)
TO_PRETTIER_NAME = {x: y for y, x in TO_DATABASE_NAME.items()}
TO_DATABASE_NAME.update({x:x for x in TO_DATABASE_NAME.keys()})
TO_PRETTIER_NAME.update({x:x for x in TO_PRETTIER_NAME.keys()})


def single_item(data, key):
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

def list_item(data, key):
    try:
        thing = data[key]
    except KeyError:
        return []
    if isinstance(thing, list):
        return thing
    else:
        assert(isinstance(thing, str))
        return [thing]
