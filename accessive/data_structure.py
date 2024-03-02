
ENTITY_TABLE_COLS = ['taxon INTEGER', 'gene_index INTEGER', 'mrna_index INTEGER', 'prot_index INTEGER']
IDENTIFIER_TABLE_COLS = ['entity_index INTEGER', 'identifier TEXT', 'taxon INTEGER', 'is_canonical INTEGER'] # Whether the index is _gene or etc is determined in metadata table
METADATA_COLS = ['identifier_type TEXT', 'entity_type TEXT']
DIRECTORY_COLS = ['identifier TEXT', 'identifier_type TEXT']
SPECIES_COLS = ['taxon INTEGER', 'name TEXT', 'common_name TEXT']

GENE_COLS = [('ensembl_gene', 'id'), ('gene_description', 'description'), ('gene_name', 'name'), ('arrayexpress', 'ArrayExpress'), ('biogrid', 'BioGRID'), 
                     ('ens_lrg_gene', 'ENS_LRG_gene'), ('entrez_gene', 'EntrezGene'), ('genecards', 'GeneCards'), ('hgnc', 'HGNC'), 
                     ('mim_gene', 'MIM_GENE'), ('pfam', 'Pfam'), ('uniprot_gene', 'Uniprot_gn'), ('wikigene', 'WikiGene')]
ISOFORM_COLS = [('ensembl_mrna', 'id'), ('ccds', 'CCDS'), ('ens_lrg_transcript', 'ENS_LRG_transcript'), ('refseq_mrna', 'RefSeq_mRNA'), ('refseq_ncrna', 'RefSeq_ncRNA'),
                        ('ucsc', 'UCSC'), ('isoform_biotype', 'biotype')]
PROTEOFORM_COLS = [('ensembl_prot', 'id'), ('uniparc', 'UniParc'), ('alphafold', 'alphafold'), ('uniprot_swissprot', 'Uniprot/SWISSPROT'), ('uniprot_trembl', 'Uniprot/SPTREMBL'),
                           ('uniprot_isoform', 'Uniprot_isoform'), ('refseq_peptide', 'RefSeq_peptide'), ('embl', 'EMBL'), ('pdb', 'PDB')]


TO_DATABASE_NAME = dict(GENE_COLS + ISOFORM_COLS + PROTEOFORM_COLS)
TO_PRETTIER_NAME = {x: y for y, x in TO_DATABASE_NAME.items()}
TO_DATABASE_NAME.update({x:x for x in TO_DATABASE_NAME.keys()})
TO_PRETTIER_NAME.update({x:x for x in TO_PRETTIER_NAME.keys()})


# AMBIGUOUS_IDENTIFIERS = ((set(GENE_COLUMNS) & set(ISOFORM_COLUMNS) | set(GENE_COLUMNS) & set(PROTEOFORM_COLUMNS) | set(ISOFORM_COLUMNS) & set(PROTEOFORM_COLUMNS)) 
#                          - set(['parent_gene', 'proteoform_list']))
# assert(len(AMBIGUOUS_IDENTIFIERS) == 1), AMBIGUOUS_IDENTIFIERS
# Protein and proteoform levels aren't represented directly in the JSON strucutre; each gene has a list of transcripts, each with
# precisely one translation listed. However, multiple transcripts can be the same protein, as going by uniprot_swissprot accession. (How? Why??)
# Group transcripts/translations by swissprot identifier to get the set of "protein", and then take each translation as a separate "proteoform."
# (Noting that not all transcripts have translations.)
# NB each translation has its own ENSP, so I'ms ure that ensembl_prot goes to proteoform. Ditto for refseq_peptide.
# NB NB Not all translations have uniprot_isoforms? What does that mean?
# ........and alphafold can show up without an isoform identifier, although they all have refseq_peptide identifiers. (They're also translation-unique.)

# PDB is *probably* at gene level, because of super high nonspecificity of how PDB works apparently? But this may not always be true! Have a check.
# Pfam is inconsistent; seems mostly same across gene, but also differs by isoform!


def column_definitions(columns):
    return ', '.join([f"{columns[0]} TEXT", f"{columns[1]} INTEGER"] + [f"{x} TEXT[]" for x in columns[2:]])


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
