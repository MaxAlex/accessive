
SPECIES_COLUMNS = ['name TEXT UNIQUE', 'display_name TEXT UNIQUE', 'taxon TEXT UNIQUE', 'valid BOOLEAN']
IDENTIFIER_COLUMNS = ['identifier TEXT', 'identifier_type TEXT', 'referent TEXT', 'referent_type TEXT', 'taxon INTEGER']

GENE_COLUMNS = ['ensembl_gene', 'taxon', 'gene_description', 'gene_name', 'arrayexpress', 'biogrid', 'ens_lrg_gene', 'entrez_gene', 'genecards', 'hgnc', 'mim_gene', 'pfam', 'uniprot_gene', 'wikigene']
ISOFORM_COLUMNS = ['ensembl_mrna', 'taxon', 'ccds', 'ens_lrg_transcript', 'refseq_mrna', 'refseq_ncrna', 'ucsc', 'isoform_biotype']
PROTEOFORM_COLUMNS = ['ensembl_prot', 'taxon', 'uniparc', 'alphafold', 'uniprot_swissprot', 'uniprot_trembl', 'uniprot_isoform', 'refseq_peptide', 'embl', 'pdb']

AMBIGUOUS_IDENTIFIERS = ((set(GENE_COLUMNS) & set(ISOFORM_COLUMNS) | set(GENE_COLUMNS) & set(PROTEOFORM_COLUMNS) | set(ISOFORM_COLUMNS) & set(PROTEOFORM_COLUMNS)) 
                         - set(['parent_gene', 'proteoform_list']))
assert(len(AMBIGUOUS_IDENTIFIERS) == 1), AMBIGUOUS_IDENTIFIERS
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
        return None
    if isinstance(thing, list):
        return thing
    else:
        assert(isinstance(thing, str))
        return [thing]
