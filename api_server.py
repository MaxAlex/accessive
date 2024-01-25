from flask import Flask, request, jsonify
from accessor import Accessive


app = Flask(__name__)
database = Accessive('accessive', 'postgres', 'postgres')

@app.route('/hello', methods=['GET'])
def hello():
    return "This whole thing is sort of accessive, isn't it?"


@app.route('/map_identifiers', methods=['GET'])
def map_identifiers():
    """
    Map a list of identifiers from one type to another.

    Parameters:
        ids: A list of identifiers to map.
        from_type: The type of the identifiers to map.
        to_types: A list of types to map the identifiers to.
        taxon: The taxon to map the identifiers in.

    Returns:
        A JSON object containing the mapping results.

    
    Available identifier types:
        Gene level:
            ensembl_gene : Ensembl Gene (Ensembl gene identifier)
            gene_description : Short plaintext description of the gene
            gene_name : Gene Name (Name of the gene)
            arrayexpress : ArrayExpress (Microarray data repository)
            biogrid : BioGRID (Biological interaction database)
            ens_lrg_gene : ENS_LRG_gene (Locus Reference Genomic gene ID)
            entrez_gene : EntrezGene (NCBI gene identifier)
            genecards : GeneCards (Human gene compendium)
            hgnc : HGNC (HUGO Gene Nomenclature Committee ID)
            mim_gene : MIM_GENE (Mendelian Inheritance in Man gene ID)
            pfam : Pfam (Protein family database)
            uniprot_gene : Uniprot_gn (UniProt gene name)
            wikigene : WikiGene (Collaborative gene wiki)
        Transcript level:
            ensembl_mrna : Ensembl mRNA (Ensembl mRNA identifier)
            ccds : CCDS (Consensus CDS project ID)
            ens_lrg_transcript : ENS_LRG_transcript (LRG transcript identifier)
            refseq_mrna : RefSeq mRNA (NCBI mRNA reference sequence)
            refseq_ncrna : RefSeq ncRNA (NCBI non-coding RNA reference)
            ucsc : UCSC (UCSC Genome Browser ID)
            isoform_biotype : Isoform Biotype (Gene isoform biotype classification)
        Protein level:
            ensembl_prot : Ensembl Protein (Ensembl protein identifier)
            uniparc : UniParc (UniProt Archive identifier)
            alphafold : AlphaFold (Protein structure prediction database)
            uniprot_swissprot : Uniprot/SWISSPROT (UniProtKB/Swiss-Prot ID)
            uniprot_trembl : Uniprot/SPTREMBL (UniProtKB/TrEMBL ID)
            uniprot_isoform : Uniprot Isoform (UniProt isoform identifier)
            refseq_peptide : RefSeq Peptide (NCBI peptide reference sequence)
            embl : EMBL (European Molecular Biology Lab sequence ID)
            pdb : PDB (Protein Data Bank identifier)
        Additional identifiers:
            taxon : Taxon (NCBI Taxonomy identifier)

    """

    try:
        # Extract parameters from the request
        ids = request.args.getlist('ids')  # Expects a list of identifiers
        from_type = request.args.get('from_type')
        to_types = request.args.getlist('to_types') if 'to_types' in request.args else None
        taxon = request.args.get('taxon')

        # Call the map function
        result = database.map(ids, from_type, to_types, taxon, return_query_info=True)

        # Return the result
        return jsonify(result)
    except Exception as e:
        # Handle exceptions and return an error message
        return jsonify({'error': str(e)}), 400
