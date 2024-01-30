import os
from flask import Flask, request, jsonify, render_template
from interface import Accessive


app = Flask(__name__)
database = Accessive()

RESULT_DIR = os.path.join(os.path.dirname(__file__), 'map_result_cache')
if not os.path.exists(RESULT_DIR):
    os.makedirs(RESULT_DIR)


@app.route('/hello', methods=['GET'])
def hello():
    return "This whole thing is sort of accessive, isn't it?"


@app.route('/map_identifiers', methods=['GET', 'POST'])
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
        format = request.args.get('format', 'json')

        if not ids:
            raise Exception('No identifiers provided.')
        if not taxon:
            raise Exception('No taxon provided; this is required. (NB: Human is 9606.)')
    
        # Call the map function
        result = database.map(ids, from_type, to_types, taxon, return_query_info=True, return_format=format)

        # Return the result
        return jsonify(result)
    except Exception as e:
        # Handle exceptions and return an error message
        return jsonify({'error': str(e)}), 400


@app.route('/map_identifiers_page', methods=['POST'])
def map_identifiers_page():
    try:
        ids = request.args.getlist('ids')  # Expects a list of identifiers
        from_type = request.args.get('from_type')
        to_types = request.args.getlist('to_types') if 'to_types' in request.args else None
        taxon = request.args.get('taxon')

        if not ids:
            return jsonify({'error': 'No identifiers provided.'}), 400
        if not taxon:
            return jsonify({'error': 'No taxon provided; this is required. (NB: Human is 9606.)'}), 400
    
        # Call the map function
        result : dict = database.map(ids, from_type, to_types, taxon, return_query_info=True, return_format='pandas') # type: ignore

        # Return the result
        result_page = render_template('result_page.html', 
                                      table_html=result['result'].to_html(classes='dataframe'), # type: ignore
                                      from_type=result['from_type'], 
                                      to_types=result['to_types'],
                                      taxon=result['taxon']
                                   ) 
        result_id = hash(result['result'].to_html()) # type: ignore
        result_file = os.path.join(RESULT_DIR, str(result_id))
        with open(result_file, 'w') as f:
            f.write(result_page)

        return jsonify({'result_id': result_id})
    except Exception as e:
        # Handle exceptions and return an error message
        return jsonify({'error': str(e)}), 400

@app.route('/result_page/<result_id>', methods=['GET', 'POST'])
def result_page(result_id):
    result_file = os.path.join(RESULT_DIR, result_id)
    if os.path.exists(result_file):
        with open(result_file, 'r') as f:
            return f.read()
    else:
        return jsonify({'error': 'No result found with that ID.'}), 400


