# Accessive

Accessive is a Python library designed to facilitate the conversion between various bioinformatic accession types, streamlining the process of working with biological data identifiers across different databases and formats.
Features

    Convert between a wide array of accession types.
    Support for major bioinformatics databases and identifier formats.
    Flexible query options with support for filtering by taxon.

## Supported Accession Types

Accessive supports the following accession types:

    Ensembl (Gene, mRNA, and Protein)
    Uniprot (SwissProt, TrEMBL, and isoforms)
    RefSeq (Gene, mRNA, ncRNA, and peptide)
    Nextprot (Gene and proteoform)
    Alphafold
    ArrayExpress
    BioGRID
    CCDS
    EMBL
    ENS_LRG_gene
    ENS_LRG_transcript
    EntrezGene
    GeneCards
    HGNC
    MIM_GENE
    PDB
    Pfam
    UCSC
    UniParc
    WikiGene

## Installation

To install Accessive, use pip:

```bash
pip install accessive
```

## Usage Example

```python

from accessive import Accessive

accessive = Accessive()

result = accessive.map(ids=['ENSG00000139618'], from_type='ensembl_gene', to_types=['uniprot_swissprot', 'refseq_peptide'])
```

