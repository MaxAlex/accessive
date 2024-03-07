# Accessive

Accessive is a Python library designed to facilitate the conversion between various bioinformatic accession types, streamlining the process of working with biological data identifiers across different databases and formats.
Features

- Convert between a wide array of accession types.
- Support for major bioinformatics databases and identifier formats.
- Flexible query options with support for filtering by taxon.

## Supported Accession Types

Accessive supports the following accession types:

- Ensembl (Gene, mRNA, and Protein)
  - Gene (ENSG00000096717)
  - mRNA (ENST00000361390)
  - Protein (ENSP00000354689)
- Uniprot 
    - Swissprot (P00750)
    - TrEMBL (A0A024R161)
    - Isoform (P00750-1)
- RefSeq 
    - Gene (NM_001278)
    - mRNA (NM_001278.2)
    - ncRNA (NR_001278)
    - Peptide (NP_001265)
- Nextprot
    - Gene (NX_P00750)
    - Proteoform (NX_P00750-1)
- Alphafold id (AF-B0QZ35-F1)
- ArrayExpress (E-GEOD-3307)
- BioGRID index
- CCDS (CCDS81469)
- EMBL (AK074805)
- Entrez Gene id
- GeneCards index
- HGNC (HGNC:142929)
- MIM_GENE index (601739)
- PDB (4IG9)
- Pfam (PF02146)
- UCSC (uc057tpn.1)
- UniParc (UPI000015D95A)
- WikiGene index

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

