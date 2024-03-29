# Accessive

Accessive is a Python library designed to facilitate the conversion between various bioinformatic accession types, streamlining the process of working with biological data identifiers across different databases and formats.
Features

- Convert between a wide array of accession types.
- Support for major bioinformatics databases and identifier formats.
- Flexible query options with support for filtering by taxon.

## Supported Accession Types

Accessive supports the following accession types:

- Ensembl
    - Gene (ENSG00000096717)
    - mRNA (ENST00000361390)
    - Proteoform (ENSP00000354689)
- Uniprot 
    - Swissprot (P00750)
    - TrEMBL (A0A024R161)
    - Isoform (P00750-1)
- RefSeq 
    - mRNA (NM_001278)
    - ncRNA (NR_001278)
    - Proteoform (NP_001265)
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
After installation, install the database via the included utiltiy:

```bash
python -m accessive.database_ops --download
```

Note that the Accessive database requires approximately 500MB of disk space. 

## Usage Example

```
>>> from accessive import Accessive
>>> acc = Accessive()
>>> acc.map(ids=['SIRT1', 'UCP1', 'VGF'], to_types=['uniprot_swissprot', 'uniprot_trembl', 'refseq_mrna'], taxon = 9606)

          uniprot_swissprot uniprot_trembl   refseq_mrna
gene_name                                               
SIRT1                  None         B0QZ35          None
SIRT1                  None         E9PC49  NM_001142498
SIRT1                Q96EB6     A0A024QZQ1     NM_012238
SIRT1                  None         B0QZ35  NM_001314049
VGF                  O15240           None          None
VGF                  O15240           None     NM_003378
UCP1                 P25874           None     NM_021833

```

## Documentation

For more information, please refer to the [documentation](https://accessive.readthedocs.io/en/latest/).
