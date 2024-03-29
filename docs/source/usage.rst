Usage
=====

.. _installation:

Installation
------------

To use Accessive, first install it using pip:

.. code-block:: console

   $ pip install accessive 

Installing the database
----------------

Accessive uses a local SQLite database to store accession mapping information. To install the database, run:

.. code-block:: console

    $ python -m accessive.database_ops --download

The database will download automatically and immediately be usable by Accessive. Note that the Accessive database 
is about 500MB in size, make sure you have sufficient disk space beforehand.


.. _accessions:

Accessions
----------

Accessive supports the following accession types, denoted by the indicated keywords:

====================== ==================== ===============
Format                 Keyword              Example Value
====================== ==================== ===============
Ensembl Gene           ensembl_gene         ENSG00000096717
Ensembl mRNA           ensembl_mrna         ENST00000361390
Ensembl Proteoform     ensembl_prot         ENSP00000354689
Uniprot Swissprot      uniprot_swissprot    P00750
Uniprot TrEMBL         uniprot_trembl       A0A024R161
Uniprot Isoform        uniprot_isoform      P00750-1
RefSeq mRNA            refseq_mrna          NM_001278
RefSeq ncRNA           refseq_ncrna         NR_001278
RefSeq Proteoform      refseq_peptide       NP_001265
Nextprot Gene          nexprot              NX_P00750
Nextprot Proteoform    nextprot_isoform     NX_P00750-1
Alphafold id           alphafold            AF-B0QZ35-F1
ArrayExpress           arrayexpress         E-GEOD-3307
BioGRID index          biogrid              113010
CCDS                   ccds                 CCDS81469
EMBL                   embl                 AK074805
Entrez Gene id         entrez_gene          100009676
GeneCards index        genecards            GC01M000000
HGNC                   hgnc                 HGNC:142929
MIM_GENE index         mim_gene             601739
PDB                    pdb                  4IG9
Pfam                   pfam                 PF02146
UCSC                   ucsc                 uc057tpn.1
UniParc                uniparc              UPI000015D95A
WikiGene index         wikigene             100287
====================== ==================== ===============



In calls to accessive.map() or accessive.get(), the accession type should be specified as a keyword argument. For
example, to get the PDB and alphafold entry accessions for a given gene by it's Ensembl Gene ID, use:

.. code-block:: python

   import accessive
   accessive.map('ENSG00000096717', from_type='ensembl_gene', to_types=['pdb', 'alphafold'])


