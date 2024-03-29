API
===

.. autosummary::
   :toctree: generated


accessive.Accessive() 
---------------------

.. autoclass:: Accessive
   :members:
   :undoc-members:
   :show-inheritance:

The primary functionality of the Accessive API is provided by the Accessive class. This class
provides methods for mapping between accession formats.


Constructor
^^^^^^^^^^^

The ``Accessive`` class can be set to use default values for later calls to map(); this can be used, for instance, to
initialize an object instance that always queries mouse genes, or which maps between Ensembl and UniProt identifiers,
without specifying the respective arguments on each separate call to map(). 

.. code-block:: python

    from accessive import Accessive
    acc = Accessive(sqlite_file=None,
                    default_from_type=None,
                    default_to_types=None,
                    default_format='pandas',
                    default_taxon=None,
                    default_require_canonical=False)

Parameters:

- ``sqlite_file``: Use a specified SQLite database; defaults to the database installed by database_ops.
- ``default_from_type``: Sets the default source identifier type. If not specified, the source type must be provided in each call to the ``map`` or ``get`` methods.
- ``default_to_types``: Defines a list of default target identifier types for mapping operations. This default can be overridden by specifying ``to_types`` in the ``map`` method.
- ``default_format``: Determines the default format for query results. Supported formats include 'pandas' (returns a Pandas DataFrame), 'json', and 'txt'. The default value is 'pandas'.
- ``default_taxon``: Sets the default taxonomic identifier to narrow down queries. If not specified, the taxon must be provided in each call to the ``map`` or ``get`` methods.
- ``default_require_canonical``: When ``True``, only canonical or 'recommended' identifiers will be returned, which helps avoid less-common gene names or outdated identifier versions.


map() 
^^^^^^^^^^

The ``map`` method allows for the conversion of a set of biological identifiers from one type to another.

.. code-block:: python

    result = acc.map(ids=['ENSG00000139618'],
                     from_type='ensembl_gene',
                     to_types=['uniprot_swissprot', 'refseq_peptide'],
                     taxon=9606)

Parameters:

- ``ids``: A list of identifiers to be converted.
- ``from_type``: The type of the input identifiers. See :ref:`the usage page <accessions>` for a list of supported types.
- ``to_types``: A list of types to convert the identifiers to. :ref:`the usage page <accessions>` for a list of supported types.
- ``taxon``: The taxonomic species identifier (optional).

The method returns a table or dict structure containing the requested identifiers.

get() 
^^^^^^^^^^

The ``get`` method is a convenience method for converting a single identifier from one type to another.

.. code-block:: python

    result = acc.get(accession='ENSG00000139618',
                     from_type='ensembl_gene',
                     to_type='uniprot_swissprot',
                     taxon=9606)

Parameters:

- ``accession``: The accession identifier to be converted.
- ``from_type``: The type of the input identifier.
- ``to_type``: The target identifier type to convert to.
- ``taxon``: The taxonomic species identifier (optional).

The method returns the requested identifier.


identify() 
^^^^^^^^^^^^^^^^

The ``identify`` method is used to identify the type of a given identifier.

.. code-block:: python

    result = acc.identify(accession='ENSG00000139618')

Parameters:

- ``accession``: The identifier to be identified.

The method returns a list of potential types of the provided identifier.


available_taxons() 
^^^^^^^^^^^^^^^^^^^^^^^^^

The ``available_taxons`` method returns a list of available species in the current database and their taxon numbers.

.. code-block:: python

    result = acc.available_taxons()

Parameters:

- None


The method returns a list of species name and taxon numbers.
