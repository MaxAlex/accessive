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

