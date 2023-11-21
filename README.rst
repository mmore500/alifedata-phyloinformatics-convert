==================================
alifedata-phyloinformatics-convert
==================================


.. image:: https://img.shields.io/pypi/v/alifedata-phyloinformatics-convert.svg
        :target: https://pypi.python.org/pypi/alifedata-phyloinformatics-convert
        :alt: PyPI Status

.. image:: https://github.com/mmore500/alifedata-phyloinformatics-convert/actions/workflows/CI.yml/badge.svg
        :target: https://github.com/mmore500/alifedata-phyloinformatics-convert/actions/workflows/CI.yml
        :alt: CI Status

.. image:: https://readthedocs.org/projects/alifedata-phyloinformatics-convert/badge/?version=latest
        :target: https://alifedata-phyloinformatics-convert.readthedocs.io/en/latest/?badge=latest
        :alt: Documentation Status




alifedata-phyloinformatics-convert helps apply traditional phyloinformatics software to alife standardized data


* Free software: MIT license
* Documentation: https://alifedata-phyloinformatics-convert.readthedocs.io.

Usage
----

Use :code:`apc`'s functional interface to convert between alife format other libraries' tree objects.

.. code-block:: python3

  import alifedata_phyloinformatics_convert as apc

  alife_df = pd.read_csv('alifedata.csv')

  # biopython
  tree = apc.alife_dataframe_tobiopython_tree(alife_df)
  frame = apc.biopython_tree_to_alife_dataframe(tree)

  # dendropy
  tree = apc.alife_dataframe_to_dendropy_tree(alife_df)
  frame = apc.dendropy_tree_to_alife_dataframe(tree)

  # ete
  ete_tree = apc.alife_dataframe_to_ete_tree(alife_df)
  frame = apc.ete_tree_to_alife_dataframe(tree)

  # networkx
  digraph = apc.alife_dataframe_to_networkx_digraph(alife_df)
  frame = apc.networkx_digraph_to_alife_dataframe(digraph)

  # phylotrackpy
  systematics = apc.alife_dataframe_to_phylotrack_systematics(alife_df)
  frame = apc.phylotrack_systematics_to_alife_dataframe(systematics)

  # partial support is also included for,
  # - adjacency lists
  # - anytree trees
  # - scipy linkage matrices
  # ... see API documentation for details

Command Line Interface
----------------------

Use :code:`apc`'s CLI :code:`toalifedata` command to convert newick, nexml, and nexus data to alife standard phylogenetics data

.. code-block:: bash

  Usage: alifedata-phyloinformatics-convert toalifedata [OPTIONS]

    convert standard alife phylogeny data to phloinformatics format

  Options:
    --input-file FILENAME           phyloinformatics data file path; default
                                    stdin
    --input-schema TEXT             phyloinformatics data format schema; options
                                    include newick, nexml, and nexus  [required]
    --output-file FILENAME          alife data file path; default stdout
    --output-format TEXT            alife data file format; default csv
    --suppress-unifurcations / --keep-unifurcations
                                    Compress sequences of nodes with single
                                    descendants
    --help                          Show this message and exit.



Use the :code:`fromalifedata` command to convert to other formats from alife standard phylogenetics data

.. code-block:: bash

  Usage: alifedata-phyloinformatics-convert fromalifedata [OPTIONS]

    convert phloinformatics data to standard alife phylogeny format

  Options:
    --input-file FILENAME           alife data file path; default stdin
    --input-format TEXT             alife data file format; default csv
    --output-file FILENAME          phyloinformatics data file path; default
                                    stdout
    --output-schema TEXT            phyloinformatics data format schema; options
                                    include newick, nexml, and nexus  [required]
    --suppress-unifurcations / --keep-unifurcations
                                    Compress sequences of nodes with single
                                    descendants
    --help                          Show this message and exit.

Installation
------------

Install from PyPi

.. code-block:: bash

  pip3 install alifedata-phyloinformatics-convert


Credits
-------

Built using the `DendroPy`_ library.
This package was created with Cookiecutter_ and the `audreyr/cookiecutter-pypackage`_ project template.

.. _DendroPy: https://github.com/jeetsukuruman/dendropy
.. _Cookiecutter: https://github.com/audreyr/cookiecutter
.. _`audreyr/cookiecutter-pypackage`: https://github.com/audreyr/cookiecutter-pypackage
