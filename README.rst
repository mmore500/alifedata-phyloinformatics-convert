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

.. image:: https://zenodo.org/badge/466241441.svg
  :target: https://zenodo.org/doi/10.5281/zenodo.10701178


alifedata-phyloinformatics-convert helps apply traditional phyloinformatics software to alife standardized data


* Free software: MIT license
* Documentation: https://alifedata-phyloinformatics-convert.readthedocs.io.

Usage
-----

Use :code:`apc`'s :code:`RosettaTree` interface for flexible conversion between phylogenetic data structures and schemas.
First, create a :code:`RosettaTree` object from any supported structure/schema

.. code-block:: python3

  import io
  import pathlib

  import alifedata_phyloinformatics_convert as apc
  import anytree
  import Bio
  import dendropy
  import ete3 as ete
  import networkx
  import pandas
  import phylotrackpy
  import treeswift

  newickstr = "((A,B),(C,D));"

  for obj in [
    anytree.AnyNode(),
    Bio.Phylo.read(io.StringIO(newickstr), "newick"),
    dendropy.Tree.get(data=newickstr, schema="newick"),
    ete.Tree(newickstr),
    treeswift.Tree.read_tree_newick(newickstr),
    networkx.DiGraph(),
    pandas.DataFrame({"id": [0], "ancestor_list": "[None]"}),  # alife standard
    phylotrackpy.systematics.Systematics(lambda x: x),
  ]:
    converter = apc.RosettaTree(obj)

  # from phyloinformatics schema
  # ... nexml and nexus also supported!
  converter = apc.RosettaTree.from_newick(newickstr)
  converter = apc.RosettaTree.from_newick(pathlib.Path("read.newick"))
  with open("read.newick", "r") as fp:
    converter = apc.RosettaTree.from_newick(fp)

  # from alife standard data via Pandas
  converter = apc.RosettaTree(pandas.read_csv("read-alifestd.csv"))

Then, convert or serialize data

.. code-block:: python3

  # ... converter created as above
  converter.as_alife  # pandas DataFrame
  converter.as_biopython
  converter.as_dendropy
  converter.as_ete
  converter.as_networkx
  converter.as_phylotrack
  converter.as_treeswift

  # serialization, nexml and nexus schemata also supported
  converter.to_newick()  # returns newick string
  converter.to_newick(pathlib.Path("write.newick"))  # writes to path
  with open("write.newick", "w") as fp:  # writes to file object
    converter.to_newick(fp)

  # alifestd serialization
  converter.as_alife.to_csv("write-alifestd.csv", index=False)

Use :code:`apc`'s functional interface to convert between alife format other libraries' tree objects

.. code-block:: python3

  import alifedata_phyloinformatics_convert as apc
  import pandas

  alife_df = pandas.read_csv('alifedata.csv')

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

  # treeswift
  treeswift_tree = apc.alife_dataframe_to_treeswift_tree(alife_df)
  frame = apc.treeswift_tree_to_alife_dataframe(tree)

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

Citing
------

If alifedata-phyloinformatics-convert is used in scientific publication, please cite it as

    Matthew Andres Moreno and Santiago Rodriguez Papa. (2024). mmore500/alifedata-phyloinformatics-convert. Zenodo. https://doi.org/10.5281/zenodo.10701178

.. code:: bibtex

    @software{moreno2024apc,
      author = {Matthew Andres Moreno AND Santiago {Rodriguez Papa}},
      title = {mmore500/alifedata-phyloinformatics-convert},
      month = feb,
      year = 2024,
      publisher = {Zenodo},
      doi = {10.5281/zenodo.10701178},
      url = {https://doi.org/10.5281/zenodo.10701178}
    }

And don't forget to leave a `star on GitHub <https://github.com/mmore500/alifedata-phyloinformatics-convert/stargazers>`__!

Credits
-------

Built using the `DendroPy`_ library.
This package was created with Cookiecutter_ and the `audreyr/cookiecutter-pypackage`_ project template.

.. _DendroPy: https://github.com/jeetsukuruman/dendropy
.. _Cookiecutter: https://github.com/audreyr/cookiecutter
.. _`audreyr/cookiecutter-pypackage`: https://github.com/audreyr/cookiecutter-pypackage
