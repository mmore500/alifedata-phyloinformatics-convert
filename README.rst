=====================================
alifedata-phyloinformatics-convert
=====================================


.. image:: https://img.shields.io/pypi/v/alifedata-phyloinformatics-convert.svg
        :target: https://pypi.python.org/pypi/alifedata-phyloinformatics-convert

.. image:: https://img.shields.io/travis/mmore500/alifedata-phyloinformatics-convert.svg
        :target: https://travis-ci.com/mmore500/alifedata-phyloinformatics-convert

.. image:: https://readthedocs.org/projects/alifedata-phyloinformatics-convert/badge/?version=latest
        :target: https://alifedata-phyloinformatics-convert.readthedocs.io/en/latest/?badge=latest
        :alt: Documentation Status




alifedata-phyloinformatics-convert helps apply traditional phyloinformatics software to alife standardized data


* Free software: MIT license
* Documentation: https://alifedata-phyloinformatics-convert.readthedocs.io.


Built using the :code:`dendropy` library.

Installation

.. code-block:: bash

  pip3 install alifedata-phyloinformatics-converters


Use it as a Python module:

.. code-block:: python3

  import alifedata-phyloinformatics-convert as apc
  import dendropy
  import pandas as pd

  alife_df = pd.read_csv('alifedata.csv')

  # get a dendropy Tree from alife-standardized phylogeny dataframe
  tree = apc.alife_dataframe_to_dendropy_tree(converted_df)

  # get a alife-standardized phylogeny dataframe from a dendropy Tree
  reconverted_alife_df = apc.dendropy_tree_to_alife_dataframe(tree)



Credits
-------

This package was created with Cookiecutter_ and the `audreyr/cookiecutter-pypackage`_ project template.

.. _Cookiecutter: https://github.com/audreyr/cookiecutter
.. _`audreyr/cookiecutter-pypackage`: https://github.com/audreyr/cookiecutter-pypackage
