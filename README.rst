=====================================
alifedata-phyloinformatics-conversion
=====================================


.. image:: https://img.shields.io/pypi/v/alifedata-phyloinformatics-conversion.svg
        :target: https://pypi.python.org/pypi/alifedata-phyloinformatics-conversion

.. image:: https://img.shields.io/travis/mmore500/alifedata-phyloinformatics-conversion.svg
        :target: https://travis-ci.com/mmore500/alifedata-phyloinformatics-conversion

.. image:: https://readthedocs.org/projects/alifedata-phyloinformatics-conversion/badge/?version=latest
        :target: https://alifedata-phyloinformatics-conversion.readthedocs.io/en/latest/?badge=latest
        :alt: Documentation Status




alifedata-phyloinformatics-conversion helps apply traditional phyloinformatics software to alife standardized data


* Free software: MIT license
* Documentation: https://alifedata-phyloinformatics-conversion.readthedocs.io.


Built using the :code:`dendropy` library.

Installation

.. code-block:: bash

  pip3 install alifedata-phyloinformatics-converters


Use it as a Python module:

.. code-block:: python3

  import alifedata-phyloinformatics-conversion as apc
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
