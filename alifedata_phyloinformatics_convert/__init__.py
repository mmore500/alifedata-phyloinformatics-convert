"""Top-level package for alifedata-phyloinformatics-convert."""

__author__ = """Matthew Andres Moreno"""
__email__ = 'm.more500@gmail.com'
__version__ = '0.1.2'

from .alife_dataframe_to_dendropy_tree import alife_dataframe_to_dendropy_tree
from .alife_dataframe_to_dendropy_trees \
    import alife_dataframe_to_dendropy_trees
from .dendropy_tree_to_alife_dataframe import dendropy_tree_to_alife_dataframe

# adapted from https://stackoverflow.com/a/31079085
__all__ = [
    'alife_dataframe_to_dendropy_tree',
    'alife_dataframe_to_dendropy_trees',
    'dendropy_tree_to_alife_dataframe',
]
