from .alife_dataframe_to_biopython_tree \
    import alife_dataframe_to_biopython_tree
from .alife_dataframe_to_dendropy_tree \
    import alife_dataframe_to_dendropy_tree
from .biopython_tree_to_alife_dataframe \
    import biopython_tree_to_alife_dataframe
from .dendropy_tree_to_alife_dataframe \
    import dendropy_tree_to_alife_dataframe

import pandas
import dendropy
import Bio

# auxiliary tree class to allow for inter-format conversions
class ALifeTree():
    def __init__(self, tree):
        # internal tree representation is an alife-formatted tree
        self._tree = self._alife_dispatcher(tree)

    @property
    def biopython(self):
        return alife_dataframe_to_biopython_tree(self._tree, setup_edge_lengths=True)

    @property
    def dendropy(self):
        return alife_dataframe_to_dendropy_tree(self._tree, setup_edge_lengths=True)

    @property
    def alife(self):
        return self._tree

    def _alife_dispatcher(self, tree):
        """
        Convert any supported tree format to ALife format
        """
        if isinstance(tree, dendropy.Tree):
            # is a Dendropy Tree
            return dendropy_tree_to_alife_dataframe(tree) #, {'name': 'taxon_label'})
        elif isinstance(tree, pandas.DataFrame) and self._is_valid_alife_tree(tree):
            # is an Alife Dataframe
            return tree
        elif isinstance(tree, Bio.Phylo.BaseTree.Tree):
            # is a biopython tree
            return biopython_tree_to_alife_dataframe(tree, {'name': 'taxon_label'})
        else:
            raise ValueError("Unsupported tree format")

    def _is_valid_alife_tree(self, tree):
        return 'id' in tree and 'ancestor_list' in tree
