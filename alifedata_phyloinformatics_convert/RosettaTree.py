import Bio
import dendropy
from functools import cached_property
import pandas

from .alife_dataframe_to_biopython_tree \
    import alife_dataframe_to_biopython_tree
from .alife_dataframe_to_dendropy_tree \
    import alife_dataframe_to_dendropy_tree
from .biopython_tree_to_alife_dataframe \
    import biopython_tree_to_alife_dataframe
from .dendropy_tree_to_alife_dataframe \
    import dendropy_tree_to_alife_dataframe


class RosettaTree:
    """Simple class for seamless interoperability between phylogenetic libraries.
    Enables user to translate across tree implementations without repeatedly calling
    the individual conversion functions.

    Attributes
    ----------
    tree:
        Tree to convert from. Can be any of the three supported by this library.
    """
    def __init__(self, tree):
        """Construct the library-agnostic tree.

        Parameters
        ----------
            tree:
                Tree to convert from. Can be any of the three supported by this library.
        """
        # convert any supported tree format to ALife format,
        # as this is our interal representation
        if isinstance(tree, dendropy.Tree):
            # is a Dendropy Tree
            self._tree = dendropy_tree_to_alife_dataframe(tree) #, {'name': 'taxon_label'})
        elif isinstance(tree, pandas.DataFrame) and self._is_valid_alife_tree(tree):
            # is an Alife Dataframe
            self._tree = tree
        elif isinstance(tree, Bio.Phylo.BaseTree.Tree):
            # is a biopython tree
            self._tree = biopython_tree_to_alife_dataframe(tree, {'name': 'taxon_label'})
        else:
            raise ValueError("Unsupported tree format")

    @cached_property
    def biopython(self):
        """Return stored tree as a biopython tree.
        """
        return alife_dataframe_to_biopython_tree(self._tree, setup_edge_lengths=True)

    @cached_property
    def dendropy(self):
        """Return stored tree as a dendropy tree.
        """
        return alife_dataframe_to_dendropy_tree(self._tree, setup_edge_lengths=True)

    @cached_property
    def alife(self):
        """Return stored tree as an alife-standardized phylogeny pandas dataframe.
        """
        return self._tree

    def _is_valid_alife_tree(self, tree):
        return 'id' in tree and 'ancestor_list' in tree
