import Bio
import dendropy
import ete3
from functools import lru_cache
import networkx as nx
import pandas
import typing

from ._impl import _try_alifestd_validate as _try_alifestd_validate

from .alife_dataframe_to_biopython_tree \
    import alife_dataframe_to_biopython_tree
from .alife_dataframe_to_dendropy_tree \
    import alife_dataframe_to_dendropy_tree
from .alife_dataframe_to_ete_tree \
    import alife_dataframe_to_ete_tree
from .alife_dataframe_to_networkx_digraph \
    import alife_dataframe_to_networkx_digraph
from .biopython_tree_to_alife_dataframe \
    import biopython_tree_to_alife_dataframe
from .dendropy_tree_to_alife_dataframe \
    import dendropy_tree_to_alife_dataframe
from .ete_tree_to_alife_dataframe \
    import ete_tree_to_alife_dataframe
from .networkx_digraph_to_alife_dataframe \
    import networkx_digraph_to_alife_dataframe


class RosettaTree:
    """Adapter class for implicit conversion between tree representations
    across phylogenetic libraries.

    Enables user to translate across tree implementations without repeatedly
    calling the individual conversion functions. This allows users to accept
    any such tree as an argument in a library-agnostic way.
    """

    _tree: pandas.DataFrame

    def __init__(self, tree: typing.Union[
        dendropy.Tree,
        ete3.Tree,
        ete3.TreeNode,
        nx.DiGraph,
        pandas.DataFrame,
        Bio.Phylo.BaseTree.Tree
    ]) -> None:
        """Set up underlying data to be viewed, from any supported format."""
        # convert any supported tree format to ALife format,
        # as this is our interal representation
        if isinstance(tree, dendropy.Tree):
            # is a Dendropy Tree
            self._tree = dendropy_tree_to_alife_dataframe(tree)
        elif isinstance(tree, (ete3.Tree, ete3.TreeNode)):
            # is a Dendropy Tree
            self._tree = ete_tree_to_alife_dataframe(tree)
        elif isinstance(tree, Bio.Phylo.BaseTree.Tree):
            # is a biopython tree
            self._tree = biopython_tree_to_alife_dataframe(
                tree, {'name': 'taxon_label'}
            )
        elif isinstance(tree, nx.DiGraph):
            # is a networkx digraph
            self._tree = networkx_digraph_to_alife_dataframe(tree)
        elif (
            isinstance(tree, pandas.DataFrame)
            and len(tree)
            and _try_alifestd_validate(tree)
        ):
            # is an Alife Dataframe
            self._tree = tree
        else:
            raise ValueError("Unsupported tree format")

    @property
    @lru_cache(maxsize=None)
    def as_biopython(self: "RosettaTree") -> Bio.Phylo.BaseTree.Tree:
        """Return stored tree as a BioPython tree."""
        return alife_dataframe_to_biopython_tree(
            self._tree, setup_branch_lengths=True
        )

    @property
    @lru_cache(maxsize=None)
    def as_dendropy(self: "RosettaTree") -> dendropy.Tree:
        """Return stored tree as a DendroPy tree."""
        return alife_dataframe_to_dendropy_tree(
            self._tree, setup_edge_lengths=True
        )

    @property
    @lru_cache(maxsize=None)
    def as_ete(self: "RosettaTree") -> dendropy.Tree:
        """Return stored tree as an ete tree."""
        return alife_dataframe_to_ete_tree(
            self._tree, setup_dists=True
        )

    @property
    @lru_cache(maxsize=None)
    def as_networkx(self: "RosettaTree") -> nx.DiGraph:
        """Return stored tree as a NetworkX DiGraph tree."""
        return alife_dataframe_to_networkx_digraph(
            self._tree, setup_edge_lengths=True
        )

    @property
    @lru_cache(maxsize=None)
    def as_newick(self: "RosettaTree") -> str:
        """Return stored tree as a NetworkX DiGraph tree."""
        return self.as_dendropy.as_string(schema="newick")

    @property
    @lru_cache(maxsize=None)
    def as_alife(self: "RosettaTree") -> pandas.DataFrame:
        """Return stored tree as a dataframe in alife standard format."""
        return self._tree
