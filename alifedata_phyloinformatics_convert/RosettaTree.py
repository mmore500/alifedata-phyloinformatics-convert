import anytree
import Bio
import dendropy
import contextlib
from deprecated.sphinx import deprecated
from functools import lru_cache
import networkx as nx
import pandas
import pathlib
import treeswift
import typing
import typing_extensions
import validators
import warnings
import yarl

from ._impl import alifestd_validate as alifestd_validate
from ._impl import ete3
from ._impl import phytrack_Systematics

from .alife_dataframe_to_biopython_tree \
    import alife_dataframe_to_biopython_tree
from .alife_dataframe_to_dendropy_tree \
    import alife_dataframe_to_dendropy_tree
from .alife_dataframe_to_ete_tree \
    import alife_dataframe_to_ete_tree
from .alife_dataframe_to_treeswift_tree \
    import alife_dataframe_to_treeswift_tree
from .alife_dataframe_to_networkx_digraph \
    import alife_dataframe_to_networkx_digraph
from .alife_dataframe_to_phylotrack_systematics \
    import alife_dataframe_to_phylotrack_systematics
from .anytree_tree_to_alife_dataframe \
    import anytree_tree_to_alife_dataframe
from .biopython_tree_to_alife_dataframe \
    import biopython_tree_to_alife_dataframe
from .dendropy_tree_to_alife_dataframe \
    import dendropy_tree_to_alife_dataframe
from .ete_tree_to_alife_dataframe \
    import ete_tree_to_alife_dataframe
from .networkx_digraph_to_alife_dataframe \
    import networkx_digraph_to_alife_dataframe
from .phylotrack_systematics_to_alife_dataframe \
    import phylotrack_systematics_to_alife_dataframe
from .treeswift_tree_to_alife_dataframe \
    import treeswift_tree_to_alife_dataframe


class RosettaTree:
    """Adapter class for implicit conversion between tree representations
    across phylogenetic libraries.

    Enables user to translate across tree implementations without repeatedly
    calling the individual conversion functions. This allows users to accept
    any such tree as an argument in a library-agnostic way.
    """

    _tree: pandas.DataFrame

    def __init__(
        self,
        tree: typing.Union[
            anytree.NodeMixin,
            dendropy.Tree,
            ete3.Tree,
            ete3.TreeNode,
            nx.DiGraph,
            pandas.DataFrame,
            Bio.Phylo.BaseTree.Tree,
            phytrack_Systematics,
            treeswift.Tree
        ],
        validate: typing.Literal["warn", "error", "ignore"] = "warn",
    ) -> None:
        """Load phylogeny from any supported data structure.

        Supported data structures include:
        - `anytree.NodeMixin`
        - `Bio.Phylo.BaseTree.Tree` (biopython)
        - `dendropy.Tree`
        - `ete3.Tree`
        - `networkx.Digraph`
        - `pandas.DataFrame` (alife standard format)
        - `phylotrackpy.systematics.Systematics`
        """
        # convert any supported tree format to ALife format,
        # as this is our interal representation
        if isinstance(tree, anytree.node.NodeMixin):
            # is an AnyTree tree
            self._tree = anytree_tree_to_alife_dataframe(tree)
        elif isinstance(tree, dendropy.Tree):
            # is a Dendropy Tree
            self._tree = dendropy_tree_to_alife_dataframe(tree)
        elif isinstance(tree, (ete3.Tree, ete3.TreeNode)):
            # is a ete Tree
            self._tree = ete_tree_to_alife_dataframe(tree)
        elif isinstance(tree, Bio.Phylo.BaseTree.Tree):
            # is a biopython tree
            self._tree = biopython_tree_to_alife_dataframe(
                tree, {'name': 'taxon_label'}
            )
        elif isinstance(tree, nx.DiGraph):
            # is a networkx digraph
            self._tree = networkx_digraph_to_alife_dataframe(tree)
        elif isinstance(tree, phytrack_Systematics):
            # is a phylotrack Systematics object
            self._tree = phylotrack_systematics_to_alife_dataframe(tree)
        elif isinstance(tree, treeswift.Tree):
            # is a phylotrack Systematics object
            self._tree = treeswift_tree_to_alife_dataframe(tree)
        elif (
            isinstance(tree, pandas.DataFrame)
            and "id" in tree.columns
            # i.e., ancestor_id or ancestor_list
            and tree.columns.str.startswith("ancestor_").any()
        ):
            if validate == "ignore":
                pass
            elif not alifestd_validate(tree):
                if validate == "error":
                    raise ValueError(
                        "Tree does not comply with alife data standards.",
                    )
                elif validate == "warn":
                    warnings.warn(
                        "Tree does not comply with alife data standards.",
                    )
            self._tree = tree
        else:
            raise ValueError(
                f"Unsupported tree format tree={tree} of type {type(tree)}",
            )

    @property
    @lru_cache(maxsize=None)
    def as_biopython(
        self: "RosettaTree",
    ) -> typing.Optional[Bio.Phylo.BaseTree.Tree]:
        """Return stored tree as a BioPython tree."""
        return alife_dataframe_to_biopython_tree(
            self._tree, setup_branch_lengths=True
        )

    @property
    @lru_cache(maxsize=None)
    def as_dendropy(self: "RosettaTree") -> typing.Optional[dendropy.Tree]:
        """Return stored tree as a DendroPy tree."""
        return alife_dataframe_to_dendropy_tree(
            self._tree, setup_edge_lengths=True
        )

    @property
    @lru_cache(maxsize=None)
    def as_ete(self: "RosettaTree") -> typing.Optional[ete3.Tree]:
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
    def as_phylotrack(self: "RosettaTree") -> phytrack_Systematics:
        """Return stored tree as a phylotrack Systematics object."""
        return alife_dataframe_to_phylotrack_systematics(self._tree)

    @property
    @lru_cache(maxsize=None)
    def as_treeswift(self: "RosettaTree") -> phytrack_Systematics:
        """Return stored tree as a treeswift object."""
        return alife_dataframe_to_treeswift_tree(
            self._tree, setup_edge_lengths=True
        )

    @property
    @deprecated(version="0.15.0", reason="Use to_newick instead.")
    @lru_cache(maxsize=None)
    def as_newick(self: "RosettaTree") -> str:
        """Return stored tree as a Newick string."""
        return self.to_newick()

    @property
    @lru_cache(maxsize=None)
    def as_alife(self: "RosettaTree") -> pandas.DataFrame:
        """Return stored tree as a dataframe in alife standard format."""
        return self._tree

    def to_schema(
        self: "RosettaTree",
        schema: typing_extensions.Literal[
            "newick",
            "nexus",
            "nexml",
        ],
        file: typing.Union[None, str, pathlib.Path, typing.IO] = None,
    ) -> typing.Optional[str]:
        """Serialize the stored tree to `schema` format."""
        if len(self._tree) == 0:
            if file is None:
                return None
            else:
                raise ValueError(
                    f"Schema {schema} cannot represent an empty tree.",
                )
        try:
            if file is None:
                return self.as_dendropy.as_string(schema=schema)
            elif isinstance(file, (str, pathlib.Path)):
                self.as_dendropy.write_to_path(dest=file, schema=schema)
            else:
                self.as_dendropy.write_to_stream(dest=file, schema=schema)
        except Exception as e:
            raise ValueError(
                f"Exception '{e}' ocurred. If provided, argument file={file} "
                "must be file path or stream handle.",
            )

    def to_newick(
        self: "RosettaTree",
        file: typing.Union[None, str, pathlib.Path, typing.IO] = None,
    ) -> typing.Optional[str]:
        """Convert the stored tree to Newick format."""
        return self.to_schema(schema="newick", file=file)

    def to_nexus(
        self: "RosettaTree",
        file: typing.Union[None, str, pathlib.Path, typing.IO] = None,
    ) -> typing.Optional[str]:
        """Convert the stored tree to Nexus format."""
        return self.to_schema(schema="nexus", file=file)

    def to_nexml(
        self: "RosettaTree",
        file: typing.Union[None, str, pathlib.Path, typing.IO] = None,
    ) -> typing.Optional[str]:
        """Convert the stored tree to Nexml format."""
        return self.to_schema(schema="nexml", file=file)

    @classmethod
    def from_schema(
        self: "RosettaTree",
        schema: typing_extensions.Literal[
            "newick",
            "nexus",
            "nexml",
        ],
        source: typing.Union[None, str, pathlib.Path, yarl.URL, typing.IO],
    ) -> "RosettaTree":
        """Serialize the stored tree to `schema` format."""

        def safe_is_file() -> bool:
            with contextlib.suppress(Exception):
                return pathlib.Path(source).is_file()
            return False

        if isinstance(source, pathlib.Path):
            tree = dendropy.Tree.get(path=str(source), schema=schema)
        elif isinstance(source, yarl.URL):
            tree = dendropy.Tree.get(url=str(source), schema=schema)
        elif isinstance(source, str) and safe_is_file():
            warnings.warn(
                f"String source={source} is ambiguous, interpreting as path. "
                "Pass argument as pathlib.Path object to suppress warning."
            )
            tree = dendropy.Tree.get(path=source, schema=schema)
        elif isinstance(source, str) and validators.url(source):
            warnings.warn(
                f"String source={source} is ambiguous, interpreting as url. "
                "Pass argument as yarl.URL object to suppress warning."
            )
            tree = dendropy.Tree.get(url=source, schema=schema)
        elif isinstance(source, str):
            tree = dendropy.Tree.get(data=source, schema=schema)
        else:
            tree = dendropy.Tree.get(file=source, schema=schema)

        return RosettaTree(tree)

    @classmethod
    def from_newick(
        cls: typing.Type,
        source: typing.Union[None, str, pathlib.Path, yarl.URL, typing.IO],
    ) -> "RosettaTree":
        """Open data in Newick format."""
        return cls.from_schema(schema="newick", source=source)

    @classmethod
    def from_nexus(
        cls: typing.Type,
        source: typing.Union[None, str, pathlib.Path, yarl.URL, typing.IO],
    ) -> "RosettaTree":
        """Open data in Nexus format."""
        return cls.from_schema(schema="nexus", source=source)

    @classmethod
    def from_nexml(
        cls: typing.Type,
        source: typing.Union[None, str, pathlib.Path, yarl.URL, typing.IO],
    ) -> "RosettaTree":
        """Open data in Nexml format."""
        return cls.from_schema(schema="nexml", source=source)
