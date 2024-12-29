from iterpop import iterpop as ip
from nanto import isanan, nantonone, nantozero
import pandas as pd
import string
import typing

from ._impl import ete3
from ._impl import keydefaultdict
from ._impl import parse_ancestor_list as _parse_ancestor_list


def alife_dataframe_to_ete_trees(
    df: pd.DataFrame,
    setattrs: typing.Optional[typing.Union[
        typing.Iterable[str],
        typing.Mapping[str, str],
    ]] = None,
    setup_dists: bool = False,
) -> typing.List[ete3.TreeNode]:
    """Open a phylogeny dataframe formatted to the artificial life community
    data format standards as zero or more ete trees, depending on the number
    of clades with no common ancestor.

    The following columns will automatically be applied as attributes to
    generated TreeNode objects:
        * dist,
        * id,
        * origin_time, and
        * name.

    Parameters
    ----------
    df:
        Pandas DataFrame to convert.
    setattrs: optional
        Dataframe columns that should be attached as attributes to TreeNode
        objects within the trees. If a map is provided, values at columns in
        keys will be attached with the corresponding value as the attr name.
    setup_dists: bool, optional
        Should we try to set up edge lengths using the origin_time column?
        Will not override if dist is provided as a column of df.
    """

    df = df.copy()

    df['parsed_ancestor_list'] \
        = df['ancestor_list'].apply(_parse_ancestor_list)

    # maps id to node
    def setup_node(id: int) -> ete3.TreeNode:
        res = ete3.TreeNode()
        res.add_features(id=id)
        return res
    nodes = keydefaultdict(setup_node)

    root_nodes = []

    for __, row in df.iterrows():
        node = nodes[row['id']]

        if setattrs is not None:
            for col_name in setattrs:
                attr_name: str
                try:
                    attr_name = setattrs[col_name]
                except TypeError:
                    attr_name = col_name
                assert attr_name not in (
                    'dist',
                    'origin_time',
                    'name',
                )
                assert not hasattr(node, attr_name)
                node.add_features(**{attr_name: row[col_name]})

        node.add_features(
            origin_time=nantonone(row.get('origin_time', default=None)),
        )
        name = row.get('name', default=None)
        if name not in ('None', None):
            node.name = name
        node.dist = nantozero(row.get('dist', 1.0))

        if row['parsed_ancestor_list']:
            ancestor_id = ip.popsingleton(row['parsed_ancestor_list'])
            nodes[ancestor_id].add_child(node)
        else:
            root_nodes.append(node)

    # set up edge lengths
    if setup_dists:
        for id, parent_node in nodes.items():
            for child_node in parent_node.children:
                if child_node.dist == 1.0 and None not in (
                    getattr(child_node, 'origin_time', None),
                    getattr(parent_node, 'origin_time', None),
                ):
                    assert child_node.origin_time >= parent_node.origin_time
                    child_node.dist \
                        = child_node.origin_time - parent_node.origin_time

        for root_node in root_nodes:
            if (
                root_node.dist == 1.0
                and getattr(root_node, 'origin_time', None) is not None
            ):
                assert not isanan(root_node.origin_time)
                root_node.dist = root_node.origin_time

    return root_nodes  # Tree is an alias of TreeNode
