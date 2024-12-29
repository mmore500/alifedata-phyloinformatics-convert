import treeswift
from iterpop import iterpop as ip
from nanto import isanan, nantonone
import opytional as opyt
import pandas as pd
import typing

from ._impl import keydefaultdict
from ._impl import parse_ancestor_list as _parse_ancestor_list

def _treeswift_Tree_with_root(root: treeswift.Node) -> treeswift.Tree:
    res = treeswift.Tree(is_rooted=True)
    res.is_rooted = True
    res.root = root
    return res

def alife_dataframe_to_treeswift_trees(
    df: pd.DataFrame,
    setattrs: typing.Optional[typing.Union[
        typing.Iterable[str],
        typing.Mapping[str, str],
    ]] = None,
    setup_edge_lengths: bool = False,
) -> typing.List[treeswift.Tree]:
    """Open a phylogeny dataframe formatted to the artificial life community
    data format standards as zero or more treeswift trees, depending on the
    number of clades with no common ancestor.

    The following columns will automatically be applied as attributes to
    generated Node objects:
        * edge_length,
        * id,
        * label,
        * origin_time, and
        * taxon_label.

    Parameters
    ----------
    df:
        Pandas DataFrame to convert.
    setattrs: optional
        Dataframe columns that should be attached as attributes to Node
        objects within the trees. If a map is provided, values at columns in
        keys will be attached with the corresponding value as the attr name.
    setup_edge_lengths: bool, optional
        Should we try to set up edge lengths using the origin_time column?
        Will not override if edge_length is provided as a column of df.
    """

    df = df.copy()

    df['parsed_ancestor_list'] \
        = df['ancestor_list'].apply(_parse_ancestor_list)

    # maps id to node
    def setup_node(id: int) -> treeswift.Node:
        res = treeswift.Node()
        res.id = id
        return res
    nodes = keydefaultdict(setup_node)

    # maps id to origin time
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
                    'edge_length',
                    'origin_time',
                    'taxon',
                    'taxon_label',
                )
                assert not hasattr(node, attr_name)
                setattr(node, attr_name, row[col_name])

        node.origin_time = nantonone(row.get('origin_time', default=None))
        node.label = row.get('label', default=None)
        taxon_label = row.get('taxon_label', default=None)
        if taxon_label not in ('None', None):
            node.taxon = treeswift.Taxon(label=taxon_label)
        opyt.apply_if(
            opyt.apply_if(nantonone(row.get('edge_length', None)), float),
            node.set_edge_length,
        )

        if row['parsed_ancestor_list']:
            ancestor_id = ip.popsingleton(row['parsed_ancestor_list'])
            nodes[ancestor_id].add_child(node)
        else:
            root_nodes.append(node)

    # set up edge lengths
    if setup_edge_lengths:
        for id, parent_node in nodes.items():
            for child_node in parent_node.children:
                if child_node.edge_length is None and None not in (
                    getattr(child_node, 'origin_time', None),
                    getattr(parent_node, 'origin_time', None),
                ):
                    assert child_node.origin_time >= parent_node.origin_time
                    child_node.set_edge_length(
                        float(child_node.origin_time - parent_node.origin_time),
                    )

        for root_node in root_nodes:
            if (
                root_node.edge_length is None
                and getattr(root_node, 'origin_time', None) is not None
            ):
                assert not isanan(root_node.origin_time)
                root_node.set_edge_length(float(root_node.origin_time))

    res = [
        _treeswift_Tree_with_root(root_node)
        for root_node in root_nodes
    ]
    return res
