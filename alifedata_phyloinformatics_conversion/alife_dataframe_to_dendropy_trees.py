import dendropy
from iterpop import iterpop as ip
from lyncs_utils import keydefaultdict
import math
import pandas as pd
import string
import typing


def alife_dataframe_to_dendropy_trees(
    df: pd.DataFrame,
    setup_edge_lengths: bool = False,
) -> typing.List[dendropy.Tree]:

    # maps id to node
    def setup_node(id: int) -> dendropy.Node:
        res = dendropy.Node()
        res.id = id
        return res
    nodes = keydefaultdict(setup_node)

    # maps id to origin time
    root_nodes = []

    for __, row in df.iterrows():
        node = nodes[row['id']]

        node.origin_time = row.get('origin_time', default=None)
        if (
            isinstance(node.origin_time, float)
            and math.isnan(node.origin_time)
        ):
            node.origin_time = None
        node.label = row.get('label', default=None)
        taxon_label = row.get('taxon_label', default=None)
        if taxon_label not in ('None', None):
            node.taxon = dendropy.Taxon(label=taxon_label)
        if 'edge_length' in row and not (
            isinstance(row['edge_length'], float)
            and math.isnan(row['edge_length'])
        ):
            node.edge_length = row['edge_length']
        try:
            ancestor_list = eval(row['ancestor_list'])
            assert len(ancestor_list) < 2, "Sexual lineages not supported."
            assert not any(c in ancestor_list for c in string.whitespace), \
                "Whitespace separated ancestor list not supported."
            if ip.poursingleton(ancestor_list) is not None:
                ancestor_id = ip.popsingleton(ancestor_list)
                nodes[ancestor_id].add_child(node)
            else:
                # ancestor_list is empty... trigger exeption handler below
                # to execute root node case logic
                raise NameError
        except NameError:
            # if ancestor list contains a placeholder string like NONE,
            # it will trigger a NameError when eval'ed
            root_nodes.append(node)

    # set up edge lengths
    if setup_edge_lengths:
        for id, parent_node in nodes.items():
            for child_node in parent_node.child_node_iter():
                if child_node.edge_length is None and None not in (
                    getattr(child_node, 'origin_time', None),
                    getattr(parent_node, 'origin_time', None),
                ):
                    assert child_node.origin_time >= parent_node.origin_time
                    child_node.edge_length \
                        = child_node.origin_time - parent_node.origin_time

        for root_node in root_nodes:
            if (
                root_node.edge_length is None
                and getattr(root_node, 'origin_time', None) is not None
            ):
                assert not (
                    isinstance(root_node.origin_time, float)
                    and math.isnan(root_node.origin_time)
                )
                root_node.edge_length = root_node.origin_time

    return([
        dendropy.Tree(seed_node=root_node)
        for root_node in root_nodes
    ])
