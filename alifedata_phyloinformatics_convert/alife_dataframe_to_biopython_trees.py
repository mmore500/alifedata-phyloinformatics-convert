from Bio import Phylo
from iterpop import iterpop as ip
from lyncs_utils import keydefaultdict
from nanto import isanan, nantonone
import pandas as pd
import string
import typing


def _parse_ancestor_list(raw: str) -> typing.List[int]:
    assert not any(
        c in raw for c in string.whitespace
    ), "Whitespace separated ancestor list not supported."
    try:
        ancestor_list = eval(raw)
        if ancestor_list == [None]:
            return []
        else:
            assert None not in ancestor_list
            return ancestor_list
    except NameError:
        # if ancestor list contains a placeholder string like NONE,
        # it will trigger a NameError when eval'ed
        return []


def alife_dataframe_to_biopython_trees(
    df: pd.DataFrame,
    setattrs: typing.Optional[typing.Union[
        typing.Iterable[str],
        typing.Mapping[str, str],
    ]] = None,
    setup_branch_lengths: bool = False,
):
    """Open a phylogeny dataframe formatted to the artificial life community
    data format standards as zero or more biopython trees, depending on the
    number of clades with no common ancestor.

    The following columns will automatically be applied as attributes to
    generated Clade objects:
        * branch_length,
        * id,
        * name, and
        * origin_time.

    Parameters
    ----------
    df:
        Pandas DataFrame to convert.
    setattrs: optional
        Dataframe columns that should be attached as attributes to Clade
        objects within the trees. If a map is provided, values at columns in
        keys will be attached with the corresponding value as the attr name.
    setup_branch_lengths: bool, optional
        Should we try to set up branch lengths using the origin_time column?
        Will not override if branch_length is provided as a column of df.
    """

    df = df.copy()

    df['parsed_ancestor_list'] \
        = df['ancestor_list'].apply(_parse_ancestor_list)

    # maps id to node
    def setup_node(id: int) -> Phylo.BaseTree.Clade:
        res = Phylo.BaseTree.Clade()
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
                    'name',
                    'branch_length',
                    'origin_time',
                )
                assert not hasattr(node, attr_name)
                setattr(node, attr_name, row[col_name])

        node.origin_time = nantonone(row.get('origin_time', default=None))
        node.name = row.get('name', default=None)
        node.branch_length = nantonone(row.get('branch_length', None))

        if row['parsed_ancestor_list']:
            ancestor_id = ip.popsingleton(row['parsed_ancestor_list'])
            nodes[ancestor_id].clades.append(node)
        else:
            root_nodes.append(node)

    # set up branch lengths
    if setup_branch_lengths:
        for id, parent_node in nodes.items():
            for child_node in parent_node.clades:
                if child_node.branch_length is None and None not in (
                    getattr(child_node, 'origin_time', None),
                    getattr(parent_node, 'origin_time', None),
                ):
                    assert child_node.origin_time >= parent_node.origin_time
                    child_node.branch_length \
                        = child_node.origin_time - parent_node.origin_time

        for root_node in root_nodes:
            if (
                root_node.branch_length is None
                and getattr(root_node, 'origin_time', None) is not None
            ):
                assert not isanan(root_node.origin_time)
                root_node.branch_length = root_node.origin_time

    return([
        Phylo.BaseTree.Tree(root=root_node)
        for root_node in root_nodes
    ])
