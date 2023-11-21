#!/usr/bin/env python

"""
`RosettaTree` tests for `alifedata-phyloinformatics-convert` package.
"""

from contextlib import redirect_stdout
import io
from os.path import dirname, realpath
import pathlib
import tempfile

import anytree
from Bio import Phylo as BioPhylo
import dendropy as dp
import networkx as nx
import pandas as pd
import pytest

import alifedata_phyloinformatics_convert as apc


@pytest.mark.parametrize(
    "original_df",
    [
        pd.read_csv(f"{dirname(realpath(__file__))}/assets/alifedata.csv"),
        pd.read_csv(
            f"{dirname(realpath(__file__))}/assets/alifedata_minimal.csv",
        ),
    ],
)
def test_alife_to_biopython(original_df):
    expected_tree = apc.alife_dataframe_to_biopython_tree(
        original_df, setup_branch_lengths=True
    )
    rosetta_tree = apc.RosettaTree(original_df)

    # twice to test caching
    for __ in range(2):
        converted_tree = rosetta_tree.as_biopython
        assert str(converted_tree) == str(expected_tree)


@pytest.mark.parametrize(
    "original_df",
    [
        pd.read_csv(f"{dirname(realpath(__file__))}/assets/alifedata.csv"),
        pd.read_csv(
            f"{dirname(realpath(__file__))}/assets/alifedata_minimal.csv",
        ),
    ],
)
def test_alife_to_dendropy(original_df):
    expected_tree = apc.alife_dataframe_to_dendropy_tree(
        original_df, setup_edge_lengths=True
    )
    rosetta_tree = apc.RosettaTree(original_df)

    # twice to test caching
    for __ in range(2):
        converted_tree = rosetta_tree.as_dendropy
        assert str(converted_tree) == str(expected_tree)


@pytest.mark.parametrize(
    "original_df",
    [
        pd.read_csv(f"{dirname(realpath(__file__))}/assets/alifedata.csv"),
        pd.read_csv(
            f"{dirname(realpath(__file__))}/assets/alifedata_minimal.csv",
        ),
    ],
)
def test_alife_to_ete(original_df):
    expected_tree = apc.alife_dataframe_to_ete_tree(
        original_df, setup_dists=True
    )
    rosetta_tree = apc.RosettaTree(original_df)

    # twice to test caching
    for __ in range(2):
        converted_tree = rosetta_tree.as_ete
        assert converted_tree.get_ascii() == expected_tree.get_ascii()


@pytest.mark.parametrize(
    "original_df",
    [
        pd.read_csv(f"{dirname(realpath(__file__))}/assets/alifedata.csv"),
        pd.read_csv(
            f"{dirname(realpath(__file__))}/assets/alifedata_minimal.csv",
        ),
    ],
)
def test_alife_to_networkx(original_df):
    expected_tree = apc.alife_dataframe_to_networkx_digraph(
        original_df, setup_edge_lengths=True
    )
    rosetta_tree = apc.RosettaTree(original_df)

    # twice to test caching
    for __ in range(2):
        converted_tree = rosetta_tree.as_networkx
        assert "\n".join(
            nx.generate_gml(converted_tree)
        ) == "\n".join(
            nx.generate_gml(expected_tree)
        )


@pytest.mark.parametrize(
    "original_df",
    [
        pd.read_csv(f"{dirname(realpath(__file__))}/assets/alifedata.csv"),
        pd.read_csv(
            f"{dirname(realpath(__file__))}/assets/alifedata_minimal.csv",
        ),
    ],
)
def test_alife_to_phylotrack(original_df):
    expected_tree = apc.alife_dataframe_to_phylotrack_systematics(original_df)
    with io.StringIO() as buf, redirect_stdout(buf):
        expected_tree.print_status()
        expected_string = buf.getvalue()

    rosetta_tree = apc.RosettaTree(original_df)

    # twice to test caching
    for __ in range(2):
        converted_tree = rosetta_tree.as_phylotrack

        with io.StringIO() as buf, redirect_stdout(buf):
            converted_tree.print_status()
            converted_string = buf.getvalue()

        assert converted_string == expected_string


@pytest.mark.parametrize(
    "original_df",
    [
        pd.read_csv(f"{dirname(realpath(__file__))}/assets/alifedata.csv"),
        pd.read_csv(
            f"{dirname(realpath(__file__))}/assets/alifedata_minimal.csv",
        ),
    ],
)
def test_biopython_to_alife(original_df):
    original_tree = apc.alife_dataframe_to_biopython_tree(
        original_df, setup_branch_lengths=True
    )
    expected_df = apc.biopython_tree_to_alife_dataframe(
        original_tree, {"name": "taxon_label"}
    )
    rosetta_tree = apc.RosettaTree(original_tree)

    # twice to test caching
    for __ in range(2):
        converted_df = rosetta_tree.as_alife
        assert converted_df.equals(expected_df)


@pytest.mark.parametrize(
    "original_df",
    [
        pd.read_csv(f"{dirname(realpath(__file__))}/assets/alifedata.csv"),
        pd.read_csv(
            f"{dirname(realpath(__file__))}/assets/alifedata_minimal.csv",
        ),
    ],
)
def test_biopython_to_dendropy(original_df):
    original_tree = apc.alife_dataframe_to_biopython_tree(
        original_df, setup_branch_lengths=True
    )
    expected_tree = apc.alife_dataframe_to_dendropy_tree(
        original_df, setup_edge_lengths=True
    )
    rosetta_tree = apc.RosettaTree(original_tree)

    # twice to test caching
    for __ in range(2):
        converted_tree = rosetta_tree.as_dendropy
        assert str(converted_tree) == str(expected_tree)


@pytest.mark.parametrize(
    "original_df",
    [
        pd.read_csv(f"{dirname(realpath(__file__))}/assets/alifedata.csv"),
        pd.read_csv(
            f"{dirname(realpath(__file__))}/assets/alifedata_minimal.csv",
        ),
    ],
)
def test_dendropy_to_alife(original_df):
    original_tree = apc.alife_dataframe_to_dendropy_tree(
        original_df, setup_edge_lengths=True
    )
    expected_df = apc.dendropy_tree_to_alife_dataframe(original_tree)
    rosetta_tree = apc.RosettaTree(original_tree)

    # twice to test caching
    for __ in range(2):
        converted_df = rosetta_tree.as_alife
        assert converted_df.equals(expected_df)


@pytest.mark.parametrize(
    "original_df",
    [
        pd.read_csv(f"{dirname(realpath(__file__))}/assets/alifedata.csv"),
        pd.read_csv(
            f"{dirname(realpath(__file__))}/assets/alifedata_minimal.csv",
        ),
    ],
)
def test_ete_to_alife(original_df):
    original_tree = apc.alife_dataframe_to_ete_tree(
        original_df, setup_dists=True
    )
    expected_df = apc.ete_tree_to_alife_dataframe(original_tree)
    rosetta_tree = apc.RosettaTree(original_tree)

    # twice to test caching
    for __ in range(2):
        converted_df = rosetta_tree.as_alife
        assert converted_df.equals(expected_df)


@pytest.mark.parametrize(
    "original_df",
    [
        pd.read_csv(f"{dirname(realpath(__file__))}/assets/alifedata.csv"),
        pd.read_csv(
            f"{dirname(realpath(__file__))}/assets/alifedata_minimal.csv",
        ),
    ],
)
def test_dendropy_to_biopython(original_df):
    original_tree = apc.alife_dataframe_to_dendropy_tree(
        original_df, setup_edge_lengths=True
    )
    expected_tree = apc.alife_dataframe_to_biopython_tree(
        original_df, setup_branch_lengths=True
    )
    rosetta_tree = apc.RosettaTree(original_tree)

    # twice to test caching
    for __ in range(2):
        converted_tree = rosetta_tree.as_biopython
        assert str(converted_tree) == str(expected_tree)



@pytest.mark.parametrize(
    "original_df",
    [
        pd.read_csv(f"{dirname(realpath(__file__))}/assets/alifedata.csv"),
        pd.read_csv(
            f"{dirname(realpath(__file__))}/assets/alifedata_minimal.csv",
        ),
    ],
)
def test_networkx_to_alife(original_df):
    original_tree = apc.alife_dataframe_to_networkx_digraph(
        original_df, setup_edge_lengths=True
    )
    expected_df = apc.networkx_digraph_to_alife_dataframe(original_tree)
    rosetta_tree = apc.RosettaTree(original_tree)

    # twice to test caching
    for __ in range(2):
        converted_df = rosetta_tree.as_alife
        assert converted_df.equals(expected_df)


@pytest.mark.parametrize(
    "original_df",
    [
        pd.read_csv(f"{dirname(realpath(__file__))}/assets/alifedata.csv"),
        pd.read_csv(
            f"{dirname(realpath(__file__))}/assets/alifedata_minimal.csv",
        ),
    ],
)
def test_networkx_to_phylotrack(original_df):
    original_tree = apc.alife_dataframe_to_phylotrack_systematics(original_df)
    expected_df = apc.phylotrack_systematics_to_alife_dataframe(original_tree)
    rosetta_tree = apc.RosettaTree(original_tree)

    # twice to test caching
    for __ in range(2):
        converted_df = rosetta_tree.as_alife
        assert converted_df.equals(expected_df)


@pytest.mark.parametrize(
    "original_df",
    [
        pd.read_csv(f"{dirname(realpath(__file__))}/assets/alifedata.csv"),
        pd.read_csv(
            f"{dirname(realpath(__file__))}/assets/alifedata_minimal.csv",
        ),
    ],
)
def test_networkx_to_biopython(original_df):
    original_tree = apc.alife_dataframe_to_networkx_digraph(
        original_df, setup_edge_lengths=True
    )
    expected_tree = apc.alife_dataframe_to_biopython_tree(
        original_df, setup_branch_lengths=True
    )
    rosetta_tree = apc.RosettaTree(original_tree)

    # twice to test caching
    for __ in range(2):
        converted_tree = rosetta_tree.as_biopython
        assert str(converted_tree) == str(expected_tree)


@pytest.mark.parametrize(
    "original_df",
    [
        pd.read_csv(f"{dirname(realpath(__file__))}/assets/alifedata.csv"),
        pd.read_csv(
            f"{dirname(realpath(__file__))}/assets/alifedata_minimal.csv",
        ),
    ],
)
def test_networkx_to_dendropy(original_df):
    original_tree = apc.alife_dataframe_to_networkx_digraph(
        original_df, setup_edge_lengths=True
    )
    expected_tree = apc.alife_dataframe_to_dendropy_tree(
        original_df, setup_edge_lengths=True
    )
    rosetta_tree = apc.RosettaTree(original_tree)

    # twice to test caching
    for __ in range(2):
        converted_tree = rosetta_tree.as_dendropy
        assert str(converted_tree) == str(expected_tree)


@pytest.mark.parametrize(
    "original_df",
    [
        pd.read_csv(f"{dirname(realpath(__file__))}/assets/alifedata.csv"),
        pd.read_csv(
            f"{dirname(realpath(__file__))}/assets/alifedata_minimal.csv",
        ),
    ],
)
def test_dendropy_to_newick(original_df):
    expected_tree = apc.alife_dataframe_to_dendropy_tree(
        original_df, setup_edge_lengths=True
    ).as_string(schema="newick")
    rosetta_tree = apc.RosettaTree(original_df)

    # twice to test caching
    for __ in range(2):
        converted_tree = rosetta_tree.as_newick
        assert converted_tree == expected_tree


@pytest.mark.parametrize(
    "original_df",
    [
        pd.read_csv(f"{dirname(realpath(__file__))}/assets/alifedata.csv"),
        pd.read_csv(
            f"{dirname(realpath(__file__))}/assets/alifedata_minimal.csv",
        ),
        pd.read_csv(
            f"{dirname(realpath(__file__))}/assets/alifedata_empty.csv",
        ),
    ],
)
def test_dendropy_as_newick(original_df):
    source_tree = apc.alife_dataframe_to_dendropy_tree(
        original_df, setup_edge_lengths=True
    )
    expected_tree = (
        source_tree.as_string(schema="newick")
        if source_tree is not None
        else None
    )
    rosetta_tree = apc.RosettaTree(original_df)

    # twice to test caching
    for __ in range(2):
        converted_tree = rosetta_tree.as_newick
        assert converted_tree == expected_tree


@pytest.mark.parametrize(
    "original_df",
    [
        pd.read_csv(f"{dirname(realpath(__file__))}/assets/alifedata.csv"),
        pd.read_csv(
            f"{dirname(realpath(__file__))}/assets/alifedata_minimal.csv",
        ),
        pd.read_csv(
            f"{dirname(realpath(__file__))}/assets/alifedata_empty.csv",
        ),
    ],
)
def test_to_newick(original_df):
    source_tree = apc.alife_dataframe_to_dendropy_tree(
        original_df, setup_edge_lengths=True
    )
    expected_tree = (
        source_tree.as_string(schema="newick")
        if source_tree is not None
        else None
    )
    rosetta_tree = apc.RosettaTree(original_df)

    converted_tree = rosetta_tree.to_newick()
    assert converted_tree == expected_tree

    with tempfile.NamedTemporaryFile("w") as fp:
        if expected_tree is not None:
            rosetta_tree.to_newick(fp)
            fp.flush()
            assert (
                pathlib.Path(fp.name).read_text().strip()
                == expected_tree.strip()
            )
        else:
            with pytest.raises(ValueError):
                rosetta_tree.to_newick(fp)

    with tempfile.NamedTemporaryFile("w") as fp:
        if expected_tree is not None:
            rosetta_tree.to_newick(fp.name)
            assert (
                pathlib.Path(fp.name).read_text().strip()
                == expected_tree.strip()
            )
        else:
            with pytest.raises(ValueError):
                rosetta_tree.to_newick(fp.name)


@pytest.mark.parametrize(
    "original_df",
    [
        pd.read_csv(f"{dirname(realpath(__file__))}/assets/alifedata.csv"),
        pd.read_csv(
            f"{dirname(realpath(__file__))}/assets/alifedata_minimal.csv",
        ),
        pd.read_csv(
            f"{dirname(realpath(__file__))}/assets/alifedata_empty.csv",
        ),
    ],
)
def test_to_nexus(original_df):
    source_tree = apc.alife_dataframe_to_dendropy_tree(
        original_df, setup_edge_lengths=True
    )
    expected_tree = (
        source_tree.as_string(schema="nexus")
        if source_tree is not None
        else None
    )
    rosetta_tree = apc.RosettaTree(original_df)

    converted_tree = rosetta_tree.to_nexus()
    assert converted_tree == expected_tree

    with tempfile.NamedTemporaryFile("w") as fp:
        if expected_tree is not None:
            rosetta_tree.to_nexus(fp)
            fp.flush()
            assert (
                pathlib.Path(fp.name).read_text().strip()
                == expected_tree.strip()
            )
        else:
            with pytest.raises(ValueError):
                rosetta_tree.to_nexus(fp)

    with tempfile.NamedTemporaryFile("w") as fp:
        if expected_tree is not None:
            rosetta_tree.to_nexus(fp.name)
            assert (
                pathlib.Path(fp.name).read_text().strip()
                == expected_tree.strip()
            )
        else:
            with pytest.raises(ValueError):
                rosetta_tree.to_nexus(fp.name)


@pytest.mark.parametrize(
    "original_df",
    [
        pd.read_csv(f"{dirname(realpath(__file__))}/assets/alifedata.csv"),
        pd.read_csv(
            f"{dirname(realpath(__file__))}/assets/alifedata_minimal.csv",
        ),
        pd.read_csv(
            f"{dirname(realpath(__file__))}/assets/alifedata_empty.csv",
        ),
    ],
)
def test_to_nexml(original_df):
    source_tree = apc.alife_dataframe_to_dendropy_tree(
        original_df, setup_edge_lengths=True
    )
    expected_tree = (
        source_tree.as_string(schema="nexml")
        if source_tree is not None
        else None
    )
    rosetta_tree = apc.RosettaTree(original_df)

    converted_tree = rosetta_tree.to_nexml()
    assert converted_tree == expected_tree

    with tempfile.NamedTemporaryFile("w") as fp:
        if expected_tree is not None:
            rosetta_tree.to_nexml(fp)
            fp.flush()
            assert (
                pathlib.Path(fp.name).read_text().strip()
                == expected_tree.strip()
            )
        else:
            with pytest.raises(ValueError):
                rosetta_tree.to_nexml(fp)

    with tempfile.NamedTemporaryFile("w") as fp:
        if expected_tree is not None:
            rosetta_tree.to_nexml(fp.name)
            assert (
                pathlib.Path(fp.name).read_text().strip()
                == expected_tree.strip()
            )
        else:
            with pytest.raises(ValueError):
                rosetta_tree.to_nexml(fp.name)
