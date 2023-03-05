#!/usr/bin/env python

'''
`alife_dataframe_to_dendropy_tree` tests for
`alifedata-phyloinformatics-convert` package.
'''

from click.testing import CliRunner
import filecmp
from os.path import dirname, getsize, realpath
import tempfile

from alifedata_phyloinformatics_convert import cli

import pytest

def test_toalifedata_nexml_csv():
    runner = CliRunner()
    scriptdir = dirname(realpath(__file__))
    with tempfile.TemporaryDirectory() as tempdir:
        result = runner.invoke(
            cli.toalifedata,
            f'--input-file {scriptdir}/assets/pythonidae.annotated.nexml '
            '--input-schema nexml '
            f'--output-file {tempdir}/pythonidae.annotated.csv '
            '--output-format csv '
        )
        assert result.exit_code == 0
        assert filecmp.cmp(
            f'{scriptdir}/converted_toalifedata/pythonidae.annotated.csv',
            f'{tempdir}/pythonidae.annotated.csv',
        )


def test_toalifedata_nexml_json():
    runner = CliRunner()
    scriptdir = dirname(realpath(__file__))
    with tempfile.TemporaryDirectory() as tempdir:
        result = runner.invoke(
            cli.toalifedata,
            f'--input-file {scriptdir}/assets/pythonidae.annotated.nexml '
            '--input-schema nexml '
            f'--output-file {tempdir}/pythonidae.annotated.json '
            '--output-format json '
        )
        assert result.exit_code == 0
        assert filecmp.cmp(
            f'{scriptdir}/converted_toalifedata/pythonidae.annotated.json',
            f'{tempdir}/pythonidae.annotated.json',
        )


def test_toalifedata_newick_csv():
    runner = CliRunner()
    scriptdir = dirname(realpath(__file__))
    with tempfile.TemporaryDirectory() as tempdir:
        result = runner.invoke(
            cli.toalifedata,
            f'--input-file {scriptdir}/assets/APG_Angiosperms.newick '
            '--input-schema newick '
            f'--output-file {tempdir}/APG_Angiosperms.csv '
            '--output-format csv '
        )
        assert result.exit_code == 0
        assert filecmp.cmp(
            f'{scriptdir}/converted_toalifedata/APG_Angiosperms.csv',
            f'{tempdir}/APG_Angiosperms.csv',
        )


def test_toalifedata_newick_json():
    runner = CliRunner()
    scriptdir = dirname(realpath(__file__))
    with tempfile.TemporaryDirectory() as tempdir:
        result = runner.invoke(
            cli.toalifedata,
            f'--input-file {scriptdir}/assets/APG_Angiosperms.newick '
            '--input-schema newick '
            f'--output-file {tempdir}/APG_Angiosperms.json '
            '--output-format json '
        )
        assert result.exit_code == 0
        assert filecmp.cmp(
            f'{scriptdir}/converted_toalifedata/APG_Angiosperms.json',
            f'{tempdir}/APG_Angiosperms.json',
        )


def test_fromalifedata_newick():
    runner = CliRunner()
    scriptdir = dirname(realpath(__file__))
    with tempfile.TemporaryDirectory() as tempdir:
        result = runner.invoke(
            cli.fromalifedata,
            f'--input-file {scriptdir}/assets/alifedata.csv '
            '--input-format csv '
            f'--output-file {tempdir}/alifedata.newick '
            '--output-schema newick '
        )
        assert result.exit_code == 0
        assert filecmp.cmp(
            f'{scriptdir}/converted_fromalifedata/keep-unifurcations/alifedata.newick',
            f'{tempdir}/alifedata.newick',
        )


def test_fromalifedata_minimal_newick():
    runner = CliRunner()
    scriptdir = dirname(realpath(__file__))
    with tempfile.TemporaryDirectory() as tempdir:
        result = runner.invoke(
            cli.fromalifedata,
            f'--input-file {scriptdir}/assets/alifedata_minimal.csv '
            '--input-format csv '
            f'--output-file {tempdir}/alifedata_minimal.newick '
            '--output-schema newick '
        )
        assert result.exit_code == 0
        assert filecmp.cmp(
            f'{scriptdir}/converted_fromalifedata/keep-unifurcations/alifedata_minimal.newick',
            f'{tempdir}/alifedata_minimal.newick',
        )


def test_fromalifedata_nexml():
    runner = CliRunner()
    scriptdir = dirname(realpath(__file__))
    with tempfile.TemporaryDirectory() as tempdir:
        result = runner.invoke(
            cli.fromalifedata,
            f'--input-file {scriptdir}/assets/alifedata.csv '
            '--input-format csv '
            f'--output-file {tempdir}/alifedata.nexml '
            '--output-schema nexml '
        )
        assert result.exit_code == 0
        assert filecmp.cmp(
            f'{scriptdir}/converted_fromalifedata/keep-unifurcations/alifedata.nexml',
            f'{tempdir}/alifedata.nexml',
        )


def test_fromalifedata_minimal_nexml():
    runner = CliRunner()
    scriptdir = dirname(realpath(__file__))
    with tempfile.TemporaryDirectory() as tempdir:
        result = runner.invoke(
            cli.fromalifedata,
            f'--input-file {scriptdir}/assets/alifedata_minimal.csv '
            '--input-format csv '
            f'--output-file {tempdir}/alifedata_minimal.nexml '
            '--output-schema nexml '
        )
        assert result.exit_code == 0
        assert filecmp.cmp(
            f'{scriptdir}/converted_fromalifedata/keep-unifurcations/alifedata_minimal.nexml',
            f'{tempdir}/alifedata_minimal.nexml',
        )


def test_fromalifedata_nexus():
    runner = CliRunner()
    scriptdir = dirname(realpath(__file__))
    with tempfile.TemporaryDirectory() as tempdir:
        result = runner.invoke(
            cli.fromalifedata,
            f'--input-file {scriptdir}/assets/alifedata.csv '
            '--input-format csv '
            f'--output-file {tempdir}/alifedata.nexus '
            '--output-schema nexus '
        )
        assert result.exit_code == 0
        assert filecmp.cmp(
            f'{scriptdir}/converted_fromalifedata/keep-unifurcations/alifedata.nexus',
            f'{tempdir}/alifedata.nexus',
        )


def test_fromalifedata_minimal_nexus():
    runner = CliRunner()
    scriptdir = dirname(realpath(__file__))
    with tempfile.TemporaryDirectory() as tempdir:
        result = runner.invoke(
            cli.fromalifedata,
            f'--input-file {scriptdir}/assets/alifedata_minimal.csv '
            '--input-format csv '
            f'--output-file {tempdir}/alifedata_minimal.nexus '
            '--output-schema nexus '
        )
        assert result.exit_code == 0
        assert filecmp.cmp(
            f'{scriptdir}/converted_fromalifedata/keep-unifurcations/alifedata_minimal.nexus',
            f'{tempdir}/alifedata_minimal.nexus',
        )


@pytest.mark.parametrize(
    "output_schema",
    [
        "nexml",
        "nexus",
        "newick",
    ],
)
def test_fromalifedata_keepsuppress_unifurcations(output_schema):
    runner = CliRunner()
    scriptdir = dirname(realpath(__file__))
    with tempfile.TemporaryDirectory() as tempdir:
        result = runner.invoke(
            cli.fromalifedata,
            f'--input-file {scriptdir}/assets/alifedata.csv '
            '--input-format csv '
            f'--output-file {tempdir}/keep_unifurcations '
            f'--output-schema {output_schema} '
            f'--keep-unifurcations'
        )
        assert result.exit_code == 0
        result = runner.invoke(
            cli.fromalifedata,
            f'--input-file {scriptdir}/assets/alifedata.csv '
            '--input-format csv '
            f'--output-file {tempdir}/suppress_unifurcations '
            f'--output-schema {output_schema} '
            f'--suppress-unifurcations'
        )
        assert result.exit_code == 0

        assert getsize(
            f'{tempdir}/keep_unifurcations'
        ) > getsize(
            f'{tempdir}/suppress_unifurcations'
        )
        for option in "keep", "suppress":
            assert filecmp.cmp(
                f"{scriptdir}/"
                "converted_fromalifedata/"
                f"{option}-unifurcations/"
                f"alifedata.{output_schema}",
                f'{tempdir}/{option}_unifurcations',
            )


@pytest.mark.parametrize(
    "input_schema",
    [
        "nexus",
        "newick",
    ],
)
@pytest.mark.parametrize(
    "output_format",
    [
        "csv",
        "json",
    ],
)
def test_toalifedata_keepsuppress_unifurcations(input_schema, output_format):
    runner = CliRunner()
    scriptdir = dirname(realpath(__file__))
    with tempfile.TemporaryDirectory() as tempdir:
        result = runner.invoke(
            cli.toalifedata,
            f'--input-file {scriptdir}/assets/alifedata.{input_schema} '
            f'--input-schema {input_schema} '
            f'--output-file {tempdir}/keep_unifurcations '
            f'--output-format {output_format} '
            f'--keep-unifurcations'
        )
        assert result.exit_code == 0
        result = runner.invoke(
            cli.toalifedata,
            f'--input-file {scriptdir}/assets/alifedata.{input_schema} '
            f'--input-schema {input_schema} '
            f'--output-file {tempdir}/suppress_unifurcations '
            f'--output-format {output_format} '
            f'--suppress-unifurcations'
        )
        assert result.exit_code == 0

        assert getsize(
            f'{tempdir}/keep_unifurcations'
        ) > getsize(
            f'{tempdir}/suppress_unifurcations'
        )
        for option in "keep", "suppress":
            assert filecmp.cmp(
                f"{scriptdir}/"
                "converted_toalifedata/"
                f"{option}-unifurcations/"
                f"alifedata.{output_format}",
                f'{tempdir}/{option}_unifurcations',
            )
