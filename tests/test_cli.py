#!/usr/bin/env python

'''
`alife_dataframe_to_dendropy_tree` tests for
`alifedata-phyloinformatics-convert` package.
'''

from click.testing import CliRunner
import filecmp
from os.path import dirname, realpath
import tempfile

from alifedata_phyloinformatics_convert import cli


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
            f'{scriptdir}/converted_fromalifedata/alifedata.newick',
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
            f'{scriptdir}/converted_fromalifedata/alifedata_minimal.newick',
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
            f'{scriptdir}/converted_fromalifedata/alifedata.nexml',
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
            f'{scriptdir}/converted_fromalifedata/alifedata_minimal.nexml',
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
            f'{scriptdir}/converted_fromalifedata/alifedata.nexus',
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
            f'{scriptdir}/converted_fromalifedata/alifedata_minimal.nexus',
            f'{tempdir}/alifedata_minimal.nexus',
        )
