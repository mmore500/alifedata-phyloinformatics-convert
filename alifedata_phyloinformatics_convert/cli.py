import click
import dendropy
import pandas as pd

from .alife_dataframe_to_dendropy_tree import alife_dataframe_to_dendropy_tree
from .dendropy_tree_to_alife_dataframe import dendropy_tree_to_alife_dataframe


@click.group()
def cli():
    pass


@cli.command(
    help='convert standard alife phylogeny data to phloinformatics format',
)
@click.option(
    '--input-file',
    default='-',
    help='phyloinformatics data file path; default stdin',
    type=click.File('r'),
)
@click.option(
    '--input-schema',
    help=(
        'phyloinformatics data format schema; '
        'options include newick, nexml, and nexus'
    ),
    required=True,
)
@click.option(
    '--output-file',
    help='alife data file path; default stdout',
    default='-',
    type=click.File('w'),
)
@click.option(
    '--output-format',
    help='alife data file format; default csv',
    default='csv',
)
def toalifedata(
    input_file,
    input_schema,
    output_file,
    output_format,
):

    tree = dendropy.Tree.get(
        file=input_file,
        schema=input_schema,
    )

    converted_df = dendropy_tree_to_alife_dataframe(tree)

    {
        'csv': lambda file: converted_df.to_csv(file, index=False),
        'json': converted_df.to_json,
        'html': converted_df.to_html,
        'excel': converted_df.to_excel,
        'hdf': converted_df.to_hdf,
        'feather': converted_df.to_feather,
        'parquet': converted_df.to_parquet,
        'stata': converted_df.to_stata,
        'pickle': converted_df.to_pickle,
        'sql': converted_df.to_sql,
        'gbq': converted_df.to_gbq,
    }[output_format](
        output_file,
    )


@cli.command(
    help='convert phloinformatics data to standard alife phylogeny format',
)
@click.option(
    '--input-file',
    default='-',
    help='alife data file path; default stdin',
    type=click.File('r'),
)
@click.option(
    '--input-format',
    help='alife data file format; default csv',
    default='csv',
)
@click.option(
    '--output-file',
    help='phyloinformatics data file path; default stdout',
    default='-',
    type=click.File('w'),
)
@click.option(
    '--output-schema',
    help=(
        'phyloinformatics data format schema; '
        'options include newick, nexml, and nexus'
    ),
    required=True,
)
def fromalifedata(
    input_file,
    input_format,
    output_file,
    output_schema,
):
    df = {
        'csv': pd.read_csv,
        'fwf': pd.read_fwf,
        'json': pd.read_json,
        'html': pd.read_html,
        'excel': pd.read_excel,
        'hdf': pd.read_hdf,
        'feather': pd.read_feather,
        'parquet': pd.read_parquet,
        'orc': pd.read_orc,
        'stata': pd.read_stata,
        'sass': pd.read_sas,
        'spss': pd.read_spss,
        'pickle': pd.read_pickle,
        'sql': pd.read_sql,
        'gbq': pd.read_gbq,
    }[input_format](
        input_file,
    )

    converted_tree = alife_dataframe_to_dendropy_tree(
        df,
        setup_edge_lengths=True,
    )

    converted_tree.write(
        file=output_file,
        schema=output_schema,
    )
