import logging
import sys
from typing import Optional
import click
import pandas as pd
from data_processing.name_resolver import NameResolutionMethod, NameResolver


@click.command()
@click.option(
    "--input_path",
    help="csv with query names to apply name resolution on",
    type=click.Path(exists=True, file_okay=True, readable=True),
    required=True,
)
@click.option(
    "--name_col",
    help="name of the field with names to resolve",
    type=str,
    required=False,
    default="species_name"
)
@click.option(
    "--output_path",
    help="csv with name resolution results on query names",
    type=click.Path(exists=False),
    required=True,
)
@click.option(
    "--log_path",
    help="path to log file",
    type=click.Path(exists=False),
    required=True,
)
@click.option(
    "--method",
    help="name resolution method to use",
    type=click.Choice(["gnr", "tnrs", "wfo"]),
    required=False,
    default="gnr",
)
@click.option(
    "--gnr_data_source",
    help="data source to use for GNR",
    type=str,
    required=False,
    default=None, #"World Flora Online consortium",
)
@click.option(
    "--batch_size",
    help="size of batches to run name resolution on",
    type=int,
    required=False,
    default=20,
)
def resolve_names(input_path: str, name_col: str, output_path: str, log_path: str, method: NameResolutionMethod, gnr_data_source: Optional[str], batch_size: int):
    # initialize the logger
    logging.basicConfig(
        level=logging.INFO,
        format="%(asctime)s module: %(module)s function: %(funcName)s line %(lineno)d: %(message)s",
        handlers=[logging.StreamHandler(sys.stdout), logging.FileHandler(str(log_path)), ],
        force=True,  # run over root logger settings to enable simultaneous writing to both stdout and file handler
    )

    names = pd.read_csv(input_path)[name_col].unique().tolist()
    name_resolver = NameResolver(method=method)
    name_resolver.resolve(names=names, output_path=output_path, batch_size=batch_size, gnr_data_source=gnr_data_source)


if __name__ == '__main__':
    resolve_names()
