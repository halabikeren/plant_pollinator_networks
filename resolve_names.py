import logging
import sys
from typing import Optional
import click
import pandas as pd
from data_processing.name_resolver import NameResolutionMethod, NameResolver

logger = logging.getLogger(__name__)

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
    "--gnr_data_sources",
    help="data source to use for GNR",
    type=str,
    required=False,
    default=None, #"World Flora Online consortium",
)
@click.option(
    "--tnrs_data_source",
    help="data source to use for TNRS",
    type=click.Choice(["tpl","gcc","ildis","tropicos","usda","ncbi"], case_sensitive=True),
    required=False,
    default=None,
)
@click.option(
    "--batch_size",
    help="size of batches to run name resolution on",
    type=int,
    required=False,
    default=20,
)
@click.option(
    "--max_jobs_in_parallel",
    help="maximal number of jobs that should be executed in parallel",
    type=int,
    required=False,
    default=30, # for API request, this limitation will prevent too many simultaneous requests issue
)
def resolve_names(input_path: str, name_col: str, output_path: str, log_path: str, method: NameResolutionMethod, gnr_data_sources: Optional[str], tnrs_data_source: Optional[str], batch_size: int, max_jobs_in_parallel: int):

    logging.basicConfig(
        level=logging.INFO,
        format="%(asctime)s module: %(module)s function: %(funcName)s line %(lineno)d: %(message)s",
        handlers=[logging.StreamHandler(sys.stdout), logging.FileHandler(str(log_path)), ],
        force=True,  # run over root logger settings to enable simultaneous writing to both stdout and file handler
    )

    names = pd.read_csv(input_path)[name_col].str.capitalize().unique().tolist()
    name_resolver = NameResolver(method=method)
    name_resolver.resolve(names=names, output_path=output_path, batch_size=batch_size, gnr_data_sources=gnr_data_sources, tnrs_data_source=tnrs_data_source, max_jobs_in_parallel=max_jobs_in_parallel)


if __name__ == '__main__':
    resolve_names()
