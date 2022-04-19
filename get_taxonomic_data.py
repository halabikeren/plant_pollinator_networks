import os
import sys
from typing import List
import click
import numpy as np
import pandas as pd
import sqlite3
import logging

logger = logging.getLogger(__name__)

def get_db(db_link: str, db_dir: str):
    if not os.path.exists(f"{db_dir}/itisSqlite.zip"):
        res = os.system(f"cd {db_dir}; wget {db_link}; unzip itisSqlite.zip")
    db_path = f"{db_dir}/itisSqlite032822/ITIS.sqlite"
    return db_path

def get_plant_kingdom_id(db_connection: sqlite3.Connection) -> int:
    return pd.read_sql_query("SELECT kingdom_id from kingdoms WHERE kingdom_name LIKE 'Plant%' --case-insensitive;", db_connection).values[0][0]

def get_rank_data(db_connection: sqlite3.Connection, kingdom_id: int) -> pd.DataFrame:
    rank_data_query = f"SELECT rank_id, rank_name FROM taxon_unit_types WHERE kingdom_id IS {kingdom_id};"
    ranks_data = pd.read_sql_query(rank_data_query, db_connection)
    ranks_data.rank_name = ranks_data.rank_name.str.lower()
    return ranks_data


def get_taxonomic_data(db_connection: sqlite3.Connection, kingdom_id: int, names: List[str]):
    taxonomic_units_quey = f"SELECT * FROM taxonomic_units WHERE kingdom_id IS {kingdom_id}"
    taxonomic_units_data = pd.read_sql_query(taxonomic_units_quey, db_connection)
    taxonomic_units_data.complete_name = taxonomic_units_data.complete_name.str.lower()
    relevant_taxonomic_data = taxonomic_units_data.loc[taxonomic_units_data.complete_name.isin(names)]
    logger.info(f"% names covered by db taxonomic data = {np.round(relevant_taxonomic_data.shape[0]/len(names)*100, 2)}%")
    return relevant_taxonomic_data


def get_rank_name(item: pd.Series, rank_name: str, ranks_data: pd.DataFrame, taxonomic_data: pd.DataFrame) -> str:
    ordered_ranks = ranks_data.sort_values("rank_id")["rank_name"].to_list()
    tax_order_diff = ordered_ranks.index(rank_name) - ordered_ranks.index(item.rank_name)
    if tax_order_diff < 0: # the rank of item is higher than that of the required rank name, so it cannot be a subset of rank_name
        return np.nan
    while tax_order_diff > 0: # go up by one rank
        matches = taxonomic_data.loc[taxonomic_data.tsn == item.parent_tsn]
        if matches.shape[0] == 0:
            return np.nan
        item = matches.iloc[0]
        tax_order_diff = ordered_ranks.index(rank_name) - ordered_ranks.index(item.rank_name) # re-compute the taxonomic diff
    if tax_order_diff < 0:
        return np.nan
    return item.complete_name




@click.command()
@click.option(
    "--input_path",
    help="csv with query names to find their taxonomic data for",
    type=click.Path(exists=True, file_okay=True, readable=True),
    required=True,
)
@click.option(
    "--input_col",
    help="name of column consisting of names to cross with the db",
    type=str,
    required=False,
    default="matched_name",
)
@click.option(
    "--output_path",
    help="csv with query names and their their taxonomic data",
    type=click.Path(exists=False),
    required=True,
)
@click.option(
    "--log_path",
    help="path to the script's logger",
    type=click.Path(exists=False),
    required=True,
)
@click.option(
    "--itis_db_link",
    help="link to download the sqlite db from",
    type=str,
    required=False,
    default="https://www.itis.gov/downloads/itisSqlite.zip"
)
def get_missing_taxonomic_data(input_path: str, input_col:str, output_path: str, log_path: str, itis_db_link: str):

    logging.basicConfig(
        level=logging.INFO,
        format="%(asctime)s module: %(module)s function: %(funcName)s line %(lineno)d: %(message)s",
        handlers=[logging.StreamHandler(sys.stdout), logging.FileHandler(str(log_path)), ],
        force=True,  # run over root logger settings to enable simultaneous writing to both stdout and file handler
    )

    input_df = pd.read_csv(input_path)
    input_df[f"{input_col}_lowercase"] = input_df[input_col].str.lower()
    input_names = input_df[f"{input_col}_lowercase"].to_list()

    db_dir = os.path.dirname(output_path)
    db_path = get_db(db_link=itis_db_link, db_dir=db_dir)
    connection = sqlite3.connect(db_path)

    plant_kingdom_id = get_plant_kingdom_id(db_connection=connection)

    ranks_data = get_rank_data(db_connection=connection, kingdom_id=plant_kingdom_id)
    rank_id_to_name = ranks_data.set_index("rank_id")["rank_name"].to_dict()

    taxonomic_data = get_taxonomic_data(db_connection=connection, kingdom_id=plant_kingdom_id, names=input_names)
    taxonomic_data["rank_name"] = taxonomic_data["rank_id"].apply(lambda rank_id: rank_id_to_name[rank_id])
    taxonomic_data["genus_name"] = taxonomic_data.apply(lambda row: get_rank_name(item=row, rank_name="genus", ranks_data=ranks_data, taxonomic_data=taxonomic_data), axis=1)
    taxonomic_data["family_name"] = taxonomic_data.apply(lambda row: get_rank_name(item=row, rank_name="family", ranks_data=ranks_data, taxonomic_data=taxonomic_data), axis=1)

    input_df.set_index(f"{input_col}_lowercase", inplace=True)
    input_df["taxon_rank"].fillna(value=taxonomic_data.set_index("complete_name")["rank_name"].to_dict(), inplace=True)
    input_df["genus"].fillna(value=taxonomic_data.set_index("complete_name")["genus_name"].to_dict(), inplace=True)
    input_df["family"].fillna(value=taxonomic_data.set_index("complete_name")["family_name"].to_dict(), inplace=True)
    input_df.to_csv(output_path, index=False)

if __name__ == '__main__':
    get_missing_taxonomic_data()




