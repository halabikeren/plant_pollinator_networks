import logging
import os
import re
from enum import Enum
from typing import List

import numpy as np
import pandas as pd
from taxon_names_resolver import Resolver
from opentree import OTWebServiceWrapper
from flatten_json import flatten

from services import PBSService

logger = logging.getLogger(__name__)


class NameResolutionMethod(str, Enum):
    TNRS = "tnrs"  # https://tnrs.biendata.org/
    GNR = "gnr"  # https://resolver.globalnames.org/
    WFO = "wfo"  # https://github.com/rizqirizqi/scientific-name-scraper


name_resolution_result_fields = ["query", "matched_name", "has_mismatches", "score", "taxon_rank", "name_sources"]


class NameResolver:

    def __init__(self, method: NameResolutionMethod = "gnr", aux_dir: str = os.getcwd()):
        self.method = method
        self.aux_dir = aux_dir

    def resolve(self, names: List[str], output_path: str, batch_size: int = 20):

        input_dir = f"{self.aux_dir}/query_names/"
        output_dir = f"{self.aux_dir}/resolved_names/"
        work_dir = f"{self.aux_dir}/name_resolution_jobs/"
        os.makedirs(input_dir, exist_ok=True)
        os.makedirs(output_dir, exist_ok=True)
        os.makedirs(work_dir, exist_ok=True)
        logger.info(f"set up working environment in {work_dir} and output environment in {output_dir}")

        input_batches = [names[i: i + batch_size] for i in range(0, len(names), batch_size)]
        for i in range(len(input_batches)):
            input_path = f"{input_dir}/batch_{i}.csv"
            pd.Series(input_batches[i]).to_csv(input_path, header=["species_name"])
        logger.info(f"generated {len(input_batches)} input batches of size {batch_size}")

        job_commands = []
        for path in os.listdir(input_dir):
            input_path = f"'{input_dir}{path}'"
            batch_output_dir = f"{output_dir}{path.split('.')[0]}"
            os.makedirs(batch_output_dir, exist_ok=True)
            output_path = f"'{batch_output_dir}/output.csv'"
            parent_path = f"'{os.path.abspath(os.curdir)}'"
            if not os.path.exists(f"{batch_output_dir}/output.csv"):
                command = f'''python -c "import sys;sys.path.append({parent_path});from data_processing' 
                          import NameResolver;NameResolver.exec_{self.method}(' 
                          input_path={input_path}, output_path={output_path})"'''
                commands = [
                    f"cd {batch_output_dir}",
                    command
                ]
                job_commands.append(commands)
        PBSService.execute_job_array(work_dir=work_dir, output_dir=output_dir, jobs_commands=job_commands)
        logger.info(f"completed name resolution on {len(input_batches)} batches, will now merge results")

        output_dfs = [pd.read_csv(f"{output_dir}/{path}") for path in os.listdir(output_dir)]
        df = pd.concat(output_dfs)[name_resolution_result_fields]
        df.to_csv(output_path)

    @staticmethod
    def exec_tnrs(input_path: str, output_path: str):
        names = pd.read_csv(input_path).iloc[:, 1].to_list()
        ot = OTWebServiceWrapper(api_endpoint='production')
        response = ot.tnrs_match_names(names=names, do_approximate_matching=True, include_suppressed=True)
        dfs = []
        for result in response.response_dict['results']:
            if len(result['matches']) > 0:
                sub_df = [pd.json_normalize(match) for match in result['matches']]
                dfs.append(sub_df)
        df = pd.concat(dfs)
        df.rename(columns={
            "search_string": "query",
            "is_approximate_match": "has_mismatches",
            "taxon_tax_sources": "name_sources"
        }, inplace=True)
        df.sort_values(["query", "score"], ascending=[True, False], inplace=True)  # select matches with the best score
        df.drop_duplicates(subset=["query"], keep='first', inplace=True)
        df.to_csv(output_path)

    @staticmethod
    def exec_wfo(input_path: str, output_path: str):
        res = os.system(
            f"Rscript --vanilla {os.path.dirname(__file__)}/wfo_name_resolution.R --input_path={input_path} --output_path={output_path}")
        df = pd.read_csv(output_path).rename(
            columns={"spec.name.ORIG": "query",
                     "Fuzzy": "has_mismatches",
                     "scientificName": "matched_name",
                     "taxonRank": "taxon_rank",
                     "references": "name_sources"})
        df["score"] = 1 - (df["Fuzzy.dist"] / 100)
        df.dropna(subset=["query"], inplace=True)
        df.to_csv()

    @staticmethod
    def exec_gnr(input_path: str, output_path: str):
        resolver = Resolver(input_path)
        resolver.main()  # to run the search
        resolver.write()  # to output the csv file
        df = pd.read_csv(output_path)
        df.rename(
            columns={"query_name": "query", "data_source_title": "name_sources", "name_string": "matched_name"},
            inplace=True)
        df["taxon_rank"] = np.nan
        df["has_mismatches"] = np.where(["query"] != df["matched_name"], True, False)
        df.to_csv(output_path)
