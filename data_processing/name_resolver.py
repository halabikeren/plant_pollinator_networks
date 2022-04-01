import shutil
import sys
import os
sys.path.append(os.path.dirname(os.path.dirname(__file__)))

import glob
import logging
from enum import Enum
from typing import List
import numpy as np
import pandas as pd
from taxon_names_resolver import Resolver
from opentree import OTWebServiceWrapper

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

    def _set_env(self) -> tuple[str, str, str]:
        input_dir = f"{self.aux_dir}/{self.method}_query_names/"
        output_dir = f"{self.aux_dir}/{self.method}_resolved_names/"
        work_dir = f"{self.aux_dir}/{self.method}_name_resolution_jobs/"
        os.makedirs(input_dir, exist_ok=True)
        os.makedirs(output_dir, exist_ok=True)
        os.makedirs(work_dir, exist_ok=True)
        logger.info(f"set up working environment in {work_dir} and output environment in {output_dir}")
        return input_dir, output_dir, work_dir

    @staticmethod
    def _get_input_batches(items: List[str], write_dir: str, batch_size: int = 200):
        input_batches = [items[i: i + batch_size] for i in range(0, len(items), batch_size)]
        for i in range(len(input_batches)):
            input_path = f"{write_dir}/batch_{i}.csv"
            pd.Series(input_batches[i]).to_csv(input_path, header=["species_name"])
        logger.info(f"generated {len(input_batches)} input batches of size {batch_size}")

    def _exec_name_resolution(self, input_dir: str, output_dir: str, work_dir: str):
        job_commands = []
        num_batches = 0
        for path in os.listdir(input_dir):
            num_batches += 1
            input_path = f"'{input_dir}{path}'"
            batch_output_dir = f"{output_dir}{path.split('.')[0]}"
            os.makedirs(batch_output_dir, exist_ok=True)
            output_path = f"'{batch_output_dir}/output.csv'"
            if not os.path.exists(f"{batch_output_dir}/output.csv"):
                func_name = f"exec_{self.method}"
                command = f"python {os.path.abspath(__file__)} {func_name} {input_path} {output_path}"
                commands = [
                    f"cd {batch_output_dir}",
                    command
                ]
                job_commands.append(commands)
        PBSService.execute_job_array(work_dir=work_dir, output_dir=output_dir, jobs_commands=job_commands)
        logger.info(f"completed name resolution on {num_batches} batches, will now merge results")

    @staticmethod
    def _unite_batch_data(read_dir: str, output_path: str):
        dfs_paths = glob.glob(read_dir + f"**/*csv", recursive=True)
        dfs = [pd.read_csv(path) for path in dfs_paths]
        df = pd.concat(dfs)[name_resolution_result_fields]
        df.to_csv(output_path)

    def resolve(self, names: List[str], output_path: str, batch_size: int = 200):
        input_dir, output_dir, work_dir = self._set_env()
        self._get_input_batches(items=names, write_dir=input_dir, batch_size=batch_size)
        self._exec_name_resolution(input_dir=input_dir, output_dir=output_dir, work_dir=work_dir)
        self._unite_batch_data(read_dir=output_dir, output_path=output_path)
        shutil.rmtree(input_dir, ignore_errors=True)
        shutil.rmtree(output_dir, ignore_errors=True)

    @staticmethod
    def exec_tnrs(input_path: str, output_path: str):
        names = pd.read_csv(input_path).iloc[:, 1].to_list()
        ot = OTWebServiceWrapper(api_endpoint='production')
        response = ot.tnrs_match_names(names=names, do_approximate_matching=True, include_suppressed=True)
        dfs = []
        for result in response.response_dict['results']:
            if len(result['matches']) > 0:
                sub_df = pd.concat([pd.json_normalize(match) for match in result['matches']])
                dfs.append(sub_df)
        df = pd.concat(dfs)
        df.rename(columns={
            "search_string": "query",
            "is_approximate_match": "has_mismatches",
            "taxon.tax_sources": "name_sources",
            "taxon.rank": "taxon_rank"
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


if __name__ == '__main__':
    func_to_call = sys.argv[1]
    func = getattr(NameResolver, func_to_call)
    func(*sys.argv[2:])
