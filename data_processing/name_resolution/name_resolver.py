import shutil
import sys
import os
import time

sys.path.append((os.path.dirname(os.path.dirname(os.path.dirname(__file__)))))
from services.pbs_service import PBSService

import glob
import logging
from enum import Enum
from typing import List, Optional
import numpy as np
import pandas as pd
from taxon_names_resolver import Resolver
from opentree import OTWebServiceWrapper

logger = logging.getLogger(__name__)

from dotenv import load_dotenv, find_dotenv
load_dotenv(find_dotenv())


class NameResolutionMethod(str, Enum):
    TNRS = "tnrs"  # https://tnrs.biendata.org/
    GNR = "gnr"  # https://resolver.globalnames.org/
    WFO = "wfo"  # https://github.com/rizqirizqi/scientific-name-scraper

class NameResolver:
    def __init__(self, method: NameResolutionMethod = "gnr", aux_dir: str = os.getcwd()):
        self.method = method
        self.aux_dir = aux_dir

    def _set_env(self) -> tuple[str, str, str]:
        timestamp = time.time()
        input_dir = f"{self.aux_dir}/{self.method}_query_names_{timestamp}/"
        output_dir = f"{self.aux_dir}/{self.method}_resolved_names_{timestamp}/"
        work_dir = f"{self.aux_dir}/{self.method}_name_resolution_jobs_{timestamp}/"
        os.makedirs(input_dir, exist_ok=True)
        os.makedirs(output_dir, exist_ok=True)
        os.makedirs(work_dir, exist_ok=True)
        logger.info(f"set up working environment in {work_dir} and output environment in {output_dir} with timestamp {timestamp}")
        return input_dir, output_dir, work_dir

    @staticmethod
    def _get_input_batches(items: List[str], write_dir: str, batch_size: int = 20):
        input_batches = [items[i: i + batch_size] for i in range(0, len(items), batch_size)]
        for i in range(len(input_batches)):
            input_path = f"{write_dir}/batch_{i}.csv"
            pd.Series(input_batches[i]).to_csv(input_path, header=["species_name"], index=False)
        logger.info(f"generated {len(input_batches)} input batches of size {batch_size}")

    def _exec_name_resolution(self, input_dir: str, output_dir: str, work_dir: str, gnr_data_sources: Optional[str] = None, tnrs_data_source: Optional[str] = None, max_jobs_in_parallel: int = 30, max_retries_num: int = 10):
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
                if self.method == NameResolutionMethod.GNR and gnr_data_sources is not None:
                    command += f" '{gnr_data_sources}'"
                elif self.method == NameResolutionMethod.TNRS and tnrs_data_source is not None:
                    command += f" '{tnrs_data_source}'"
                commands = [
                    f"cd {batch_output_dir}",
                    os.environ.get("CONDA_ACT_CMD", ""),
                    command
                ]
                job_commands.append(commands)
        PBSService.execute_job_array(work_dir=work_dir, output_dir=output_dir, jobs_commands=job_commands, max_parallel_jobs=max_jobs_in_parallel)
        logger.info(f"completed name resolution on {num_batches} batches, will now check for failed requests")

        # now check which batches failed and retry on them, until none of the batches failed
        failed_datasets_ids = NameResolver.get_failed_dataset_ids(log_paths=glob.glob(output_dir + f"*out", recursive=True))
        retry_num = 1
        while len(failed_datasets_ids) > 0 and retry_num < max_retries_num:
            retry_num += 1
            logger.info(f"retrying applying requests for {len(failed_datasets_ids)} batches")
            jobs_paths = [f"{work_dir}{dataset_id}.sh" for dataset_id in failed_datasets_ids]
            PBSService.execute_jobs(jobs_paths=jobs_paths, max_parallel_jobs=max_jobs_in_parallel)
            failed_datasets_ids = NameResolver.get_failed_dataset_ids(
                log_paths=glob.glob(output_dir + f"*out", recursive=True))

        logger.info(f"completed name resolution on retried batches, will now merge results")


    @staticmethod
    def _unite_batch_data(read_dir: str, output_path: str):
        dfs_paths =  glob.glob(read_dir + f"**/**/output.csv", recursive=True)
        dfs = [pd.read_csv(path) for path in dfs_paths]
        df = pd.concat(dfs)
        os.makedirs(os.path.dirname(output_path), exist_ok=True)
        df.to_csv(output_path)

    def resolve(self, names: List[str], output_path: str, batch_size: int = 20, gnr_data_sources: Optional[str] = None, tnrs_data_source: Optional[str] = None, max_jobs_in_parallel: int = 30):
        input_dir, output_dir, work_dir = self._set_env()
        self._get_input_batches(items=names, write_dir=input_dir, batch_size=batch_size)
        self._exec_name_resolution(input_dir=input_dir, output_dir=output_dir, work_dir=work_dir, gnr_data_sources=gnr_data_sources, tnrs_data_source=tnrs_data_source, max_jobs_in_parallel=max_jobs_in_parallel)
        self._unite_batch_data(read_dir=output_dir, output_path=output_path)
        # shutil.rmtree(input_dir, ignore_errors=True)
        # shutil.rmtree(output_dir, ignore_errors=True)

    @staticmethod
    def exec_tnrs(input_path: str, output_path: str, data_source: Optional[str] = None):
        if not os.path.exists(output_path):
            names = pd.read_csv(input_path).iloc[:, 0].to_list()
            ot = OTWebServiceWrapper(api_endpoint='production')
            response = ot.tnrs_match_names(names=names, do_approximate_matching=True, include_suppressed=True, data_source=data_source)
            dfs = []
            for result in response.response_dict['results']:
                if len(result['matches']) > 0:
                    sub_df = pd.concat([pd.json_normalize(match) for match in result['matches']])
                    dfs.append(sub_df)
            df = pd.concat(dfs)
            df.rename(columns={
                "search_string": "original_name",
                "is_approximate_match": "has_mismatches",
                "taxon.tax_sources": "name_sources",
                "taxon.rank": "taxon_rank",

            }, inplace=True)
            df["author"] = np.nan
            df["orig_status"] = df[["taxon.synonyms", "original_name"]].apply(lambda row: "SYNONYM" if row["original_name"].capitalize() in row["taxon.synonyms"] else "ACCEPTED", axis=1)
            df.sort_values(["original_name", "score"], ascending=[True, False], inplace=True)  # select matches with the best score
            df.drop_duplicates(subset=["original_name"], keep='first', inplace=True)
            df.to_csv(output_path)

    @staticmethod
    def exec_wfo(input_path: str, output_path: str):
        res = os.system(
            f"Rscript --vanilla {os.path.dirname(__file__)}/wfo_name_resolution.R --input_path={input_path} --output_path={output_path}")
        df = pd.read_csv(output_path).rename(
            columns={"spec.name.ORIG": "original_name",
                     "Fuzzy": "has_mismatches",
                     "scientificName": "resolved_name",
                     "taxonRank": "taxon_rank",
                     "references": "name_sources",
                     "scientificNameAuthorship": "author",
                     "Old.status": "orig_status"})
        df["score"] = 1 - (df["Fuzzy.dist"] / 100)
        df.dropna(subset=["original_name"], inplace=True)
        df.to_csv(output_path)

    @staticmethod
    def exec_gnr(input_path: str, output_path: str, gnr_data_sources: Optional[str] = None):
        base_output_path = f"{os.path.dirname(output_path)}/resolved_names/search_results.csv"
        if not os.path.exists(base_output_path):
            orig_dir = os.getcwd()
            os.chdir(os.path.dirname(output_path))
            if gnr_data_sources is not None:
                resolver = Resolver(input_path, datasource=gnr_data_sources.split(","))
            else:
                resolver = Resolver(input_path)
            resolver.main()  # to run the search
            resolver.write()  # to output the csv file
            os.chdir(orig_dir)
        df = pd.read_csv(base_output_path)
        df.rename(
            columns={"query_name": "original_name",
                     "name_string": "matched_name",
                     "current_name_string": "resolved_name",
                     "data_source_title": "name_sources",},
            inplace=True)
        df["taxon_rank"] = df["classification_path_ranks"].apply(lambda rank: rank.split("|")[-1] if pd.notna(rank) else np.nan)
        df.taxon_rank = df.taxon_rank.replace(
            {'spec.': 'species', 'gen.': 'genus', 'f.': 'form', 'var.': 'variety', 'subsp.': 'subspecies'})
        df["orig_status"] = df["resolved_name"].apply(lambda name: "SYNONYM" if pd.notna(name) else "ACCEPTED")
        df["resolved_name"].fillna(value=df["matched_name"].to_dict(), inplace=True)
        df["has_mismatches"] = np.where(df["original_name"] != df["matched_name"], True, False)
        df["author"] = np.nan
        df.to_csv(output_path)

    @staticmethod
    def get_failed_dataset_ids(log_paths: list[str]) -> list[str]:
        failed_str = "----- Returning nothing : GNR server may be down -----"
        failed_datasets_ids = []
        for path in log_paths:
            dataset_id = os.path.basename(path).replace(".out", "")
            with open(path, "r") as log_file:
                log_content = log_file.read()
            if failed_str in log_content:
                failed_datasets_ids.append(dataset_id)
        return failed_datasets_ids


if __name__ == '__main__':
    func_to_call = sys.argv[1]
    func = getattr(NameResolver, func_to_call)
    func(*sys.argv[2:])
